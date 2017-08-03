function [fwhm, fi, fx, fy, fz] = resestsmooth(res, mmpv, opts)
% resestsmooth  - estimate smoothness from regression residual
%
% FORMAT:       [fwhm, fi] = resestsmooth(res, mmpv [, opts])
%
% Input fields:
%
%       res         regression residual
%       mmpv        mm per voxel (voxel size)
%       opts        optional settings struct
%        .tdim      temporal dimension (default: last)
%
% Output fields:
%
%       fwhm        [X, Y, Z] overall FWHM estimate
%       fi          smoothness image
%
% Note: this function uses ideas taken from a FMRIB webpage at
%       http://www.fmrib.ox.ac.uk/analysis/techrep/tr00df1/tr00df1/node6.html

% Version:  v1.0
% Build:    14110322
% Date:     Nov-03 2014, 10:44 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% argument check
if nargin < 2 || ...
   ~isnumeric(res) || ...
    ndims(res) ~= 4 || ...
   ~isa(mmpv, 'double') || ...
   (numel(mmpv) ~= 1 && ...
    numel(mmpv) ~= 3) || ...
    any(isinf(mmpv) | isnan(mmpv) | mmpv <= 0)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if numel(mmpv) == 1
    mmpv = mmpv([1, 1, 1]);
end
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'tdim') || ...
   ~isa(opts.tdim, 'double') || ...
    numel(opts.tdim) ~= 1 || ...
    isinf(opts.tdim) || ...
    isnan(opts.tdim) || ...
   ~any(opts.tdim == 1:4)
    opts.tdim = 4;
end

% make sure residuals are good
if opts.tdim == 4
    res = double(res);
else
    res = double(permute(res(:, :, :, :), [setdiff(1:4, opts.tdim), opts.tdim]));
end
rsz = size(res);
res = reshape(res, prod(rsz(1:3)), size(res, 4));
aii = any(isinf(res), 2);
aii = aii | any(isnan(res), 2);
res(aii, :) = NaN;
res = reshape(res, rsz);

% perform computation
[cv, rx] = cov_nd(res(1:end-1, :, :, :), res(2:end, :, :, :));
[cv, ry] = cov_nd(res(:, 1:end-1, :, :), res(:, 2:end, :, :));
[cv, rz] = cov_nd(res(:, :, 1:end-1, :), res(:, :, 2:end, :));
rx(isinf(rx) | isnan(rx)) = 1 - eps;
ry(isinf(ry) | isnan(ry)) = 1 - eps;
rz(isinf(rz) | isnan(rz)) = 1 - eps;

% fisher-transform (for linear computation of average)
fii = [Inf; 1; 1] * [1, 1, 1];
fin = [0; 1; 0];
fia = [0; 1; 1; 0; 0];
rx = fisherr2z(cat(1, flexinterpn(fisherr2z(rx), ...
    [fii; size(rx)], {fia, fin, fin}, {1, 1, 1}), rx(end, :, :)), true);
ry = fisherr2z(cat(2, flexinterpn(fisherr2z(ry), ...
    [fii; size(ry)], {fin, fia, fin}, {1, 1, 1}), ry(:, end, :)), true);
rz = fisherr2z(cat(3, flexinterpn(fisherr2z(rz), ...
    [fii; size(rz)], {fin, fin, fia}, {1, 1, 1}), rz(:, :, end)), true);

% replace bad values
seps = 4 * sqrt(eps);
rx(isinf(rx) | isnan(rx) | rx < seps) = seps;
ry(isinf(ry) | isnan(ry) | ry < seps) = seps;
rz(isinf(rz) | isnan(rz) | rz < seps) = seps;

% compute FWHM estimate
ft = sqrt(-1 ./ (4 .* log(seps)));
sef = sqrt(8 * log(2));
fx = (mmpv(1) * sef) .* sqrt(-1 ./ (4 .* log(rx)));
fy = (mmpv(2) * sef) .* sqrt(-1 ./ (4 .* log(ry)));
fz = (mmpv(3) * sef) .* sqrt(-1 ./ (4 .* log(rz)));
fx(fx <= ft | fx > 1e8) = 0;
fy(fy <= ft | fy > 1e8) = 0;
fz(fz <= ft | fz > 1e8) = 0;

% discard boundary
fx = fx .* erode3d(fx > 0);
fy = fy .* erode3d(fy > 0);
fz = fz .* erode3d(fz > 0);

% compute outputs
fwhm = [ ...
    median(fx(fx > 0)), ...
    median(fy(fy > 0)), ...
    median(fz(fz > 0))];
if nargout > 1
    fi = harmmean(cat(4, fx, fy, fz), 4);
end
