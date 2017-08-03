function [slindex, slidist] = sltargets(mask, opts)
%SLTARGETS  Create indices of searchlight targets within a 3D mask.
%   SLINDEX = SLTARGETS(MASK) creates a NUMEL(MASK)x1 cell array with
%   indices to elements for a searchlight analysis. MASK must be a logical
%   or numeric array (in which case the expression MASK > 0 is used).
%
%   The default is to use all within-mask elements that fit within a
%   regular 5x5 cube where the central 3x3 cube is all set and then the
%   most central voxel in each outside plane is set as well (rotated
%   diamond shaped).
%
%   SLINDEX = SLTARGETS(MASK, OPTS) allows to use optional settings for a
%   1x1 struct with fields:
%    .dist     1x1 double, maximum distance of target voxels/vertices (6)
%    .res      1x3 double, resolution of mask, default: FLOOR(2^8/LENGTH(MASK))
%    .samenum  1x1 logical, if true force the number of targets to its median
%
%   [SLINDEX, SLIDIST] = SLTARGETS(MASK) in addition returns the distance
%   of each target to the searchlight (center)

% Version:  v1.1
% Build:    16040421
% Date:     Apr-04 2016, 9:25 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2016, Jochen Weber
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

% mask argument
if nargin < 1 || (~isnumeric(mask) && ~islogical(mask)) || ndims(mask) > 3 || any(size(mask) < 5)
    error('neuroelf:general:badArgument', 'Bad or missing argument');
end

% create distances
cdist = (nargout > 1);

% for numeric mask
if isnumeric(mask)

    % convert to logical
    mask = (mask > 0);
end

% number of elements, size, and total search space
nmask = numel(mask);
smask = size(mask);
summask = sum(mask(:));

% options
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end

% distance (in mm)
if ~isfield(opts, 'dist') || ~isa(opts.dist, 'double') || numel(opts.dist) ~= 1 || ...
    isinf(opts.dist) || isnan(opts.dist) || opts.dist <= 0
    opts.dist = 6;
else
    opts.dist = min(24, opts.dist);
end

% resolution (per dim)
if ~isfield(opts, 'res') || ~isa(opts.res, 'double') || numel(opts.res) ~= 3 || ...
    any(isinf(opts.res) | isnan(opts.res) | opts.res <= 0)
    opts.res = floor(256 / max(size(mask))) .* [1, 1, 1];
else
    opts.res = opts.res(:)';
end

% same number of targets
if ~isfield(opts, 'samenum') || ~islogical(opts.samenum) || numel(opts.samenum) ~= 1
    opts.samenum = false;
else
    cdist = (cdist || opts.samenum);
end

% outputs
slindex = cell(nmask, 1);
if cdist
    slidist = cell(nmask, 1);
end

% create numeric mask
nummask = reshape(uint32(1:nmask), smask);
nummask(~mask) = 0;

% create required shape
xt = floor(opts.dist ./ opts.res + sqrt(eps));
if numel(smask) < 3
    [sh1, sh2] = ndgrid(-xt(1):xt(1), -xt(2):xt(2));
    sh3 = zeros(size(sh1));
    s3 = false;
else
    [sh1, sh2, sh3] = ndgrid(-xt(1):xt(1), -xt(2):xt(2), -xt(3):xt(3));
    s3 = true;
end
sh = [sh1(:), sh2(:), sh3(:)] .* (ones(numel(sh1), 1) * opts.res);

% restrict shape to distance
gd = sqrt(sum(sh .* sh, 2));
[gd, gdi] = sort(gd);
gdi(gd > opts.dist) = [];
ngdi = numel(gdi);

% create required indices
if s3
    sh = [sh1(:), sh2(:), sh3(:)];
else
    sh = [sh1(:), sh2(:)];
end
sh = sh(gdi, :);

% create numeric table
nslindex = uint32(0);
nslindex(summask, ngdi) = 0;

% first element is the mask!
nslindex(:, 1) = nummask(mask);

% iterate over shape elements
i2 = [Inf, Inf];
i3 = [Inf, Inf, Inf];
o2 = [1, 1];
o3 = [1, 1, 1];
for ec = 2:ngdi

    % access mask at shifted positions
    if s3
        si = o3 + sh(ec, :);
        ti = smask + sh(ec, :);
        shindex = indexarraynb(nummask, [i3; si; o3; ti]);
    else
        si = o2 + sh(ec, :);
        ti = smask + sh(ec, :);
        shindex = indexarraynb(nummask, [i2; si; o2; ti]);
    end
    nslindex(:, ec) = shindex(mask);
end

% translate to indices
nsltarget = find(mask(:));

% first pass
for ec = 1:summask
    usei = (nslindex(ec, :) > 0);
    slindex{nsltarget(ec)} = nslindex(ec, usei)';
    if cdist
        slidist{nsltarget(ec)} = gd(usei);
    end
end
