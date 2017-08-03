function [fa, evec, svec] = samplefaevec(x, y, z, fam, fem, evdir, opts)
% samplefaevec  - sample FA and eigenvector maps at X,Y,Z coordinates
%
% FORMAT:       [fa, evec] = samplefaevec(x, y, z, fa, evec [, evdir [, opts]])
%
% Input fields:
%
%       x, y, z     voxel coordinates into fa and evec
%       fa          3D FA map
%       evec        4D (3D-x-3) (first) eigenvector map
%       evdir       either of {1}, 1, or numel(x)-by-3 vector
%       opts        optional settings
%        .faweight  weight eigenvector interpolation by FA value (false)
%        .svec      second eigenvector (returned as third optional output)
%
% Output fields:
%
%       fa          numel(x)-by-1 interpolated FA values
%       evec        numel(x)-by-3 interpolated first eigenvectors
%       svec        numel(x)-by-3 interpolated second eigenvectors

% Version:  v0.9d
% Build:    14071611
% Date:     Jul-16 2014, 11:07 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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
if nargin < 5 || ...
   ~isa(x, 'double') || ...
   ~isa(y, 'double') || ...
   ~isa(z, 'double') || ...
    numel(x) ~= numel(y) || ...
    numel(x) ~= numel(z) || ...
   (~isa(fam, 'single') && ...
    ~isa(fam, 'double')) || ...
    ndims(fam) ~= 3 || ...
   (~isa(fem, 'single') && ...
    ~isa(fem, 'double')) || ...
    ndims(fem) ~= 4 || ...
    size(fam, 1) ~= size(fem, 1) || ...
    size(fam, 2) ~= size(fem, 2) || ...
    size(fam, 3) ~= size(fem, 3) || ...
    size(fem, 4) ~= 3
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
nx = numel(x);
if nargin < 6 || ...
   ~isa(evdir, 'double') || ...
   (numel(evdir) ~= 1 && ...
    ~isequal(size(evdir), [nx, 3]))
    evdir = 1;
elseif numel(evdir) == 1
    evdir = sign(evdir);
    if isnan(evdir) || ...
        evdir == 0
        evdir = 1;
    end
end
if nargin < 7 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'faweight') || ...
   ~islogical(opts.faweight) || ...
    numel(opts.faweight) ~= 1
    opts.faweight = false;
end
if ~isfield(opts, 'svec') || ...
   (~isa(opts.svec, 'double') && ...
    ~isa(opts.svec, 'single')) || ...
    ~isequal(size(opts.svec), size(fem))
    opts.svec = [];
    dosvec = false;
else
    dosvec = true;
end
if nargout < 3
    dosvec = false;
else
    svec = [];
end

% which voxels are well defined
sz = size(fam);
sx = sz(1);
sy = sz(2);
sz = sz(3);
sxy = sx * sy;
x = x(:);
y = y(:);
z = z(:);
welldef = (x >= 1 & x <= sx & y >= 1 & y <= sy & z >= 1 & z <= sz);
wdnum = numel(welldef);
wdsum = sum(welldef);
wdmax = 0.4 * nx;
if wdsum < wdmax
    welldef = find(welldef);
end
if wdsum < wdnum
    x = x(welldef);
    y = y(welldef);
    z = z(welldef);
    if numel(evdir) > 1
        evdir = evdir(welldef, :);
    end

    % prepare outputs
    fa = zeros(nx, 1);
    evec = zeros(nx, 3);
    if dosvec
        svec = zeros(nx, 3);
    end
end

% reshape arguments
fam = fam(:);
fem = reshape(fem, sxy * sz, 3);
if dosvec
    sem = reshape(opts.svec, size(fem));
end

% get base coordinate, interpolation weights, and base index
ix = floor(x);
iy = floor(y);
iz = floor(z);
cx = min(ix + 1, sx) - ix;
cy = sx .* (min(iy + 1, sy) - iy);
cxy = cx + cy;
cz = sxy .* (min(iz + 1, sz) - iz);
x = x - ix;
dx = 1 - x;
y = y - iy;
dy = 1 - y;
z = z - iz;
dz = 1 - z;
bxyz = ix + sx .* (iy - 1) + sxy .* (iz - 1);
zxyz = bxyz + cz;

% access data for FA
v111 = fam(bxyz);
v112 = fam(bxyz + cx);
v121 = fam(bxyz + cy);
v122 = fam(bxyz + cxy);
v211 = fam(zxyz);
v212 = fam(zxyz + cx);
v221 = fam(zxyz + cy);
v222 = fam(zxyz + cxy);

% remove NaNs
v111(isnan(v111)) = 0;
v112(isnan(v112)) = 0;
v121(isnan(v121)) = 0;
v122(isnan(v122)) = 0;
v211(isnan(v211)) = 0;
v212(isnan(v212)) = 0;
v221(isnan(v221)) = 0;
v222(isnan(v222)) = 0;

% weights
w111 = dx .* dy .* dz;
w112 =  x .* dy .* dz;
w121 = dx .*  y .* dz;
w122 =  x .*  y .* dz;
w211 = dx .* dy .*  z;
w212 =  x .* dy .*  z;
w221 = dx .*  y .*  z;
w222 =  x .*  y .*  z;

% reweight
if opts.faweight

    % compute FA sum
    fasum = v111 + v112 + v121 + v122 + v211 + v212 + v221 + v222;
    fasum(fasum == 0) = 1;

    % divide by ~= 0 elements
    fasum = fasum ./ sum(cat(2, v111 ~= 0, v112 ~= 0, v121 ~= 0, v122 ~= 0, ...
        v211 ~= 0, v212 ~= 0, v221 ~= 0, v222 ~= 0), 2);

    % reweight
    w111 = (v111 ./ fasum) .* w111;
    w112 = (v112 ./ fasum) .* w112;
    w121 = (v121 ./ fasum) .* w121;
    w122 = (v122 ./ fasum) .* w122;
    w211 = (v211 ./ fasum) .* w211;
    w212 = (v212 ./ fasum) .* w212;
    w221 = (v221 ./ fasum) .* w221;
    w222 = (v222 ./ fasum) .* w222;
    wsum = w111 + w112 + w121 + w122 + w211 + w212 + w221 + w222;
    w111 = w111 ./ wsum;
    w112 = w112 ./ wsum;
    w121 = w121 ./ wsum;
    w122 = w122 ./ wsum;
    w211 = w211 ./ wsum;
    w212 = w212 ./ wsum;
    w221 = w221 ./ wsum;
    w222 = w222 ./ wsum;
end

% interpolate
if wdsum < wdnum
    fa(welldef) = w111 .* v111 + w112 .* v112 + w121 .* v121 + w122 .* v122 + ...
        w211 .* v211 + w212 .* v212 + w221 .* v221 + w222 .* v222;
else
    fa = w111 .* v111 + w112 .* v112 + w121 .* v121 + w122 .* v122 + ...
        w211 .* v211 + w212 .* v212 + w221 .* v221 + w222 .* v222;
end

% access first eigenvector
v111 = vecorient(fem(bxyz, :), evdir);
v112 = vecorient(fem(bxyz + cx, :), evdir);
v121 = vecorient(fem(bxyz + cy, :), evdir);
v122 = vecorient(fem(bxyz + cxy, :), evdir);
v211 = vecorient(fem(zxyz, :), evdir);
v212 = vecorient(fem(zxyz + cx, :), evdir);
v221 = vecorient(fem(zxyz + cy, :), evdir);
v222 = vecorient(fem(zxyz + cxy, :), evdir);

% remove NaNs
v111(isnan(v111)) = 0;
v112(isnan(v112)) = 0;
v121(isnan(v121)) = 0;
v122(isnan(v122)) = 0;
v211(isnan(v211)) = 0;
v212(isnan(v212)) = 0;
v221(isnan(v221)) = 0;
v222(isnan(v222)) = 0;

% interpolate
o3 = [1, 1, 1];
vecs = w111(:, o3) .* v111 + w112(:, o3) .* v112 + ...
    w121(:, o3) .* v121 + w122(:, o3) .* v122 + ...
    w211(:, o3) .* v211 + w212(:, o3) .* v212 + ...
    w221(:, o3) .* v221 + w222(:, o3) .* v222;

% scale
vecs = vecs ./ repmat(sqrt(sum(vecs .* vecs, 2)), 1, 3);
if wdsum < wdnum
    evec(welldef, :) = vecs;
else
    evec = vecs;
end

% second vector as well
if dosvec
    v111 = vecorient(sem(bxyz, :), evdir);
    v112 = vecorient(sem(bxyz + cx, :), evdir);
    v121 = vecorient(sem(bxyz + cy, :), evdir);
    v122 = vecorient(sem(bxyz + cxy, :), evdir);
    v211 = vecorient(sem(zxyz, :), evdir);
    v212 = vecorient(sem(zxyz + cx, :), evdir);
    v221 = vecorient(sem(zxyz + cy, :), evdir);
    v222 = vecorient(sem(zxyz + cxy, :), evdir);

    % remove NaNs
    v111(isnan(v111)) = 0;
    v112(isnan(v112)) = 0;
    v121(isnan(v121)) = 0;
    v122(isnan(v122)) = 0;
    v211(isnan(v211)) = 0;
    v212(isnan(v212)) = 0;
    v221(isnan(v221)) = 0;
    v222(isnan(v222)) = 0;

    % interpolate
    vecs = w111(:, o3) .* v111 + w112(:, o3) .* v112 + ...
        w121(:, o3) .* v121 + w122(:, o3) .* v122 + ...
        w211(:, o3) .* v211 + w212(:, o3) .* v212 + ...
        w221(:, o3) .* v221 + w222(:, o3) .* v222;

    % scale
    vecs = vecs ./ repmat(sqrt(sum(vecs .* vecs, 2)), 1, 3);
    if wdsum < wdnum
        svec(welldef, :) = vecs;
    else
        svec = vecs;
    end
end



% sub function for flipping vectors in right direction
function vec = vecorient(vec, vdir)
if numel(vdir) == 1
    mp = maxpos(abs(vec), [], 2) - 1;
    mp = (1:size(vec, 1)) + (size(vec, 1) .* mp');
    if vdir == 1
        mp = (vec(mp) < 0);
    else
        mp = (vec(mp) > 0);
    end
else
    mp = (sum(vec .* vdir, 2) < 0);
end
vec(mp, :) = -vec(mp, :);
