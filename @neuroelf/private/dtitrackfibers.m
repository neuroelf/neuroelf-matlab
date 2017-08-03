function [fibers, flookup] = dtitrackfibers(fa, evec1, opts)
% dtitrackfibers  - track DTI-based fibers (using FA and first eigenvector)
%
% FORMAT:       [fibers, flookup] = dtitrackfibers(fa, evec1 [, opts])
%
% Input fields:
%
%       fa          XxYxZ FA map (scaled between 0 and 1)
%       evec1       XxYxZx3 first eigenvector array
%       opts        optional settings
%        .angthresh angle threshold, greater stops tracking (default: 70)
%        .fathresh  FA threshold (default: 0.2)
%        .faweight  weight vectors by FA when interpolating (true)
%        .flmin     minimum fiber length (in mm, default: 30)
%        .maxtiter  maximum tracking iterations per direction (default: 480)
%        .oversmp   seedmask oversampling factor (default: 1, max: 11)
%        .pbar      progress bar object
%        .res       1x3 voxel dimension resolution (default: [2, 2, 2])
%        .seedfa    FA threshold for seed mask (default: fathresh)
%        .seedmask  XxYxZ boolean mask (default: FA >= seedfa)
%        .stepsize  stepsize (mm, default: 1)
%
% Output fields:
%
%       fibers      Fx1 cell array with fibers (empty if below threshold)
%       flookup     XxYxZ cell array with indices into fibers array
%
% Note: this function has been developed based on DTISearch/Streamline
%       available at http://www.mathworks.com/matlabcentral/fileexchange/34008

% Version:  v0.9d
% Build:    14072016
% Date:     Jul-20 2014, 4:09 PM EST
% Author:   Chang Chia-Hao <oh75420 (at) gmail.com>
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
% Source:   http://www.mathworks.com/matlabcentral/fileexchange/34008

% Copyright (c) 2011 - 2012, Chang Chia-Hao; (c) 2014, Jochen Weber
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

% main argument check
if nargin < 2 || ...
    isempty(fa) || ...
   (~isa(fa, 'single') && ...
    ~isa(fa, 'double')) || ...
    ndims(fa) ~= 3 || ...
   (~isa(evec1, 'single') && ...
    ~isa(evec1, 'double')) || ...
   ~isequal(size(evec1), [size(fa), 3])
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% work on double arrays
if ~isa(fa, 'double')
    fa = double(fa);
end
vxs = size(fa);
xs = vxs(1);
ys = vxs(2);
zs = vxs(3);
if ~isa(evec1, 'double')
    evec1 = double(evec1);
end

% options
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end

% angle threshold
if ~isfield(opts, 'angthresh') || ...
   ~isa(opts.angthresh, 'double') || ...
    numel(opts.angthresh) ~= 1 || ...
    isinf(opts.angthresh) || ...
    isnan(opts.angthresh) || ...
    opts.angthresh <= 0 || ...
    opts.angthresh >= 90
    tangle = 70;
else
    tangle = min(85, max(30, opts.angthresh));
end

% FA threshold
if ~isfield(opts, 'fathresh') || ...
   ~isa(opts.fathresh, 'double') || ...
    numel(opts.fathresh) ~= 1 || ...
    isinf(opts.fathresh) || ...
    isnan(opts.fathresh) || ...
    opts.fathresh <= 0 || ...
    opts.fathresh >= 1
    tfa = 0.2;
else
    tfa = max(0.05, min(0.95, opts.fathresh));
end

% weighing by FA
if ~isfield(opts, 'faweight') || ...
   ~islogical(opts.faweight) || ...
    numel(opts.faweight) ~= 1
    faweight = true;
else
    faweight = opts.faweight;
end

% minimum fiber length
if ~isfield(opts, 'flmin') || ...
   ~isa(opts.flmin, 'double') || ...
    numel(opts.flmin) ~= 1 || ...
    isinf(opts.flmin) || ...
    isnan(opts.flmin) || ...
    opts.flmin <= 1
    flmin = 30;
else
    flmin = min(225, max(2, opts.flmin));
end

% maximum number of tracking iterations (in both directions)
if ~isfield(opts, 'maxtiter') || ...
   ~isa(opts.maxtiter, 'double') || ...
    numel(opts.maxtiter) ~= 1 || ...
    isinf(opts.maxtiter) || ...
    isnan(opts.maxtiter) || ...
    opts.maxtiter < 3
    maxtiter = 480;
else
    maxtiter = max(2 * flmin, min(1000, round(opts.maxtiter)));
end

% over-sampling
if ~isfield(opts, 'oversmp') || ...
   ~isa(opts.oversmp, 'double') || ...
    numel(opts.oversmp) ~= 1 || ...
    isinf(opts.oversmp) || ...
    isnan(opts.oversmp) || ...
   ~any((1:11) == opts.oversmp)
    oversmp = 1;
else
    oversmp = opts.oversmp;
end

% progress bar
if ~isfield(opts, 'pbar') || ...
    numel(opts.pbar) ~= 1 || ...
   (~isa(opts.pbar, 'xprogress') && ...
    ~isxfigure(opts.pbar, true))
    pbar = [];
else
    pbar = opts.pbar;
    pbarvis = pbar.Visible;
    pbar.Progress(0, 'Tracking DTI fibers...');
    pbar.Visible = 'on';
    drawnow;
end

% resolution
if ~isfield(opts, 'res') || ...
   ~isa(opts.res, 'double') || ...
    numel(opts.res) ~= 3 || ...
    any(isinf(opts.res) | isnan(opts.res) | opts.res <= 0 | opts.res > 8)
    vxr = [2, 2, 2];
else
    vxr = opts.res;
end

% seed-FA
if ~isfield(opts, 'seedfa') || ...
   ~isa(opts.seedfa, 'double') || ...
    numel(opts.seedfa) ~= 1 || ...
    isinf(opts.seedfa) || ...
    isnan(opts.seedfa) || ...
    opts.seedfa <= 0 || ...
    opts.seedfa > 1
    seedfa = tfa;
else
    seedfa = opts.seedfa;
end

% seedmask
if ~isfield(opts, 'seedmask') || ...
   ~islogical(opts.seedmask) || ...
   ~isequal(size(opts.seedmask), vxs)
    seedmsk = [];
else
    seedmsk = opts.seedmask;
end

% stepsize
if ~isfield(opts, 'stepsize') || ...
   ~isa(opts.stepsize, 'double') || ...
    numel(opts.stepsize) ~= 1 || ...
    isinf(opts.stepsize) || ...
    isnan(opts.stepsize) || ...
    opts.stepsize <= 0
    stepsize = 1;
else
    stepsize = max(0.1, min(4, opts.stepsize));
end
stepx = stepsize ./ vxr(1);
stepy = stepsize ./ vxr(2);
stepz = stepsize ./ vxr(3);

% create seed mask
if isempty(seedmsk)
    seedmsk = (fa >= seedfa);
end

% over-sampling
if oversmp > 1
    osstep = 1 / oversmp;
    osfrom = 0.5 + 0.5 * osstep;
    osx = round(osfrom:osstep:(xs+0.5));
    osy = round(osfrom:osstep:(ys+0.5));
    osz = round(osfrom:osstep:(zs+0.5));

    % still manageable (20MB logical array)
    if prod(oversmp .* vxs) < 2.1e7

        % proceed as normal
        seedmsk = seedmsk(osx, osy, osz);
        seedmsk = find(seedmsk(:));
        [sx, sy, sz] = ind2sub(vxs, seedmsk(:));
        sx = sx(:);
        sy = sy(:);
        sz = sz(:);

    % different tactic
    else

        % compute and allocate required voxels first
        numseed = sum(seedmsk(:)) * oversmp * oversmp * oversmp;
        sx = zeros(numseed, 1);
        sy = zeros(numseed, 1);
        sz = zeros(numseed, 1);
        ti = 1;

        % then apply along Z dimension
        osxy = numel(osx) * numel(osy);
        ovxs = oversmp .* vxs;
        for zc = 1:numel(osz)
            seedmski = seedmsk(osx, osy, osz(zc));
            seedmski = find(seedmski(:)) + (zc - 1) .* osxy;
            if isempty(seedmski)
                continue;
            end
            tt = ti + numel(seedmski) - 1;
            [sx(ti:tt), sy(ti:tt), sz(ti:tt)] = ind2sub(ovxs, seedmski);
            ti = tt + 1;
        end
    end
    sx = osfrom + osstep .* (sx - 1);
    sy = osfrom + osstep .* (sy - 1);
    sz = osfrom + osstep .* (sz - 1);
    
% no over-sampling
else
    seedmsk = find(seedmsk(:));
    [sx, sy, sz] = ind2sub(vxs, seedmsk(:));
    sx = sx(:);
    sy = sy(:);
    sz = sz(:);
end

% number of voxels in mask
numseed = numel(sx);
seedmski = 1:numel(sx);

% progress bar update
if ~isempty(pbar)
    pbar.Progress(0, sprintf('Tracking %d DTI fibers (direction 1)...', numseed));
    maxprogress = 2.5 * maxtiter;
    stepprogress = 1 / 86400;
    nextprogress = now + stepprogress;
end

% create output arrays
fibers = repmat({zeros(ceil(maxtiter / oversmp), 8)}, numseed, 1);
flookup = cell(vxs);

% temporary counting array
fcount = uint32(0);
fcount(xs, ys, zs) = 0;

% initialize counters for two directions!
s = ones(numseed, 1);
s2 = zeros(numseed, 1);

% rad-to-degree factor
r2d = 180 / pi;

% starting fiber tracking
for sdir = 1:2

    % at beginning (or either direction) we start with seed voxels
    x = sx;
    y = sy;
    z = sz;
    
    % sample FA and first eigenvector at coordinates (bilinear interpolation)
    [faxyz, ev1] = samplefaevec(x, y, z, fa, evec1, 3 - 2 * sdir, faweight);
    theta = zeros(numseed, 1);

    % keep track of voxels (within mask) that are still being tracked
    fmatch = 1:numseed;

    % while voxels remain
    while ~isempty(fmatch)

        % for first (positive) direction
        if sdir == 1

            % get next index into tracked fiber arrays (same for all!)
            si = s(fmatch(1));

            % iterate over fibers
            for mc = 1:numel(fmatch)

                % get target index (within mask!)
                ti = seedmski(fmatch(mc));

                % add to fiber (at that voxel and target index)
                % added are: 1x3 coordinate, 1x1 FA, 1x3 vector, 1x1 angle
                fibers{ti}(si, :) = ...
                    [x(mc), y(mc), z(mc), faxyz(mc), ev1(mc, :), theta(mc)];
            end

            % increase counter
            s(fmatch) = s(fmatch) + 1;

            % if iterations are reached
            if si >= maxtiter
                break;
            end

            % progress
            if ~isempty(pbar) && ...
                now > nextprogress
                pbar.Progress(si / maxprogress, ...
                    sprintf('Tracking %d DTI fibers (direction 1, step %d)...', ...
                    numel(fmatch), si));
                nextprogress = now + stepprogress;
            end
            
        % for second (negative) direction
        else

            % get the smallest *additional* index at which to add data
            s2i = s2(fmatch(1));

            % iterate over fibers
            for mc = 1:numel(fmatch)

                % get target index (this time it is used twice!)
                ti = fmatch(mc);

                % then add to fibers
                fibers{seedmski(ti)}(s2i + s(ti), :) = ...
                    [x(mc), y(mc), z(mc), faxyz(mc), ev1(mc, :), theta(mc)];
            end

            % increase counter
            s2(fmatch) = s2(fmatch) + 1;

            % if iterations in this direction are exceeded
            if s2i >= maxtiter
                break;
            end

            % progress
            if ~isempty(pbar) && ...
                now > nextprogress
                pbar.Progress((maxtiter + s2i) / maxprogress, ...
                    sprintf('Tracking %d DTI fibers (direction 2, step %d)...', ...
                    numel(fmatch), s2i));
                nextprogress = now + stepprogress;
            end
        end

        % keep copy of previous (first) eigenvector
        oev1 = ev1;

        % increase coordinates
        x = x + stepx .* ev1(:, 1);
        y = y + stepy .* ev1(:, 2);
        z = z + stepz .* ev1(:, 3);

        % sample FA and eigenvector at current coordinates
        [faxyz, ev1] = samplefaevec(x, y, z, fa, evec1, oev1, faweight);

        % compute angle between previous and new vectors
        theta = r2d .* real(acos(sum(oev1 .* ev1, 2)));

        % figure out which ones to keep (remaining within fmatch)
        remmatch = (x >= 1 & x <= xs & y >= 1 & y <= ys & z >= 1 & z <= zs & ...
            faxyz >= tfa & theta <= tangle);

        % mask matching array
        fmatch = fmatch(remmatch);

        % done?
        if isempty(fmatch)
            break;
        end

        % mask coordinates and vectors
        x = x(remmatch);
        y = y(remmatch);
        z = z(remmatch);
        ev1 = ev1(remmatch, :);
        theta = theta(remmatch);
    end

    % reduce number of elements by one
    if sdir == 1
        s = s - 1;
    else
        s2 = s2 - 1;
    end

    % compute total number of elements in fibers (after first direction,
    % this is only the number of elements in one direction!)
    s12 = s + s2;

    % iterate over ALL tracked fibers
    for fc = 1:numseed

        % get the number of elements for this fiber (used three times)
        si = s12(fc);

        % if more than one element
        if si > 1

            % then reverse the order of elements (between and after directions!)
            fibers{seedmski(fc)}(1:si, :) = fibers{seedmski(fc)}(si:-1:1, :);
        end
    end
end

% compute the length (in millimeters)
ss = s12 .* stepsize;

% and find those that make the threshold
t = (ss >= flmin);

% remove the others (cell(1) is better than {[]}!)
fibers(seedmski(~t)) = [];

% get the lengths of remaining fibers
s = s12(t);

% progress
if ~isempty(pbar)
    pbar.Progress(0.8, 'Creating voxel-to-fiber lookup matrix...');
end

% iterate over those
for fc = 1:numel(s)

    % cut each remaining fiber to the actual length (remove trailing 0s)
    fibers{fc} = single(fibers{fc}(1:s(fc), :));

    % get the (rounded) coordinates this fiber passes through
    tci = double(round(fibers{fc}(:, 1:3)));
    vxi = unique(sub2ind(vxs, tci(:, 1), tci(:, 2), tci(:, 3)));

    % increate the counter
    ti = fcount(vxi) + 1;
    fcount(vxi) = ti;

    % and iterate over those coordinates
    for ic = 1:numel(vxi)

        % and add to the list of fibers in that voxel
        flookup{vxi(ic)}(ti(ic)) = fc;
    end

    % progress
    if ~isempty(pbar) && ...
        now > nextprogress
        pbar.Progress(0.8 + 0.2 * fc / numel(s));
        nextprogress = now + stepprogress;
    end
end

% re-set progress bar
if ~isempty(pbar)
    pbar.Visible = pbarvis;
end



% sub-function samplefaevec, apply linear interpolation to FA an evec maps
function [fa, evec] = samplefaevec(x, y, z, fam, fem, evdir, faweight)

% number of samples
nx = numel(x);

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

% filter voxels (if needed)
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

    % prepare outputs (only needed if not all!)
    fa = zeros(nx, 1);
    evec = zeros(nx, 3);
    
    % nothing to do
    if isempty(x)
        return;
    end
end

% reshape arguments
fam = fam(:);
fem = reshape(fem, sxy * sz, 3);

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

% weigh with FA?
if faweight

    % compute FA sum
    fasum = v111 + v112 + v121 + v122 + v211 + v212 + v221 + v222;
    fasum(fasum == 0) = 1;

    % divide by ~= 0 elements
    fasum = fasum ./ sum(cat(2, v111 ~= 0, v112 ~= 0, v121 ~= 0, v122 ~= 0, ...
        v211 ~= 0, v212 ~= 0, v221 ~= 0, v222 ~= 0), 2);

    % reweigh
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

% access first eigenvector (in desired direction!)
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



% sub function for flipping vectors in right direction
function vec = vecorient(vec, vdir)
if numel(vdir) == 1
    [maxval, revpos] = max(abs(vec), [], 2);
    revpos = revpos(:)' - 1;
    revpos = (1:size(vec, 1)) + (size(vec, 1) .* revpos);
    if vdir == 1
        revpos = (vec(revpos) < 0);
    else
        revpos = (vec(revpos) > 0);
    end
else
    revpos = (sum(vec .* vdir, 2) < 0);
end
vec(revpos, :) = -vec(revpos, :);
