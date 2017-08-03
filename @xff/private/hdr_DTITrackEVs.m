function srf = hdr_DTITrackEVs(xo, opts)
% HDR::DTITrackFEVs  - track DTI-based directions (using FA and EVs)
%
% FORMAT:       srf = hdr.DTITrackEVs([opts])
%
% Input fields:
%
%       opts        optional settings
%        .angthresh angle threshold, greater stops tracking (default: 20)
%        .atlas     xff object to sample endpoints for label lookup
%        .atlssph   sampling sphere (highest count is selected, default: 5)
%        .fathresh  FA threshold (default: 0.2)
%        .faweight  weight vectors by FA when interpolating (true)
%        .flmin     minimum fiber length (in mm, default: 30)
%        .maxtiter  maximum tracking iterations per direction (default: 480)
%        .oversmp   seedmask oversampling factor (default: 1, max: 11)
%        .pbar      progress bar object
%        .seedfa    FA threshold for seed mask (default: fathresh)
%        .seedmask  XxYxZ boolean mask (default: FA >= seedfa)
%        .stepsize  stepsize (mm, default: 1)
%        .title     title for tracking (default: auto)
%
% Output fields:
%
%       srf         SRF object with fibers in VertexCoordinate and
%                   additional information set in VertexColor,
%                   VertexNormal and RunTimeVars fields
%
% Note: this function has been developed based on DTISearch/Streamline
%       available at http://www.mathworks.com/matlabcentral/fileexchange/34008
%
% Using: bvcoordconv, ddeblank, findfirst, histcount, icbm2tal, lsqueeze,
%        maxpos, splittocellc, tdclient.

% Version:  v1.1
% Build:    16020515
% Date:     Feb-05 2016, 3:50 PM EST
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

% using many library functions
global ne_methods;
maxpos = ne_methods.maxpos;

% main argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'hdr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
rtv = bc.RunTimeVars;
cfr = hdr_CoordinateFrame(xo);
trf = cfr.Trf;
itrf = inv(trf);
trft = trf';
if ~isfield(rtv, 'Map') || ~isstruct(rtv.Map) || numel(rtv.Map) < 7 || ...
   ~isfield(rtv, 'BVals') || ~isfield(rtv, 'BVecs')
    error('neuroelf:xff:badArgument', 'Not a DTI-HDR/DTI-NII object.');
end
mnames = {rtv.Map.Name};
if ~any(strcmp(mnames, 'FA')) || ~any(strcmp(mnames, 'EVec1(1)')) || ...
   ~any(strcmp(mnames, 'EVec1(2)')) || ~any(strcmp(mnames, 'EVec1(3)')) || ...
   ~any(strcmp(mnames, 'EVec2(1)')) || ~any(strcmp(mnames, 'EVec2(2)')) || ...
   ~any(strcmp(mnames, 'EVec2(3)'))
    error('neuroelf:xff:badArgument', 'Not a DTI-HDR/DTI-NII object.');
end

% get FA and EVec1 maps
findfirst = ne_methods.findfirst;
fa = double(bc.VoxelData(:, :, :, findfirst(strcmp(mnames, 'FA'))));
evec1 = double(bc.VoxelData(:, :, :, [findfirst(strcmp(mnames, 'EVec1(1)')), ...
    findfirst(strcmp(mnames, 'EVec1(2)')), findfirst(strcmp(mnames, 'EVec1(3)'))]));
evec2 = double(bc.VoxelData(:, :, :, [findfirst(strcmp(mnames, 'EVec2(1)')), ...
    findfirst(strcmp(mnames, 'EVec2(2)')), findfirst(strcmp(mnames, 'EVec2(3)'))]));
vxr = sqrt(sum(trf(1:3, 1:3) .^ 2, 1));
vxs = size(fa);
xs = vxs(1);
ys = vxs(2);
zs = vxs(3);

% check options
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'angthresh') || ~isa(opts.angthresh, 'double') || numel(opts.angthresh) ~= 1 || ...
    isinf(opts.angthresh) || isnan(opts.angthresh) || opts.angthresh <= 0 || opts.angthresh >= 90
    tangle = 20;
else
    tangle = max(0.1, min(89, opts.angthresh));
end
tangle = cos((pi / 180) * tangle);
atlasfile = '';
labstrings = {};
if ~isfield(opts, 'atlas') || numel(opts.atlas) ~= 1 || ...
   ~xffisobject(opts.atlas, true, {'hdr', 'head', 'vmr'})
   	atlas = [];
else
    atlas = opts.atlas;
    atlasvol = aft_GetVolume(atlas, 1);
    atlasmax = double(max(atlasvol(:)));
    if isfield(atlas.C.RunTimeVars, 'AtlasLabels') && ...
        iscell(atlas.C.RunTimeVars.AtlasLabels) && ...
        all(cellfun(@ischar, atlas.C.RunTimeVars.AtlasLabels(:))) && ...
       ~any(cellfun('isempty', atlas.C.RunTimeVars.AtlasLabels(:))) && ...
        numel(atlas.C.RunTimeVars.AtlasLabels) >= atlasmax
        atlasfile = atlas.F;
        labstrings = atlas.C.RunTimeVars.AtlasLabels(:);
        asv = size(atlasvol);
        asx = asv(1);
        asy = asv(2);
        asz = asv(3);
        switch lower(atlas.S.Extensions{1})
            case 'hdr'
                atlastrf = hdr_CoordinateFrame(atlas);
                atlastrf = inv(atlastrf.Trf)';
            case 'head'
                atlastrf = head_CoordinateFrame(atlas);
                atlastrf = inv(atlastrf.Trf)';
            case 'vmr'
                atlastrf = ne_methods.bvcoordconv(zeros(0, 3), 'tal2bvx', ...
                    aft_BoundingBox(atlas))';
        end
    else
        atlas = [];
    end
end
if ~isfield(opts, 'atlssph') || ~isa(opts.atlssph, 'double') || numel(opts.atlssph) ~= 1 || ...
    isinf(opts.atlssph) || isnan(opts.atlssph) || opts.atlssph < 0
    opts.atlssph = 5;
else
    opts.atlssph = min(20, opts.atlssph);
end
if ~isempty(atlas) && opts.atlssph >= 1
    sx = floor(opts.atlssph);
    [sx, sy, sz] = ndgrid(-sx:sx, -sx:sx, -sx:sx);
    sx = sx(:);
    sy = sy(:);
    sz = sz(:);
    sxd = (sqrt(sx .* sx + sy .* sy + sz .* sz) <= opts.atlssph);
    atlssph = [sx(sxd), sy(sxd), sz(sxd)];
    atlssph(:, 4) = 0;
    atlsones = ones(size(atlssph, 1), 1);
elseif ~isempty(atlas)
    atlssph = [0, 0, 0, 0];
    atlsones = 1;
end
if ~isfield(opts, 'fathresh') || ~isa(opts.fathresh, 'double') || numel(opts.fathresh) ~= 1 || ...
    isinf(opts.fathresh) || isnan(opts.fathresh) || opts.fathresh <= 0 || opts.fathresh >= 1
    tfa = 0.2;
else
    tfa = max(0.05, min(0.95, opts.fathresh));
end
if ~isfield(opts, 'faweight') || ~islogical(opts.faweight) || numel(opts.faweight) ~= 1
    faweight = true;
else
    faweight = opts.faweight;
end
if ~isfield(opts, 'flmin') || ~isa(opts.flmin, 'double') || numel(opts.flmin) ~= 1 || ...
    isinf(opts.flmin) || isnan(opts.flmin) || opts.flmin <= 1
    flmin = 30;
else
    flmin = min(225, max(2, opts.flmin));
end
if ~isfield(opts, 'maxtiter') || ~isa(opts.maxtiter, 'double') || numel(opts.maxtiter) ~= 1 || ...
    isinf(opts.maxtiter) || isnan(opts.maxtiter) || opts.maxtiter < 3
    maxtiter = 480;
else
    maxtiter = max(2 * flmin, min(1000, round(opts.maxtiter)));
end
if ~isfield(opts, 'oversmp') || ~isa(opts.oversmp, 'double') || numel(opts.oversmp) ~= 1 || ...
    isinf(opts.oversmp) || isnan(opts.oversmp) || ~any((1:11) == opts.oversmp)
    oversmp = 1;
else
    oversmp = opts.oversmp;
end
if ~isfield(opts, 'pbar') || numel(opts.pbar) ~= 1 || ...
   (~isa(opts.pbar, 'xprogress') && ~isxfigure(opts.pbar, true))
    pbar = [];
else
    pbar = opts.pbar;
    pbarvis = pbar.Visible;
end
if ~isfield(opts, 'seedfa') || ~isa(opts.seedfa, 'double') || numel(opts.seedfa) ~= 1 || ...
    isinf(opts.seedfa) || isnan(opts.seedfa) || opts.seedfa <= 0 || opts.seedfa > 1
    seedfa = tfa;
else
    seedfa = opts.seedfa;
end
if ~isfield(opts, 'seedmask') || ~islogical(opts.seedmask) || ~isequal(size(opts.seedmask), vxs)
    seedmsk = [];
else
    seedmsk = opts.seedmask;
end
if isempty(seedmsk)
    seedmsk = (fa >= seedfa);
end
if ~isfield(opts, 'stepsize') || ~isa(opts.stepsize, 'double') || numel(opts.stepsize) ~= 1 || ...
    isinf(opts.stepsize) || isnan(opts.stepsize) || opts.stepsize <= 0
    stepsize = 1;
else
    stepsize = max(0.1, min(4, opts.stepsize));
end
stepx = stepsize ./ vxr(1);
stepy = stepsize ./ vxr(2);
stepz = stepsize ./ vxr(3);
if ~isfield(opts, 'title') || ~ischar(opts.title) || isempty(opts.title)
    opts.title = '';
end
if isempty(opts.title)
    if sum(seedmsk(:)) * prod(vxr) >= 1e5
        opts.title = 'whole-brain';
    else
        [sx, sy, sz] = ind2sub(vxs, find(seedmsk(:)));
        nsx = numel(sx);
        sx = mean(sx);
        sy = mean(sy);
        sz = mean(sz);
        mnicoord = trf * [sx; sy; sz; 1];
        talcoord = ne_methods.icbm2tal(mnicoord(1:3, 1)');
        tdcube = ne_methods.ddeblank(ne_methods.splittocellc( ...
            ne_methods.tdclient(talcoord, 'cube', 9), char(10)));
        tdcube(cellfun('isempty', regexpi(tdcube, 'white\s+matter'))) = [];
        if isempty(tdcube)
            tdlabel = 'white matter';
        else
            if tdcube{1}(end) == '*'
                tdlabel = regexprep(regexprep(lower(tdcube{1}), ...
                    '\,\s*white\s+matter.*$', ''), '^.*\,([^\,]+)$', '$1');
            else
                tdlabel = regexprep(lower(tdcube{1}), '^.*\,([^\,]+)$', '$1');
            end
        end
        opts.title = sprintf('%s ROI (%d voxels around [%d, %d, %d])', ...
            tdlabel, nsx, round(sx), round(sy), round(sz));
    end
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
        [sx, sy, sz] = ind2sub(oversmp .* vxs, seedmsk(:));
        sx = sx(:);
        sy = sy(:);
        sz = sz(:);

    % different tactic
    else

        % compute and allocate required voxels first
        seedmskz = ne_methods.lsqueeze(sum(sum(seedmsk, 1), 2));
        numseed = sum(seedmskz) * oversmp * oversmp * oversmp;
        seedmskz = find(seedmskz);
        sx = zeros(numseed, 1);
        sy = zeros(numseed, 1);
        sz = zeros(numseed, 1);
        ti = 1;

        % then apply along Z dimension
        osxy = numel(osx) * numel(osy);
        ovxs = oversmp .* vxs;
        for zc = seedmskz(:)'
            seedmski = seedmsk(osx, osy, zc);
            seedmski = find(seedmski(:)) + (zc - 1) .* (osxy * oversmp);
            for oc = 1:oversmp
                tt = ti + numel(seedmski) - 1;
                [sx(ti:tt), sy(ti:tt), sz(ti:tt)] = ind2sub(ovxs, seedmski);
                ti = tt + 1;
                seedmski = seedmski + osxy;
            end
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
    pbar.Progress(0, sprintf('Tracking %d DTI fibers (direction 1, step 1)...', numseed));
    pbar.Visible = 'on';
    maxprogress = 2.5 * maxtiter;
    stepprogress = 1 / 86400;
    nextprogress = now + stepprogress;
    drawnow;
end

% create tracking arrays
fibers = repmat({zeros(ceil(maxtiter / oversmp), 11)}, numseed, 1);

% initialize counters for two directions!
s = ones(numseed, 1);
s2 = zeros(numseed, 1);

% starting fiber tracking
for sdir = 1:2

    % at beginning (or either direction) we start with seed voxels
    x = sx;
    y = sy;
    z = sz;
    
    % sample FA and first eigenvector at coordinates (bilinear interpolation)
    [faxyz, ev1, ev2, theta] = ...
        samplefaevec(x, y, z, fa, evec1, evec2, 3 - 2 * sdir, faweight);

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
                    [x(mc), y(mc), z(mc), faxyz(mc), ev1(mc, :), ev2(mc, :), theta(mc)];
            end

            % increase counter
            s(fmatch) = s(fmatch) + 1;

            % if iterations are reached
            if si >= maxtiter
                break;
            end

            % progress
            if ~isempty(pbar) && now > nextprogress
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
                    [x(mc), y(mc), z(mc), faxyz(mc), ev1(mc, :), ev2(mc, :), theta(mc)];
            end

            % increase counter
            s2(fmatch) = s2(fmatch) + 1;

            % if iterations in this direction are exceeded
            if s2i >= maxtiter
                break;
            end

            % progress
            if ~isempty(pbar) && now > nextprogress
                pbar.Progress((maxtiter + s2i) / maxprogress, ...
                    sprintf('Tracking %d DTI fibers (direction 2, step %d)...', ...
                    numel(fmatch), s2i));
                nextprogress = now + stepprogress;
            end
        end

        % increase coordinates
        x = x + stepx .* ev1(:, 1);
        y = y + stepy .* ev1(:, 2);
        z = z + stepz .* ev1(:, 3);

        % sample FA and eigenvector at current coordinates
        [faxyz, ev1, ev2, theta] = samplefaevec(x, y, z, fa, evec1, evec2, ev1, faweight);

        % figure out which ones to keep (remaining within fmatch)
        remmatch = (x >= 1 & x <= xs & y >= 1 & y <= ys & z >= 1 & z <= zs & ...
            faxyz >= tfa & theta >= tangle);

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

% remove the others
fibers(seedmski(~t)) = [];

% get the lengths of remaining fibers
s = s12(t);
tlength = sum(s) + numel(fibers);

% arrays
srffibers = nan(tlength, 3);
srfcolors = nan(tlength, 3);
srfnormal = nan(tlength, 3);

% progress
if ~isempty(pbar)
    pbar.Progress(0.8, 'Creating voxel-to-fiber lookup matrix...');
end

% prepare index, lookup, and count arrays
sti = 1;
flookup = cell(vxs);
fcount = uint32(0);
fcount(xs, ys, zs) = 0;
if ~isempty(atlas)
    alabels = zeros(numel(s), 2);
else
    alabels = [];
end

% iterate over those
for fc = 1:numel(s)

    % cut each remaining fiber to the actual length (remove trailing 0s)
    fiber = single(fibers{fc}(1:s(fc), :));
    fl = s(fc) - 1;

    % get the (rounded) coordinates this fiber passes through
    tci = double(round(fiber(:, 1:3)));
    vxi = unique(sub2ind(vxs, tci(:, 1), tci(:, 2), tci(:, 3)));

    % increate the counter
    ti = fcount(vxi) + 1;
    fcount(vxi) = ti;

    % and iterate over those coordinates
    for ic = 1:numel(vxi)

        % and add to the list of fibers in that voxel
        flookup{vxi(ic)}(ti(ic), 1) = fc;
    end

    % store if SRF fibers
    stit = sti + fl;
    srffibers(sti:stit, :) = fiber(:, 1:3);
    srfcolors(sti:stit, :) = fiber(:, 5:7);
    srfnormal(sti:stit, :) = fiber(:, 8:10);

    % atlas labels
    if ~isempty(atlas)
        scoord1 = [srffibers(sti, :), 1] * trft;
        scoord2 = [srffibers(stit, :), 1] * trft;
        scoord1 = round((scoord1(atlsones, :) + atlssph) * atlastrf);
        scoord2 = round((scoord2(atlsones, :) + atlssph) * atlastrf);
        scoord1(any(scoord1 < 1, 2) | scoord1(:, 1) > asx | scoord1(:, 2) > asy | scoord1(:, 3) > asz, :) = [];
        scoord2(any(scoord2 < 1, 2) | scoord2(:, 1) > asx | scoord2(:, 2) > asy | scoord2(:, 3) > asz, :) = [];
        if ~isempty(scoord1)
            lab1 = atlasvol(sub2ind(asv, scoord1(:, 1), scoord1(:, 2), scoord1(:, 3)));
            lab1(lab1 == 0) = [];
            if ~isempty(lab1)
                lab1 = maxpos(ne_methods.histcount(lab1, 1, atlasmax, 1));
            else
                lab1 = 0;
            end
        else
            lab1 = 0;
        end
        if ~isempty(scoord2)
            lab2 = atlasvol(sub2ind(asv, scoord2(:, 1), scoord2(:, 2), scoord2(:, 3)));
            lab2(lab2 == 0) = [];
            if ~isempty(lab2)
                lab2 = maxpos(ne_methods.histcount(lab2, 1, atlasmax, 1));
            else
                lab2 = 0;
            end
        else
            lab2 = 0;
        end
        if lab1 <= lab2
            alabels(fc, :) = [lab1, lab2];
        else
            alabels(fc, :) = [lab2, lab1];
        end
    end

    % progress
    sti = sti + fl + 2;
    if ~isempty(pbar) && now > nextprogress
        pbar.Progress(0.8 + 0.2 * fc / numel(s));
        nextprogress = now + stepprogress;
    end
end

% convert count
mcount = max(fcount(:));
if mcount < 256
    fcount = uint8(fcount);
elseif mcount < 65536
    fcount = uint16(fcount);
else
    fcount = uint32(fcount);
end

% convert to BV coordinate system
srffibers(:, 4) = 1;
srffibers = srffibers * trft;

% work on colors
srfcolors = [nan(tlength, 1), round(255 .* srfcolors)];

% produce output object
srf = ne_methods.emptysrf();
srfc = srf.C;

% set in object
srfc.NrOfVertices = tlength;
srfc.NrOfTriangles = 0;
srfc.MeshCenter = [128, 128, 128];
srfc.VertexCoordinate = 128 - srffibers(:, [2, 3, 1]);
srfc.VertexNormal = srfnormal(:, [2, 3, 1]);
srfc.VertexColor = abs(srfcolors);
srfc.Neighbors = repmat({0, []}, tlength, 1);
srfc.RunTimeVars.AutoSave = true;
srfc.RunTimeVars.BVals = rtv.BVals;
srfc.RunTimeVars.BVecs = rtv.BVecs;
srfc.RunTimeVars.FiberAtlasFile = atlasfile;
srfc.RunTimeVars.FiberAtlasLabels = alabels;
srfc.RunTimeVars.FiberLookupMask = fcount;
srfc.RunTimeVars.FiberLookupCells = flookup(srfc.RunTimeVars.FiberLookupMask > 0);
srfc.RunTimeVars.FiberStarts = uint32(1 + [0; cumsum(1 + s(:))]);
srfc.RunTimeVars.LabelStrings = labstrings;
srfc.RunTimeVars.Oversmp = oversmp;
srfc.RunTimeVars.SourceDTIxffID = rtv.xffID;
srfc.RunTimeVars.SourceDTIFilename = xo.F;
if isfield(rtv, 'SPMsn')
    srfc.RunTimeVars.SPMsn = rtv.SPMsn;
end
srfc.RunTimeVars.TrackingStepsize = stepsize;
srfc.RunTimeVars.TrackingTitle = opts.title;
srfc.RunTimeVars.TransTalToVoxel = itrf;
srfc.RunTimeVars.TransVoxelToTal = trf;
srfc.RunTimeVars.VoxelResolution = vxr;
srf.C = srfc;

% re-set progress bar
if ~isempty(pbar)
    pbar.Visible = pbarvis;
end



% sub-function samplefaevec, apply linear interpolation to FA an evec maps
function [fa, evec, evec2, theta, theta2] = samplefaevec(x, y, z, fam, fem, sem, evdir, faweight)

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
    evec2 = zeros(nx, 3);
    theta = zeros(nx, 1);
    theta2 = zeros(nx, 1);
    
    % nothing to do
    if isempty(x)
        return;
    end
end

% make sure direction vector scales to 1
if numel(evdir) > 1
    evdir = evdir ./ repmat(sqrt(sum(evdir .* evdir, 2)), 1, 3);
end

% reshape arguments
fam = fam(:);
fem = reshape(fem, sxy * sz, 3);
sem = reshape(sem, sxy * sz, 3);

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
o3 = [1, 1, 1];
fav = w111 .* v111 + w112 .* v112 + w121 .* v121 + w122 .* v122 + ...
        w211 .* v211 + w212 .* v212 + w221 .* v221 + w222 .* v222;
if wdsum < wdnum
    fa(welldef) = fav;
else
    fa = fav;
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
vec1 = w111(:, o3) .* v111 + w112(:, o3) .* v112 + w121(:, o3) .* v121 + w122(:, o3) .* v122 + ...
    w211(:, o3) .* v211 + w212(:, o3) .* v212 + w221(:, o3) .* v221 + w222(:, o3) .* v222;

% access first eigenvector (in desired direction!)
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
vec2 = w111(:, o3) .* v111 + w112(:, o3) .* v112 + w121(:, o3) .* v121 + w122(:, o3) .* v122 + ...
    w211(:, o3) .* v211 + w212(:, o3) .* v212 + w221(:, o3) .* v221 + w222(:, o3) .* v222;

% scale
vec1 = vec1 ./ repmat(sqrt(sum(vec1 .* vec1, 2)), 1, 3);
vec2 = vec2 ./ repmat(sqrt(sum(vec2 .* vec2, 2)), 1, 3);

% theta
if numel(evdir) == 1
    tval1 = zeros(wdsum, 1);
    tval2 = zeros(wdsum, 1);
else
    tval1 = sum(vec1 .* evdir, 2);
    tval2 = sum(vec2 .* evdir, 2);
    swapvecs = (((tval2 - tval1) ./ (tval1 .* fav)) >= 2);
    nswaps = sum(swapvecs);
    if nswaps > 0 && nswaps < (0.4 * numel(swapvecs))
        swapvecs = find(swapvecs);
    end
    if nswaps > 0
        v1o = vec1(swapvecs, :);
        vec1(swapvecs, :) = vec2(swapvecs, :);
        vec2(swapvecs, :) = v1o;
        v1o = tval1(swapvecs);
        tval1(swapvecs) = tval2(swapvecs);
        tval2(swapvecs) = v1o;
    end
end
if wdsum < wdnum
    evec(welldef, :) = vec1;
    evec2(welldef, :) = vec1;
    theta(welldef) = tval1;
    theta2(welldef) = tval2;
else
    evec = vec1;
    evec2 = vec2;
    theta = tval1;
    theta2 = tval2;
end



% sub function for flipping vectors in right direction
function vec = vecorient(vec, vdir)
if numel(vdir) == 1
    [revmaxval, revpos] = max(abs(vec), [], 2);
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
