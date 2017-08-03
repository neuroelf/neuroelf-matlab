function [newfibers, colors] = srf_SelectFibers(xo, sel, opts)
% SRF::SelectFibers  - sub-select from a DTI tracking surface
%
% FORMAT:       subsrf = srf.SelectFibers(sel [, opts])
% FORMAT:       [fibers, colors] = srf.SelectFibers(sel [, opts])
%
% Input fields:
%
%       sel         either Cx3 coordinates or a 1x1 xff VOI object
%       opts        optional settings
%        .segdist   maximum number of segments from end to be considered (1)
%        .segsel    number of segments to select (default: Inf)
%        .svoisel   source VOI selection index
%        .tvoi      target (second set of regions that must be crossed)
%        .tvoisel   target VOI selection index (only if tvoi is a VOI)
%
% Output fields:
%
%       subsrf      new SRF object with sub-selection of fibers
%          -or-
%       fibers      Px3 fiber coordinates (with NaNs between fibers)
%       colors      Px3 colors (between 0 and 255 per RGB component)
%
% Note: the .segsel field can be non-integer, in which case the ending
%       vertex is interpolated between the preceding and last points
%
% Using: lsqueezecells.

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:50 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, Jochen Weber
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

% neuroelf library
global ne_methods;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get content and check it's a fiber-tracked file
bc = xo.C;
rtv = bc.RunTimeVars;
if ~isfield(rtv, 'FiberLookupMask') || ~isfield(rtv, 'FiberLookupCells') || ...
   ~isfield(rtv, 'FiberAtlasLabels') || ~isfield(rtv, 'LabelStrings') || ...
   ~isinteger(rtv.FiberLookupMask) || ~iscell(rtv.FiberLookupCells) || ...
    numel(rtv.FiberLookupCells) ~= sum(rtv.FiberLookupMask(:) > 0)
    error('neuroelf:xff:badObject', 'Invalid fiber tracks SRF object.');
end

% character selection
if ischar(sel) && ~isempty(sel)
    if isempty(rtv.FiberAtlasLabels) || isempty(rtv.LabelStrings)
        error('neuroelf:xff:badArgument', 'Tracking was done without Atlas labeling.');
    end
    sel = sel(:)';
    if any(sel == '*' | sel == '+' | sel == '?' | sel == '|')
        sel = find(~cellfun('isempty', regexpi(rtv.LabelStrings(:), sel)));
    else
        sel = find(strcmpi(rtv.LabelStrings(:), sel));
    end
    sel = sel(:);
end

% check selection argument
if ((numel(sel) ~= 1 || ~xffisobject(sel, true, 'voi')) && ...
    (~isa(sel, 'double') || isempty(sel) || ndims(sel) > 2 || ...
     ~any(size(sel, 2) == [1, 3, numel(sel)])))
    error('neuroelf:xff:badArgument', 'Bad selection argument.');
end

% VOI as selection
if numel(sel) == 1 && xffisobject(sel, true, 'voi')
    selc = sel.C;

% coordinate selection
elseif size(sel, 2) == 3 && (size(sel, 1) ~= 1 || all(abs(sel) <= 128))

    % create "fake" VOI content
    selc = xffnewcont('voi');
    selc.VOI.Voxels = sel;
    selc.VOI.NrOfVoxels = size(sel, 1);
    selc.VOI.Name = 'voxel selection';

% label selection
else
    if isempty(rtv.FiberAtlasLabels) || ...
        any(isinf(sel) | isnan(sel) | sel < 1 | sel ~= round(sel))
        error('neuroelf:xff:badArgument', ...
            'Tracking was done without Atlas labeling or bad label selection.');
    end
    fal = rtv.FiberAtlasLabels;
    selc = [];
    sel = ones(size(fal, 1), 1) * sel(:)';
    tidx = find(any((fal(:, 1) * ones(1, size(sel, 2))) == sel | ...
        (fal(:, 2) * ones(1, size(sel, 2))) == sel, 2));
end

% options
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'segdist') || ~isa(opts.segdist, 'double') || numel(opts.segdist) ~= 1 || ...
    isnan(opts.segdist) || opts.segdist < 0
    segdist = 2;
else
    segdist = 1 + opts.segdist;
end
if ~isfield(opts, 'segsel') || ~isa(opts.segsel, 'double') || numel(opts.segsel) ~= 1 || ...
    isnan(opts.segsel) || opts.segsel < 0
    opts.segsel = Inf;
end
if ~isempty(selc)
    if ~isfield(opts, 'svoisel') || ~isa(opts.svoisel, 'double') || isempty(opts.svoisel) || ...
        any(isinf(opts.svoisel(:)) | isnan(opts.svoisel(:)) | opts.svoisel(:) < 1)
        opts.svoisel = 1:numel(selc.VOI);
    else
        opts.svoisel = unique(round(min(numel(selc.VOI), opts.svoisel(:))))';
    end
end
if ~isfield(opts, 'tvoi') || ((numel(opts.tvoi) ~= 1 || ~xffisobject(opts.tvoi, true, 'voi')) && ...
    (~isa(opts.tvoi, 'double') || isempty(opts.tvoi) || ndims(opts.tvoi) > 2 || size(opts.tvoi, 2) ~= 3))
    opts.tvoi = [];
    eselc = [];
elseif numel(opts.tvoi) == 1 && xffisobject(opts.tvoi, true, 'voi')
    eselc = opts.tvoi.C;
else
    eselc = xffnewcont('voi');
    eselc.VOI.Voxels = opts.tvoi;
    eselc.VOI.NrOfVoxels = size(opts.tvoi, 1);
    eselc.VOI.Name = 'secondary voxel selection';
end
if ~isempty(eselc)
    if ~isfield(opts, 'tvoisel') || ~isa(opts.tvoisel, 'double') || isempty(opts.tvoisel) || ...
        any(isinf(opts.tvoisel(:)) | isnan(opts.tvoisel(:)) | opts.tvoisel(:) < 1)
        opts.tvoisel = 1:numel(eselc.VOI);
    else
        opts.tvoisel = unique(round(min(numel(eselc.VOI), opts.tvoisel(:))))';
    end
end

% get fibers and other info for sub-selection
fibers = bc.VertexCoordinate;
fevec2 = bc.VertexNormal;
fcolor = bc.VertexColor;

% reconstruct lookup cell array
flookup = cell(size(rtv.FiberLookupMask));
vxs = size(flookup);
vxx = vxs(1);
vxy = vxs(2);
vxz = vxs(3);
flookup(rtv.FiberLookupMask > 0) = rtv.FiberLookupCells;

% get transformation matrices
itrf = rtv.TransTalToVoxel';

% get information required for last segment to be drawn
drawstepsp = opts.segsel + 1;
drawstepst = ceil(drawstepsp);
drawstepw1 = drawstepst - drawstepsp;
drawstepw2 = 1 - drawstepw1;
drawstepsf = (drawstepsp ~= round(drawstepsp));

% get VOI content
if ~isempty(selc)
    voi = selc.VOI(opts.svoisel);
    tvox = cat(1, voi.Voxels);

    % center voxel (for title)
    cvox = round(mean(tvox));

    % apply transformation
    tvox(:, 4) = 1;
    vox = round(tvox * itrf);

    % and grab indices of fibers implicated in voxels from lookup
    vox(any(vox < 1, 2) | vox(:, 1) > vxx | vox(:, 2) > vxy | vox(:, 3) > vxz, :) = [];
    vox = unique(sub2ind(vxs, vox(:, 1), vox(:, 2), vox(:, 3)));
    numvox = numel(vox);
    tidx = ne_methods.lsqueezecells(flookup(vox));
    tidx = unique(cat(1, tidx{:}))';
else
    vox = [];
    numvox = numel(tidx);
    cvox = [0, 0, 0];
end    

% no fibers?
if isempty(tidx)

    % return the appropriate output(s)
    if nargout < 2
        newfibers = emptycopy(xo);
    else
        newfibers = zeros(0, 3);
        colors = zeros(0, 3);
    end
    return;
end

% then reconstruct indices
fromidx = rtv.FiberStarts(tidx);
toidx = rtv.FiberStarts(tidx+1) - 1;
flength = double(toidx - fromidx);

% sub-select
if ~isinf(segdist) && segdist < max(flength) && ~isempty(vox)

    % match fiber ends to mask
    endmask = false(vxs);
    endmask(vox) = true;
    tfibers = [128 - fibers(:, [3, 1, 2]), ones(size(fibers, 1), 1)] * itrf;
    fmatch = true(numel(tidx), 1);
    for fc = 1:numel(fmatch)
        f = tfibers(fromidx(fc):toidx(fc)-1, 1:3);
        if size(f, 1) <= segdist
            continue;
        end
        tcrd = round(double(f(1:segdist, :)));
        tcrd(any(tcrd < 1, 2) | tcrd(:, 1) > vxx | tcrd(:, 2) > vxy | tcrd(:, 3) > vxz, :) = [];
        if any(endmask(sub2ind(vxs, tcrd(:, 1), tcrd(:, 2), tcrd(:, 3))))
            continue;
        end
        tcrd = round(double(f(end+1-segdist:end, :)));
        tcrd(any(tcrd < 1, 2) | tcrd(:, 1) > vxx | tcrd(:, 2) > vxy | tcrd(:, 3) > vxz, :) = [];
        if any(endmask(sub2ind(vxs, tcrd(:, 1), tcrd(:, 2), tcrd(:, 3))))
            fibers(fromidx(fc):toidx(fc)-1, :) = fibers(toidx(fc)-1:-1:fromidx(fc), :);
            fevec2(fromidx(fc):toidx(fc)-1, :) = fevec2(toidx(fc)-1:-1:fromidx(fc), :);
            fcolor(fromidx(fc):toidx(fc)-1, :) = fcolor(toidx(fc)-1:-1:fromidx(fc), :);
        else
            fmatch(fc) = false;
        end
    end

    % mask
    tidx = tidx(fmatch);
    fromidx = fromidx(fmatch);
    toidx = toidx(fmatch);
    flength = flength(fmatch);
    if isempty(fromidx)
        if nargout < 2
            newfibers = emptycopy(xo);
        else
            newfibers = zeros(0, 3);
            colors = zeros(0, 3);
        end
        return;
    end
end

% additional mask (endpoints)
if ~isinf(segdist) && ~isempty(opts.tvoi)
    voi = eselc.VOI(opts.tvoisel);
    vox = cat(1, voi.Voxels);
    vox(:, 4) = 1;
    vox = round(vox * itrf);
    vox(any(vox < 1, 2) | vox(:, 1) > vxx | vox(:, 2) > vxy | vox(:, 3) > vxz, :) = [];
    vox = unique(sub2ind(vxs, vox(:, 1), vox(:, 2), vox(:, 3)));
    endmask = false(vxs);
    endmask(vox) = true;
    tfibers = [128 - fibers(:, [3, 1, 2]), ones(size(fibers, 1), 1)] * itrf;
    fmatch = true(numel(fromidx), 1);
    for fc = 1:numel(fmatch)
        f = tfibers(fromidx(fc):toidx(fc)-1, :);
        tcrd = round(double(f(max(1, end+1-segdist):end, :)));
        tcrd(any(tcrd < 1, 2) | tcrd(:, 1) > vxx | tcrd(:, 2) > vxy | tcrd(:, 3) > vxz, :) = [];
        if ~any(endmask(sub2ind(vxs, tcrd(:, 1), tcrd(:, 2), tcrd(:, 3))))
            fmatch(fc) = false;
        end
    end
    tidx = tidx(fmatch);
    fromidx = fromidx(fmatch);
    flength = flength(fmatch);
    if isempty(fromidx)
        if nargout < 2
            newfibers = emptycopy(xo);
        else
            newfibers = zeros(0, 3);
            colors = zeros(0, 3);
        end
        return;
    end
end

% compile new SRF data
if ~isinf(opts.segsel)
    flength = min(flength, opts.segsel);
end
tlength = sum(flength) + numel(flength);
flength = flength - 1;
toidx = fromidx + uint32(flength);
newfibers = nan(tlength, 3);
colors = nan(tlength, 4);
tfrom = 1;
if nargout == 1
    newfevecs = nan(tlength, 3);
    flookup = cell(vxs);
    fcount = uint32(0);
    fcount(vxx, vxy, vxz) = 0;
    for fc = 1:numel(flength)
        fl = flength(fc);
        tto = tfrom + fl;
        fidx = fromidx(fc);
        eidx = toidx(fc);
        newfibers(tfrom:tto, :) = fibers(fidx:eidx, :);
        newfevecs(tfrom:tto, :) = fevec2(fidx:eidx, :);
        colors(tfrom:tto, :) = fcolor(fidx:eidx, :);
        if drawstepsf
            newfibers(tto, :) = drawstepw1 .* newfibers(tto-1, :) + drawstepw2 .* newfibers(tto, :);
            newfevecs(tto, :) = drawstepw1 .* newfevecs(tto-1, :) + drawstepw2 .* newfevecs(tto, :);
            newfevecs(tto, :) = newfevecs(tto, :) ./ sqrt(sum(newfevecs(tto, :) .* newfevecs(tto, :)));
            colors(tto, :) = round(drawstepw1 .* colors(tto-1, :) + drawstepw2 .* colors(tto, :));
        end
        tfibers = round([128 - newfibers(tfrom:tto, [3, 1, 2]), ones(eidx + 1 - fidx, 1)] * itrf);
        tfibers(any(tfibers < 1, 2) | tfibers(:, 1) > vxx | tfibers(:, 2) > vxy | tfibers(:, 3) > vxz, :) = [];
        tfibers = unique(sub2ind(vxs, tfibers(:, 1), tfibers(:, 2), tfibers(:, 3)));
        tcount = fcount(tfibers) + 1;
        for lc = 1:numel(tfibers)
            flookup{tfibers(lc)}(tcount(lc), 1) = fc;
        end
        fcount(tfibers) = tcount;
        tfrom = tto + 2;
    end
else
    for fc = 1:numel(flength)
        fl = flength(fc);
        tto = tfrom + fl;
        fidx = fromidx(fc);
        eidx = toidx(fc);
        fibers(tfrom:tto, :) = fibers(fidx:eidx, :);
        colors(tfrom:tto, :) = fcolor(fidx:eidx, :);
        if drawstepsf
            newfibers(tto, :) = drawstepw1 .* newfibers(tto-1, :) + drawstepw2 .* newfibers(tto, :);
            colors(tto, :) = round(drawstepw1 .* colors(tto-1, :) + drawstepw2 .* colors(tto, :));
        end
    end
    return;
end

% create new SRF (copy)
fibers = newfibers;
newfibers = aft_CopyObject(xo);
bc.NrOfVertices = tlength;
bc.VertexCoordinate = fibers;
bc.VertexNormal = newfevecs;
bc.VertexColor = colors;
bc.Neighbors = repmat({0, []}, tlength, 1);
flmask = cellfun('prodofsize', flookup);
mflmask = max(flmask(:));
if ~isempty(bc.RunTimeVars.FiberAtlasLabels)
    bc.RunTimeVars.FiberAtlasLabels = bc.RunTimeVars.FiberAtlasLabels(tidx, :);
end
if mflmask < 256
    bc.RunTimeVars.FiberLookupMask = uint8(flmask);
elseif mflmask < 65536
    bc.RunTimeVars.FiberLookupMask = uint16(flmask);
else
    bc.RunTimeVars.FiberLookupMask = uint32(flmask);
end
bc.RunTimeVars.FiberLookupCells = flookup(flmask(:) > 0);
bc.RunTimeVars.FiberStarts = uint32(1 + [0; cumsum(2 + flength(:))]);
bc.RunTimeVars.TrackingTitle = ...
    sprintf('%s (masked with %d voxels around [%d, %d, %d])', ...
    rtv.TrackingTitle, numvox, cvox(1), cvox(2), cvox(3));
newfibers.C = bc;



% sub-function for an empty copy
function e = emptycopy(h)

% create empty copy
e = aft_CopyObject(h);
ec = e.C;
ec.NrOfVertices = 0;
ec.NrOfTriangles = 0;
ec.VertexCoordinate = zeros(0, 3);
ec.VertexNormal = zeros(0, 3);
ec.VertexColor = zeros(0, 4);
ec.Neighbors = cell(0, 2);
ec.TriangleVertex = zeros(0, 3);
ec.NrOfTriangleStrips = 0;
ec.TriangleStripSequence = zeros(0, 1);
ec.RunTimeVars.FiberAtlasLabels = [];
ec.RunTimeVars.FiberLookupMask = uint8(zeros(size(ec.RunTimeVars.FiberLookupMask)));
ec.RunTimeVars.FiberLookupCells = cell(0, 1);
ec.RunTimeVars.FiberStarts = 1;
ec.RunTimeVars.LabelStrings = cell(0, 1);
e.C = ec;
