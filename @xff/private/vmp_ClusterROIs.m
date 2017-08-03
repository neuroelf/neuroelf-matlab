function voi = vmp_ClusterROIs(xo, mapno, coords, opts)
% VMP::ClusterVOIs  - create ROIs based on clusters in a map
%
% FORMAT:       voi = vmp.ClusterROIs(mapsel, coords [, opts])
%
% Input fields:
%
%       mapsel      map selection, default: 1
%       coords      coordinate(s) as a Cx3 array
%       opts        optional settings
%        .cfunc     cost function (handle) @(s, d), default:
%                   ((prod(abs(s), 2) .^ (1/size(s, 2))) ./ (1 + log(1+d)))
%        .mincval   minimum cost-function value, default: 0
%        .mmres     create ROIs in 1mm resolution, default: false
%        .negative  boolean flag, also investigate negative clusters
%        .overlap   allow clusters to re-use voxels, default: false
%        .peakdist  maximum peak-to-coord distance, default: 8mm
%        .size      1x1 or Cx1 size in mm^3, default: 1000
%        .space     space in which coordinates are given, either of
%                   'bvi', 'bvs', {'tal'}, 'voxel'
%        .thresh    clustering threshold(s), default: from selected maps
%
% Output fields:
%
%       voi         VOI object with discovered ROI clusters
%
% Note: the cost function handle must accept (exactly) two arguments in
%       this order, other potential examples for fcost are:
%       @(s, d) ((prod(abs(s), 2) .^ (1/size(s, 2))) ./ sqrt(1 + d))
%  and  @(s, d) ((prod(abs(s), 2) .^ (1/size(s, 2))) ./ (1 + d))
%       which stronger penalize distance than log (forcing rounder shapes,
%       more important if clusters are likely to overlap!)
%  or   @(s, d) ((sum(abs(s), 2) ./ size(s, 2)) ./ sqrt(1 + d))
%       which rather averages the statistics of several maps
%       the map selection can be either 1x1, Cx1, 1xN, or CxN
%       in the last case (CxN) a 0 value indicates that for a particular
%       coordinate, not all N maps are to be used
%
% Using: bvcoordconv, clustercoordsc, lsqueeze.

% Version:  v1.1
% Build:    16021315
% Date:     Feb-13 2016, 3:30 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
bvcoordconv    = ne_methods.bvcoordconv;
clustercoordsc = ne_methods.clustercoordsc;

% check arguments
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp') || ...
   ~isa(mapno, 'double') || isempty(mapno) || ndims(mapno) > 2 || ...
    any(isinf(mapno(:)) | isnan(mapno(:)) | mapno(:) < 0 | mapno(:) ~= fix(mapno)) || ...
   ~isa(coords, 'double') || isempty(coords) || ndims(coords) > 2 || size(coords, 2) ~= 3 || ...
    any(isinf(coords(:)) | isnan(coords(:)) | coords(:) < -256 | coords(:) > 512) || ...
   ~any(size(mapno, 1) == [1, size(coords, 1)]) || any(all(mapno == 0, 2))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
mapres = bc.Resolution;
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'cfunc') || ~isa(opts.cfunc, 'function_handle') || numel(opts.cfunc) ~= 1
    opts.cfunc = @(s, d) (prod(abs(s), 2) .^ (1 / size(s, 2)) ./ (1 + log(1 + d)));
end
cfunc = opts.cfunc;
try
    v = cfunc(3 + 2 .* rand(10, 2), 4 .* rand(10, 1));
    if ~isequal(size(v), [10, 1]) || any(isinf(v) | isnan(v))
        error('Invalid output size or values.');
    end
catch xfferror
    error('neuroelf:xff:badArgument', 'Invalid cost function: %s.', xfferror.message);
end

% for mm-res clustering of non-mm VMPs
if isfield(opts, 'mmres') && islogical(opts.mmres) && numel(opts.mmres) == 1 && opts.mmres && mapres > 1

    % re-sample required maps
    rmap = unique(mapno(mapno > 0));
    imap = zeros(1, max(rmap));
    imap(rmap) = 1:numel(rmap);
    to = vmp_MakeHiResVMP(xo, rmap(:)');

    % re-set input arguments
    newmapno = mapno;
    newmapno(newmapno > 0) = imap(newmapno > 0);
    opts.mmres = false;

    % create VOI
    voi = vmp_ClusterROIs(to, newmapno, coords, opts);

    % clear temp object
    delete(to);
    return;
end

% otherwise
if ~isfield(opts, 'mincval') || ~isa(opts.mincval, 'double') || ...
    numel(opts.mincval) ~= 1 || isnan(opts.mincval)
    opts.mincval = 0;
end
if ~isfield(opts, 'negative') || ~islogical(opts.negative) || numel(opts.negative) ~= 1
    opts.negative = false;
end
if ~isfield(opts, 'overlap') || ~islogical(opts.overlap) || numel(opts.overlap) ~= 1
    opts.overlap = false;
end
if ~isfield(opts, 'peakdist') || ~isa(opts.peakdist, 'double') || numel(opts.peakdist) ~= 1 || ...
    isinf(opts.peakdist) || isnan(opts.peakdist) || opts.peakdist < 0
    opts.peakdist = 8;
end
if ~isfield(opts, 'size') || ~isa(opts.size, 'double') || numel(opts.size) ~= 1 || ...
    isinf(opts.size) || isnan(opts.size) || opts.size <= 0
    opts.size = 1000;
end
if ~isfield(opts, 'space') || ~ischar(opts.space) || ...
   ~any(strcmpi(opts.space(:)', {'bvi', 'bvs', 'tal', 'voxel'}))
    opts.space = 'tal';
else
    opts.space = lower(opts.space(:)');
end
if ~isfield(opts, 'thresh') || ~isa(opts.thresh, 'double') || ~isequal(size(mapno), size(opts.thresh))
    th = nan(size(mapno));
else
    th = opts.thresh;
end
if any(mapno(:) > numel(bc.Map))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
for mc = 1:numel(mapno)
    if isnan(th(mc))
        th(mc) = abs(bc.Map(mapno(mc)).LowerThreshold);
    else
        th(mc) = abs(th(mc));
    end
end
if size(mapno, 1) == 1 && size(coords, 1) > 1
    mapno = ones(size(coords, 1), 1) * mapno;
    th = ones(size(coords, 1), 1) * th;
end
msz = size(bc.Map(1).VMPData);

% get bounding box
bbox = aft_BoundingBox(xo);

% convert coordinates to voxel coordinates if required
switch (opts.space)
    case 'tal'
        voxcoords = bvcoordconv(coords, 'tal2bvc', bbox);
    case 'bvi'
        voxcoords = bvcoordconv(coords, 'bvi2bvc', bbox);
    case 'bvs'
        voxcoords = bvcoordconv(coords, 'bvs2bvc', bbox);
    case 'voxel'
        voxcoords = coords;
end

% find peaks according to maximum distance from coords
numcoords = size(coords, 1);
peaks = voxcoords;
peakdist = opts.peakdist / mapres;
pdgrid = -ceil(peakdist + 0.5):ceil(peakdist + 0.5);
[pdgrid, y, z] = ndgrid(pdgrid, pdgrid, pdgrid);
pdgrid = [pdgrid(:), y(:), z(:)];
pdgrid(sqrt(sum(pdgrid .* pdgrid, 2)) > (peakdist + 0.5), :) = [];
npdg = size(pdgrid, 1);

% iterate over coordinates
for cc = 1:numcoords

    % get maps to be considered
    maps = mapno(cc, mapno(cc, :) > 0);

    % get coordinates to include in peak search
    tcrds = round(ones(npdg, 1) * peaks(cc, :) + pdgrid);
    ucrds = all(tcrds > 0, 2) & ...
        tcrds(:, 1) <= msz(1) & tcrds(:, 2) <= msz(2) & tcrds(:, 3) <= msz(3) & ...
        sqrt(sum((tcrds - ones(npdg, 1) * peaks(cc, :)) .^ 2, 2)) <= peakdist;
    tcrds = sub2ind(msz, tcrds(ucrds, 1), tcrds(ucrds, 2), tcrds(ucrds, 3));
    tvals = zeros(numel(tcrds), numel(maps));

    % iterate over maps and get values
    for mc = 1:numel(maps)
        if opts.negative
            tvals(:, mc) = abs(bc.Map(maps(mc)).VMPData(tcrds));
        else
            tvals(:, mc) = bc.Map(maps(mc)).VMPData(tcrds);
        end
    end

    % compute cost function with 0 distance (all voxels within radius)
    [tvals, tvali] = sort(cfunc(tvals, zeros(size(tvals, 1), 1)), 'descend');
    if tvals(1) <= 0
        error('neuroelf:xff:badArgument', 'No peak detected for coordinate %d.', cc);
    end

    % find peak (first in descending sorted values) coordinate
    [peakx, peaky, peakz] = ind2sub(msz, tcrds(tvali(1)));
    peaks(cc, :) = [peakx, peaky, peakz];
end

% check that all peaks are unique (unless overlap is allowed)
if size(unique(peaks, 'rows'), 1) ~= numcoords && ~opts.overlap

    % which peaks are the same?
    samep = corrcoef([-100 .* ones(1, numcoords); peaks'; 100 .* ones(1, numcoords)]);
    samep(tril(ones(numcoords, numcoords)) > 0) = 0;
    [samep, samey] = ind2sub(size(samep), find(samep(:) == 1));

    % issue warning
    warning('neuroelf:xff:badArgument', ...
        '%d peak(s) not uniquely defined and overlap not allowed.', numel(samey));

    % for those peaks, find a slightly different solution!
    for cc = 1:numel(samey)

        % compute difference from original coordinates
        px = samep(cc);
        py = samey(cc);
        diffx = voxcoords(px, :) - peaks(px, :);
        diffy = voxcoords(py, :) - peaks(py, :);

        % if directions are "opposite"
        if sum(diffx .* diffy) <= 0

            % go into the direction with both peaks (about 2 voxel units)
            peaks(px, :) = peaks(px, :) + round((2 / sqrt(sum(diffx .* diffx))) .* diffx);
            peaks(py, :) = peaks(py, :) + round((2 / sqrt(sum(diffy .* diffy))) .* diffy);

        % otherwise
        else

            % go into the orthogonal direction (between X and Y, 1.5 units)
            diffxy = diffx - diffy;
            diffxy = round((1.5 / sqrt(sum(diffxy .* diffxy))) .* diffxy);
            peaks(px, :) = peaks(px, :) + diffxy;
            peaks(py, :) = peaks(py, :) - diffxy;
        end
    end
end

% create VOI object
voi = xff('new:voi');

% get VOI content
vc = voi.L;

% extend VOI
vc.VOI(numcoords).Name = '';

% maximum total radius is twice that needed for a sphere with size
maxrad = (2 / mapres) * ((0.75 * opts.size / pi) .^ (1 / 3));

% maximum number of voxels
maxvox = ceil(opts.size / (mapres * mapres * mapres));

% create maximally required voxels for search
pdgrid = -ceil(maxrad):ceil(maxrad);
[pdgrid, y, z] = ndgrid(pdgrid, pdgrid, pdgrid);
pdgrid = [pdgrid(:), y(:), z(:)];
pdgrid(sqrt(sum(pdgrid .* pdgrid, 2)) > maxrad, :) = [];
[y, z] = sort(sqrt(sum(pdgrid .* pdgrid, 2)));
pdgrid = pdgrid(z, :);
npdg = size(pdgrid, 1);

% overlap
if opts.overlap

    % iterate over coordinates
    for cc = 1:numcoords
        pk = peaks(cc, :);
        maps = mapno(cc, mapno(cc, :) > 0);
        mth = th(mapno(cc, :) > 0);
        tcrds = round(ones(npdg, 1) * pk + pdgrid);
        ucrds = all(tcrds > 0, 2) & ...
            tcrds(:, 1) <= msz(1) & tcrds(:, 2) <= msz(2) & tcrds(:, 3) <= msz(3);
        tcrds = sub2ind(msz, tcrds(ucrds, 1), tcrds(ucrds, 2), tcrds(ucrds, 3));
        dcrds = mapres .* sqrt(sum(pdgrid(ucrds, :) .* pdgrid(ucrds, :), 2));
        tvals = zeros(numel(tcrds), numel(maps));

        % extract map values
        for mc = 1:numel(maps)

            % cluster map
            mapv = bc.Map(maps(mc)).VMPData;
            if mapv(pk(1), pk(2), pk(3)) < 0
                mapv = -mapv;
            end
            [cs, cv] = clustercoordsc(mapv >= mth(mc));
            mapv(cv ~= cv(pk(1), pk(2), pk(3))) = 0;
            tvals(:, mc) = (1 / mapv(tcrds(1))) .* mapv(tcrds);
        end

        % compute cost function
        [tvals, tvali] = sort(cfunc(tvals, dcrds), 'descend');

        % sort coordinates and get first up to required size
        [peakx, peaky, peakz] = ind2sub(msz, tcrds(tvali(1:min(maxvox, sum(tvals > 0)))));

        % set VOI field
        vc.VOI(cc).Voxels = bvcoordconv([peakx, peaky, peakz], 'bvc2tal', bbox);
    end

% no overlap
else

    % create final tvals and coordinates arrays
    ftval = zeros(npdg, numcoords);
    ftcrd = zeros(npdg, numcoords);
    tcnum = ne_methods.lsqueeze(ones(npdg, 1) * (1:numcoords));

    % iterate over coordinates
    for cc = 1:numcoords
        pk = peaks(cc, :);
        maps = mapno(cc, mapno(cc, :) > 0);
        mth = th(mapno(cc, :) > 0);
        tcrds = round(ones(npdg, 1) * pk + pdgrid);
        ucrds = all(tcrds > 0, 2) & ...
            tcrds(:, 1) <= msz(1) & tcrds(:, 2) <= msz(2) & tcrds(:, 3) <= msz(3);
        tcrds = sub2ind(msz, tcrds(ucrds, 1), tcrds(ucrds, 2), tcrds(ucrds, 3));
        dcrds = mapres .* sqrt(sum(pdgrid(ucrds, :) .* pdgrid(ucrds, :), 2));
        tvals = zeros(numel(tcrds), numel(maps));

        % extract map values
        for mc = 1:numel(maps)

            % cluster map
            mapv = bc.Map(maps(mc)).VMPData;
            if mapv(pk(1), pk(2), pk(3)) < 0
                mapv = -mapv;
            end
            [cs, cv] = clustercoordsc(mapv >= mth(mc));
            mapv(cv ~= cv(pk(1), pk(2), pk(3))) = 0;
            tvals(:, mc) = (1 / mapv(tcrds(1))) .* mapv(tcrds);
        end

        % compute cost function
        [tvals, tvali] = sort(cfunc(tvals, dcrds), 'descend');

        % store in final arrays
        ftval(1:numel(tvals), cc) = tvals;
        ftcrd(1:numel(tvals), cc) = tcrds(tvali);
    end

    % sort linearized outcome
    [ftval, ftidx] = sort(ftval(:), 'descend');
    ftcrd = ftcrd(ftidx);
    tcnum = tcnum(ftidx);

    % select from the top
    tcidx = uint8(0);
    tcidx(max(ftcrd)) = 0;
    ccnum = zeros(1, numcoords);
    nidx = 1;
    while any(ccnum < maxvox) && nidx < numel(tcnum)

        % if not already used
        tc = tcnum(nidx);
        if tcidx(ftcrd(nidx)) == 0 && ccnum(tc) < maxvox

            % select next index
            tcidx(ftcrd(nidx)) = tc;
            ccnum(tc) = ccnum(tc) + 1;

        % otherwise
        else

            % invalidate ftcrd entry
            ftcrd(nidx) = 0;
        end

        % increase index
        nidx = nidx + 1;

        % break anyway
        if ftval(nidx) <= 0
            break;
        end
    end

    % iterate again
    for cc = 1:numcoords

        % select coordinates
        [tc, tcx, tci] = intersect(find(tcidx == cc), ftcrd);
        tc = ftcrd(tci);

        % order according to value
        [tvals, tvali] = sort(ftval(tci), 'descend');
        tc = tc(tvali);

        % set VOI field
        vc.VOI(cc).Voxels = bvcoordconv(tc, 'bvx2tal', bbox);
    end
end

% set other VOI fields
for cc = 1:numel(vc.VOI)
    vc.VOI(cc).NrOfVoxels = size(vc.VOI(cc).Voxels, 1);
    vc.VOI(cc).Name = sprintf('%dvoxel_cluster_at_%d_%d_%d', ...
        vc.VOI(cc).NrOfVoxels, vc.VOI(cc).Voxels(1, :));
    vc.VOI(cc).Color = ceil(256 .* rand(1, 3));
end

% set in content
vc.NrOfVOIs = numel(vc.VOI);
voi.C = vc;
