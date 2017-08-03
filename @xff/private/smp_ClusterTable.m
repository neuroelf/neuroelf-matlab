function [cs, ctab, smp, poi] = smp_ClusterTable(xo, srf, mapno, thresh, opts)
% SMP::ClusterTable  - generate a cluster table for a surface map
%
% FORMAT:       [c, t, m, po] = smp.ClusterTable(srf, mapno, [, thresh [, opts]])
%
% Input fields:
%
%       srf         surface file needed for neighborhood and area
%       mapno       map number (1 .. NrOfMaps)
%       thresh      either p/r values (0 .. 1) or t/F value (1 .. Inf)
%                   if not given or 0, uses the LowerThreshold of map
%       opts        optional settings
%        .altmaps   alternative maps to extract values from (default: [])
%        .cclag     flag, interpret the threshold as lag number (false)
%        .clconn    cluster connectivity ('edge', {'vertex'})
%        .icbm2tal  flag, SMP coords are passed to icbm2tal (default: false)
%        .localmax  break down larger clusters neighbor order (default: Inf)
%        .localmsz  print sub-cluster sizes (default: false)
%        .minsize   minimum cluster size (mm, default by map)
%        .mni2tal   flag, SMP coords are passed to mni2tal (default: false)
%        .showneg   flag, negative values are considered (default: false)
%        .showpos   flag, positive values are considered (default: true)
%        .sorting   either of 'maxstat', {'maxstats'}, 'size', 'x', 'y', 'z'
%        .tdclient  flag, lookup closest talairach label (default false)
%
% Output fields:
%
%       c           Cx1 struct with properties of clusters
%       t           text output
%       m           thresholded map (data for further use)
%       po          if requested, POI structure with vertices
%
% Note: if only one output is requested, the table is text table is
%       returned!!
%
% Note: icbm2tal overrides mni2tal!
%
% Using: clustermeshmap, correlinvtstat, correlpvalue, sdist.

% Version:  v1.1
% Build:    16031221
% Date:     Mar-12 2016, 9:03 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
correlinvtstat = ne_methods.correlinvtstat;
correlpvalue   = ne_methods.correlpvalue;
sdist          = ne_methods.sdist;

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'smp') || ...
    numel(srf) ~= 1 || ~xffisobject(srf, true, 'srf') || ...
   ~isa(mapno, 'double') || numel(mapno) ~= 1 || isinf(mapno) || isnan(mapno) || mapno < 1
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
sc = srf.C;
if mapno > numel(bc.Map)
    error('neuroelf:xff:badArgument', 'Map number out of bounds.');
end

% get bounding box
bb = aft_BoundingBox(srf);

% get map structure and data
mapno = fix(mapno);
map = bc.Map(mapno);
smp = map.SMPData(:);
if isfield(srf.H, 'VertexCoordinateOrig') && ...
    isequal(size(sc.VertexCoordinate), size(srf.H.VertexCoordinateOrig))
    crd = srf.H.VertexCoordinateOrig(:, :);
else
    crd = sc.VertexCoordinate(:, :);
end
tri = sc.TriangleVertex(:, :);

% continue argument check
if nargin < 4 || ~isa(thresh, 'double') || numel(thresh) ~= 1 || ...
    isinf(thresh) || isnan(thresh) || thresh < 0
    thresh = 0;
end
tfrommap = false;
if thresh == 0
    thresh = map.LowerThreshold;
    tfrommap = true;
end
if nargin < 5 || numel(opts) ~= 1 || ~isstruct(opts)
    opts = struct;
end
if ~isfield(opts, 'altmaps') || ~isa(opts.altmaps, 'double')
    opts.altmaps = [];
else
    opts.altmaps = opts.altmaps(:)';
    opts.altmaps(isinf(opts.altmaps) | isnan(opts.altmaps) | ...
        opts.altmaps < 1 | opts.altmaps > numel(bc.Map)) = [];
    altmaps = unique(round(opts.altmaps(:)));
    opts.altmapsp = true(1, numel(altmaps));
    opts.altmaps = double(zeros([numel(smp), numel(altmaps)]));
    for amc = 1:numel(altmaps)
        opts.altmaps(:, amc) = bc.Map(altmaps(amc)).SMPData(:);
    end
end
if ~isfield(opts, 'cclag') || ~islogical(opts.cclag) || numel(opts.cclag) ~=1
    opts.cclag = false;
end
if ~isfield(opts, 'clconn') || ~ischar(opts.clconn) || ...
   ~any(strcmpi(opts.clconn, {'edge', 'vertex'}))
    opts.clconn = 'vertex';
else
    opts.clconn = lower(opts.clconn(:)');
end
if ~isfield(opts, 'icbm2tal') || numel(opts.icbm2tal) ~= 1 || ...
   (~isa(opts.icbm2tal, 'double') && ~islogical(opts.icbm2tal))
    opts.icbm2tal = false;
else
    opts.icbm2tal = true && opts.icbm2tal;
end
if ~isfield(opts, 'localmax') || ~isa(opts.localmax, 'double') || ...
    numel(opts.localmax) ~= 1 || isnan(opts.localmax) || opts.localmax < 2
    opts.localmax = Inf;
else
    opts.localmax = round(opts.localmax);
end
if ~isfield(opts, 'localmsz') || ~islogical(opts.localmsz) || numel(opts.localmsz) ~= 1
    opts.localmsz = false;
end
if ~isfield(opts, 'minsize') || ~isa(opts.minsize, 'double') || numel(opts.minsize) ~= 1 || ...
    isinf(opts.minsize) || isnan(opts.minsize) || opts.minsize < 0 || opts.minsize > 1e5
    if map.EnableClusterCheck
        opts.minsize = map.ClusterSize;
    else
        opts.minsize = 0;
    end
end
if ~isfield(opts, 'mni2tal') || numel(opts.mni2tal) ~= 1 || ...
   (~isa(opts.mni2tal, 'double') && ~islogical(opts.mni2tal)) || opts.icbm2tal
    opts.mni2tal = false;
else
    opts.mni2tal = true && opts.mni2tal;
end
if ~isfield(opts, 'showneg') || ~islogical(opts.showneg) || numel(opts.showneg) ~= 1
    opts.negative = (map.ShowPositiveNegativeFlag > 1);
else
    opts.negative = opts.showneg;
end
if ~isfield(opts, 'showpos') || ~islogical(opts.showpos) || numel(opts.showpos) ~= 1
    opts.positive = (mod(map.ShowPositiveNegativeFlag, 2) == 1);
else
    opts.positive = opts.showpos;
end
if ~isfield(opts, 'sorting') || ~ischar(opts.sorting) || isempty(opts.sorting) || ...
   ~any(strcmpi(opts.sorting, {'maxstat', 'maxstats', 'size', 'x', 'y', 'z'}))
    opts.sorting = 'maxstats';
else
    opts.sorting = lower(opts.sorting(:)');
end
if ~isfield(opts, 'tdclient') || ~islogical(opts.tdclient) || numel(opts.tdclient) ~= 1
    opts.tdclient = false;
end

% some initial checks on threshold
if thresh < 0
    if nargin < 5
        opts.negative = true;
        opts.positive = false;
    end
    thresh = -thresh;
end

% remove invalid entries first
smp(isinf(smp) | isnan(smp)) = 0;

% put some additional fields in opts
opts.mat = bb.QuatB2T;
opts.mat(1:3, 4) = 128;
switch (map.Type)
    case 1
        thtype = 't-Map';
        thdegf = sprintf('%d', map.DF1);
        opts.tptype = 't';
        opts.tptypedf = map.DF1;
    case 2
        thtype = 'correlation Map';
        opts.tptype = 'r';
        if bc.FileVersion < 4
            thdegf = sprintf('%d (:= t[%d])', map.DF1, map.DF1 - 2);
            opts.tptypedf = map.DF1 - 2;
        else
            thdegf = sprintf('%d', map.DF1);
            opts.tptypedf = map.DF1;
        end
    case 3
        thtype = 'CC Map';
        thdegf = sprintf('%d', map.DF1);

        % in this case also patch the VMP data, and set fields
        if ~opts.cclag
            smp = mod(smp + 2, 1000) - 2;
            opts.tptype = 'r';
            if bc.FileVersion < 4
                thdegf = sprintf('%d (:= t[%d])', map.DF1, map.DF1 - 2);
                opts.tptypedf = map.DF1 - 2;
            else
                thdegf = sprintf('%d', map.DF1);
                opts.tptypedf = map.DF1;
            end
        else
            smp = floor(0.001 * (smp + 2));
        end
    case 4
        thtype = 'F-Map';
        thdegf = sprintf('%d, %d', map.DF1, map.DF1);
        opts.tptype = 'f';
        opts.tptypedf = [map.DF1, map.DF2];
    case 11
        thtype = 'PSC Map';
        thdegf = 'n/a';
    case 12
        thtype = 'z(ica) Map';
        thdegf = sprintf('%d', map.DF1);
    case 13
        thtype = 'Thickness Map';
        thdegf = 'n/a';
    case 16
        thtype = 'Probability Map';
        thdegf = sprintf('%d', map.DF1);
    case 20
        thtype = 'Mean Diffusivity Map';
        thdegf = 'n/a';
    case 21
        thtype = 'Fractional Anisotropy Map';
        thdegf = 'n/a';
    otherwise
        thtype = 'unknown';
        thdegf = 'n/a';
end
mthreshp = nan;
if ~tfrommap && any([1, 2, 4] == map.Type) && ...
    (thresh <= 0.05 || (thresh < 0.2 && opts.tptype ~= 'r'))
    mthreshp = thresh;
    switch (opts.tptype)
        case 'f'
            thresh = sdist('finv', 1 - thresh, opts.tptypedf(1), opts.tptypedf(2));
        case 'r'
            if opts.negative && opts.positive
                thresh = correlinvtstat(-sdist('tinv', 0.5 * thresh, ...
                    opts.tptypedf), opts.tptypedf + 2);
            else
                thresh = correlinvtstat(-sdist('tinv', thresh, ...
                    opts.tptypedf), opts.tptypedf + 2);
            end
        case 't'
            if opts.negative && opts.positive
                thresh = -sdist('tinv', 0.5 * thresh, opts.tptypedf);
            else
                thresh = -sdist('tinv', thresh, opts.tptypedf);
            end
    end
elseif any([1, 2, 4] == map.Type)
    switch (opts.tptype)
        case 'f'
            mthreshp = 1 - sdist('fcdf', thresh, opts.tptypedf(1), opts.tptypedf(2));
        case 'r'
            if opts.negative && opts.positive
                mthreshp = correlpvalue(thresh, opts.tptypedf + 2);
            else
                mthreshp = 0.5 * correlpvalue(thresh, opts.tptypedf + 2);
            end
        case 't'
            if opts.negative && opts.positive
                mthreshp = 2 - 2 * sdist('tcdf', thresh, opts.tptypedf);
            else
                mthreshp = 1 - sdist('tcdf', thresh, opts.tptypedf);
            end
    end
end
thead = cell(1, 5);
thead{1} = sprintf('   Clustertable of map:   "%s"', map.Name);
thead{2} = sprintf('           Type of map:   %s', thtype);
thead{3} = sprintf('    Degrees of freedom:   %s', thdegf);
thead{4} = sprintf('   Cluster k-threshold:   %d u^2', opts.minsize);
thead{5} = sprintf(' Applied map threshold:   %.5f (p < %.5f)', thresh, mthreshp);

% perform clustering
[ctab, cs, smp] = ne_methods.clustermeshmap(smp, tri, crd, thresh, opts.minsize, opts);

% store back in struct if cluster-threshold enabled
if opts.minsize > 0
    bc.Map(mapno).SMPDataCT = (smp ~= 0);
    xo.C = bc;
end

% extend ctab
ctab = sprintf('%s\n%s\n%s\n%s\n%s\n\n%s', thead{:}, ctab);

% create VOI structure ?
if nargout > 3
    poi = xff('new:poi');
    poic = poi.C;
    poic.FromMeshFile = srf.F;
    poic.NrOfMeshVertices = size(crd, 1);
    if ~isempty(cs)
        poic.POI(numel(cs)).Name = '';
    else
        poic.POI(:) = [];
    end
    for cc = 1:numel(cs)
        poic.POI(cc).Name = sprintf('Cluster%04d_%6.1f_%6.1f_%6.1f_%s', cc, ...
            cs(cc).peak(1, :), strrep(regexprep(char(cs(cc).talout), '\(.*$', ''), ' ', '_'));
        poic.POI(cc).InfoTextFile = '<none>';
        if cs(cc).values(1) >= 0
            poic.POI(cc).Color = floor([255, 127.999 * rand(1, 2)]);
        else
            poic.POI(cc).Color = floor([127.999 * rand(1, 2), 255]);
        end
        poic.POI(cc).LabelVertex = cs(cc).vertices(1);
        poic.POI(cc).NrOfVertices = numel(cs(cc).vertices);
        poic.POI(cc).Vertices = cs(cc).vertices;
        poic.POI(cc).VertexValues = cs(cc).values;
    end
    poic.NrOfPOIs = numel(poic.POI);
    poi.C = poic;
end

% return table if only one output is requested
if nargout < 2
    cs = ctab;
end
