function [table, list, tm] = clustermeshmap(m, tri, crd, t, st, opts)
% clustermeshmap  - cluster a SMP/SRF based map
%
% FORMAT:       [tab, list, tm] = clustermeshmap(m, tri, c, t, st, [, opts])
%
% Input fields:
%
%       m           double-value (Vx1)
%       tri         Tx3 1-based triangle list
%       t           threshold value
%       st          size threshold (in coordinate units^2)
%       crd         Vx3 coordinates (for area computation)
%       opts        optional settings
%        .altmaps   alternative maps (2D data, along 2nd dim)
%        .altmapsp  print alt maps values into table
%        .clconn    connectivity (either 'edge' or {'vertex'})
%        .icbm2tal  apply icbm2tal transform for tdclient (default: false)
%        .localmax  break down large clusters neighbor order (default: Inf)
%        .localmsz  print sub-cluster sizes (default: false)
%        .mat       4x4 transformation matrix applied to crd
%        .mni2tal   apply mni2tal transform for tdclient (default: false)
%        .negative  cluster negative end of distribution (default: false)
%        .pbar      either xprogress or xfigure:XProgress object
%        .pbarrange default: [0, 1]
%        .positive  cluster positive end of distribution (default: true)
%        .sorting   sorting of output clusters ({'maxstat'}, 'size', 'x'}
%        .tdclient  run tdclient on coords (requires .mat, default: false)
%        .tptype    convert threshold (if not raw: 'F', 'r', 't')
%        .tptypedf  required D.F. for conversion
%
% Output fields:
%
%       tab         character table output
%       list        Cx1 structure of clusters with fields
%        .coords    vertex coordinates
%        .localmax  flag set to 'L' if a local max/min within a cluster
%        .max       peak value
%        .mean      mean value
%        .peak      peak coordinate
%        .peakalt   values of peak in alternative maps (empty if not given)
%        .peakvtx   peak vertex index
%        .size      cluster size in coordinate units^2
%        .talcoords rounded TAL coords (only if .mni2tal is true)
%        .talout    tdclient output for peak
%        .talpeak   rounded TAL peak coordinate
%        .values    values of vertices
%        .vertices  vertex indices
%       tm          map with clustered mask
%
% Note: currently only the vertex connectivity is implemented!
%
% Note: icbm2tal overrides mni2tal!

% Version:  v0.9b
% Build:    13071714
% Date:     Apr-09 2011, 1:55 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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
   ~isnumeric(m) || ...
    isempty(m) || ...
    numel(m) ~= size(m, 1) || ...
   ~isa(tri, 'double') || ...
    ndims(tri) ~= 2 || ...
    size(tri, 2) ~= 3 || ...
    size(tri, 1) < numel(m) || ...
    any(isinf(tri(:)) | isnan(tri(:)) | tri(:) < 1 | tri(:) > numel(m)) || ...
   ~isa(crd, 'double') || ...
    size(crd, 1) ~= size(m, 1) || ...
    numel(crd) ~= (3 * size(m, 1)) || ...
    any(isinf(crd(:)) | isnan(crd(:))) || ...
   ~isnumeric(t) || ...
    numel(t) ~= 1 || ...
    isinf(t) || ...
    isnan(t) || ...
    t <= 0 || ...
   ~isnumeric(st) || ...
    numel(st) ~= 1 || ...
    isinf(st) || ...
    isnan(st) || ...
    st < 0
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument.' ...
    );
end
m = double(m);
mn = numel(m);
t = double(t);
st = double(st);
if nargin < 6 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'altmaps') || ...
   ~isnumeric(opts.altmaps) || ...
    ndims(opts.altmaps) ~= 2 || ...
    size(opts.altmaps, 1) ~= mn
    opts.altmaps = zeros([mn, 0]);
else
    opts.altmaps = double(opts.altmaps(:, :));
end
ea = numel(opts.altmaps) - 1;
na = size(opts.altmaps, 2);
if ~isfield(opts, 'altmapsp') || ...
   ~islogical(opts.altmapsp) || ...
    isempty(opts.altmapsp) || ...
    size(opts.altmapsp, 1) ~= 1 || ...
   ~any([1, na] == size(opts.altmapsp, 1))
    opts.altmapsp = false(1, na);
elseif numel(opts.altmapsp) ~= na
    opts.altmapsp = opts.altmapsp(1, ones(1, na));
end
altmapsh = cell(1, sum(opts.altmapsp));
altmapsd = '-';
altmapsd = altmapsd(1, ones(1, 12 * numel(altmapsh)));
altmapsf = repmat('| %9.6f ', 1, numel(altmapsh));
amt = 0;
for amc = 1:numel(opts.altmapsp)
    if opts.altmapsp(amc)
        amt = amt + 1;
        altmapsh{amt} = sprintf('|   am#%02d   ', amc);
    end
end
altmapsh = sprintf('%s', altmapsh{:});
if ~isfield(opts, 'clconn') || ...
   ~ischar(opts.clconn) || ...
   ~any(strcmpi(opts.clconn(:)', {'edge', 'vertex'}))
    opts.clconn = 'vertex';
else
    opts.clconn = lower(opts.clconn(:))';
end
if opts.clconn(1) == 'v'
    clconn = 1;
else
    clconn = 1;
end
if ~isfield(opts, 'icbm2tal') || ...
   ~islogical(opts.icbm2tal) || ...
    numel(opts.icbm2tal) ~= 1
    opts.icbm2tal = false;
end
if ~isfield(opts, 'localmax') || ...
   ~isa(opts.localmax, 'double') || ...
    numel(opts.localmax) ~= 1 || ...
   ~any((2:6) == opts.localmax)
    opts.localmax = Inf;
end
if ~isfield(opts, 'localmsz') || ...
   ~islogical(opts.localmsz) || ...
    numel(opts.localmsz) ~= 1
    opts.localmsz = false;
end
if ~isfield(opts, 'mat') || ...
   ~isa(opts.mat, 'double') || ...
   ~isequal(size(opts.mat), [4, 4]) || ...
    any(isinf(opts.mat(:)) | isnan(opts.mat(:))) || ...
    any(opts.mat(4, :) ~= [0, 0, 0, 1])
    opts.mat = [];
else
    crd(:, 4) = 1;
    crd = crd * opts.mat';
    crd(:, 4) = [];
end
if ~isfield(opts, 'mni2tal') || ...
   ~islogical(opts.mni2tal) || ...
    numel(opts.mni2tal) ~= 1 || ...
    opts.icbm2tal
    opts.mni2tal = false;
end
if ~isfield(opts, 'negative') || ...
   ~islogical(opts.negative) || ...
    numel(opts.negative) ~= 1
    opts.negative = false;
end
if ~isfield(opts, 'pbar') || ...
    numel(opts.pbar) ~= 1 || ...
   ~any(strcmpi(class(opts.pbar), {'xfigure', 'xprogress'}))
    opts.pbar = [];
end
if ~isfield(opts, 'pbarrange') || ...
   ~isa(opts.pbarrange, 'double') || ...
    numel(opts.pbarrange) ~= 2 || ...
    any(isinf(opts.pbarrange) | isnan(opts.pbarrange) | ...
        opts.pbarrange < 0 | opts.pbarrange > 1)
    opts.pbarrange = [0, 1];
else
    opts.pbarrange = sort(opts.pbarrange(:))';
end
if ~isfield(opts, 'positive') || ...
   ~islogical(opts.positive) || ...
    numel(opts.positive) ~= 1
    opts.positive = true;
end
if ~isfield(opts, 'sorting') || ...
   ~ischar(opts.sorting) || ...
   ~any(strcmpi(opts.sorting(:)', {'maxstat', 'maxstats', 'size', 'x', 'y', 'z'}))
    opts.sorting = 'maxstat';
else
    opts.sorting = lower(opts.sorting(:)');
end
if ~isfield(opts, 'tdclient') || ...
   ~islogical(opts.tdclient) || ...
    numel(opts.tdclient) ~= 1
    opts.tdclient = false;
end
if ~isfield(opts, 'tptype') || ...
   ~ischar(opts.tptype) || ...
   ~any(strcmpi(opts.tptype(:)', {'f', 'r', 'raw', 't'}))
    opts.tptype = 'raw';
else
    opts.tptype = lower(opts.tptype(:)');
end
if ~isfield(opts, 'tptypedf') || ...
   ~isa(opts.tptypedf, 'double') || ...
    isempty(opts.tptypedf) || ...
    numel(opts.tptypedf) > 2 || ...
    any(isinf(opts.tptypedf) | isnan(opts.tptypedf) | opts.tptypedf < 1)
    if opts.tptype(1) ~= 'f'
        opts.tptypedf = 1e7;
    else
        opts.tptypedf = [1, 1e7];
    end
end

% convert threshold
if opts.tptype(1) ~= 'r' && ...
    t <= 0.2
    if strcmp(opts.tptype, 'f') && ...
        numel(opts.tptypedf) ~= 2
        error( ...
            'neuroelf:BadArgument', ...
            'F-statistic needs 1x2 d.f. option in tptypedf field.' ...
        );
    end
    switch (opts.tptype)
        case {'f'}
            t = sdist('finv', 1 - t, opts.tptypedf(1), opts.tptypedf(2));
        case {'t'}
            if opts.negative && ...
                opts.positive
                t = -sdist('tinv', 0.5 * t, opts.tptypedf);
            else
                t = -sdist(t, opts.tptypedf);
            end
    end
end

% create output volume
if nargout > 2
    tm = zeros(mn, 1);
end

% get neighbors and triangle back-reference
try
    [nei, bn, trb] = mesh_trianglestoneighbors(mn, tri);
    if ~isempty(bn)
        warning( ...
            'neuroelf:BadSurface', ...
            'Cluster sizes potentially flawed. %d bad neighborhoods!', ...
            numel(bn) ...
        );
    end
    if isempty(nei{end}) || ...
        isempty(trb{end})
        error('BAD_SURFACE');
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:BadSurface', ...
        'Invalid surface, neighborhood references invalid.' ...
    );
end

% compute areas of triangles
tsa = sqrt(sum((crd(tri(:, 1), :) - crd(tri(:, 2), :)) .^ 2, 2));
tsb = sqrt(sum((crd(tri(:, 1), :) - crd(tri(:, 3), :)) .^ 2, 2));
tsc = sqrt(sum((crd(tri(:, 2), :) - crd(tri(:, 3), :)) .^ 2, 2));
tss = 0.5 * (tsa + tsb + tsc);
tra = sqrt(tss .* (tss - tsa) .* (tss - tsb) .* (tss - tsc));

% cluster positive and negative clusters
if opts.positive
    [pcs, pcv, pcl, pcc, pst] = ...
        clustermeshmapbin(m >= t, nei(:, 2), crd, tra, trb, st, clconn);
    pmv = zeros(size(pcs));
    if nargout > 2
        tm(pcv > 0) = m(pcv > 0);
    end
    pcc = cell(numel(pcs), 1);
    pcv = cell(numel(pcs), 1);
    for cc = 1:numel(pcs)
        mc = pst(cc).vertices;
        mv = m(mc);
        [sc, si] = sort(mv, 'descend');
        pcc{cc} = mc(si);
        pst(cc).vertices = pcc{cc};
        pst(cc).coords = pst(cc).coords(si, :);
        pcv{cc} = sc;
        pmv(cc) = sc(1);
    end
else
    pcc = cell(0, 1);
    pcs = zeros(0, 1);
    pcv = pcc;
    pmv = zeros(0, 1);
    pst = emptystruct({'coords', 'vertices'}, [0, 1]);
end
if opts.negative
    [ncs, ncv, ncl, ncc, nst] = ...
        clustermeshmapbin(m <= -t, nei(:, 2), crd, tra, trb, st, clconn);
    nmv = zeros(size(ncs));
    if nargout > 2
        tm(ncv > 0) = m(ncv > 0);
    end
    ncc = cell(numel(ncs), 1);
    ncv = cell(numel(ncs), 1);
    for cc = 1:numel(ncs)
        mc = nst(cc).vertices;
        mv = m(mc);
        [sc, si] = sort(mv, 'ascend');
        ncc{cc} = mc(si);
        nst(cc).vertices = ncc{cc};
        nst(cc).coords = nst(cc).coords(si, :);
        ncv{cc} = sc;
        nmv(cc) = sc(1);
    end
else
    ncc = cell(0, 1);
    ncs = zeros(0, 1);
    ncv = ncc;
    nmv = zeros(0, 1);
    nst = emptystruct({'coords', 'vertices'}, [0, 1]);
end

% join results
cc = [pcc(:); ncc(:)];
cs = [pcs(:); ncs(:)];
cv = [pcv(:); ncv(:)];
mv = [pmv(:); nmv(:)];
str = [pst(:); nst(:)];
st = true(1, numel(cc));

% sorting
switch (opts.sorting)
    case {'maxstat'}
        [sortv, sorti] = sort(abs(mv), 'descend');
    case {'maxstats'}
        [sortv, sorti] = sort(mv, 'descend');
        sortn = sorti(sortv < 0);
        sorti = [sorti(sortv >= 0); sortn(end:-1:1)];
    case {'size'}
        [sortv, sorti] = sort(cs, 'descend');
    case {'x'}
        xc = zeros(1, numel(cc));
        for xcc = 1:numel(xc)
            xc(xcc) = cc{xcc}(1, 1);
        end
        [sortv, sorti] = sort(xc);
    case {'y'}
        xc = zeros(1, numel(cc));
        for xcc = 1:numel(xc)
            xc(xcc) = cc{xcc}(1, 2);
        end
        [sortv, sorti] = sort(xc);
    case {'z'}
        xc = zeros(1, numel(cc));
        for xcc = 1:numel(xc)
            xc(xcc) = cc{xcc}(1, 3);
        end
        [sortv, sorti] = sort(xc);
end
cc = cc(sorti);
cs = cs(sorti);
cv = cv(sorti);
mv = mv(sorti);
str = str(sorti);

% split local max
if ~isinf(opts.localmax)
    c = 1;
    while c <= numel(cc)
        if opts.localmax <= size(cc{c}, 1)
%             [scs, scv, scl, scc] = splitclustercoords(cc{c}, abs(cv{c}), k);
%             if numel(scs) > 1
%                 smv = zeros(numel(scs), 1);
%                 if cv{c}(1) > 0
%                     for suc = 2:numel(scv)
%                         smv(suc) = scv{suc}(1);
%                     end
%                 else
%                     for suc = 2:numel(scv)
%                         scv{suc} = -scv{suc};
%                         smv(suc) = scv{suc}(1);
%                     end
%                 end
%                 st = [st(1:c), false(1, numel(scs) -1), st((c+1):end)];
%                 cc = [cc(1:c); scc(2:end); cc((c+1):end)];
%                 cv = [cv(1:c); scv(2:end); cv((c+1):end)];
%                 mv = [mv(1:c); smv(2:end); mv((c+1):end)];
%                 c = c + numel(scs) - 1;
%             end
        end
        c = c + 1;
    end
end
nc = numel(cc);

% create output structure
list = emptystruct({ ...
    'coords', ...
    'localmax', 'max', 'mean', ...
    'peak', 'peakalt', 'peakvtx', ...
    'size', ...
    'talcoords', 'talout', 'talpeak', ...
    'values', 'vertices'}, [nc, 1]);
table = cell(1, nc);
taltext = '';

% process results
% process results
if opts.tdclient && ...
    nc > 7
    pbn = opts.pbarrange(1);
    pbx = opts.pbarrange(2) - pbn + eps;
    if isempty(opts.pbar)
        try
            pbar = xprogress;
            xprogress(pbar, 'setposition', [80, 200, 640, 36]);
            xprogress(pbar, 'settitle', 'Retrieving TAL labels...');
            xprogress(pbar, pbn, ...
                sprintf('Getting label 1/%d...', nc), 'visible', 0, 1);
            pbarn = '';
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            pbar = [];
        end
    else
        pbar = opts.pbar;
        pbar.Progress(pbn, sprintf('tdclient: Getting label 1/%d...', nc));
        pbarn = 'tdclient: ';
    end
else
    pbar = [];
end
for c = 1:nc
    crd = str(c).coords(1, :);
    list(c).coords = str(c).coords;
    list(c).vertices = str(c).vertices;
    if st(c)
        list(c).localmax = ' ';
    else
        list(c).localmax = 'L';
    end
    list(c).max = mv(c);
    list(c).mean = sum(cv{c}) / numel(cv{c});
    list(c).peak = crd;
    list(c).peakalt = opts.altmaps((0:mn:ea) + list(c).vertices(1));
    if opts.icbm2tal
        list(c).talcoords = icbm2tal(str(c).coords);
        list(c).talpeak = list(c).talcoords(1, :);
        crd = list(c).talpeak;
    elseif opts.mni2tal
        list(c).talcoords = mni2tal(str(c).coords);
        list(c).talpeak = list(c).talcoords(1, :);
        crd = list(c).talpeak;
    end
    if opts.tdclient
        if isempty(list(c).talpeak)
            list(c).talpeak = crd;
        end
        talta = tdlabel(round(crd));
        if ~isempty(pbar)
            pbar.Progress(pbn + pbx * c / nc, ...
                sprintf('%sGetting label %d/%d...', pbarn, c + 1, nc));
            pbar.Visible = 'on';
        end
        taltext = talta{1};
    end
    list(c).size = cs(c);
    list(c).talout = taltext;
    list(c).values = cv{c};
    if ~isempty(altmapsf)
        if opts.localmsz || ...
            list(c).localmax == ' '
            table{c} = sprintf(['%6.1f %6.1f %6.1f %s|  %7.1f  | %9.6f | %9.6f  ' altmapsf '|%s'], ...
                        crd(1), crd(2), crd(3), list(c).localmax, ...
                        list(c).size, list(c).max, list(c).mean, ...
                        list(c).peakalt(opts.altmapsp), taltext);
        else
            table{c} = sprintf(['%6.1f %6.1f %6.1f %s|           | %9.6f | %9.6f  ' altmapsf '|%s'], ...
                        crd(1), crd(2), crd(3), list(c).localmax, ...
                        list(c).max, list(c).mean, ...
                        list(c).peakalt(opts.altmapsp), taltext);
        end
    else
        if opts.localmsz || ...
            list(c).localmax == ' '
            table{c} = sprintf('%6.1f %6.1f %6.1f %s|  %7.1f  | %9.6f | %9.6f  |%s', ...
                        crd(1), crd(2), crd(3), list(c).localmax, ...
                        list(c).size, list(c).max, list(c).mean, taltext);
        else
            table{c} = sprintf('%6.1f %6.1f %6.1f %s|           | %9.6f | %9.6f  |%s', ...
                        crd(1), crd(2), crd(3), list(c).localmax, ...
                        list(c).max, list(c).mean, taltext);
        end
    end
end
if ~isempty(pbar) && ...
    isempty(opts.pbar);
    closebar(pbar);
end

% join table
table = sprintf( ...
    ['     x      y      z  |      u^2  |    max    |    mean    %s| tdclient\n' ...
     '----------------------------------------------------------%s---------\n%s\n'], ...
     altmapsh, altmapsd, gluetostringc(table, char(10)));
