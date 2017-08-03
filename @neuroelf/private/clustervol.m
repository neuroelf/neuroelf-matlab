function [table, list, vol] = clustervol(v, t, k, opts)
% clustervol  - cluster volume with threshold
%
% FORMAT:       [tab, list, vol] = clustervol(vol, thresh, kthresh [, opts])
%
% Input fields:
%
%       vol         3d numeric volume
%       thresh      numeric threshold
%       kthresh     cluster threshold k (minimum cluster size)
%       opts        optional setting
%        .altmaps   alternative maps (4D data, along 4th dim)
%        .altmapsp  print alt maps values into table (either 1x1 or 1xM)
%        .altstat   either of 'mean' or {'peak'}
%        .clconn    cluster connectivity ('face', {'edge'}, 'vertex')
%        .icbm2tal  apply icbm2tal transform for tdclient (default: false)
%        .localmax  look for local maxima if size > (default: Inf)
%        .localmaxi iterate on sub-clusters (default: false)
%        .localmsz  print size of local maxima (default: false)
%        .lupcrd    look-up coordinate, either or 'center', 'cog', {'peak'}
%        .mat       4x4 transformation matrix to compute real-world coords
%        .mni2tal   apply mni2tal transform for tdclient (default: false)
%        .negative  cluster negative end of distribution (default: false)
%        .pbar      either xprogress or xfigure:XProgress object
%        .pbarrange default: [0, 1]
%        .positive  cluster positive end of distribution (default: true)
%        .sorting   sorting of output clusters ({'maxstat'}, 'size', 'x'}
%        .tdclient  run tdclient on coords (requires .mat, default: false)
%        .tpnoconv  flag, do not use threshold conversion (default: false)
%        .tptype    convert threshold (if not raw: 'F', 'r', 't')
%        .tptypedf  required D.F. for conversion
%
% Output fields:
%
%       tab         character table output
%       list        Cx1 structure of clusters with fields
%        .center    central coordinate (mean)
%        .cog       weighted central coordinate (center-of-gravity)
%        .coords    voxel coordinates
%        .localmax  flag set to 'L' if a local max/min within a cluster
%        .max       peak value
%        .mean      mean value
%        .meanalt   mean values of coordinates in alternative maps
%        .peak      peak coordinate
%        .peakalt   values of peak in alternative maps (empty if not given)
%        .rwcenter  real-world central coordinate (mean)
%        .rwcog     real-world center-of-gravity
%        .rwcoords  real-world coordinates (only filled if .mat is given)
%        .rwpeak    real-world peak (only set if .mat is given)
%        .rwsize    real-world cluster size
%        .size      cluster size
%        .talcenter central TAL coordinate (mean)
%        .talcog    weighted central TAL coordinate (center-of-gravity)
%        .talcoords rounded TAL coords (only if .mat and .mni2tal given)
%        .talout    tdclient output for peak
%        .talpeak   rounded TAL peak coordinate
%        .values    values of voxels (in their order)
%       vol         vol with clustered mask

% Version:  v0.9b
% Build:    12081711
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

% check arguments
if nargin < 3 || ...
   ~isnumeric(v) || ...
    isempty(v) || ...
    ndims(v) < 3 || ...
    any(size(v) < 3) || ...
   ~isnumeric(t) || ...
    numel(t) ~= 1 || ...
    isinf(t) || ...
    isnan(t) || ...
    t <= 0 || ...
   ~isnumeric(k) || ...
    numel(k) ~= 1 || ...
    isinf(k) || ...
    isnan(k) || ...
    k < 1
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument.' ...
    );
end
v = double(v);
vs = size(v);
ns = prod(vs);
t = double(t);
k = round(double(k));
if nargin < 4 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'altmaps') || ...
   ~isnumeric(opts.altmaps) || ...
    ndims(opts.altmaps) < 3 || ...
    size(opts.altmaps, 1) ~= vs(1) || ...
    size(opts.altmaps, 2) ~= vs(2) || ...
    size(opts.altmaps, 3) ~= vs(3)
    opts.altmaps = zeros([vs, 0]);
else
    opts.altmaps = double(opts.altmaps);
end
ea = numel(opts.altmaps) - 1;
nsea = 0:ns:ea;
nnsea = numel(nsea);
na = size(opts.altmaps, 4);
if ~isfield(opts, 'altmapsp') || ...
   ~islogical(opts.altmapsp) || ...
    isempty(opts.altmapsp) || ...
    size(opts.altmapsp, 1) ~= 1 || ...
   ~any([1, na] == numel(opts.altmapsp))
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
if ~isfield(opts, 'altstat') || ...
   ~ischar(opts.altstat) || ...
    isempty(opts.altstat) || ...
   ~any(lower(opts.altstat(1)) == 'mp')
    opts.altstat = 'peak';
else
    opts.altstat = lower(opts.altstat(1));
end
if ~isfield(opts, 'clconn') || ...
   ~ischar(opts.clconn) || ...
   ~any(strcmpi(opts.clconn(:)', {'e', 'edge', 'f', 'face', 'v', 'vertex'}))
    opts.clconn = 'edge';
else
    opts.clconn = lower(opts.clconn(:)');
end
switch (opts.clconn(1))
    case {'e'}
        opts.clconn = 2;
    case {'f'}
        opts.clconn = 1;
    case {'v'}
        opts.clconn = 3;
end
if ~isfield(opts, 'icbm2tal') || ...
   ~islogical(opts.icbm2tal) || ...
    numel(opts.icbm2tal) ~= 1
    opts.icbm2tal = false;
end
if ~isfield(opts, 'localmax') || ...
   ~isa(opts.localmax, 'double') || ...
    numel(opts.localmax) ~= 1 || ...
    isinf(opts.localmax) || ...
    isnan(opts.localmax) || ...
    opts.localmax < 3
    opts.localmax = Inf;
else
    opts.localmax = max(4, round(opts.localmax));
end
if ~isfield(opts, 'localmaxi') || ...
   ~islogical(opts.localmaxi) || ...
    numel(opts.localmaxi) ~= 1
    opts.localmaxi = false;
end
if ~isfield(opts, 'localmin') || ...
   ~isa(opts.localmin, 'double') || ...
    numel(opts.localmin) ~= 1 || ...
    isinf(opts.localmin) || ...
    isnan(opts.localmin) || ...
    opts.localmin < 1
    opts.localmin = [];
else
    opts.localmin = round(real(opts.localmin));
end
if ~isfield(opts, 'localmsz') || ...
   ~islogical(opts.localmsz) || ...
    numel(opts.localmsz) ~= 1
    opts.localmsz = false;
end
if ~isfield(opts, 'lupcrd') || ...
   ~ischar(opts.lupcrd) || ...
   ~any(strcmpi(opts.lupcrd(:)', {'center', 'cog', 'peak'}))
    opts.lupcrd = 'peak';
else
    opts.lupcrd = lower(opts.lupcrd(:)');
end
if ~isfield(opts, 'mat') || ...
   ~isa(opts.mat, 'double') || ...
   ~isequal(size(opts.mat), [4, 4]) || ...
    any(isinf(opts.mat(:)) | isnan(opts.mat(:))) || ...
    any(opts.mat(4, :) ~= [0, 0, 0, 1])
    opts.mat = [];
else
    vxsf = prod(sqrt(sum(opts.mat(1:3, 1:3) .^ 2, 1)));
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
if ~isfield(opts, 'tpnoconv') || ...
   ~islogical(opts.tpnoconv) || ...
    numel(opts.tpnoconv) ~= 1
    opts.tpnoconv = false;
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
if ~strcmp(opts.tptype, 'raw') && ...
    t <= 0.1 && ...
   ~opts.tpnoconv
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
        case {'r'}
            if opts.negative && ...
                opts.positive
                t = correlinvtstat(-sdist('tinv', 0.5 * t, opts.tptypedf), ...
                    opts.tptypedf + 2);
            else
                t = correlinvtstat(-sdist('tinv', t, opts.tptypedf), ...
                    opts.tptypedf + 2);
            end
        case {'t'}
            if opts.negative && ...
                opts.positive
                t = -sdist('tinv', 0.5 * t, opts.tptypedf);
            else
                t = -sdist('tinv', t, opts.tptypedf);
            end
    end
end

% create output volume
if nargout > 2
    vol = zeros(vs);
end

% cluster negative and positive clusters
if opts.positive
    [pcs, pcv, pcl, pcc] = clustercoordsc(v >= t, opts.clconn, k);
    pmv = zeros(size(pcs));
    if nargout > 2
        vol(pcv > 0) = v(pcv > 0);
    end
    pcv = cell(numel(pcs), 1);
    for cc = 1:numel(pcs)
        mc = sub2ind(vs, pcc{cc}(:, 1), pcc{cc}(:, 2), pcc{cc}(:, 3));
        mv = v(mc);
        [sc, si] = sort(mv, 'descend');
        pcc{cc} = pcc{cc}(si, :);
        pcv{cc} = sc;
        pmv(cc) = sc(1);
    end
else
    pcc = cell(0, 1);
    pcs = zeros(0, 1);
    pcv = pcc;
    pmv = zeros(0, 1);
end
if opts.negative
    [ncs, ncv, ncl, ncc] = clustercoordsc(-v >= t, opts.clconn, k);
    nmv = zeros(size(ncs));
    if nargout > 2
        vol(ncv > 0) = v(ncv > 0);
    end
    ncv = cell(numel(ncs), 1);
    for cc = 1:numel(ncs)
        mc = sub2ind(vs, ncc{cc}(:, 1), ncc{cc}(:, 2), ncc{cc}(:, 3));
        mv = v(mc);
        [sc, si] = sort(mv, 'ascend');
        ncc{cc} = ncc{cc}(si, :);
        ncv{cc} = sc;
        nmv(cc) = sc(1);
    end
else
    ncc = cell(0, 1);
    ncs = zeros(0, 1);
    ncv = ncc;
    nmv = zeros(0, 1);
end

% join results
cc = [pcc(:); ncc(:)];
cs = [pcs(:); ncs(:)];
cv = [pcv(:); ncv(:)];
mv = [pmv(:); nmv(:)];
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
            xc(xcc) = cc{xcc}(1, 1) + 0.0001 * cc{xcc}(1, 2);
        end
        [sortv, sorti] = sort(xc);
    case {'y'}
        xc = zeros(1, numel(cc));
        for xcc = 1:numel(xc)
            xc(xcc) = cc{xcc}(1, 2) + 0.0001 * cc{xcc}(1, 1);
        end
        [sortv, sorti] = sort(xc);
    case {'z'}
        xc = zeros(1, numel(cc));
        for xcc = 1:numel(xc)
            xc(xcc) = cc{xcc}(1, 3) + 0.0001 * cc{xcc}(1, 1);
        end
        [sortv, sorti] = sort(xc);
end
cc = cc(sorti);
cv = cv(sorti);
mv = mv(sorti);

% split local max
if ~isinf(opts.localmax)
    if ~isempty(opts.localmin)
        k = opts.localmin;
    else
        k = 2;
    end
    c = 1;
    while c <= numel(cc)
        if opts.localmax <= size(cc{c}, 1)
            [scs, scv, scl, scc] = splitclustercoords(cc{c}, abs(cv{c}), k);
            if numel(scs) > 1
                smv = zeros(numel(scs), 1);
                if cv{c}(1) > 0
                    for suc = 1:numel(scv)
                        smv(suc) = scv{suc}(1);
                    end
                else
                    for suc = 1:numel(scv)
                        scv{suc} = -scv{suc};
                        smv(suc) = scv{suc}(1);
                    end
                end
                st = [st(1:c), false(1, numel(scs)), st((c+1):end)];
                cc = [cc(1:c); scc(:); cc((c+1):end)];
                cv = [cv(1:c); scv(:); cv((c+1):end)];
                mv = [mv(1:c); smv(:); mv((c+1):end)];
                if ~opts.localmaxi
                    c = c + numel(scs);
                end
            end
        end
        c = c + 1;
    end
end
nc = numel(cc);

% create output structure
list = emptystruct({ ...
    'center', 'cog', 'coords', ...
    'localmax', 'max', 'mean', 'meanalt', ...
    'peak', 'peakalt', ...
    'rwcenter', 'rwcog', 'rwcoords', 'rwpeak', 'rwsize', ...
    'size', ...
    'talcenter', 'talcog', 'talcoords', 'talout', 'talpeak', ...
    'values'}, [nc, 1]);
table = cell(1, nc);
taltext = '';

% process results
if opts.tdclient && ...
    nc > 7
    pbn = opts.pbarrange(1);
    pbx = opts.pbarrange(2) - pbn + eps;
    pbd = 1 / 172800;
    pbl = now + pbd;
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
    crd = cc{c}(1, :);
    numv = size(cc{c}, 1);
    list(c).center = sum(cc{c}, 1) ./ numv;
    list(c).cog = sum(abs(cv{c}(:, [1, 1, 1])) .* cc{c}, 1) ./ sum(abs(cv{c}));
    list(c).coords = cc{c};
    if st(c)
        list(c).localmax = ' ';
    else
        list(c).localmax = 'L';
    end
    list(c).max = mv(c);
    list(c).mean = sum(cv{c}) / numv;
    if ~isempty(opts.altmaps)
        crdx = sub2ind(vs, cc{c}(:, 1), cc{c}(:, 2), cc{c}(:, 3));
        list(c).meanalt = meannoinfnan( ...
            opts.altmaps(crdx(:) * ones(1, nnsea) + ones(numel(crdx), 1) * nsea), 1, true);
    end
    list(c).peak = crd;
    list(c).peakalt = opts.altmaps(nsea + ...
        sub2ind(vs, crd(1), crd(2), crd(3)));
    if ~isempty(opts.mat)
        list(c).rwcoords = (opts.mat * [cc{c}, ones(numv, 1)]')';
        list(c).rwcoords(:, end) = [];
        list(c).rwcenter = sum(list(c).rwcoords, 1) ./ numv;
        list(c).rwcog = sum(abs(cv{c}(:, [1, 1, 1])) .* list(c).rwcoords, 1) ./ sum(abs(cv{c}));
        list(c).rwpeak = list(c).rwcoords(1, :);
        list(c).rwsize = numv * vxsf;
        if opts.icbm2tal
            list(c).talcoords = icbm2tal(list(c).rwcoords);
        elseif opts.mni2tal
            list(c).talcoords = mni2tal(list(c).rwcoords);
        else
            list(c).talcoords = list(c).rwcoords;
        end
        list(c).talcenter = sum(list(c).talcoords, 1) ./ numv;
        list(c).talcog = ...
            sum(abs(cv{c}(:, [1, 1, 1])) .* list(c).talcoords, 1) ./ sum(abs(cv{c}));
        list(c).talpeak = list(c).talcoords(1, :);
        switch (opts.lupcrd)
            case {'peak'}
                crd = round(list(c).talpeak);
            case {'center'}
                crd = round(list(c).talcenter);
            case {'cog'}
                crd = round(list(c).talcog);
        end
        if opts.tdclient
            talta = tdlabel(round(crd));
            if ~isempty(pbar) && ...
                now >= pbl
                pbl = now + pbd;
                pbar.Progress(pbn + pbx * c / nc, ...
                    sprintf('%sGetting label %d/%d...', pbarn, c + 1, nc));
                pbar.Visible = 'on';
            end
            taltext = talta{1};
        end
    else
        switch (opts.lupcrd)
            case {'peak'}
                crd = round(list(c).peak);
            case {'center'}
                crd = round(list(c).center);
            case {'cog'}
                crd = round(list(c).cog);
        end
    end
    list(c).size = numv;
    list(c).talout = taltext;
    list(c).values = cv{c};
    if ~isempty(altmapsf)
        if opts.altstat == 'm'
            if opts.localmsz || ...
                list(c).localmax == ' '
                table{c} = sprintf(['%4d %4d %4d %s| %5d  | %9.6f | %9.6f  ' altmapsf '|%s'], ...
                            crd(1), crd(2), crd(3), list(c).localmax, ...
                            list(c).size, list(c).max, list(c).mean, ...
                            list(c).meanalt(opts.altmapsp), taltext);
            else
                table{c} = sprintf(['%4d %4d %4d %s|        | %9.6f | %9.6f  ' altmapsf '|%s'], ...
                            crd(1), crd(2), crd(3), list(c).localmax, ...
                            list(c).max, list(c).mean, ...
                            list(c).meanalt(opts.altmapsp), taltext);
            end
        else
            if opts.localmsz || ...
                list(c).localmax == ' '
                table{c} = sprintf(['%4d %4d %4d %s| %5d  | %9.6f | %9.6f  ' altmapsf '|%s'], ...
                            crd(1), crd(2), crd(3), list(c).localmax, ...
                            list(c).size, list(c).max, list(c).mean, ...
                            list(c).peakalt(opts.altmapsp), taltext);
            else
                table{c} = sprintf(['%4d %4d %4d %s|        | %9.6f | %9.6f  ' altmapsf '|%s'], ...
                            crd(1), crd(2), crd(3), list(c).localmax, ...
                            list(c).max, list(c).mean, ...
                            list(c).peakalt(opts.altmapsp), taltext);
            end
        end
    else
        if opts.localmsz || ...
            list(c).localmax == ' '
            table{c} = sprintf('%4d %4d %4d %s| %5d  | %9.6f | %9.6f  |%s', ...
                        crd(1), crd(2), crd(3), list(c).localmax, ...
                        list(c).size, list(c).max, list(c).mean, taltext);
        else
            table{c} = sprintf('%4d %4d %4d %s|        | %9.6f | %9.6f  |%s', ...
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
    ['   x    y    z  |     k  |    max    |    mean    %s| tdclient\n' ...
     '--------------------------------------------------%s---------\n%s\n'], ...
     altmapsh, altmapsd, gluetostringc(table, char(10)));
