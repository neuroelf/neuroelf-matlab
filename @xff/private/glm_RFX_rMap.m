function map = glm_RFX_rMap(xo, c, r, opts)
% GLM::RFX_rMap  - calculate a second-level r contrast map
%
% FORMAT:       map = glm.RFX_rMap(c, r [, opts])
%
% Input fields:
%
%       c           NxC contrast vector
%       r           SxR regression data (or VOI object)
%       opts        structure with optional fields
%        .allrs     boolean flag, treat multiple Rs as one model (false)
%        .alphasim  call alphasim with data to assess CDT-based FPs (false)
%        .autobal   auto-balance contrast if pos and abs(neg) == 1 (true)
%        .bbox      bounding box for VMP if necessary (default: full MNI)
%        .bvcomp    BV-compatible map names (length restriction, true)
%        .cnames    1xC cell array with contrast names
%        .const     also create (t-) map of constant term
%        .estfwhm   estimate smoothness and store in Map.RunTimeVars (true)
%        .estsmap   store estimated smoothness in separate map (false)
%        .groups    Gx2 cell array with group names and subject IDs
%        .imeth     interpolation method if necessary (default: 'cubic')
%        .meanr     boolean flag, remove mean from map (added as cov)
%        .meanrmsk  mask to get mean from (object or XxYxZ logical)
%        .minnum    minumum number of subjects to compute (2 * sqrt(N))
%        .names     1xN cell array with map names
%        .rank      flag, rank-transform data before regression
%        .resvtc    create and open residual VTC in GUI (default: false)
%        .rnames    1xR cell array with regressor names
%        .robust    flag, use robust regression in addition to OLS
%        .robwmaps  create robust-weight summary maps (default: false)
%        .smk       smoothing kernel for data (prior to regression, 0)
%        .srf       SRF object (required for multi-study MTC-based stats)
%        .subsel    subject selection (otherwise all subjects)
%        .swmaps    weight map per subject and regression (default: false)
%        .thresh    1x2 threshold (lower, upper), as p-values!
%        .voiidx    index into VOI list (only used if r is a VOI object)
%
% Output fields:
%
%       map         MAP/VMP/SMP object with maps
%
% Using: alphasim, correlinvtstat, findfirst, fitrobustbisquare_img,
%        limitrangec, lsqueeze, newnatresvmp, ranktrans, resestsmooth,
%        robustt, sdist, smoothdata3, ztrans.

% Version:  v1.1
% Build:    16101016
% Date:     Oct-10 2016, 4:31 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2015, 2016, Jochen Weber
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
alphasim       = ne_methods.alphasim;
correlinvtstat = ne_methods.correlinvtstat;
findfirst      = ne_methods.findfirst;
limitrangec    = ne_methods.limitrangec;
ranktrans      = ne_methods.ranktrans;
sdist          = ne_methods.sdist;
ztrans         = ne_methods.ztrans;

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm') || ...
   ~isa(c, 'double') || isempty(c) || any(isnan(c(:)) | isinf(c(:))) || ...
   ((~isa(r, 'double') || isempty(r) || any(isinf(r(:)))) && ...
    (numel(r) ~= 1 || ~xffisobject(r, true, 'voi')))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
bcrtv = bc.RunTimeVars;
glmfile = xo.F;
glmid = xo.L;
if isempty(glmfile)
    glmfile = glmid;
end
if bc.ProjectTypeRFX == 0 && bc.SeparatePredictors ~= 2
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
isrfx = (bc.ProjectTypeRFX > 0);
ffxspred = glm_SubjectPredictors(xo);
ffxsubs = glm_Subjects(xo);
if isrfx
    numsubs = numel(bc.GLMData.Subject);
    numspred = size(bc.GLMData.Subject(1).BetaMaps, ndims(bc.GLMData.Subject(1).BetaMaps));
else
    ffxpred = bc.Predictor;
    ffxpred = {ffxpred(:).Name2};
    ffxpred = ffxpred(:);
    numsubs = numel(ffxsubs);
    numspred = numel(ffxspred) + 1;
end
if ~any(bc.ProjectType == [1, 2])
    error('neuroelf:xff:unsupported', ...
        'RFX correlation maps of FMRs are not yet supported.');
end
if numsubs < 3
    error('neuroelf:xff:badArgument', 'Invalid RFX GLM object.');
end
if bc.ProjectType == 1
    if isrfx
        msz = size(bc.GLMData.RFXGlobalMap);
    else
        msz = size(bc.GLMData.MCorrSS);
    end
else
    if isrfx
        msz = numel(bc.GLMData.RFXGlobalMap);
    else
        msz = numel(bc.GLMData.MCorrSS);
    end
end
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'allrs') || ~islogical(opts.allrs) || numel(opts.allrs) ~= 1
    opts.allrs = false;
end
if ~isfield(opts, 'alphasim') || ~islogical(opts.alphasim) || numel(opts.alphasim) ~= 1
    opts.alphasim = false;
else
    athrs = (0.001 .* [50, 20, 10, 5, 2, 1, .5, .1, .05, .01, .005, .001])';
end
if ~isfield(opts, 'autobal') || ~islogical(opts.autobal) || numel(opts.autobal) ~= 1
    opts.autobal = true;
end
if ~isfield(opts, 'bbox') || ~isa(opts.bbox, 'double') || ...
   ~isequal(size(opts.bbox), [2, 3]) || ...
    any(isinf(opts.bbox(:)) | isnan(opts.bbox(:)) | ...
        opts.bbox(:) < 0 | opts.bbox(:) > 256 | opts.bbox(:) ~= fix(opts.bbox(:))) || ...
    any(diff(opts.bbox) < 0)
    opts.bbox = [44, 38, 44; 242, 194, 212];
end
if ~isfield(opts, 'bvcomp') || ~islogical(opts.bvcomp) || numel(opts.bvcomp) ~= 1
    opts.bvcomp = true;
end
if ~isfield(opts, 'cnames') || ~iscell(opts.cnames)
    opts.cnames = {};
end
if ~isfield(opts, 'const') || ~islogical(opts.const) || numel(opts.const) ~= 1
    opts.const = false;
end
if ~isfield(opts, 'estfwhm') || ~islogical(opts.estfwhm) || numel(opts.estfwhm) ~= 1
    opts.estfwhm = true;
end
if ~isfield(opts, 'estsmap') || ~islogical(opts.estsmap) || numel(opts.estsmap) ~= 1
    opts.estsmap = false;
elseif opts.estsmap
    opts.estfwhm = true;
end
if ~isfield(opts, 'groups') || ~iscell(opts.groups) || size(opts.groups, 2) ~= 2
    ngrp = 0;
    if ~isfield(opts, 'subsel') || ~isa(opts.subsel, 'double') || ...
        isempty(opts.subsel) || any(isinf(opts.subsel(:)) | isnan(opts.subsel(:)))
        ga = ones(numsubs, 1);
    else
        opts.subsel = opts.subsel(:)';
        opts.subsel(opts.subsel < 1 | opts.subsel > numsubs) = [];
        ga = zeros(numsubs, 1);
        ga(unique(round(opts.subsel))) = 1;
    end
    opts.groups = [];
    gamx = 1;
else
    ga = zeros(numsubs, 1);
    for gc = 1:size(opts.groups, 1)
        if ~ischar(opts.groups{gc, 1}) || isempty(opts.groups{gc, 1}) || ...
            isempty(opts.groups{gc, 2}) || ...
            any(isinf(opts.groups{gc, 2}(:)) | isnan(opts.groups{gc, 2}(:))) || ...
            any(opts.groups{gc, 2}(:) < 1 | opts.groups{gc, 2}(:) > numsubs) || ...
            any(ga(round(opts.groups{gc, 2}(:))) > 0)
            error('neuroelf:xff:badArgument', 'Invalid group assignment.');
        end
        opts.groups{gc, 2} = unique(round(opts.groups{gc, 2}(:)));
        ga(opts.groups{gc, 2}) = gc;
    end
    ngrp = size(opts.groups, 1);
    nval = sum(ga > 0) - ngrp;
    if nval < 3
        error('neuroelf:xff:badArgument', 'Too few subjects for grouping.');
    end
    gas = ga(ga > 0);
    gag = tril(ones(ngrp));
    gag(1:(ngrp + 1):end) = 0;
    gamx = sum(gag(:));
    gag(gag > 0) = 1:gamx;
end
gax = find(ga > 0);
if ~isfield(opts, 'imeth') || ~ischar(opts.imeth) || isempty(opts.imeth) || ...
    isempty(regexpi(opts.imeth(:)', '^(cubic|lanczos\d|linear|nearest)$'))
    opts.imeth = 'cubic';
else
    opts.imeth = lower(opts.imeth(:)');
end
if ~isfield(opts, 'meanr') || ~islogical(opts.meanr) || numel(opts.meanr) ~= 1
    opts.meanr = false;
end
if isfield(opts, 'meanrmsk') && numel(opts.meanrmsk) == 1 && ...
    xffisobject(opts.meanrmsk, true, 'msk')
    mbc = opts.meanrmsk.C;
    if numel(mbc.Mask) == prod(msz)
        opts.meanrmsk = ne_methods.lsqueeze(mbc.Mask > 0);
    else
        opts.meanrmsk = [];
    end
elseif isfield(opts, 'meanrmsk') && islogical(opts.meanrmsk) && numel(opts.meanrmsk) == prod(msz)
    opts.meanrmsk = ne_methods.lsqueeze(opts.meanrmsk);
else
    opts.meanrmsk = [];
end
if ~isfield(opts, 'minnum') || ~isa(opts.minnum, 'double') || numel(opts.minnum) ~= 1 || ...
    isinf(opts.minnum) || isnan(opts.minnum) || opts.minnum < 0
    opts.minnum = 0;
end
if ~isfield(opts, 'names') || ~iscell(opts.names) || isempty(opts.names)
    opts.names = {};
end
if ~isfield(opts, 'rank') || ~islogical(opts.rank) || numel(opts.rank) ~= 1 || ...
   ~opts.rank
    opts.rank = false;
    ranktxt = '';
else
    ranktxt = 'Rank-';
end
if ~isfield(opts, 'resvtc') || ~islogical(opts.resvtc) || numel(opts.resvtc) ~= 1 || bc.ProjectType == 2
    opts.resvtc = false;
end
resvtc = [];
if ~isfield(opts, 'rnames') || ~iscell(opts.rnames)
    opts.rnames = {};
end
if ~isfield(opts, 'robust') || ~islogical(opts.robust) || numel(opts.robust) ~= 1
    opts.robust = false;
end
if ~isfield(opts, 'robwmaps') || ~islogical(opts.robwmaps) || numel(opts.robwmaps) ~= 1
    opts.robwmaps = false;
end
if ~isfield(opts, 'smk') || numel(opts.smk) ~= 1 || ~isa(opts.smk, 'double') || ...
    isinf(opts.smk) || isnan(opts.smk) || opts.smk <= (0.5 * bc.Resolution)
    opts.smk = 0;
else
    opts.smk = min(opts.smk, 8 * bc.Resolution);
end
if ~isfield(opts, 'srf') || numel(opts.srf) ~= 1 || ~xffisobject(opts.srf, true, 'srf')
    opts.srf = [];
end
if ~isfield(opts, 'subsel') || ~isa(opts.subsel, 'double') || isempty(opts.subsel) || ...
    any(isinf(opts.subsel(:)) | isnan(opts.subsel(:))) || ...
    numel(unique(round(opts.subsel(:)))) ~= numel(opts.subsel) || ...
    any(opts.subsel(:) < 1 | opts.subsel(:) > numsubs)
    opts.subsel = 1:numsubs;
else
    opts.subsel = round(opts.subsel(:)');
end
if ~isfield(opts, 'swmaps') || ~islogical(opts.swmaps) || numel(opts.swmaps) ~= 1
    opts.swmaps = false;
end
if ~isfield(opts, 'voiidx') || ~isa(opts.voiidx, 'double') || numel(opts.voiidx) ~= 1 || ...
    isinf(opts.voiidx) || isnan(opts.voiidx) || opts.voiidx < 1
    opts.voiidx = 1;
else
    opts.voiidx = floor(opts.voiidx);
end
if numel(r) == 1
    try
        rbc = r.C;
        rbc = numel(rbc.VOI);
        if size(c, 2) ~= numspred && size(c, 2) ~= (numspred - 1)
            ct = c';
        else
            ct = c;
        end
        r = glm_VOIBetas(xo, r, struct('c', ct, 'vl', min(opts.voiidx, rbc)));
    catch xfferror
        rethrow(xfferror);
    end
end
if any(size(r) == numsubs) && numel(opts.subsel) ~= numsubs
    rsr = repmat({':'}, 1, ndims(r));
    rsr{findfirst(size(r) == numsubs)} = opts.subsel;
    r = r(rsr{:});
end
if ~any(size(r) == numel(opts.subsel))
    error('neuroelf:xff:badArgument', ...
        'Correlation regressors must match in size with number of subjects.');
end
rsdim = findfirst(size(r) == numel(opts.subsel));
nanr = isnan(r);
for dc = 1:ndims(r)
    if dc == rsdim
        continue;
    end
    nanr = any(nanr, dc);
end
if any(nanr)
    rsr = repmat({':'}, 1, ndims(r));
    rsr{rsdim} = find(~nanr);
    r = r(rsr{:});
    opts.subsel = opts.subsel(~nanr);
end
if ~isfield(opts, 'thresh') || ~isa(opts.thresh, 'double') || numel(opts.thresh) ~= 2 || ...
    any(isinf(opts.thresh) | isnan(opts.thresh) | opts.thresh <= 0 | opts.thresh >= 0.5)
    opts.thresh = [0.005, 0.0001];
else
    opts.thresh = -sort(-opts.thresh(:)');
end
subsel = opts.subsel;
numsubs = numel(subsel);
if opts.minnum == 0
    opts.minnum = ceil(2 * sqrt(numsubs));
elseif opts.minnum < 1
    opts.minnum = ceil(opts.minnum * numsubs);
end
opts.minnum = min(numsubs, opts.minnum);
if opts.meanr
    mvm = zeros(numsubs, 1);
end
thresh = opts.thresh;
if size(r, 2) == numel(opts.subsel) && size(r, 1) ~= numel(opts.subsel)
    r = r';
end
if size(c, 1) == 1 && (size(c, 2) == (numspred - 1) || size(c, 2) == numspred)
    c = c';
end
if size(c, 1) == (numspred - 1)
    c(end+1,:) = 0;
end
if size(c, 1) ~= numspred
    error('neuroelf:xff:badArgument', 'Contrast vector must span all conditions.');
end
nummaps = size(c, 2);
numrs = size(r, 2);
if opts.allrs
    nval = numsubs - (1 + numrs);
else
    nval = numsubs - 2;
end
if numel(opts.cnames) ~= nummaps
    sprednames = glm_SubjectPredictors(xo);
    opts.cnames = cell(1, nummaps);
    for cc = 1:nummaps
        con = find(c(:, cc) > 0);
        coff = find(c(:, cc) < 0);
        conn = cell(1, numel(con));
        coffn = cell(1, numel(coff));
        for occ = 1:numel(con)
            if c(con(occ), cc) ~= 1
                conn{occ} = sprintf('%g * %s', sprednames{con(occ)});
            else
                conn{occ} = sprednames{con(occ)};
            end
        end
        for occ = 1:numel(coff)
            if c(coff(occ), cc) ~= -1
                coffn{occ} = sprintf('%g * %s', sprednames{coff(occ)});
            else
                coffn{occ} = sprednames{coff(occ)};
            end
        end
        if ~isempty(con)
            if numel(con) > 1
                connc = sprintf('%s + ', conn{:});
                connc = sprintf('(%s)', connc(1:end-3));
            else
                connc = conn{1};
            end
        else
            connc = 'Baseline';
        end
        if ~isempty(coff)
            if numel(coff) > 1
                coffnc = sprintf('%s + ', coffn{:});
                coffnc = sprintf('(%s)', coffnc(1:end-3));
            else
                coffnc = coffn{1};
            end
        else
            coffnc = '';
        end
        if ~isempty(coffnc)
            opts.cnames{cc} = sprintf('%s > %s', connc, coffnc);
        else
            opts.cnames{cc} = connc;
        end
    end
end
if numel(opts.rnames) ~= numrs
    opts.rnames = cell(1, numrs);
    for cc = 1:numrs
        opts.rnames = sprintf('reg%02d', cc);
    end
end
nummapst = nummaps * (numrs + opts.const + opts.estsmap) * (1 + opts.meanr) * (1 + opts.robust) * ...
    (1 + double(opts.robust) .* (2 * double(opts.robwmaps) + numsubs * double(opts.swmaps)));
if bc.ProjectType == 1
    subsa = {':', ':', ':', []};
    subsr = {':', ':', ':', []};
else
    subsa = {':', []};
    subsr = {':', []};
end
rpma = [msz, 1];
if numel(opts.names) ~= (nummaps * numrs)
    opts.names = cell(1, nummaps * numrs);
end
for cc = 1:numel(opts.names)
    if ~ischar(opts.names{cc})
        ccr = 1 + mod(cc - 1, numrs);
        ccc = 1 + round((cc - ccr) / numrs);
        opts.names{cc} = sprintf('%sCorr: (%s, %s)', ranktxt, opts.cnames{ccc}, opts.rnames{ccr});
    end
    if opts.bvcomp && numel(opts.names{cc}) > 96
        opts.names{cc} = [opts.names{cc}(1:46) ' ... ' opts.names{cc}(end-45:end)];
    else
        opts.names{cc} = opts.names{cc}(:)';
    end
end
if isrfx
    szmap = size(bc.GLMData.RFXGlobalMap);
else
    szmap = size(bc.GLMData.MCorrSS);
end
szmapo = [szmap, 1];

% create map container
copymaps = true;
if bc.ProjectType == 0
    error('neuroelf:xff:notYetSupported', 'FMR->MAP generation not yet supported.');
end

% VTC-based
if bc.ProjectType == 1
    mdim = 4;
    mones = [1, 1, 1];
    map = ne_methods.newnatresvmp();
    mapc = map.C;
    mapc.Resolution = bc.Resolution;
    if ~isfield(bcrtv, 'SubjectSPMsn') || ~isstruct(bcrtv.SubjectSPMsn) || ...
        isempty(fieldnames(bcrtv.SubjectSPMsn))
        if isrfx
            szmap = size(bc.GLMData.RFXGlobalMap);
        else
            szmap = size(bc.GLMData.MCorrSS);
        end
        mapc.XStart = bc.XStart;
        mapc.XEnd = bc.XEnd;
        mapc.YStart = bc.YStart;
        mapc.YEnd = bc.YEnd;
        mapc.ZStart = bc.ZStart;
        mapc.ZEnd = bc.ZEnd;
    else
        opts.bbox(2, :) = opts.bbox(1, :) + mapc.Resolution .* ...
            ceil((opts.bbox(2, :) - opts.bbox(1, :)) ./ mapc.Resolution - 0.01);
        mapc.XStart = opts.bbox(1, 1);
        mapc.XEnd = opts.bbox(2, 1);
        mapc.YStart = opts.bbox(1, 2);
        mapc.YEnd = opts.bbox(2, 2);
        mapc.ZStart = opts.bbox(1, 3);
        mapc.ZEnd = opts.bbox(2, 3);
        szmap = round((opts.bbox(2, :) - opts.bbox(1, :)) ./ mapc.Resolution);
        copymaps = false;
        sbbox = struct('BBox', opts.bbox, 'ResXYZ', mapc.Resolution);
        rpma = round(diff(opts.bbox) ./ mapc.Resolution);
        msz = rpma;
    end
    mapc.RunTimeVars.AutoSave = true;
    mapc.RunTimeVars.TrfPlus = bcrtv.TrfPlus;
    mapf = 'VMPData';

    % residuals VTC container
    if opts.resvtc
        resvtc = xff('new:vtc');
        resvtcc = resvtc.C;
        resvtcc.NameOfSourceFMR = glmfile;
        resvtcc.DataType = 2;
        resvtcc.Resolution = mapc.Resolution;
        resvtcc.XStart = mapc.XStart;
        resvtcc.XEnd = mapc.XEnd;
        resvtcc.YStart = mapc.YStart;
        resvtcc.YEnd = mapc.YEnd;
        resvtcc.ZStart = mapc.ZStart;
        resvtcc.ZEnd = mapc.ZEnd;
        resvtcc.TR = 29554; % hex('sr')
    end
    subsur = {':', ':', ':', []};

% MTC/SMP
elseif bc.ProjectType == 2
    mdim = 2;
    mones = 1;
    map = xff('new:smp');
    mapc = map.C;
    mapc.NrOfVertices = bc.NrOfVertices;
    mapf = 'SMPData';
    szmap = bc.NrOfVertices;
    if isrfx
        szmapo = size(bc.GLMData.RFXGlobalMap);
    else
        szmapo = size(bc.GLMData.MCorrSS);
    end
    rpma = szmapo;
    if opts.estfwhm && (isempty(opts.srf) || opts.srf.C.NrOfVertices ~= szmap)
        opts.estfwhm = false;
        opts.estsmap = false;
    end
    subsur = {':', []};
end

% set some common fields
mapc.NrOfMaps = nummapst;
mapc.Map.Type = 2;
mapc.Map.LowerThreshold = correlinvtstat(-sdist('tinv', thresh(1), nval), numsubs);
mapc.Map.UpperThreshold = correlinvtstat(-sdist('tinv', thresh(2), nval), numsubs);
mapc.Map.DF1 = nval;
mapc.Map.DF2 = 0;
mapc.Map.NrOfFDRThresholds = 0;
mapc.Map.FDRThresholds = zeros(0, 3);
mapc.Map.(mapf) = single(zeros(rpma));

% replicate
mapc.Map = mapc.Map(1, ones(1, nummapst));

% rank-transform
if opts.rank
    r = ranktrans(r, 1);
end

% what models?
if opts.allrs
    micc = 1:numrs:(nummaps * numrs);
    numrsi = numrs;
else
    micc = 1:(nummaps * numrs);
    numrsi = 1;
end

% computation
conmaps = zeros([msz, numsubs]);
tmc = 1;
lcc = 0;
for icc = 1:numel(micc)

    % which contrast and regressor
    cr = 1 + mod(micc(icc) - 1, numrs);
    cc = 1 + round((micc(icc) - cr) / numrs);

    % zero out and fill conmaps if necessary
    if lcc ~= cc
        conmaps(:) = 0;
        keepsubs = true(1, numsubs);

        % auto-balancing
        if abs(sum(c(:, cc))) < sqrt(eps) && abs(sum(c(c(:, cc) > 0, cc)) - 1) < sqrt(eps) && ...
            abs(1 + sum(c(c(:, cc) < 0, cc))) < sqrt(eps) && opts.autobal
            autobal = true;
            possum = zeros(numsubs, 1);
            negsum = zeros(numsubs, 1);
            nconmaps = conmaps;
        else
            autobal = false;
        end

        % fill contrast maps
        for pc = 1:numspred
            if c(pc, cc) ~= 0
                subsr{end} = pc;
                for sc = 1:numsubs
                    if ~keepsubs(sc)
                        continue;
                    end
                    subsa{end} = sc;
                    subsur{end} = sc;
                    if isrfx
                        if copymaps
                            cpmap = bc.GLMData.Subject(subsel(sc)).BetaMaps(subsr{:}); 
                        else
                            cpmap = aft_SampleBVBox(xo, sbbox, ...
                                (sc - 1) * numspred + pc, opts.imeth);
                        end
                        if autobal
                            if all(isinf(cpmap(:)) | isnan(cpmap(:)) | cpmap(:) == 0)
                                continue;
                            end
                            if c(pc, cc) > 0
                                possum(sc) = possum(sc) + c(pc, cc);
                                conmaps(subsur{:}) = conmaps(subsur{:}) + c(pc, cc) .* cpmap;
                            else
                                negsum(sc) = negsum(sc) + c(pc, cc);
                                nconmaps(subsur{:}) = nconmaps(subsur{:}) + c(pc, cc) .* cpmap;
                            end
                        else
                            conmaps(subsur{:}) = conmaps(subsur{:}) + c(pc, cc) .* cpmap;
                        end
                    else
                        keepsubi = findfirst(~cellfun('isempty', regexpi(ffxpred, ...
                            sprintf('^subject\\s+%s:\\s*%s', ...
                            ffxsubs{gax(sc)}, ffxspred{pc}))));
                        if isempty(keepsubi) && autobal
                            continue;
                        elseif ~isempty(keepsubi)
                            subsr{end} = keepsubi;
                            cpmap = bc.GLMData.BetaMaps(subsr{:});
                            if autobal
                                if all(isinf(cpmap(:)) | isnan(cpmap(:)) | cpmap(:) == 0)
                                    continue;
                                end
                                if c(pc, cc) > 0
                                    possum(sc) = possum(sc) + c(pc, cc);
                                    conmaps(subsur{:}) = conmaps(subsur{:}) + ...
                                        c(pc, cc) .* cpmap;
                                else
                                    negsum(sc) = negsum(sc) + c(pc, cc);
                                    nconmaps(subsur{:}) = nconmaps(subsur{:}) + ...
                                        c(pc, cc) .* cpmap;
                                end
                            else
                                if all(isinf(cpmap(:)) | isnan(cpmap(:)) | cpmap(:) == 0)
                                    conmaps(subsur{:}) = 0;
                                    keepsubs(sc) = false;
                                end
                                conmaps(subsur{:}) = conmaps(subsur{:}) + c(pc, cc) .* cpmap;
                            end
                        else
                            conmaps(subsur{:}) = 0;
                            keepsubs(sc) = false;
                        end
                    end
                end
            end
        end

        % auto-balancing
        if autobal
            if ~all(abs(possum - 1) < sqrt(eps))
                conmaps = reshape(reshape(conmaps, prod(szmap), numsubs) * diag(1 ./ possum), [szmap, numsubs]);
            end
            if ~all(abs(1 + negsum) < sqrt(eps))
                nconmaps = reshape(reshape(nconmaps, prod(szmap), numsubs) * diag(-1 ./ negsum), [szmap, numsubs]);
            end
            conmaps = conmaps + nconmaps;
            badsum = (possum == 0 | negsum == 0);
            if any(badsum)
                subsur{end} = badsum;
                conmaps(subsur{:}) = NaN;
            end
        end

        % check maps
        for sc = 1:numsubs
            if ~keepsubs(sc)
                continue;
            end
            subsa{end} = sc;
            subsas = struct('type', '()', 'subs', {subsa});
            cmtest = ne_methods.lsqueeze(subsref(conmaps, subsas));
            if all(isnan(cmtest) | isinf(cmtest) | cmtest == 0)
                conmaps = subsasgn(conmaps, subsas, 0);
                keepsubs(sc) = false;
            end
        end
        keepsubsn = sum(keepsubs);
        if opts.allrs
            nval = keepsubsn - (1 + numrs);
        else
            nval = keepsubsn - 2;
        end

        % minumum criterion not met? then set to 0 for all subjects
        conmapsa = (~isinf(conmaps));
        conmapsa = conmapsa & (~isnan(conmaps));
        conmapsa = conmapsa & (conmaps ~= 0);
        conmapsa = sum(conmapsa, numel(subsr));
        conmaps(repmat(conmapsa < opts.minnum, [ones(1, numel(msz)), numsubs])) = 0;

        % rank transform?
        if opts.rank
            subsa{end} = keepsubs;
            subsas = struct('type', '()', 'subs', {subsa});
            conmaps = subsasgn(conmaps, subsas, ...
                ranktrans(subsref(conmaps, subsas), numel(subsr), ...
                struct('meancenter', true, 'nozero', true)));
        else
            conmaps(isinf(conmaps) | isnan(conmaps)) = 0;
        end

        % keep track!
        lcc = cc;
    end

    % generate design matrix/ces
    if opts.allrs
        X = [ztrans(r(keepsubs, :)), ones(keepsubsn, 1)];
    else
        X = [ztrans(r(keepsubs, cr)), ones(keepsubsn, 1)];
    end
    iXX = pinv(X' * X);
    if opts.meanr
        if isempty(opts.meanrmsk)
            meanrmsk = all(bc.GLMData.Subject(findfirst(keepsubs)).BetaMaps ~= 0, ...
                ndims(bc.GLMData.Subject(findfirst(keepsubs)).BetaMaps));
            for sc = 1:numsubs
                if keepsubs(sc)
                    subsa{end} = sc;
                    subsas = struct('type', '()', 'subs', {subsa});
                    meanrmsk = (meanrmsk & subsref(conmaps, subsas) ~= 0);
                end
            end
            meanrmsk = ne_methods.lsqueeze(meanrmsk);
        else
            meanrmsk = ne_methods.lsqueeze(opts.meanrmsk);
        end
        for sc = 1:keepsubsn
            subsa{end} = sc;
            subsas = struct('type', '()', 'subs', {subsa});
            mv = subsref(conmaps, subsas);
            mvm(sc) = ne_methods.meannoinfnan(mv(meanrmsk));
        end
        if opts.rank
            Xm = [X, ztrans(ranktrans(mvm(keepsubs), 1))];
        else
            Xm = [X, ztrans(mvm(keepsubs))];
        end
        iXXm = pinv(Xm' * Xm);
    end

    % smooth data
    if opts.smk > 0 && ndims(conmaps) == 4
        resim = (isinf(conmaps) | isnan(conmaps) | conmaps == 0);
        conmaps = ne_methods.smoothdata3(conmaps, opts.smk / mapc.Resolution);
        conmaps(resim) = 0;
    end

    % set additional data
    artv = struct( ...
        'SourceGLM',   glmfile, ...
        'SourceGLMID', glmid, ...
        'Contrast',    c(:, lcc), ...
        'FWHMResEst',  [], ...
        'FWHMResImg',  [], ...
        'GlobSigMap',  [], ...
        'MeanRem',     opts.meanr, ...
        'Regressors',  r, ...
        'RFXGLM',      true, ...
        'Robust',      opts.robust, ...
        'SubPreds',    {ffxspred}, ...
        'SubSel',      subsel(:));

    % OLS computations first
    subsa{end} = keepsubs;
    subsas = struct('type', '()', 'subs', {subsa});
    conmaps = subsref(conmaps, subsas);
    betas = iXX * X' * reshape(conmaps, prod(msz), keepsubsn)';
    resim = conmaps - reshape((X * betas)', [msz, keepsubsn]);
    stder = sqrt(sum(resim .^ 2, ndims(conmaps))) .* sqrt(1 / nval);
    if opts.estfwhm && mdim == 4
        [artv.FWHMResEst, artv.FWHMResImg] = ne_methods.resestsmooth(resim, bc.Resolution);
    elseif opts.estfwhm && mdim == 2
        [artv.FWHMResEst, artv.FWHMResImg] = ne_methods.resestsmoothsrf(resim, opts.srf);
    end

    % store and visualize residuals as VTC
    if opts.resvtc
        resvtcc.VTCData = permute(single(resim), [4, 1, 2, 3]);
        resvtcc.NrOfVolumes = size(resvtcc.VTCData, 1);
        resvtcc.RunTimeVars.Subjects = ffxsubs(gax(keepsubs));
        resvtc.C = resvtcc;
        cresvtc = aft_CopyObject(resvtc);
        cresvtc.H.InstCorrTC = single(ztrans(resvtcc.VTCData));
        cresvtc.H.InstCorrTCFilt = zeros(size(resvtcc.VTCData, 1), 0);
        neuroelf_gui('openfile', cresvtc, true);
    end

    % first maps
    for irc = 1:numrsi
        tmap = betas(irc, :)' ./ (stder(:) .* sqrt(iXX(irc, irc)));
        tmap(isinf(tmap) | isnan(tmap)) = 0;
        mapc.Map(tmc).(mapf) = reshape(single(correlinvtstat(tmap, keepsubsn)), [msz, 1]);
        mapc.Map(tmc).Name = opts.names{icc+irc-1};
        if opts.allrs
            artv.Regressors = r(:, irc);
        else
            artv.Regressors = r(:, cr);
        end
        if opts.alphasim && mdim == 4
            regmodsc = zeros(1, size(X, 2));
            regmodsc(irc) = 1;
            [athr, anull, amap, amax] = alphasim(msz, struct('regmaps', conmaps, ...
                'regmodel', X, 'regmodsc', regmodsc, 'thr', athrs));
            kthrs = zeros(size(athrs));
            for atc = 1:numel(athrs)
                kthrs(atc) = findfirst(athr{atc}(:, end) < 0.05);
            end
            artv.AlphaSim12 = {0.05, [athrs, kthrs]};
            artv.AlphaSim12zMaps = amap;
            artv.AlphaSimMaxDist = amax;
        end
        mapc.Map(tmc).DF1 = nval;
        mapc.Map(tmc).RunTimeVars = artv;
        tmc = tmc + 1;
    end
    if isfield(artv, 'AlphaSim12')
        artv = rmfield(artv, 'AlphaSim12');
    end
    if isfield(artv, 'AlphaSim12zMaps')
        artv = rmfield(artv, 'AlphaSim12zMaps');
    end
    if isfield(artv, 'AlphaSimMaxDist')
        artv = rmfield(artv, 'AlphaSimMaxDist');
    end
    if opts.const
        tmap = betas(numrsi+1, :)' ./ (stder(:) .* sqrt(iXX(numrsi+1, numrsi+1)));
        tmap(isinf(tmap) | isnan(tmap)) = 0;
        mapc.Map(tmc).(mapf) = reshape(single(tmap), [msz, 1]);
        mapc.Map(tmc).Type = 1;
        mapc.Map(tmc).Name = sprintf('%s (intercept-t)', opts.names{icc});
        mapc.Map(tmc).LowerThreshold = -sdist('tinv', thresh(1), nval);
        mapc.Map(tmc).UpperThreshold = -sdist('tinv', thresh(2), nval);
        artv.Regressors = [];
        mapc.Map(tmc).DF1 = nval;
        mapc.Map(tmc).RunTimeVars = artv;
        tmc = tmc + 1;
    end

    % add smoothness estimate map
    if opts.estsmap
        mapc.Map(tmc).Name = sprintf('%s (smoothness estimate)', opts.names{icc});
        mapc.Map(tmc).Type = 83;
        mapc.Map(tmc).(mapf) = single(artv.FWHMResImg);
        if mdim == 4
            mapc.Map(tmc).LowerThreshold = 2 * bc.Resolution + 0.5 * opts.smk;
            mapc.Map(tmc).UpperThreshold = 6 * bc.Resolution + 2.5 * opts.smk;
        else
            mapc.Map(tmc).LowerThreshold = 5;
            mapc.Map(tmc).UpperThreshold = 20;
        end
        mapc.Map(tmc).ShowPositiveNegativeFlag = 1;
        tmc = tmc + 1;
    end

    % with mean removed ?
    if opts.meanr
        betas = iXXm * Xm' * reshape(conmaps, prod(msz), keepsubsn)';
        resim = conmaps - reshape((Xm * betas)', [msz, keepsubsn]);
        stder = sqrt(sum(resim .^ 2, ndims(conmaps))) .* sqrt(1 / (nval - 1));
        if opts.estfwhm && mdim == 4
            [artv.FWHMResEst, artv.FWHMResImg] = ne_methods.resestsmooth(resim, bc.Resolution);
        elseif opts.estfwhm && mdim == 2
            [artv.FWHMResEst, artv.FWHMResImg] = ne_methods.resestsmoothsrf(resim, opts.srf);
        end
        artv.GlobSigMap = single(reshape(meanrmsk, rpma));

        % first maps
        for irc = 1:numrsi
            tmap = betas(irc, :)' ./ (stder(:) .* sqrt(iXXm(irc, irc)));
            tmap(isinf(tmap) | isnan(tmap)) = 0;
            mapc.Map(tmc).(mapf) = reshape(single(correlinvtstat(tmap, keepsubsn)), [msz, 1]);
            mapc.Map(tmc).Name = sprintf('%s (mean-rem)', opts.names{icc+irc-1});
            mapc.Map(tmc).DF1 = nval - 1;
            mapc.Map(tmc).LowerThreshold = ...
                correlinvtstat(-sdist('tinv', thresh(1), nval - 1), nval + 1);
            mapc.Map(tmc).UpperThreshold = ...
                correlinvtstat(-sdist('tinv', thresh(2), nval - 1), nval + 1);
            if opts.allrs
                artv.Regressors = r(:, irc);
            else
                artv.Regressors = r(:, cr);
            end
            if opts.alphasim && mdim == 4
                regmodsc = zeros(1, size(Xm, 2));
                regmodsc(irc) = 1;
                [athr, anull, amap, amax] = alphasim(msz, struct('regmaps', conmaps, ...
                    'regmodel', Xm, 'regmodsc', regmodsc, 'thr', athrs));
                kthrs = zeros(size(athrs));
                for atc = 1:numel(athrs)
                    kthrs(atc) = findfirst(athr{atc}(:, end) < 0.05);
                end
                artv.AlphaSim12 = {0.05, [athrs, kthrs]};
                artv.AlphaSim12zMaps = amap;
                artv.AlphaSimMaxDist = amax;
            end
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        if isfield(artv, 'AlphaSim12')
            artv = rmfield(artv, 'AlphaSim12');
        end
        if isfield(artv, 'AlphaSim12zMaps')
            artv = rmfield(artv, 'AlphaSim12zMaps');
        end
        if isfield(artv, 'AlphaSimMaxDist')
            artv = rmfield(artv, 'AlphaSimMaxDist');
        end
        if opts.const
            tmap = betas(numrsi+1, :)' ./ (stder(:) .* sqrt(iXXm(numrsi+1, numrsi+1)));
            tmap(isinf(tmap) | isnan(tmap)) = 0;
            mapc.Map(tmc).(mapf) = reshape(single(tmap), [msz, 1]);
            mapc.Map(tmc).Type = 1;
            mapc.Map(tmc).DF1 = nval - 1;
            mapc.Map(tmc).Name = sprintf('%s (intercept-t, mean-rem)', opts.names{icc});
            mapc.Map(tmc).LowerThreshold = -sdist('tinv', thresh(1), nval - 1);
            mapc.Map(tmc).UpperThreshold = -sdist('tinv', thresh(2), nval - 1);
            artv.Regressors = [];
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end

        % add smoothness estimate map
        if opts.estsmap
            mapc.Map(tmc).Name = sprintf('%s (smoothness estimate, mean-rem)', opts.names{icc});
            mapc.Map(tmc).Type = 83;
            mapc.Map(tmc).(mapf) = single(artv.FWHMResImg);
            if mdim == 4
                mapc.Map(tmc).LowerThreshold = 2 * bc.Resolution + 0.5 * opts.smk;
                mapc.Map(tmc).UpperThreshold = 6 * bc.Resolution + 2.5 * opts.smk;
            else
                mapc.Map(tmc).LowerThreshold = 5;
                mapc.Map(tmc).UpperThreshold = 20;
            end
            mapc.Map(tmc).ShowPositiveNegativeFlag = 1;
            tmc = tmc + 1;
        end
        artv.GlobSigMap = [];
    end

    % compute robust stats
    if opts.robust

        % perform fit
        [b, w] = ne_methods.fitrobustbisquare_img(X, conmaps);
        if opts.estfwhm
            ptc = zeros(size(conmaps));
            for bmc = 1:size(X, 2)
                if mdim == 4
                    ptc = ptc + repmat(b(:, :, :, bmc), [mones, size(X, 1)]) .* ...
                        repmat(reshape(X(:, bmc), [mones, size(X, 1)]), szmapo);
                else
                    ptc = ptc + repmat(b(:, bmc), [mones, size(X, 1)]) .* ...
                        repmat(reshape(X(:, bmc), [mones, size(X, 1)]), szmapo);
                end
            end
            ptc = w .* ptc + (1 - w) .* conmaps;
            if mdim == 4
                [artv.FWHMResEst, artv.FWHMResImg] = ...
                    ne_methods.resestsmooth(conmaps - ptc, bc.Resolution);
            else
                [artv.FWHMResEst, artv.FWHMResImg] = ...
                    ne_methods.resestsmoothsrf(conmaps - ptc, opts.srf);
            end
        end
        rsicv = zeros(1, size(X, 2));
        for irc = 1:numrsi
            rsicv(:) = 0;
            rsicv(irc) = 1;
            rt = ne_methods.robustt(X, conmaps, b, w, rsicv);
            rt(isinf(rt) | isnan(rt)) = 0;
            corm = correlinvtstat(rt, keepsubsn);
            corm(isinf(corm) | isnan(corm)) = 0;
            mapc.Map(tmc).(mapf) = single(corm);
            mapc.Map(tmc).Name = sprintf('%s (robust)', opts.names{icc+irc-1});
            if opts.allrs
                artv.Regressors = r(:, irc);
            else
                artv.Regressors = r(:, cr);
            end
            mapc.Map(tmc).DF1 = nval;
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        if opts.const
            rsicv(:) = 0;
            rsicv(numrsi+1) = 1;
            rt = ne_methods.robustt(X, conmaps, b, w, rsicv);
            rt(isinf(rt) | isnan(rt)) = 0;
            mapc.Map(tmc).Type = 1;
            mapc.Map(tmc).(mapf) = single(rt);
            mapc.Map(tmc).Name = sprintf('%s (robust, intercept-t)', opts.names{icc});
            mapc.Map(tmc).LowerThreshold = -sdist('tinv', thresh(1), nval);
            mapc.Map(tmc).UpperThreshold = -sdist('tinv', thresh(2), nval);
            artv.Regressors = [];
            mapc.Map(tmc).DF1 = nval;
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        if opts.estsmap
            mapc.Map(tmc).Name = sprintf('%s (smoothness estimate, robust)', opts.names{icc});
            mapc.Map(tmc).Type = 83;
            mapc.Map(tmc).(mapf) = single(artv.FWHMResImg);
            if mdim == 4
                mapc.Map(tmc).LowerThreshold = 2 * bc.Resolution + 0.5 * opts.smk;
                mapc.Map(tmc).UpperThreshold = 6 * bc.Resolution + 2.5 * opts.smk;
            else
                mapc.Map(tmc).LowerThreshold = 5;
                mapc.Map(tmc).UpperThreshold = 20;
            end
            mapc.Map(tmc).ShowPositiveNegativeFlag = 1;
            tmc = tmc + 1;
        end
        if isfield(artv, 'FWHMResEst')
            artv.FWHMResEst = [];
            artv.FWHMResImg = [];
        end
        if opts.robwmaps
            mapc.Map(tmc).Type = 145;
            mapc.Map(tmc).Name = sprintf('%s (mean robust weight)', opts.names{icc});
            mapc.Map(tmc).LowerThreshold = 0.75;
            mapc.Map(tmc).UpperThreshold = 1;
            mapc.Map(tmc).(mapf) = single((1 / size(w, ndims(w))) .* sum(w, ndims(w)));
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
            mapc.Map(tmc).Name = sprintf('%s (weighted number of outliers)', opts.names{icc});
            mapc.Map(tmc).Type = 146;
            mapc.Map(tmc).LowerThreshold = 1;
            mapc.Map(tmc).UpperThreshold = 0.5 * size(w, ndims(w));
            mapc.Map(tmc).(mapf) = single(-3 .* sum( ...
                limitrangec(w - 2/3, -1/3, 0, -1/6), ndims(w)));
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        if opts.swmaps
            for irc = find(keepsubs(:)')
                mapc.Map(tmc).Name = sprintf('%s (%s outlier)', ...
                    opts.names{icc}, ffxsubs{subsel(irc)});
                mapc.Map(tmc).Type = 144;
                mapc.Map(tmc).LowerThreshold = 0.25;
                mapc.Map(tmc).UpperThreshold = 1;
                if mdim == 4
                    mapc.Map(tmc).(mapf) = single(1 - w(:, :, :, irc));
                else
                    mapc.Map(tmc).(mapf) = single(1 - w(:, irc));
                end
                mapc.Map(tmc).RunTimeVars = artv;
                tmc = tmc + 1;
            end
        end
        if opts.meanr
            [b, w] = ne_methods.fitrobustbisquare_img(Xm, conmaps);
            artv.GlobSigMap = single(reshape(meanrmsk, rpma));
            if opts.estfwhm
                ptc = zeros(size(conmaps));
                for bmc = 1:size(Xm, 2)
                    if mdim == 4
                        ptc = ptc + repmat(b(:, :, :, bmc), [mones, size(Xm, 1)]) .* ...
                            repmat(reshape(Xm(:, bmc), [mones, size(Xm, 1)]), szmapo);
                    else
                        ptc = ptc + repmat(b(:, bmc), [mones, size(Xm, 1)]) .* ...
                            repmat(reshape(Xm(:, bmc), [mones, size(Xm, 1)]), szmapo);
                    end
                end
                ptc = w .* ptc + (1 - w) .* conmaps;
                if mdim == 4
                    [artv.FWHMResEst, artv.FWHMResImg] = ...
                        ne_methods.resestsmooth(conmaps - ptc, bc.Resolution);
                else
                    [artv.FWHMResEst, artv.FWHMResImg] = ...
                        ne_methods.resestsmoothsrf(conmaps - ptc, opts.srf);
                end
            end
            rsicv = zeros(1, size(Xm, 2));
            for irc = 1:numrsi
                rsicv(:) = 0;
                rsicv(irc) = 1;
                rt = ne_methods.robustt(Xm, conmaps, b, w, rsicv);
                rt(isinf(rt) | isnan(rt)) = 0;
                corm = correlinvtstat(rt, keepsubsn);
                corm(isinf(corm) | isnan(corm)) = 0;
                mapc.Map(tmc).(mapf) = single(corm);
                mapc.Map(tmc).Name = sprintf('%s (robust, mean-rem)', opts.names{icc+irc-1});
                mapc.Map(tmc).DF1 = nval - 1;
                mapc.Map(tmc).LowerThreshold = ...
                    correlinvtstat(-sdist('tinv', thresh(1), nval - 1), nval + 1);
                mapc.Map(tmc).UpperThreshold = ...
                    correlinvtstat(-sdist('tinv', thresh(2), nval - 1), nval + 1);
                if opts.allrs
                    artv.Regressors = r(:, irc);
                else
                    artv.Regressors = r(:, cr);
                end
                mapc.Map(tmc).RunTimeVars = artv;
                tmc = tmc + 1;
            end
            if opts.const
                rsicv(:) = 0;
                rsicv(numrsi+1) = 1;
                rt = ne_methods.robustt(Xm, conmaps, b, w, rsicv);
                rt(isinf(rt) | isnan(rt)) = 0;
                mapc.Map(tmc).Type = 1;
                mapc.Map(tmc).(mapf) = single(rt);
                mapc.Map(tmc).Name = sprintf('%s (robust, intercept-t, mean-rem)', opts.names{icc});
                mapc.Map(tmc).DF1 = nval - 1;
                mapc.Map(tmc).LowerThreshold = -sdist('tinv', thresh(1), nval - 1);
                mapc.Map(tmc).UpperThreshold = -sdist('tinv', thresh(2), nval - 1);
                artv.Regressors = [];
                mapc.Map(tmc).RunTimeVars = artv;
                tmc = tmc + 1;
            end
            if opts.estsmap
                mapc.Map(tmc).Name = sprintf('%s (smoothness estimate, robust, mean-rem)', opts.names{icc});
                mapc.Map(tmc).Type = 83;
                mapc.Map(tmc).(mapf) = single(artv.FWHMResImg);
                if mdim == 4
                    mapc.Map(tmc).LowerThreshold = 2 * bc.Resolution + 0.5 * opts.smk;
                    mapc.Map(tmc).UpperThreshold = 6 * bc.Resolution + 2.5 * opts.smk;
                else
                    mapc.Map(tmc).LowerThreshold = 5;
                    mapc.Map(tmc).UpperThreshold = 20;
                end
                mapc.Map(tmc).ShowPositiveNegativeFlag = 1;
                tmc = tmc + 1;
            end
            if opts.robwmaps
                mapc.Map(tmc).Type = 145;
                mapc.Map(tmc).Name = sprintf('%s (mean robust weight, w/o mean)', opts.names{icc});
                mapc.Map(tmc).LowerThreshold = 0.75;
                mapc.Map(tmc).UpperThreshold = 1;
                mapc.Map(tmc).(mapf) = single((1 / size(w, ndims(w))) .* sum(w, ndims(w)));
                mapc.Map(tmc).RunTimeVars = artv;
                tmc = tmc + 1;
                mapc.Map(tmc).Name = sprintf('%s (weighted number of outliers, w/o mean)', opts.names{icc});
                mapc.Map(tmc).Type = 146;
                mapc.Map(tmc).LowerThreshold = 1;
                mapc.Map(tmc).UpperThreshold = 0.5 * size(w, ndims(w));
                mapc.Map(tmc).(mapf) = single(-3 .* sum( ...
                    limitrangec(w - 2/3, -1/3, 0, -1/6), ndims(w)));
                mapc.Map(tmc).RunTimeVars = artv;
                tmc = tmc + 1;
            end
            if opts.swmaps
                for irc = 1:find(keepsubs(:)')
                    mapc.Map(tmc).Name = sprintf('%s (%s outlier, w/o mean)', ...
                        opts.names{icc}, ffxsubs{subsel(irc)});
                    mapc.Map(tmc).Type = 144;
                    mapc.Map(tmc).LowerThreshold = 0.25;
                    mapc.Map(tmc).UpperThreshold = 1;
                    if mdim == 4
                        mapc.Map(tmc).(mapf) = single(1 - w(:, :, :, irc));
                    else
                        mapc.Map(tmc).(mapf) = single(1 - w(:, irc));
                    end
                    mapc.Map(tmc).RunTimeVars = artv;
                    tmc = tmc + 1;
                end
            end
            artv.GlobSigMap = [];
        end
    end
end
if tmc <= numel(mapc.Map)
    mapc.Map(tmc:end) = [];
end

% put back
map.C = mapc;

% clear residual VTC
delete(resvtc);
