function [map, mse] = glm_RFX_StatsMap(xo, c, x, c2, mapopts)
% GLM::RFX_StatsMap  - calculate a second-level (RFX) map
%
% FORMAT:       map = glm.RFX_StatsMap(c, x [, c2 [, mapopts]])
%
% Input fields:
%
%       c           1xC contrast vector (used to compute first level maps)
%       x           SxR second level design matrix (full design)
%       c2          NxC second level contrasts, either of
%                   non-NaN  (regular t-stats)
%                   some-NaN (non-NaN indices are reduced model, F-stats)
%       mapopts     structure with optional fields
%        .bvcomp    BV-compatible map names (length restriction, true)
%        .c2names   contrast names (default: try to construct from c2)
%        .cname     contrast name (default: try to construct from c)
%        .estfwhm   estimate smoothness and store in Map.RunTimeVars (true)
%        .meanr     boolean flag, add mean of first-level maps to model
%        .meanrmsk  mask to mean from (object or XxYxZ logical, def: all)
%        .minnum    minumum number of subjects to compute (2 * sqrt(N))
%        .robust    flag, use robust regression in addition to OLS
%        .rweights  if provided, either Sx1 (1xS) or XxYxZxS double data
%        .subsel    subject selection (otherwise all subjects)
%        .thresh    1x2 threshold (lower, upper), as p-values!
%
% Output fields:
%
%       map         MAP/VMP/SMP object with maps
%
% Note: for a reduced-models test, the MSE is computed from the actually
%       reduced model (not with an internal contrast logic), which will
%       require additional run-time for every model!
%
% Using: correlinvtstat, findfirst, fitrobustbisquare_img, lsqueeze,
%        newnatresvmp, ranktrans, resestsmooth, robustt, sdist, ztrans.

% Version:  v1.1
% Build:    16012916
% Date:     Jan-29 2016, 4:06 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2012, 2014, 2016, Jochen Weber
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
ranktrans      = ne_methods.ranktrans;
sdist          = ne_methods.sdist;

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm') || ...
   ~isa(c, 'double') || isempty(c) || any(isnan(c(:)) | isinf(c(:))) || ...
   ~isa(x, 'double') || isempty(x) || any(isinf(x(:))) || size(x, 2) >= size(x, 1)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if size(x, 1) ~= bc.NrOfSubjects
    error('neuroelf:xff:badArgument', 'NrOfSubjects must match with design matrix size.');
end
if nargin < 4 || ~isa(c2, 'double') || isempty(c2) || any(isinf(c2(:)) | isnan(c2(:)))
    c2 = eye(size(x, 1));
end
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
if isrfx
    numsubs = numel(bc.GLMData.Subject);
    numspred = size(bc.GLMData.Subject(1).BetaMaps, ndims(bc.GLMData.Subject(1).BetaMaps));
else
    ffxpred = bc.Predictor;
    ffxpred = {ffxpred(:).Name2};
    ffxpred = ffxpred(:);
    ffxsubs = glm_Subjects(xo);
    numsubs = numel(ffxsubs);
    numspred = numel(ffxspred) + 1;
end
if ~any(numel(c, 2) == [numspred, numspred - 1])
    error('neuroelf:xff:badArgument', 'Invalid first-level contrast spec.');
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
cmaps = glm_RFX_conmaps(xo, c);
if nargin < 5 || ~isstruct(mapopts) || numel(mapopts) ~= 1
    mapopts = struct;
end
if ~isfield(mapopts, 'bvcomp') || ~islogical(mapopts.bvcomp) || numel(mapopts.bvcomp) ~= 1
    mapopts.bvcomp = true;
end
if ~isfield(mapopts, 'c2names') || ~iscell(mapopts.c2names)
    mapopts.c2names = {};
end
if ~isfield(mapopts, 'estfwhm') || ~islogical(mapopts.estfwhm) || numel(mapopts.estfwhm) ~= 1
    mapopts.estfwhm = true;
end
if ~isfield(mapopts, 'meanr') || ~islogical(mapopts.meanr) || numel(mapopts.meanr) ~= 1
    mapopts.meanr = false;
end
if isfield(mapopts, 'meanrmsk') && numel(mapopts.meanrmsk) == 1 && ...
    xffisobject(mapopts.meanrmsk, true, 'msk')
    mbc = mapopts.meanrmsk.C;
    if numel(mbc.Mask) == prod(msz)
        mapopts.meanrmsk = ne_methods.lsqueeze(mbc.Mask > 0);
    else
        mapopts.meanrmsk = [];
    end
elseif isfield(mapopts, 'meanrmsk') && islogical(mapopts.meanrmsk) && ...
    numel(mapopts.meanrmsk) == prod(msz)
    mapopts.meanrmsk = ne_methods.lsqueeze(mapopts.meanrmsk);
else
    mapopts.meanrmsk = [];
end
if isempty(mapopts.meanrmsk) && mapopts.meanr
    mapopts.meanrmsk = all(bc.GLMData.Subject(1).BetaMaps ~= 0, ...
        ndims(bc.GLMData.Subject(1).BetaMaps));
    for sc = 1:numsubs
        mapopts.meanrmsk = (mapopts.meanrmsk & ...
            all(bc.GLMData.Subject(1).BetaMaps ~= 0, ndims(bc.GLMData.Subject(1).BetaMaps)));
    end
    mapopts.meanrmsk = ne_methods.lsqueeze(mapopts.meanrmsk);
end
if ~isfield(mapopts, 'minnum') || ~isa(mapopts.minnum, 'double') || numel(mapopts.minnum) ~= 1 || ...
    isinf(mapopts.minnum) || isnan(mapopts.minnum) || mapopts.minnum < 0
    mapopts.minnum = 0;
end
if ~isfield(mapopts, 'names') || ~iscell(mapopts.names) || isempty(mapopts.names)
    mapopts.names = {};
end
if ~isfield(mapopts, 'rank') || ~islogical(mapopts.rank) || numel(mapopts.rank) ~= 1 || ...
   ~mapopts.rank
    mapopts.rank = false;
    ranktxt = '';
else
    ranktxt = 'Rank-';
end
if ~isfield(mapopts, 'rnames') || ~iscell(mapopts.rnames)
    mapopts.rnames = {};
end
if ~isfield(mapopts, 'robust') || ~islogical(mapopts.robust) || numel(mapopts.robust) ~= 1
    mapopts.robust = false;
end
if ~isfield(mapopts, 'subsel') || ~isa(mapopts.subsel, 'double') || isempty(mapopts.subsel) || ...
    any(isinf(mapopts.subsel(:)) | isnan(mapopts.subsel(:))) || ...
    numel(unique(round(mapopts.subsel(:)))) ~= numel(mapopts.subsel) || ...
    any(mapopts.subsel(:) < 1 | mapopts.subsel(:) > numsubs)
    mapopts.subsel = 1:numsubs;
else
    mapopts.subsel = round(mapopts.subsel(:)');
end
if ~isfield(mapopts, 'voiidx') || ~isa(mapopts.voiidx, 'double') || numel(mapopts.voiidx) ~= 1 || ...
    isinf(mapopts.voiidx) || isnan(mapopts.voiidx) || mapopts.voiidx < 1
    mapopts.voiidx = 1;
else
    mapopts.voiidx = floor(mapopts.voiidx);
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
        r = glm_VOIBetas(xo, r, struct('c', ct, 'vl', min(mapopts.voiidx, rbc)));
    catch xfferror
        rethrow(xfferror);
    end
end
if any(size(r) == numsubs) && numel(mapopts.subsel) ~= numsubs
    rsr = repmat({':'}, 1, ndims(r));
    rsr{ne_methods.findfirst(size(r) == numsubs)} = mapopts.subsel;
    r = r(rsr{:});
end
if ~any(size(r) == numel(mapopts.subsel))
    error('neuroelf:xff:badArgument', ...
        'Correlation regressors must match in size with number of subjects.');
end
rsdim = ne_methods.findfirst(size(r) == numel(mapopts.subsel));
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
    mapopts.subsel = mapopts.subsel(~nanr);
end
if ~isfield(mapopts, 'thresh') || ~isa(mapopts.thresh, 'double') || numel(mapopts.thresh) ~= 2 || ...
    any(isinf(mapopts.thresh) | isnan(mapopts.thresh) | mapopts.thresh <= 0 | mapopts.thresh >= 0.5)
    mapopts.thresh = [0.005, 0.0001];
else
    mapopts.thresh = -sort(-mapopts.thresh(:)');
end
subsel = mapopts.subsel;
numsubs = numel(subsel);
if mapopts.minnum == 0
    mapopts.minnum = ceil(2 * sqrt(numsubs));
elseif mapopts.minnum < 1
    mapopts.minnum = ceil(mapopts.minnum * numsubs);
end
mapopts.minnum = min(numsubs, mapopts.minnum);
if mapopts.meanr
    mvm = zeros(numsubs, 1);
end
thresh = mapopts.thresh;
if size(r, 2) == numel(mapopts.subsel) && size(r, 1) ~= numel(mapopts.subsel)
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
if mapopts.allrs
    nval = numsubs - (1 + numrs);
else
    nval = numsubs - 2;
end
if numel(mapopts.cnames) ~= nummaps
    sprednames = glm_SubjectPredictors(xo);
    mapopts.cnames = cell(1, nummaps);
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
            mapopts.cnames{cc} = sprintf('%s > %s', connc, coffnc);
        else
            mapopts.cnames{cc} = connc;
        end
    end
end
if numel(mapopts.rnames) ~= numrs
    mapopts.rnames = cell(1, numrs);
    for cc = 1:numrs
        mapopts.rnames = sprintf('reg%02d', cc);
    end
end
nummapst = nummaps * (numrs + mapopts.const) * (1 + mapopts.meanr) * (1 + mapopts.robust);
if bc.ProjectType == 1
    subsa = {':', ':', ':', []};
    subsr = {':', ':', ':', []};
else
    subsa = {':', []};
    subsr = {':', []};
end
rpma = [msz, 1];
if numel(mapopts.names) ~= (nummaps * numrs)
    mapopts.names = cell(1, nummaps * numrs);
end
for cc = 1:numel(mapopts.names)
    if ~ischar(mapopts.names{cc})
        ccr = 1 + mod(cc - 1, numrs);
        ccc = 1 + round((cc - ccr) / numrs);
        mapopts.names{cc} = sprintf('%sCorr: (%s, %s)', ranktxt, ...
            mapopts.cnames{ccc}, mapopts.rnames{ccr});
    end
    if mapopts.bvcomp && numel(mapopts.names{cc}) > 96
        mapopts.names{cc} = [mapopts.names{cc}(1:46) ' ... ' mapopts.names{cc}(end-45:end)];
    else
        mapopts.names{cc} = mapopts.names{cc}(:)';
    end
end
if isrfx
    szmap = size(bc.GLMData.RFXGlobalMap);
else
    szmap = size(bc.GLMData.MCorrSS);
end

% create map container
switch (bc.ProjectType)

    % VTC/VMP
    case 1
        map = ne_methods.newnatresvmp();
        mapc = map.C;
        mapc.XStart = bc.XStart;
        mapc.XEnd = bc.XEnd;
        mapc.YStart = bc.YStart;
        mapc.YEnd = bc.YEnd;
        mapc.ZStart = bc.ZStart;
        mapc.ZEnd = bc.ZEnd;
        mapc.Resolution = bc.Resolution;
        mapc.RunTimeVars.TrfPlus = bc.RunTimeVars.TrfPlus;
        mapf = 'VMPData';

    % MTC/SMP
    case 2
        map = xff('new:smp');
        mapc = map.C;
        mapc.NrOfVertices = bc.NrOfVertices;
        mapf = 'SMPData';
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
if mapopts.rank
    r = ranktrans(r, 1);
end

% what models?
if mapopts.allrs
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

        % fill contrast maps
        for pc = 1:numspred
            if c(pc, cc) ~= 0
                subsr{end} = pc;
                subsrs = struct('type', '()', 'subs', {subsr});
                for sc = 1:numsubs
                    subsa{end} = sc;
                    subsas = struct('type', '()', 'subs', {subsa});
                    conmaps = subsasgn(conmaps, subsas, ...
                        subsref(conmaps, subsas) + c(pc, cc) .* ...
                        subsref(bc.GLMData.Subject(subsel(sc)).BetaMaps, subsrs));
                end
            end
        end

        % minumum criterion not met? then set to 0 for all subjects
        conmapsa = (~isinf(conmaps));
        conmapsa = conmapsa & (~isnan(conmaps));
        conmapsa = conmapsa & (conmaps ~= 0);
        conmapsa = sum(conmapsa, numel(subsr));
        conmaps(repmat(conmapsa < mapopts.minnum, [ones(1, numel(msz)), numsubs])) = 0;

        % rank transform?
        if mapopts.rank
            conmaps = ranktrans(conmaps, numel(subsr), ...
                struct('meancenter', true, 'nozero', true));
        else
            conmaps(isinf(conmaps) | isnan(conmaps)) = 0;
        end

        % keep track!
        lcc = cc;
    end

    % generate design matrix/ces
    if mapopts.allrs
        X = [ne_methods.ztrans(r), ones(numsubs, 1)];
    else
        X = [ne_methods.ztrans(r(:, cr)), ones(numsubs, 1)];
    end
    iXX = pinv(X' * X);
    if mapopts.meanr
        for sc = 1:numsubs
            subsa{end} = sc;
            subsas = struct('type', '()', 'subs', {subsa});
            mv = subsref(conmaps, subsas);
            mvm(sc) = mean(mv(mapopts.meanrmsk));
        end
        if mapopts.rank
            Xm = [X, ne_methods.ztrans(ranktrans(mvm, 1))];
        else
            Xm = [X, ne_methods.ztrans(mvm)];
        end
        iXXm = pinv(Xm' * Xm);
    end

    % set additional data
    artv = struct( ...
        'SourceGLM',   glmfile, ...
        'SourceGLMID', glmid, ...
        'Contrast',    c(:, lcc), ...
        'MeanRem',     mapopts.meanr, ...
        'Regressors',  r, ...
        'RFXGLM',      true, ...
        'Robust',      mapopts.robust, ...
        'SubPreds',    {ffxspred}, ...
        'SubSel',      subsel(:));

    % OLS computations first
    betas = iXX * X' * reshape(conmaps, prod(msz), numsubs)';
    resim = conmaps - reshape((X * betas)', [msz, numsubs]);
    stder = sqrt(sum(resim .^ 2, ndims(conmaps))) .* sqrt(1 / nval);
    if mapopts.estfwhm
        [artv.FWHMResEst, artv.FWHMResImg] = ne_methods.resestsmooth(resim, bc.Resolution);
    end

    % first maps
    for irc = 1:numrsi
        tmap = betas(irc, :)' ./ (stder(:) .* sqrt(iXX(irc, irc)));
        tmap(isinf(tmap) | isnan(tmap)) = 0;
        mapc.Map(tmc).(mapf) = reshape(single(correlinvtstat(tmap, numsubs)), [msz, 1]);
        mapc.Map(tmc).Name = mapopts.names{icc+irc-1};
        if mapopts.allrs
            artv.Regressors = r(:, irc);
        else
            artv.Regressors = r(:, cr);
        end
        mapc.Map(tmc).RunTimeVars = artv;
        tmc = tmc + 1;
    end
    if mapopts.const
        tmap = betas(numrsi+1, :)' ./ (stder(:) .* sqrt(iXX(numrsi+1, numrsi+1)));
        tmap(isinf(tmap) | isnan(tmap)) = 0;
        mapc.Map(tmc).(mapf) = reshape(single(tmap), [msz, 1]);
        mapc.Map(tmc).Type = 1;
        mapc.Map(tmc).Name = sprintf('%s (intercept-t)', mapopts.names{icc});
        mapc.Map(tmc).LowerThreshold = -sdist('tinv', thresh(1), nval);
        mapc.Map(tmc).UpperThreshold = -sdist('tinv', thresh(2), nval);
        artv.Regressors = [];
        mapc.Map(tmc).RunTimeVars = artv;
        tmc = tmc + 1;
    end

    % with mean removed ?
    if mapopts.meanr
        betas = iXXm * Xm' * reshape(conmaps, prod(msz), numsubs)';
        resim = conmaps - reshape((Xm * betas)', [msz, numsubs]);
        stder = sqrt(sum(resim .^ 2, ndims(conmaps))) .* sqrt(1 / (nval - 1));
        if mapopts.estfwhm
            [artv.FWHMResEst, artv.FWHMResImg] = ...
                ne_methods.resestsmooth(resim, bc.Resolution);
        end

        % first maps
        for irc = 1:numrsi
            tmap = betas(irc, :)' ./ (stder(:) .* sqrt(iXXm(irc, irc)));
            tmap(isinf(tmap) | isnan(tmap)) = 0;
            mapc.Map(tmc).(mapf) = reshape(single(correlinvtstat(tmap, numsubs)), [msz, 1]);
            mapc.Map(tmc).Name = sprintf('%s (mean-rem)', mapopts.names{icc+irc-1});
            mapc.Map(tmc).DF1 = nval - 1;
            mapc.Map(tmc).LowerThreshold = correlinvtstat(-sdist('tinv', thresh(1), nval - 1), nval + 1);
            mapc.Map(tmc).UpperThreshold = correlinvtstat(-sdist('tinv', thresh(2), nval - 1), nval + 1);
            if mapopts.allrs
                artv.Regressors = r(:, irc);
            else
                artv.Regressors = r(:, cr);
            end
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        if mapopts.const
            tmap = betas(numrsi+1, :)' ./ (stder(:) .* sqrt(iXXm(numrsi+1, numrsi+1)));
            tmap(isinf(tmap) | isnan(tmap)) = 0;
            mapc.Map(tmc).(mapf) = reshape(single(tmap), [msz, 1]);
            mapc.Map(tmc).Type = 1;
            mapc.Map(tmc).DF1 = nval - 1;
            mapc.Map(tmc).Name = sprintf('%s (intercept-t, mean-rem)', mapopts.names{icc});
            mapc.Map(tmc).LowerThreshold = -sdist('tinv', thresh(1), nval - 1);
            mapc.Map(tmc).UpperThreshold = -sdist('tinv', thresh(2), nval - 1);
            artv.Regressors = [];
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
    end

    % compute robust stats
    if mapopts.robust

        % perform fit
        [b, w] = ne_methods.fitrobustbisquare_img(X, conmaps);
        if mapopts.estfwhm
            ptc = zeros(size(conmaps));
            for bmc = 1:size(X, 2)
                ptc = ptc + repmat(b(:, :, :, bmc), [1, 1, 1, size(X, 1)]) .* ...
                    repmat(reshape(X(:, bmc), [1, 1, 1, size(X, 1)]), szmap);
            end
            ptc = w .* ptc + (1 - w) .* conmaps;
            [artv.FWHMResEst, artv.FWHMResImg] = ...
                ne_methods.resestsmooth(conmaps - ptc, bc.Resolution);
        end
        rsicv = zeros(1, size(X, 2));
        for irc = 1:numrsi
            rsicv(:) = 0;
            rsicv(irc) = 1;
            rt = ne_methods.robustt(X, conmaps, b, w, rsicv);
            rt(isinf(rt) | isnan(rt)) = 0;
            corm = correlinvtstat(rt, numsubs);
            corm(isinf(corm) | isnan(corm)) = 0;
            mapc.Map(tmc).(mapf) = single(corm);
            mapc.Map(tmc).Name = sprintf('%s (robust)', mapopts.names{icc+irc-1});
            if mapopts.allrs
                artv.Regressors = r(:, irc);
            else
                artv.Regressors = r(:, cr);
            end
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        if mapopts.const
            rsicv(:) = 0;
            rsicv(numrsi+1) = 1;
            rt = ne_methods.robustt(X, conmaps, b, w, rsicv);
            rt(isinf(rt) | isnan(rt)) = 0;
            mapc.Map(tmc).Type = 1;
            mapc.Map(tmc).(mapf) = single(rt);
            mapc.Map(tmc).Name = sprintf('%s (robust, intercept-t)', mapopts.names{icc});
            mapc.Map(tmc).LowerThreshold = -sdist('tinv', thresh(1), nval);
            mapc.Map(tmc).UpperThreshold = -sdist('tinv', thresh(2), nval);
            artv.Regressors = [];
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        if mapopts.meanr
            [b, w] = ne_methods.fitrobustbisquare_img(Xm, conmaps);
            if mapopts.estfwhm
                ptc = zeros(size(conmaps));
                for bmc = 1:size(Xm, 2)
                    ptc = ptc + repmat(b(:, :, :, bmc), [1, 1, 1, size(Xm, 1)]) .* ...
                        repmat(reshape(Xm(:, bmc), [1, 1, 1, size(Xm, 1)]), szmap);
                end
                ptc = w .* ptc + (1 - w) .* b;
                [artv.FWHMResEst, artv.FWHMResImg] = ...
                    ne_methods.resestsmooth(conmaps - ptc, bc.Resolution);
            end
            rsicv = zeros(1, size(Xm, 2));
            for irc = 1:numrsi
                rsicv(:) = 0;
                rsicv(irc) = 1;
                rt = ne_methods.robustt(Xm, conmaps, b, w, rsicv);
                rt(isinf(rt) | isnan(rt)) = 0;
                corm = correlinvtstat(rt, numsubs);
                corm(isinf(corm) | isnan(corm)) = 0;
                mapc.Map(tmc).(mapf) = single(corm);
                mapc.Map(tmc).Name = sprintf('%s (robust, mean-rem)', mapopts.names{icc+irc-1});
                mapc.Map(tmc).DF1 = nval - 1;
                mapc.Map(tmc).LowerThreshold = ...
                    correlinvtstat(-sdist('tinv', thresh(1), nval - 1), nval + 1);
                mapc.Map(tmc).UpperThreshold = ...
                    correlinvtstat(-sdist('tinv', thresh(2), nval - 1), nval + 1);
                if mapopts.allrs
                    artv.Regressors = r(:, irc);
                else
                    artv.Regressors = r(:, cr);
                end
                mapc.Map(tmc).RunTimeVars = artv;
                tmc = tmc + 1;
            end
            if mapopts.const
                rsicv(:) = 0;
                rsicv(numrsi+1) = 1;
                rt = ne_methods.robustt(Xm, conmaps, b, w, rsicv);
                rt(isinf(rt) | isnan(rt)) = 0;
                mapc.Map(tmc).Type = 1;
                mapc.Map(tmc).(mapf) = single(rt);
                mapc.Map(tmc).Name = sprintf('%s (robust, intercept-t, mean-rem)', mapopts.names{icc});
                mapc.Map(tmc).DF1 = nval - 1;
                mapc.Map(tmc).LowerThreshold = -sdist('tinv', thresh(1), nval - 1);
                mapc.Map(tmc).UpperThreshold = -sdist('tinv', thresh(2), nval - 1);
                artv.Regressors = [];
                mapc.Map(tmc).RunTimeVars = artv;
                tmc = tmc + 1;
            end
        end
    end
end

% put back
map.C = mapc;
