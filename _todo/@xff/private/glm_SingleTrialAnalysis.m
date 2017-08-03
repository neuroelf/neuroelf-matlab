function map = glm_SingleTrialAnalysis(xo, type, opts)
% GLM::SingleTrialAnalysis  - compute a single-trial based analysis
%
% FORMAT:       map = glm.SingleTrialAnalysis(type [, opts])
%
% Input fields:
%
%       type        one of: {'bsc'}, 'icc'
%       opts        structure with optional fields
%        .conds     conditions to run over (default: all with split trials)
%        .covs      covariates (for RFX regressions)
%        .estfwhm   estimate smoothness and store in Map.RunTimeVars (true)
%        .pbar      progress bar object
%        .robust    flag, use robust regression in addition to OLS
%        .subsel    subject selection (otherwise all subjects)
%        .voi       VOI object (used to extract betas from, required for bsc)
%        .voiidx    index into VOI list (default: all)
%
% Output fields:
%
%       map         MAP/VMP/SMP object with maps
%
% Using: correlinvtstat, cov_nd, fisherr2z, fitrobustbisquare_img, lsqueeze,
%        meannoinfnan, newnatresvmp, ranktrans, resestsmooth, robustt, sdist,
%        ztrans.

% Version:  v1.1
% Build:    16020412
% Date:     Feb-04 2016, 12:56 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2013, 2014, 2016, Jochen Weber
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
cov_nd         = ne_methods.cov_nd;
fisherr2z      = ne_methods.fisherr2z;
meannoinfnan   = ne_methods.meannoinfnan;
sdist          = ne_methods.sdist;
ztrans         = ne_methods.ztrans;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
glmfile = xo.F;
glmid = xo.L;
if isempty(glmfile)
    glmfile = glmid;
end
if bc.ProjectTypeRFX == 0 && bc.SeparatePredictors ~= 2
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
isrfx = (bc.ProjectTypeRFX > 0);
if ~isrfx
    error('neuroelf:xff:badObject', 'Only valid for RFX GLMs at this time.');
end
subids = glm_Subjects(xo);
numsubs = numel(bc.GLMData.Subject);
ffxspred = glm_SubjectPredictors(xo);
numspred = size(bc.GLMData.Subject(1).BetaMaps, ndims(bc.GLMData.Subject(1).BetaMaps));
ffxcpred = regexprep(ffxspred(~cellfun('isempty', regexpi(ffxspred, '_0+1$'))), '_0+1$', '');
numcpred = numel(ffxcpred);
if ~any(bc.ProjectType == [1, 2])
    error('neuroelf:xff:unsupported', 'RFX correlation maps of FMRs are not yet supported.');
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
macc = repmat({':'}, 1, numel(msz));
if nargin < 2 || ~ischar(type) || isempty(type) || ~any(strcmpi(type(:)', {'b', 'bsc', 'i', 'icc'}))
    type = 'b';
else
    type = lower(type(1));
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'conds') || ~iscell(opts.conds) || isempty(opts.conds)
    opts.conds = ffxcpred(:)';
else
    opts.conds = opts.conds(:)';
    for cc = 1:numel(opts.conds)
        if ~ischar(opts.conds{cc}) || isempty(opts.conds{cc})
            % MISSING
        end
    end
end
if ~isfield(opts, 'covs') || ~isa(opts.covs, 'double') || isempty(opts.covs) || ...
    any(isinf(opts.covs(:)) | isnan(opts.covs(:)) | opts.covs(:) < 1) || ...
   ~isfield(bc.RunTimeVars, 'CovariatesNames')
    opts.covs = [];
else
    opts.covs = unique(min(numel(bc.RunTimeVars.CovariatesNames), round(opts.covs(:))))';
end
if ~isfield(opts, 'estfwhm') || ~islogical(opts.estfwhm) || numel(opts.estfwhm) ~= 1
    opts.estfwhm = true;
end
if ~isfield(opts, 'meanr') || ~islogical(opts.meanr) || numel(opts.meanr) ~= 1
    opts.meanr = false;
end
if isfield(opts, 'meanrmsk') && numel(opts.meanrmsk) == 1 && xffisobject(opts.meanrmsk, true, 'msk')
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
if isempty(opts.meanrmsk) && opts.meanr
    opts.meanrmsk = all(bc.GLMData.Subject(1).BetaMaps ~= 0, ndims(bc.GLMData.Subject(1).BetaMaps));
    for sc = 1:numsubs
        opts.meanrmsk = (opts.meanrmsk & all(bc.GLMData.Subject(1).BetaMaps ~= 0, ...
            ndims(bc.GLMData.Subject(1).BetaMaps)));
    end
    opts.meanrmsk = ne_methods.lsqueeze(opts.meanrmsk);
end
if ~isfield(opts, 'pbar') || numel(opts.pbar) ~= 1 || (~isxfigure(opts.pbar) && ~isa(opts.pbar, 'xprogress'))
    opts.pbar = [];
end
pbar = opts.pbar;
if ~isfield(opts, 'robust') || ~islogical(opts.robust) || numel(opts.robust) ~= 1
    opts.robust = false;
end
if ~isfield(opts, 'subsel') || ~isa(opts.subsel, 'double') || isempty(opts.subsel) || ...
    any(isinf(opts.subsel(:)) | isnan(opts.subsel(:))) || ...
    numel(unique(round(opts.subsel(:)))) ~= numel(opts.subsel) || ...
    any(opts.subsel(:) < 1 | opts.subsel(:) > numsubs)
    opts.subsel = 1:numsubs;
else
    opts.subsel = round(opts.subsel(:)');
end
if ~isfield(opts, 'voi') || numel(opts.voi) ~= 1 || ~xffisobject(opts.voi, true, 'voi')
    opts.voi = [];
    voic = [];
else
    voic = opts.voi.C;
end
if ~isfield(opts, 'voiidx') || ~isa(opts.voiidx, 'double') || isempty(opts.voiidx) || ...
    any(isinf(opts.voiidx(:)) | isnan(opts.voiidx(:)) | opts.voiidx(:) < 1)
    if ~isempty(voic)
        opts.voiidx = 1:numel(voic.VOI);
    else
        opts.voiidx = [];
    end
else
    if ~isempty(voic)
        opts.voiidx = unique(min(numel(voic.VOI), floor(opts.voiidx(:))))';
    else
        opts.voiidx = [];
    end
end
if ~isempty(opts.voi)
    try
        vb = glm_VOIBetas(xo, opts.voi, struct('vl', opts.voiidx));
    catch xfferror
        rethrow(xfferror);
    end
end
if ~isfield(opts, 'thresh') || ~isa(opts.thresh, 'double') || numel(opts.thresh) ~= 2 || ...
    any(isinf(opts.thresh) | isnan(opts.thresh) | opts.thresh <= 0 | opts.thresh >= 0.5)
    opts.thresh = [0.005, 0.0001];
else
    opts.thresh = -sort(-opts.thresh(:)');
end
subsel = opts.subsel;
numsubs = numel(subsel);
if opts.meanr
    mvm = zeros(numsubs, size(vb, 2));
    for sc = 1:numsubs
        for mc = 1:size(mvm, 2)
            m = bc.GLMData.Subject(subsel(sc)).BetaMaps(macc{:}, mc);
            mvm(sc, mc) = meannoinfnan(m(opts.meanrmsk));
        end
    end
end
thresh = opts.thresh;
nummaps = numel(opts.conds);
numrs = numel(opts.covs);
nval = numsubs - (1 + numrs);
nummapst = nummaps * (numsubs + 1 + numrs) * (1 + opts.meanr) * (1 + opts.robust);
subsa = [macc, {[]}];
subsr = [macc, {[]}];
rpma = [msz, 1];
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
    case {2}
        map = xff('new:smp');
        mapc = map.C;
        mapc.NrOfVertices = bc.NrOfVertices;
        mapf = 'SMPData';
end

% set some common fields
mapc.NrOfMaps = nummapst;
mapc.Map.Type = 1;
mapc.Map.LowerThreshold = -sdist('tinv', thresh(1), nval);
mapc.Map.UpperThreshold = -sdist('tinv', thresh(2), nval);
mapc.Map.DF1 = nval;
mapc.Map.DF2 = 1;
mapc.Map.NrOfFDRThresholds = 0;
mapc.Map.FDRThresholds = zeros(0, 3);
mapc.Map.(mapf) = single(zeros(rpma));

% replicate
mapc.Map = mapc.Map(1, ones(1, nummapst));

% get indices for conditions
spi = cell(1, nummaps);
for cc = 1:nummaps
    spi{cc} = find(~cellfun('isempty', regexp(ffxspred, ['^' opts.conds{cc} '_\d+$'])));
end

% start with the single-subject maps first
tmi = 1;
smi = zeros(numsubs, nummaps);
for sc = 1:numsubs
    for cc = 1:nummaps
        nin = ~isnan(vb(sc, spi{cc}));
        [cv, cr] = cov_nd(double( ...
            bc.GLMData.Subject(sc).BetaMaps(macc{:}, spi{cc}(nin))), ...
            repmat(reshape(vb(sc, spi{cc}(nin)), [1, 1, 1, sum(nin)]), rpma));
        cr(isinf(cr)|isnan(cr)) = 0;
        mapc.Map(tmi).(mapf) = single(fisherr2z(cr));
        mapc.Map(tmi).Name = sprintf('%s: %s', subids{sc}, opts.conds{cc});
        mapc.Map(tmi).Type = 2;
        mapc.Map(tmi).DF1 = sum(nin) - 2;
        mapc.Map(tmi).LowerThreshold = fisherr2z( ...
            correlinvtstat(-sdist('tinv', 0.001, sum(nin) - 2), sum(nin)));
        mapc.Map(tmi).UpperThreshold = fisherr2z( ...
            correlinvtstat(-sdist('tinv', 1e-6, sum(nin) - 2), sum(nin)));
        smi(sc, cc) = tmi;
        tmi = tmi + 1;
    end
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
                    conmaps = subsasgn(conmaps, subsas, subsref(conmaps, subsas) + c(pc, cc) .* ...
                        subsref(bc.GLMData.Subject(subsel(sc)).BetaMaps, subsrs));
                end
            end
        end

        % minumum criterion not met? then set to 0 for all subjects
        conmapsa = (~isinf(conmaps));
        conmapsa = conmapsa & (~isnan(conmaps));
        conmapsa = conmapsa & (conmaps ~= 0);
        conmapsa = sum(conmapsa, numel(subsr));
        conmaps(repmat(conmapsa < opts.minnum, [ones(1, numel(msz)), numsubs])) = 0;

        % rank transform?
        if opts.rank
            conmaps = ne_methods.ranktrans(conmaps, numel(subsr), ...
                struct('meancenter', true, 'nozero', true));
        else
            conmaps(isinf(conmaps) | isnan(conmaps)) = 0;
        end

        % keep track!
        lcc = cc;
    end

    % generate design matrix/ces
    if opts.allrs
        X = [ztrans(r), ones(numsubs, 1)];
    else
        X = [ztrans(r(:, cr)), ones(numsubs, 1)];
    end
    iXX = pinv(X' * X);
    if opts.meanr
        for sc = 1:numsubs
            subsa{end} = sc;
            subsas = struct('type', '()', 'subs', {subsa});
            mv = subsref(conmaps, subsas);
            mvm(sc) = mean(mv(opts.meanrmsk));
        end
        if opts.rank
            Xm = [X, ztrans(ne_methods.ranktrans(mvm, 1))];
        else
            Xm = [X, ztrans(mvm)];
        end
        iXXm = pinv(Xm' * Xm);
    end

    % set additional data
    artv = struct( ...
        'SourceGLM',   glmfile, ...
        'SourceGLMID', glmid, ...
        'Contrast',    c(:, lcc), ...
        'MeanRem',     opts.meanr, ...
        'Regressors',  r, ...
        'RFXGLM',      true, ...
        'Robust',      opts.robust, ...
        'SubPreds',    {ffxspred}, ...
        'SubSel',      subsel(:));

    % OLS computations first
    betas = iXX * X' * reshape(conmaps, prod(msz), numsubs)';
    resim = conmaps - reshape((X * betas)', [msz, numsubs]);
    stder = sqrt(sum(resim .^ 2, ndims(conmaps))) .* sqrt(1 / nval);
    if opts.estfwhm
        [artv.FWHMResEst, artv.FWHMResImg] = ne_methods.resestsmooth(resim, bc.Resolution);
    end

    % first maps
    for irc = 1:numrsi
        tmap = betas(irc, :)' ./ (stder(:) .* sqrt(iXX(irc, irc)));
        tmap(isinf(tmap) | isnan(tmap)) = 0;
        mapc.Map(tmc).(mapf) = reshape(single(correlinvtstat(tmap, numsubs)), [msz, 1]);
        mapc.Map(tmc).Name = opts.names{icc+irc-1};
        if opts.allrs
            artv.Regressors = r(:, irc);
        else
            artv.Regressors = r(:, cr);
        end
        mapc.Map(tmc).RunTimeVars = artv;
        tmc = tmc + 1;
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
        mapc.Map(tmc).RunTimeVars = artv;
        tmc = tmc + 1;
    end

    % with mean removed ?
    if opts.meanr
        betas = iXXm * Xm' * reshape(conmaps, prod(msz), numsubs)';
        resim = conmaps - reshape((Xm * betas)', [msz, numsubs]);
        stder = sqrt(sum(resim .^ 2, ndims(conmaps))) .* sqrt(1 / (nval - 1));
        if opts.estfwhm
            [artv.FWHMResEst, artv.FWHMResImg] = ne_methods.resestsmooth(resim, bc.Resolution);
        end

        % first maps
        for irc = 1:numrsi
            tmap = betas(irc, :)' ./ (stder(:) .* sqrt(iXXm(irc, irc)));
            tmap(isinf(tmap) | isnan(tmap)) = 0;
            mapc.Map(tmc).(mapf) = reshape(single(correlinvtstat(tmap, numsubs)), [msz, 1]);
            mapc.Map(tmc).Name = sprintf('%s (mean-rem)', opts.names{icc+irc-1});
            mapc.Map(tmc).DF1 = nval - 1;
            mapc.Map(tmc).LowerThreshold = correlinvtstat(-sdist('tinv', thresh(1), nval - 1), nval + 1);
            mapc.Map(tmc).UpperThreshold = correlinvtstat(-sdist('tinv', thresh(2), nval - 1), nval + 1);
            if opts.allrs
                artv.Regressors = r(:, irc);
            else
                artv.Regressors = r(:, cr);
            end
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
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
    end

    % compute robust stats
    if opts.robust

        % perform fit
        [b, w] = ne_methods.fitrobustbisquare_img(X, conmaps);
        if opts.estfwhm
            ptc = zeros(size(conmaps));
            for bmc = 1:size(X, 2)
                ptc = ptc + repmat(b(:, :, :, bmc), [1, 1, 1, size(X, 1)]) .* ...
                    repmat(reshape(X(:, bmc), [1, 1, 1, size(X, 1)]), szmap);
            end
            ptc = w .* ptc + (1 - w) .* conmaps;
            [artv.FWHMResEst, artv.FWHMResImg] = ne_methods.resestsmooth(conmaps - ptc, bc.Resolution);
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
            mapc.Map(tmc).Name = sprintf('%s (robust)', opts.names{icc+irc-1});
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
            rt = ne_methods.robustt(X, conmaps, b, w, rsicv);
            rt(isinf(rt) | isnan(rt)) = 0;
            mapc.Map(tmc).Type = 1;
            mapc.Map(tmc).(mapf) = single(rt);
            mapc.Map(tmc).Name = sprintf('%s (robust, intercept-t)', opts.names{icc});
            mapc.Map(tmc).LowerThreshold = -sdist('tinv', thresh(1), nval);
            mapc.Map(tmc).UpperThreshold = -sdist('tinv', thresh(2), nval);
            artv.Regressors = [];
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        if opts.meanr
            [b, w] = ne_methods.fitrobustbisquare_img(Xm, conmaps);
            if opts.estfwhm
                ptc = zeros(size(conmaps));
                for bmc = 1:size(Xm, 2)
                    ptc = ptc + repmat(b(:, :, :, bmc), [1, 1, 1, size(Xm, 1)]) .* ...
                        repmat(reshape(Xm(:, bmc), [1, 1, 1, size(Xm, 1)]), szmap);
                end
                ptc = w .* ptc + (1 - w) .* conmaps;
                [artv.FWHMResEst, artv.FWHMResImg] = ne_methods.resestsmooth(conmaps - ptc, bc.Resolution);
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
                mapc.Map(tmc).Name = sprintf('%s (robust, mean-rem)', opts.names{icc+irc-1});
                mapc.Map(tmc).DF1 = nval - 1;
                mapc.Map(tmc).LowerThreshold = correlinvtstat(-sdist('tinv', thresh(1), nval - 1), nval + 1);
                mapc.Map(tmc).UpperThreshold = correlinvtstat(-sdist('tinv', thresh(2), nval - 1), nval + 1);
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
        end
    end
end

% put back
map.C = mapc;
