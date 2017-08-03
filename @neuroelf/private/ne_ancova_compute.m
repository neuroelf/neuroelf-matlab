function varargout = ne_ancova_compute(varargin)
% ne_ancova_compute  - compute AN(C)OVA for selected GLM (from UI config)
%
% FORMAT:       ne_ancova_compute([SRC, EVT, varargin])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       varargin    no inputs
%
% Example:
%
%   ne_ancova_compute(0, 0, []);
%
% Notes: the ANCOVA UI *must* be loaded at the time of the call.

% Version:  v1.1
% Build:    17061221
% Date:     Jun-12 2017, 9:25 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, 2017, Jochen Weber
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

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg.AC;
ch = ne_gcfg.h.AC.h;
cp = ne_gcfg.h.Progress;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get access to GLM
glm = ne_gcfg.h.AC.ACFig.UserData.lastglm;
if glm.ProjectTypeRFX == 0
    uiwait(warndlg('Only implemented for RFX-GLMs for now.', 'NeuroElf - info', 'modal'));
    return;
end
srf = [];
if glm.ProjectType < 2
    mtype = 'v';
else
    mtype = 'm';
    glmh = glm.Handles;
    if isfield(glmh, 'Surface') && isxff(glmh.Surface, 'srf')
        srf = glmh.Surface;
    end
end
rtv = glm.RunTimeVars;
allsubs = lsqueeze(glm.Subjects);
subpreds = glm.SubjectPredictors;
cons = rtv.Contrasts(:, 2);

% group configuration
if ch.UseGroups.Value > 0
    groupi = ch.Groups.Value;
    groups = lsqueezecells(rtv.Groups(groupi, 2));
else
    groupi = [];
    groups = {ch.Subjects.Value};
end
if any(cellfun('prodofsize', groups) < 3)
    uiwait(warndlg('Subjects -> too small group (N<3) detected.', ...
        'NeuroElf - warning', 'modal'));
    return;
end
numgroups = numel(groups);
subjnums = cat(1, groups{:});
if numel(subjnums) ~= numel(unique(subjnums))
    uiwait(warndlg('Subjects -> Groups assignment not correct.', ...
        'NeuroElf - warning', 'modal'));
    return;
end
subjects = allsubs(subjnums);
subgroups = groups;
si = 1;
for gc = 1:numgroups
    subgroups{gc} = (si:si+(numel(groups{gc})-1))';
    si = si + numel(groups{gc});
end
numsubs = numel(subjects);

% get options -> covariates
covv = ch.Covs.Value;
coviact = (ch.CovInteractions.Value > 0);
covnames = rtv.CovariatesNames(covv);
covs = rtv.CovariatesData(subjnums, covv);
if isempty(covs)
    covc = '';
else
    for cvc = 1:size(covs, 2)
        for gc = 1:numgroups
            if sum(~isnan(covs(subgroups{gc}, cvc))) < 3
                uiwait(warndlg(sprintf( ...
                    'Covariates -> Too many NaN values for %s in group %d', ...
                    covnames{cvc}, gc), 'NeuroElf - warning', 'modal'));
                return;
            end
        end
        if numel(covnames{cvc}) > 24
            covnames{cvc} = sprintf('%s....%s', ...
                covnames{cvc}(1:10), covnames{cvc}(end-9:end));
        end
    end
    covc = 'C';
end
if ch.MRCovNone.Value > 0
    mrcov = 0;
elseif ch.MRCovGroup.Value > 0
    mrcov = 2;
else
    mrcov = 1;
end

% smooth data?
smk = 0;
lbmin = -Inf;
lbmax = Inf;
if ch.SmoothData.Value > 0 && mtype == 'v'
    try
        smk = str2double(ch.SmoothDataKernel.String);
        if numel(smk) ~= 1 || ...
            isinf(smk) || ...
            isnan(smk) || ...
            smk < 0
            smk = 0;
        elseif ch.BRange.Value > 0
            lbmin = str2double(ch.BRangeMin.String);
            lbmax = str2double(ch.BRangeMax.String);
            if numel(lbmin) ~= 1 || ...
                numel(lbmax) ~= 1 || ...
                isnan(lbmin) || ...
                isnan(lbmax) || ...
                lbmin > 0 || ...
                lbmax < 0 || ...
                lbmin >= lbmax
                lbmin = -Inf;
                lbmax = Inf;
            end
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        smk = 0;
        lbmin = -Inf;
        lbmax = Inf;
    end
end

% factors
layout1 = ch.Layout1.Value;
layoutX = ch.LayoutX.Value;
switch layoutX
    case {1}
        levx = [1, 1];
    case {2}
        levx = [2, 1];
    case {3}
        levx = [3, 1];
    case {4}
        levx = [2, 2];
    case {5}
        levx = [3, 2];
end

% regular factor
numst = 0;
if layout1 <= 10
    flevels = [layout1, levx];
    flevcis = zeros(flevels);
    for lc3 = 1:flevels(3)
        for lc2 = 1:flevels(2)
            for lc1 = 1:flevels(1)
                flevcis(lc1, lc2, lc3) = ch.CellCont{lc1, lc2, lc3}.Value - 1;
            end
        end
    end
    
% FIR
elseif layout1 == 11

    % TODO: get number of FIR regressors per condition
    numfirs = ch.CellCont{1, 1, 1}.String{ch.CellCont{1, 1, 1}.Value};
    numfirs = str2double(regexprep(numfirs, '^.*_D(\d+)$', '$1')) + 1;
    if isnan(numfirs) || numfirs < 2
        uiwait(warndlg('Invalid number of FIR predictors in first factor.', ...
            'NeuroElf - warning', 'modal'));
        return;
    end
    flevels = [numfirs, levx];
    flevcis = zeros(flevels);
    for lc3 = 1:flevels(3)
        for lc2 = 1:flevels(2)
            firmain = ch.CellCont{lc2, lc3, 1}.String{ch.CellCont{lc2, lc3, 1}.Value};
            firmain = regexprep(firmain, '^(.*_)D\d+$', '$1D0');
            if ~strcmp(firmain(end-1:end), 'D0')
                uiwait(warndlg('Invalid FIR predictors selected.', ...
                    'NeuroElf - warning', 'modal'));
                return;
            end
            firmain = findfirst(strcmpi(subpreds, firmain));
            if isempty(firmain)
                uiwait(warndlg('Invalid FIR predictors selected.', ...
                    'NeuroElf - warning', 'modal'));
                return;
            end
            for lc1 = 1:flevels(1)
                flevcis(lc1, lc2, lc3) = firmain - 1 + lc1;
            end
        end
    end
    
% single-trials
else

    % TODO: get number of single trials (max.)
    numst = 1;
    flevels = [numst, levx];
end
tlevels = prod(flevels);
l1 = flevels(1) - 1;
l2 = flevels(2) - 1;
l3 = flevels(3) - 1;
l12 = l1 * l2;
l13 = l1 * l3;
l23 = l2 * l3;
l123 = l1 * l2 * l3;
tlevcs = l1 + l2 + l3 + l12 + l13 + l23 + l123;

% invalid
if any(flevcis(:) < 1)
    return;
end

% hide figure and show progress bar
ne_gcfg.h.AC.ACFig.Visible = 'off';
cprog = ne_progress(0, 0, ...
    {true, 0, sprintf('Preparing data and designs for AN%sOVAs...', covc)});
drawnow;

% covariates
covs = repmat(covs, tlevels, 1);

% leave space for global mean
covgm = (ch.AddGlobalMean.Value > 0);
if covgm
    covs(:, end + 1) = NaN;
    covnames{end+1} = 'global signal';
end

% prepare design matrix for within-subjects factors (requiring RFX b0's)
if tlevels > 1

    % generate matrix
    X = zeros(numsubs * tlevels, numsubs + (tlevcs + size(covs, 2)) * numgroups);

% or without within-subject factors (no RFX b0's)
else
    X = zeros(numsubs, (1 + size(covs, 2)) * numgroups);
end

% one of many VMP-selection
spmsns = false;
vmpbbox = [];
if mtype == 'v'
    vmpbbox = glm.BoundingBox;
    if ~isempty(cc.vmp)
        cc.vmp = cc.vmp(ch.StoreInVMP.Value);
        if ~isempty(cc.vmp) && ...
           ~isempty(cc.vmp{1}) && ...
            isxff(cc.vmp{1}, 'vmp')
            vmpbbox = cc.vmp{1}.BoundingBox;
            if isfield(rtv, 'SubjectSPMsn') && ...
               ~isempty(fieldnames(rtv.SubjectSPMsn))
                spmsns = true;
            end
        else
            cc.vmp = [];
            if isfield(rtv, 'SubjectSPMsn') && ...
               ~isempty(fieldnames(rtv.SubjectSPMsn))
                spmsns = true;
                vmpbbox = struct('BBox', [44, 38, 44; 241, 193, 211], ...
                    'ResXYZ', min(2, glm.Resolution) .* ones(1, 3));
            end
        end
    end

    % create VMP
    mext = 'vmp';
    vmpres = vmpbbox.ResXYZ(1);
    newvmp = newnatresvmp(vmpbbox.BBox, vmpres);
    newvmp.Map = newvmp.Map(1);
    mfld = 'VMPData';
    mdim = 4;
    sref = {':', ':', ':', []};
    sasg = {[], ':', ':', ':'};
else

    % check selected SMP
    if ~isempty(cc.vmp)
        cc.vmp = cc.vmp(ch.StoreInVMP.Value);
        if isempty(cc.vmp) || ...
            isempty(cc.vmp{1}) || ...
           ~isxff(cc.vmp{1}, 'smp') || ...
            cc.vmp{1}.NrOfVertices ~= glm.NrOfVertices
            cc.vmp = {[]};
        end
    end

    % create SMP
    mext = 'smp';
    numvert = glm.NrOfVertices;
    newvmp = xff('new:smp');
    newvmp.NrOfVertices = numvert;
    szmapo = [numvert, 1];
    newvmp.Map = newvmp.Map(1);
    mfld = 'SMPData';
    newvmp.Map.(mfld) = single(zeros(szmapo));
    mdim = 2;
    sref = {':', []};
    sasg = {[], ':'};
end

% get data -> predictors
msize = size(newvmp.Map.(mfld));
mnumel = numel(newvmp.Map.(mfld));
ydata = zeros([size(X, 1), msize]);
sgd = glm.GLMData.Subject;
if spmsns
    nsm = size(sgd(1).BetaMaps, mdim);
end
if ch.CellType.Value == 1
    ctype = 'predictors';
    cnames = reshape(subpreds(flevcis(:)), size(flevcis));
    yi = 0;

    % predictors or FIRs
    if numst == 0
        for lc3 = 1:flevels(3)
            for lc2 = 1:flevels(2)
                for lc1 = 1:flevels(1)
                    for sc = 1:numsubs
                        yi = yi + 1;
                        if spmsns
                            cdatap = reshape(glm.SampleBVBox(vmpbbox, ...
                                (subjnums(sc) - 1) * nsm + flevcis(lc1, lc2, lc3), 'cubic'), ...
                                msize);
                        else
                            sref{end} = flevcis(lc1, lc2, lc3);
                            cdatap = ...
                                sgd(subjnums(sc)).BetaMaps(sref{:});
                        end
                        if smk > 0 && mdim == 4
                            if ~isinf(lbmin)
                                cdatap(cdatap < lbmin) = NaN;
                            end
                            if ~isinf(lbmax)
                                cdatap(cdatap > lbmax) = NaN;
                            end
                            sasg{1} = yi;
                            ydata(sasg{:}) = reshape(smoothdata3(cdatap, ...
                                (smk / vmpres) .* ones(1, 3)), [1, msize]);
                        else
                            sasg{1} = yi;
                            ydata(sasg{:}) = reshape(cdatap, [1, msize]);
                        end
                    end
                end
            end
        end

    % single-trial based factor
    else
    end

% contrasts (precludes FIR/ST models)
else
    ctype = 'contrasts';
    cnames = reshape(rtv.Contrasts(flevcis(:), 1), size(flevcis));
    yi = 0;
    for lc3 = 1:flevels(3)
        for lc2 = 1:flevels(2)
            for lc1 = 1:flevels(1)
                if ~spmsns
                    cdata = glm.RFX_conmaps(cons{flevcis(lc1, lc2, lc3)});
                    if smk > 0 && mdim == 4
                        if ~isinf(lbmin)
                            cdata(cdata < lbmin) = NaN;
                        end
                        if ~isinf(lbmax)
                            cdata(cdata > lbmin) = NaN;
                        end
                        cdata = smoothdata3(cdata, (smk / vmpres) .* ones(1, 3));
                    end
                end
                for sc = 1:numsubs
                    yi = yi + 1;
                    if spmsns
                        cdata = zeros(msize);
                        cvec = cons{flevcis(lc1, lc2, lc3)};
                        for cvic = 1:numel(cvec)
                            if cvec(cvic) ~= 0
                                cdatap = ...
                                    glm.SampleBVBox(vmpbbox, (subjnums(sc) - 1) * nsm + cvic);
                                if ~isinf(lbmin)
                                    cdatap(cdatap < lbmin) = NaN;
                                end
                                if ~isinf(lbmax)
                                    cdatap(cdatap > lbmin) = NaN;
                                end
                                cdata = cdata + cvec(cvic) .* cdatap;
                            end
                        end
                        if smk > 0
                            cdata = smoothdata3(cdata, (smk / vmpres) .* ones(1, 3));
                        end
                        sasg{1} = yi;
                        ydata(sasg{:}) = reshape(cdata, [1, msize]);
                    else
                        sref{end} = subjnums(sc);
                        sasg{1} = yi;
                        ydata(sasg{:}) = reshape(cdata(sref{:}), [1, msize]);
                    end
                end
            end
        end
    end
end
ydata = reshape(ydata, size(X, 1), mnumel);

% global covariate
if covgm
    for sc = 1:numsubs
        sref{end} = flevcis(:);
        smsk = reshape( ...
            ~any(isinf(sgd(subjnums(sc)).BetaMaps(sref{:})) | ...
            isnan(sgd(subjnums(sc)).BetaMaps(sref{:})), mdim) & ...
            any(sgd(subjnums(sc)).BetaMaps(sref{:}) ~= 0, mdim), 1, mnumel);
        for lc1 = sc:numsubs:size(covs, 1)
            mnin = meannoinfnan(ydata(lc1, smsk), 2);
            if isempty(mnin)
                mnin = 0;
            end
            covs(lc1, end) = mnin;
        end
    end
end

% main intercept
si = 1;
seps = 4 * sqrt(eps);
X(:, si) = 1;
% Xi = reshape(1:(numsubs*tlevels), [numsubs, flevels]);

% first generate subjects-only matrix
Xs = zeros(numsubs * tlevels, numsubs);
for sc = 1:numsubs
    Xs(sc:numsubs:end, sc) = 1;
end

% add group differences for groups 2:G
grpis = si + 1;
for gc = 2:numgroups
    si = si + 1;
    X(:, si) = sum(Xs(:, subgroups{gc}), 2) - sum(Xs(:, subgroups{gc-1}), 2);
end
grpis = grpis:si;

% for covariates
if ~isempty(covs)
    
    % first remove NaN's
    for cvc = 1:size(covs, 2)
        nancovs = isnan(covs(:, cvc));
        if coviact
            for gc = 1:numgroups
                grem = 1:numsubs;
                grem(subgroups{gc}) = [];
                nancovg = nancovs;
                nancovg(grem) = false;
                if any(nancovg(1:numsubs))
                    mnin = meannoinfnan(covs(subgroups{gc}, cvc));
                    if isempty(mnin)
                        mnin = 0;
                    end
                    covs(repmat(nancovg(1:numsubs, 1), tlevels, 1), cvc) = mnin;
                end
            end
        elseif any(nancovs)
            mnin = meannoinfnan(covs(:, cvc));
            if isempty(mnin)
                mnin = 0;
            end
            covs(nancovs, cvc) = mnin;
        end
    end
    
    % remove main effect from covariates (but no rescaling)
    if mrcov == 1
        b = calcbetas(X(:, 1), covs, 1);
        covs = covs - X(:, 1) * b';
        
    % remove main+group effect from covariates
    elseif mrcov == 2
        b = calcbetas(X(:, 1:si), covs, 1);
        covs = covs - X(:, 1:si) * b';
    end
end

% add covariates
covis = si + 1;
for sc = 1:size(covs, 2)
    si = si + 1;
    X(:, si) = covs(:, sc);
end

% and keep track of indices
covis = covis:si;
if isempty(covis)
    scovis = {};
else
    scovis = cell(1, numel(covis));
    for cvc = 1:numel(covis)
        scovis{cvc} = covis(cvc);
    end
end

% then add group-by-covariate interactions
sgrpcovis = scovis;
grpcovnames = covnames;
if coviact
    for cvc = 1:numel(scovis)
        sgrpcovis{cvc} = si + 1;
        for gc = 2:numgroups
            si = si + 1;
            X(:, si) = X(:, grpis(gc-1)) .* X(:, covis(cvc));
        end
        sgrpcovis{cvc} = sgrpcovis{cvc}:si;
        grpcovnames{cvc} = sprintf('group-by-%s', covnames{cvc});
    end
end

% all between factors (and interactions)
betweenis = [{grpis}, scovis, sgrpcovis];
betweennames = [{'group'}, covnames(:)', grpcovnames(:)'];
betempty = cellfun('isempty', betweenis);
betweenis(betempty) = [];
betweennames(betempty) = [];
abetweenis = cat(2, betweenis{:});

% then add subject-differences (within groups, also works for single G)
subis = si + 1;
if tlevels > 1
    rsubc = 1;
    rsubis = false(1, numsubs);
    for gc = 1:numgroups
        rsubis(rsubc) = true;
        si = si + 1;
        rsubc = rsubc + 1;
        X(subgroups{gc}(1):numsubs:end, si) = 1;
        for sc = 2:numel(subgroups{gc})
            si = si + 1;
            rsubc = rsubc + 1;
            X(subgroups{gc}(sc):numsubs:end, si) = 1;
        end
    end
end
subis = subis:si;

% within-levels
if tlevels > 1
    
    % generate factor matrices
    Xf = repmat({zeros([numsubs, tlevels])}, flevels);
    for xc = 1:numel(Xf)
        Xf{xc}(:, xc) = 1;
        Xf{xc} = Xf{xc}(:);
    end
    
    % along factor 1
    ws1is = si + 1;
    for fc = 2:flevels(1)
        si = si + 1;
        X(:, si) = max(cat(2, Xf{fc, :, :}), [], 2) - max(cat(2, Xf{fc-1, :, :}), [], 2);
    end
    ws1is = ws1is:si;

    % and for subject levels
    subws1is = si + 1;
    for fc = 2:flevels(1)
        for sc = subis
            si = si + 1;
            X(:, si) = X(:, sc) .* X(:, ws1is(fc-1));
        end
    end
    subws1is = subws1is:si;
    
    % along factor 2
    ws2is = si + 1;
    for fc = 2:flevels(2)
        si = si + 1;
        X(:, si) = max(cat(2, Xf{:, fc, :}), [], 2) - max(cat(2, Xf{:, fc-1, :}), [], 2);
    end
    ws2is = ws2is:si;
    subws2is = si + 1;
    for fc = 2:flevels(2)
        for sc = subis
            si = si + 1;
            X(:, si) = X(:, sc) .* X(:, ws2is(fc-1));
        end
    end
    subws2is = subws2is:si;
    
    % along factor 3
    ws3is = si + 1;
    for fc = 2:flevels(3)
        si = si + 1;
        X(:, si) = max(cat(2, Xf{:, :, fc}), [], 2) - max(cat(2, Xf{:, :, fc-1}), [], 2);
    end
    ws3is = ws3is:si;
    subws3is = si + 1;
    for fc = 2:flevels(3)
        for sc = subis
            si = si + 1;
            X(:, si) = X(:, sc) .* X(:, ws3is(fc-1));
        end
    end
    subws3is = subws3is:si;
    
    % within-only interactions
    ws1ws2is = si + 1;
    for fc2 = 2:flevels(2)
        for fc1 = 2:flevels(1)
            si = si + 1;
            X(:, si) = X(:, ws1is(fc1-1)) .* X(:, ws2is(fc2-1));
        end
    end
    ws1ws2is = ws1ws2is:si;

    % and subject levels
    subws1ws2is = si + 1;
    for fc1 = ws1ws2is
        for sc = subis
            si = si + 1;
            X(:, si) = X(:, sc) .* X(:, fc1);
        end
    end
    subws1ws2is = subws1ws2is:si;
    ws1ws3is = si + 1;
    for fc3 = 2:flevels(3)
        for fc1 = 2:flevels(1)
            si = si + 1;
            X(:, si) = X(:, ws1is(fc1-1)) .* X(:, ws3is(fc3-1));
        end
    end
    ws1ws3is = ws1ws3is:si;
    subws1ws3is = si + 1;
    for fc1 = ws1ws3is
        for sc = subis
            si = si + 1;
            X(:, si) = X(:, sc) .* X(:, fc1);
        end
    end
    subws1ws3is = subws1ws3is:si;
    ws2ws3is = si + 1;
    for fc3 = 2:flevels(3)
        for fc2 = 2:flevels(2)
            si = si + 1;
            X(:, si) = X(:, ws2is(fc2-1)) .* X(:, ws3is(fc3-1));
        end
    end
    ws2ws3is = ws2ws3is:si;
    subws2ws3is = si + 1;
    for fc1 = ws2ws3is
        for sc = subis
            si = si + 1;
            X(:, si) = X(:, sc) .* X(:, fc1);
        end
    end
    subws2ws3is = subws2ws3is:si;
    ws1ws2ws3is = si + 1;
    for fc3 = 2:flevels(3)
        for fc2 = 2:flevels(2)
            for fc1 = 2:flevels(1)
                si = si + 1;
                X(:, si) = X(:, ws1is(fc1-1)) .* X(:, ws2is(fc2-1)) .* X(:, ws3is(fc3-1));
            end
        end
    end
    ws1ws2ws3is = ws1ws2ws3is:si;

    % between-within interaction (for explained levels!)
    if ~isempty(ws1is)
        betws1is = betweenis;
    else
        betws1is = {};
    end
    for gc = 1:numel(betws1is)
        betws1is{gc} = si + 1;
        for bc = 1:numel(betweenis{gc})
            for fc = 1:numel(ws1is)
                si = si + 1;
                X(:, si) = X(:, betweenis{gc}(bc)) .* X(:, ws1is(fc));
            end
        end
        betws1is{gc} = betws1is{gc}:si;
        if ~isempty(betws1is{gc})
            betws1is{gc}(std(X(:, betws1is{gc})) <= seps) = [];
        end
    end
    if (~isempty(betws1is) && ...
        all(cellfun('isempty', betws1is))) || ...
        isempty(betws1is)
        betws1is = {[]};
        abetws1is = [];
    else
        abetws1is = cat(2, betws1is{:});
    end
    if ~isempty(ws2is)
        betws2is = betweenis;
    else
        betws2is = {};
    end
    for gc = 1:numel(betws2is)
        betws2is{gc} = si + 1;
        for bc = 1:numel(betweenis{gc})
            for fc = 1:numel(ws2is)
                si = si + 1;
                X(:, si) = X(:, betweenis{gc}(bc)) .* X(:, ws2is(fc));
            end
        end
        betws2is{gc} = betws2is{gc}:si;
        if ~isempty(betws2is{gc})
            betws2is{gc}(std(X(:, betws2is{gc})) <= seps) = [];
        end
    end
    if (~isempty(betws2is) && ...
        all(cellfun('isempty', betws2is))) || ...
        isempty(betws2is)
        betws2is = {[]};
        abetws2is = [];
    else
        abetws2is = cat(2, betws2is{:});
    end
    if ~isempty(ws3is)
        betws3is = betweenis;
    else
        betws3is = {};
    end
    for gc = 1:numel(betws3is)
        betws3is{gc} = si + 1;
        for bc = 1:numel(betweenis{gc})
            for fc = 1:numel(ws3is)
                si = si + 1;
                X(:, si) = X(:, betweenis{gc}(bc)) .* X(:, ws3is(fc));
            end
        end
        betws3is{gc} = betws3is{gc}:si;
        if ~isempty(betws3is{gc})
            betws3is{gc}(std(X(:, betws3is{gc})) <= seps) = [];
        end
    end
    if (~isempty(betws3is) && ...
        all(cellfun('isempty', betws3is))) || ...
        isempty(betws3is)
        betws3is = {[]};
        abetws3is = [];
    else
        abetws3is = cat(2, betws3is{:});
    end

    % between-and-two-within interactions
    if ~isempty(ws1ws2is)
        betws1ws2is = betweenis;
    else
        betws1ws2is = {};
    end
    for gc = 1:numel(betws1ws2is)
        betws1ws2is{gc} = si + 1;
        for bc = 1:numel(betweenis{gc})
            for fc = 1:numel(ws1ws2is)
                si = si + 1;
                X(:, si) = X(:, betweenis{gc}(bc)) .* X(:, ws1ws2is(fc));
            end
        end
        betws1ws2is{gc} = betws1ws2is{gc}:si;
        if ~isempty(betws1ws2is{gc})
            betws1ws2is{gc}(std(X(:, betws1ws2is{gc})) <= seps) = [];
        end
    end
    if (~isempty(betws1ws2is) && ...
        all(cellfun('isempty', betws1ws2is))) || ...
        isempty(betws1ws2is)
        betws1ws2is = {[]};
        abetws1ws2is = [];
    else
        abetws1ws2is = cat(2, betws1ws2is{:});
    end
    if ~isempty(ws1ws3is)
        betws1ws3is = betweenis;
    else
        betws1ws3is = {};
    end
    for gc = 1:numel(betws1ws3is)
        betws1ws3is{gc} = si + 1;
        for bc = 1:numel(betweenis{gc})
            for fc = 1:numel(ws1ws3is)
                si = si + 1;
                X(:, si) = X(:, betweenis{gc}(bc)) .* X(:, ws1ws3is(fc));
            end
        end
        betws1ws3is{gc} = betws1ws3is{gc}:si;
        if ~isempty(betws1ws3is{gc})
            betws1ws3is{gc}(std(X(:, betws1ws3is{gc})) <= seps) = [];
        end
    end
    if (~isempty(betws1ws3is) && ...
        all(cellfun('isempty', betws1ws3is))) || ...
        isempty(betws1ws3is)
        betws1ws3is = {[]};
        abetws1ws3is = [];
    else
        abetws1ws3is = cat(2, betws1ws3is{:});
    end
    if ~isempty(ws2ws3is)
        betws2ws3is = betweenis;
    else
        betws2ws3is = {};
    end
    for gc = 1:numel(betws2ws3is)
        betws2ws3is{gc} = si + 1;
        for bc = 1:numel(betweenis{gc})
            for fc = 1:numel(ws2ws3is)
                si = si + 1;
                X(:, si) = X(:, betweenis{gc}(bc)) .* X(:, ws2ws3is(fc));
            end
        end
        betws2ws3is{gc} = betws2ws3is{gc}:si;
        if ~isempty(betws2ws3is{gc})
            betws2ws3is{gc}(std(X(:, betws2ws3is{gc})) <= seps) = [];
        end
    end
    if (~isempty(betws2ws3is) && ...
        all(cellfun('isempty', betws2ws3is))) || ...
        isempty(betws2ws3is)
        betws2ws3is = {[]};
        abetws2ws3is = [];
    else
        abetws2ws3is = cat(2, betws2ws3is{:});
    end

    % between-and-three-within interaction
    if ~isempty(ws1ws2ws3is)
        betws1ws2ws3is = betweenis;
    else
        betws1ws2ws3is = {};
    end
    for gc = 1:numel(betws1ws2ws3is)
        betws1ws2ws3is{gc} = si + 1;
        for bc = 1:numel(betweenis{gc})
            for fc = 1:numel(ws1ws2ws3is)
                si = si + 1;
                X(:, si) = X(:, betweenis{gc}(bc)) .* X(:, ws1ws2ws3is(fc));
            end
        end
        betws1ws2ws3is{gc} = betws1ws2ws3is{gc}:si;
        if ~isempty(betws1ws2ws3is{gc})
            betws1ws2ws3is{gc}(std(X(:, betws1ws2ws3is{gc})) <= seps) = [];
        end
    end
    if (~isempty(betws1ws2ws3is) && ...
        all(cellfun('isempty', betws1ws2ws3is))) || ...
        isempty(betws1ws2ws3is)
        betws1ws2ws3is = {[]};
        abetws1ws2ws3is = [];
    else
        abetws1ws2ws3is = cat(2, betws1ws2ws3is{:});
    end
else
    ws1is = [];
    ws2is = [];
    ws3is = [];
    ws1ws2is = [];
    ws1ws3is = [];
    ws2ws3is = [];
    ws1ws2ws3is = [];
    subws1is = [];
    subws2is = [];
    subws3is = [];
    subws1ws2is = [];
    subws1ws3is = [];
    subws2ws3is = [];
    betws1is = {[]};
    betws2is = {[]};
    betws3is = {[]};
    betws1ws2is = {[]};
    betws1ws3is = {[]};
    betws2ws3is = {[]};
    betws1ws2ws3is = {[]};
    abetws1is = [];
    abetws2is = [];
    abetws3is = [];
    abetws1ws2is = [];
    abetws1ws3is = [];
    abetws2ws3is = [];
    abetws1ws2ws3is = [];
end

% remove bad rows
gXrows = ~any(isnan(X), 2) & ~all(isnan(ydata), 2);
sX1 = sum(gXrows);
sX1sq = sqrt(sqrt(sX1));

% reduced DF
reddf = lsqueeze(sum(ydata == 0, 1));
df1os = ones(size(reddf));

% remove bad columns
bXcols = (sum(abs(X(gXrows, :)), 1) < seps);

% adjust column indices
ws1is(bXcols(ws1is)) = [];
ws2is(bXcols(ws2is)) = [];
ws3is(bXcols(ws3is)) = [];
ws1ws2is(bXcols(ws1ws2is)) = [];
ws1ws3is(bXcols(ws1ws3is)) = [];
ws2ws3is(bXcols(ws2ws3is)) = [];
ws1ws2ws3is(bXcols(ws1ws2ws3is)) = [];
subis(bXcols(subis)) = [];
subws1is(bXcols(subws1is)) = [];
subws2is(bXcols(subws2is)) = [];
subws3is(bXcols(subws3is)) = [];
subws1ws2is(bXcols(subws1ws2is)) = [];
subws1ws3is(bXcols(subws1ws3is)) = [];
subws2ws3is(bXcols(subws2ws3is)) = [];
abetweenis(bXcols(abetweenis)) = [];
abetws1is(bXcols(abetws1is)) = [];
abetws2is(bXcols(abetws2is)) = [];
abetws3is(bXcols(abetws3is)) = [];
abetws1ws2is(bXcols(abetws1ws2is)) = [];
abetws1ws3is(bXcols(abetws1ws3is)) = [];
abetws2ws3is(bXcols(abetws2ws3is)) = [];

% keep track of cost
compscost = 0;

% which models to run?
wm = struct;

% between effects (error: subject)
if ~isempty(abetweenis)
    betids = betweennames;
    for gc = 1:numel(betweennames)
        betids{gc} = sprintf('bet%d', gc);
        wm.(betids{gc}) = {[1, abetweenis], [1, setdiff(abetweenis, betweenis{gc})], ...
            sprintf('main effect of %s', betweennames{gc}), 0, [1, 1, 1]};
        wm.(betids{gc}){4} = sqrt(numel(wm.(betids{gc}){1}) * numel(wm.(betids{gc}){2})) + sX1sq;
        compscost = compscost + wm.(betids{gc}){4};
    end
else
    betids = cell(1, 0);
end

% W1 effect (error: subject-x-W1)
if ~isempty(ws1is)
    wm.ws1 = {[subis, abetws1is, ws1is], [subis, abetws1is], ...
        'main effect: within-factor 1', 0, [0, 1, 1]};
    wm.ws1{4} = sqrt(numel(wm.ws1{1}) * numel(wm.ws1{2})) + sX1sq;
    compscost = compscost + wm.ws1{4};
end
if ~isempty(abetws1is)
    betws1ids = betweennames;
    for gc = numel(betweennames):-1:1
        if isempty(betws1is{gc})
            betws1ids(gc) = [];
            continue;
        end
        betws1ids{gc} = sprintf('betws1%d', gc);
        wm.(betws1ids{gc}) = {[subis, ws1is, abetws1is], [subis, ws1is, setdiff(abetws1is, betws1is{gc})], ...
            sprintf('interaction: %s-by-within-factor 1', betweennames{gc}), 0, [0, 1, 1]};
        wm.(betws1ids{gc}){4} = sqrt(numel(wm.(betws1ids{gc}){1}) * numel(wm.(betws1ids{gc}){2})) + sX1sq;
        compscost = compscost + wm.(betws1ids{gc}){4};
    end
else
    betws1ids = cell(1, 0);
end

% W2 effect (error: subject-x-W2)
if ~isempty(ws2is)
    wm.ws2 = {[subis, abetws2is, ws2is], [subis, abetws2is], ...
        'main effect: within-factor 2', 0, [1, 0, 1]};
    wm.ws2{4} = sqrt(numel(wm.ws2{1}) * numel(wm.ws2{2})) + sX1sq;
    compscost = compscost + wm.ws2{4};
end
if ~isempty(abetws2is)
    betws2ids = betweennames;
    for gc = numel(betweennames):-1:1
        if isempty(betws2is{gc})
            betws2ids(gc) = [];
            continue;
        end
        betws2ids{gc} = sprintf('betws2%d', gc);
        wm.(betws2ids{gc}) = {[subis, ws2is, abetws2is], [subis, ws2is, setdiff(abetws2is, betws2is{gc})], ...
            sprintf('interaction: %s-by-within-factor 2', betweennames{gc}), 0, [1, 0, 1]};
        wm.(betws2ids{gc}){4} = sqrt(numel(wm.(betws2ids{gc}){1}) * numel(wm.(betws2ids{gc}){2})) + sX1sq;
        compscost = compscost + wm.(betws2ids{gc}){4};
    end
else
    betws2ids = cell(1, 0);
end

% W3 effect
if ~isempty(ws3is)
    wm.ws3 = {[subis, abetws3is, ws3is], [subis, abetws3is], ...
        'main effect: within-factor 3', 0, [1, 1, 0]};
    wm.ws3{4} = sqrt(numel(wm.ws3{1}) * numel(wm.ws3{2})) + sX1sq;
    compscost = compscost + wm.ws3{4};
end
if ~isempty(abetws3is)
    betws3ids = betweennames;
    for gc = numel(betweennames):-1:1
        if isempty(betws3is{gc})
            betws3ids(gc) = [];
            continue;
        end
        betws3ids{gc} = sprintf('betws3%d', gc);
        wm.(betws3ids{gc}) = {[subis, ws3is, abetws3is], [subis, ws3is, setdiff(abetws3is, betws3is{gc})], ...
            sprintf('interaction: %s-by-within-factor 3', betweennames{gc}), 0, [1, 1, 0]};
        wm.(betws3ids{gc}){4} = sqrt(numel(wm.(betws3ids{gc}){1}) * numel(wm.(betws3ids{gc}){2})) + sX1sq;
        compscost = compscost + wm.(betws3ids{gc}){4};
    end
else
    betws3ids = cell(1, 0);
end

% W1xW2 effect (error: subject-x-W1-x-W2)
if ~isempty(ws1ws2is)
    wm.ws1ws2 = { ...
        [subis, subws1is, subws2is, abetws1ws2is, ws1ws2is], ...
        [subis, subws1is, subws2is, abetws1ws2is], ...
        'interaction: within-factors 1 and 2', 0, [0, 0, 1]};
    wm.ws1ws2{4} = sqrt(numel(wm.ws1ws2{1}) * numel(wm.ws1ws2{2})) + sX1sq;
    compscost = compscost + wm.ws1ws2{4};
end
if ~isempty(abetws1ws2is)
    betws1ws2ids = betweennames;
    for gc = numel(betweennames):-1:1
        if isempty(betws1ws2is{gc})
            betws1ws2ids(gc) = [];
            continue;
        end
        betws1ws2ids{gc} = sprintf('betws1ws2%d', gc);
        wm.(betws1ws2ids{gc}) = { ...
            [subis, subws1is, subws2is, ws1ws2is, abetws1ws2is], ...
            [subis, subws1is, subws2is, ws1ws2is, setdiff(abetws1ws2is, betws1ws2is{gc})], ...
            sprintf('interaction: %s-by-within-factors 1+2', betweennames{gc}), 0, [0, 0, 1]};
        wm.(betws1ws2ids{gc}){4} = sqrt(numel(wm.(betws1ws2ids{gc}){1}) * numel(wm.(betws1ws2ids{gc}){2})) + sX1sq;
        compscost = compscost + wm.(betws1ws2ids{gc}){4};
    end
else
    betws1ws2ids = cell(1, 0);
end

% W1xW3 effect (error: subject-x-W1-x-W3)
if ~isempty(ws1ws3is)
    wm.ws1ws3 = { ...
        [subis, subws1is, subws3is, abetws1ws3is, ws1ws3is], ...
        [subis, subws1is, subws3is, abetws1ws3is], ...
        'interaction: within-factors 1 and 3', 0, [0, 1, 0]};
    wm.ws1ws3{4} = sqrt(numel(wm.ws1ws3{1}) * numel(wm.ws1ws3{2})) + sX1sq;
    compscost = compscost + wm.ws1ws3{4};
end
if ~isempty(abetws1ws3is)
    betws1ws3ids = betweennames;
    for gc = numel(betweennames):-1:1
        if isempty(betws1ws3is{gc})
            betws1ws3ids(gc) = [];
            continue;
        end
        betws1ws3ids{gc} = sprintf('betws1ws3%d', gc);
        wm.(betws1ws3ids{gc}) = { ...
            [subis, subws1is, subws3is, ws1ws3is, abetws1ws3is], ...
            [subis, subws1is, subws3is, ws1ws3is, setdiff(abetws1ws3is, betws1ws3is{gc})], ...
            sprintf('interaction: %s-by-within-factors 1+3', betweennames{gc}), 0, [0, 1, 0]};
        wm.(betws1ws3ids{gc}){4} = sqrt(numel(wm.(betws1ws3ids{gc}){1}) * numel(wm.(betws1ws3ids{gc}){2})) + sX1sq;
        compscost = compscost + wm.(betws1ws3ids{gc}){4};
    end
else
    betws1ws3ids = cell(1, 0);
end

% W2xW3 effect (error: subject-x-W2-x-W3)
if ~isempty(ws2ws3is)
    wm.ws2ws3 = { ...
        [subis, subws2is, subws3is, abetws2ws3is, ws2ws3is], ...
        [subis, subws2is, subws3is, abetws2ws3is], ...
        'interaction: within-factors 2 and 3', 0, [1, 0, 0]};
    wm.ws2ws3{4} = sqrt(numel(wm.ws2ws3{1}) * numel(wm.ws2ws3{2})) + sX1sq;
    compscost = compscost + wm.ws2ws3{4};
end
if ~isempty(abetws2ws3is)
    betws2ws3ids = betweennames;
    for gc = numel(betweennames):-1:1
        if isempty(betws2ws3is{gc})
            betws2ws3ids(gc) = [];
            continue;
        end
        betws2ws3ids{gc} = sprintf('betws2ws3%d', gc);
        wm.(betws2ws3ids{gc}) = { ...
            [subis, subws2is, subws3is, ws2ws3is, abetws2ws3is], ...
            [subis, subws2is, subws3is, ws2ws3is, setdiff(abetws2ws3is, betws2ws3is{gc})], ...
            sprintf('interaction: %s-by-within-factors 2+3', betweennames{gc}), 0, [1, 0, 0]};
        wm.(betws2ws3ids{gc}){4} = sqrt(numel(wm.(betws2ws3ids{gc}){1}) * numel(wm.(betws2ws3ids{gc}){2})) + sX1sq;
        compscost = compscost + wm.(betws2ws3ids{gc}){4};
    end
else
    betws2ws3ids = cell(1, 0);
end

% adapt model
xmod = [subis, subws1is, subws2is, subws3is, ...
    subws1ws2is, subws1ws3is, subws2ws3is];

% W1xW2xW3 effect (error: subject-x-W1-x-W2-x-W3)
if ~isempty(ws1ws2ws3is)
    wm.ws1ws2ws3 = {[xmod, abetws1ws2ws3is, ws1ws2ws3is], [xmod, abetws1ws2ws3is], ...
        'interaction: within-factors 1, 2, and 3', 0, [0, 0, 0]};
    wm.ws1ws2ws3{4} = sqrt(numel(wm.ws1ws2ws3{1}) * numel(wm.ws1ws2ws3{2})) + sX1sq;
    compscost = compscost + wm.ws1ws2ws3{4};
end
if ~isempty(abetws1ws2ws3is)
    betws1ws2ws3ids = betweennames;
    for gc = numel(betweennames):-1:1
        if isempty(betws1ws2ws3is{gc})
            betws1ws2ws3ids(gc) = [];
            continue;
        end
        betws1ws2ws3ids{gc} = sprintf('betws1ws2ws3%d', gc);
        wm.(betws1ws2ws3ids{gc}) = { ...
            [xmod, ws1ws2ws3is, abetws1ws2ws3is], ...
            [xmod, ws1ws2ws3is, setdiff(abetws1ws2ws3is, betws1ws2ws3is{gc})], ...
            sprintf('interaction: %s-by-within-factors 1+2+3', betweennames{gc}), 0, [0, 0, 0]};
        wm.(betws1ws2ws3ids{gc}){4} = ...
            sqrt(numel(wm.(betws1ws2ws3ids{gc}){1}) * numel(wm.(betws1ws2ws3ids{gc}){2})) + sX1sq;
        compscost = compscost + wm.(betws1ws2ws3ids{gc}){4};
    end
else
    betws1ws2ws3ids = cell(1, 0);
end

% update VMP
newvmp.Map.Type = 4;
newvmp.Map.RunTimeVars.SourceGLM = glm.FilenameOnDisk;
newvmp.Map.RunTimeVars.SourceGLMID = glm.RunTimeVars.xffID;
newvmp.Map.RunTimeVars.RFXGLM = true;
newvmp.Map.RunTimeVars.ANCOVA = struct( ...
    'BetweenLevels', numel(groupi), ...
    'CollapseLevels', [0, 0, 0], ...
    'Covariates', {covnames}, ...
    'GlobalMeanAdded', covgm, ...
    'Groups', {groups}, ...
    'Models', {{[], []}}, ...
    'Subjects', {subjects}, ...
    'SubSel', sort(subjnums(:)), ...
    'WithinLevels', flevels, ...
    'WithinLevelNames', {cnames}, ...
    'WithinType', ctype);
mi = 1;

% begin regressions
cost = 0;
for mfc = [betids, {'ws1'}, betws1ids, {'ws2'}, betws2ids, {'ws3'}, betws3ids, ...
        {'ws1ws2'}, betws1ws2ids, {'ws1ws3'}, betws1ws3ids, ...
        {'ws2ws3'}, betws2ws3ids, {'ws1ws2ws3'}, betws1ws2ws3ids]
    if isfield(wm, mfc{1}) && ...
       ~isequal(wm.(mfc{1}){1}, wm.(mfc{1}){2})
        [newvmp, mi, cost] = sfancova(newvmp, mi, cp, cost, wm.(mfc{1}), ...
            compscost, X, gXrows, ydata, df1os, reddf, covc, msize, flevels, srf);
    end
end

% new VMP requested/necessary
if isempty(cc.vmp) || ...
    isempty(cc.vmp{1}) || ...
   ~isxff(cc.vmp{1}, mext)

    % set up correctly
    cc.vmp = newvmp;
    tmapi = 1;
    created = true;

% existing VMP
else
    cc.vmp = cc.vmp{1};
    tmapi = numel(cc.vmp.Map) + 1;
    cc.vmp.Map = catstruct(cc.vmp.Map(:), newvmp.Map(:))';
    newvmp.ClearObject;
    created = false;
end

% make sure the colors in new maps are set
cc.vmp.SetColors(tmapi:numel(cc.vmp.Map), 'xauto');

% set SourceGLM handle in VMP object
cc.vmp.SetHandle('SourceGLM', glm);

% and set RTV saving to auto
cc.vmp.RunTimeVars.AutoSave = true;

% show/update correct VMP
ne_openfile(0, 0, cc.vmp, created);

% and then show the first of the newly created maps
if mext(1) == 'v'
    ne_gcfg.fcfg.StatsVarIdx = tmapi;
    ne_gcfg.h.StatsVarMaps.Value = tmapi;
    ne_setcstatmap;
else
    ne_gcfg.fcfg.SurfStatsVarIdx = tmapi;
    ne_gcfg.h.SurfStatsVarMaps.Value = tmapi;
    ne_setcsrfstatmap;
end

% finally bring up UI again
ne_progress(0, 0, cprog);
ne_gcfg.h.AC.ACFig.Visible = 'on';



% sub-function to compute a single model
function [v, mi, cost] = sfancova(v, mi, p, cost, task, cc, X, gXrows, ydata, df1os, reddf, covc, msize, flevels, srf)

% global var for error
global ne_gcfg;

% keep track of a memory error
persistent memok;
if isempty(memok)
    memok = true;
end

% models
mf = task{1};
mr = task{2};
mname = task{3};
clevels = task{5};

% task cost
tc = task{4};
taskcost = tc / cc;
p.Progress(cost, ['Computing ' mname '...']);

% collapsing
if any(clevels > 0)
    X = reshape(X, [round(size(X, 1) / prod(flevels)), flevels, size(X, 2)]);
    Xns = size(X);
    gXrows = reshape(gXrows, [Xns(1:end-1), 1]);
    ydata = reshape(ydata, [Xns(1:end-1), size(ydata, 2)]);
    if clevels(1) > 0
        X = mean(X, 2);
        gXrows = all(gXrows, 2);
        ydata = mean(ydata, 2);
        flevels(1) = 1;
    end
    if clevels(2) > 0
        X = mean(X, 3);
        gXrows = all(gXrows, 3);
        ydata = mean(ydata, 3);
        flevels(2) = 1;
    end
    if clevels(3) > 0
        X = mean(X, 4);
        gXrows = all(gXrows, 4);
        ydata = mean(ydata, 4);
        flevels(3) = 1;
    end
    X = reshape(X, [size(X, 1) * prod(flevels), Xns(end)]);
    gXrows = reshape(gXrows, size(X, 1), 1);
    ydata = reshape(ydata, [size(X, 1), size(ydata, ndims(ydata))]);
end

% remove columns with all 0
mf(all(X(gXrows, mf) == 0)) = [];
mr(all(X(gXrows, mr) == 0)) = [];

% if memory condition is good
if memok
    try
        [f, df1, df2, bfull, ptc] = modelcomp(X(gXrows, mf), X(gXrows, mr), ydata(gXrows, :), 1);
        if numel(msize) == 3
            [fwhm, fwhmi] = resestsmooth(reshape(ydata(gXrows, :) - ptc, [sum(gXrows), msize]), ...
                v.Resolution, struct('tdim', 1));
        else
            [fwhm, fwhmi] = resestsmoothsrf(reshape(ydata(gXrows, :) - ptc, [sum(gXrows), prod(msize)]), ...
                srf, struct('tdim', 1));
        end
    catch ne_eo;
        memok = false;
        ne_gcfg.c.lasterr = ne_eo;
    end
end

% failure
if ~memok
    if numel(msize) < 3
        msize(3) = 1;
    end
    ydata = reshape(ydata, [size(X, 1), msize]);
    f = zeros(msize);
    res = zeros([sum(gXrows), msize(1:2), 3]);
    fwhmi = zeros(msize);
    fx = zeros(msize);
    fy = zeros(msize);
    fz = zeros(msize);
    for slc = 1:msize(3)
        if slc > 3
            res(:, :, :, 1:2) = res(:, :, :, 2:3);
        end
        [fp, df1, df2, bfull, ptc] = modelcomp(X(gXrows, mf), X(gXrows, mr), ydata(gXrows, :, :, slc), 1);
        f(:, :, slc) = fp;
        if slc < 3
            res(:, :, :, slc) = ydata(gXrows, :, :, slc) - ptc;
        else
            res(:, :, :, 3) = ydata(gXrows, :, :, slc) - ptc;
            [fwhm, fip, fxp, fyp, fzp] = resestsmooth(res, ...
                v.Resolution, struct('tdim', 1));
            fwhmi(:, :, slc-1) = fip(:, :, 2);
            fx(:, :, slc-1) = fxp(:,  :, 2);
            fy(:, :, slc-1) = fyp(:, :, 2);
            fz(:, :, slc-1) = fzp(:, :, 2);
        end
        p.Progress(cost + taskcost * slc / msize(3));
    end
    if msize(3) > 1
        fwhm = [median(fx(fx > 0)), median(fy(fy > 0)), median(fz(fz > 0))];
    end
end

% make sure stats are valid and correct for missing DF
f(isinf(f) | isnan(f)) = 0;
f = sdist('finv', sdist('fcdf', ...
    f(:), df1 .* df1os, max(0.001, df2 - reddf), true), df1, df2, true);
f(isinf(f) | isnan(f)) = 0;

% set into VMP
v.Map(mi) = v.Map(1);
v.Map(mi).Name = ['AN' covc 'OVA: ' mname];
v.Map(mi).DF1 = df1;
v.Map(mi).DF2 = df2;
v.Map(mi).LowerThreshold = sdist('finv', 0.995, df1, df2);
v.Map(mi).UpperThreshold = sdist('finv', 0.9999, df1, df2);
if numel(msize) > 2 && msize(3) > 1
    v.Map(mi).VMPData = reshape(single(f), msize);
else
    v.Map(mi).SMPData = single(f(:));
end
v.Map(mi).RunTimeVars.FWHMResEst = fwhm;
v.Map(mi).RunTimeVars.FWHMResImg = fwhmi;
v.Map(mi).RunTimeVars.ANCOVA.CollapseLevels = clevels;
v.Map(mi).RunTimeVars.ANCOVA.Models{1} = sparse(X(:, mr));
v.Map(mi).RunTimeVars.ANCOVA.Models{2} = sparse(X(:, setdiff(mf, mr)));

% update cost
cost = cost + taskcost;
p.Progress(cost);
mi = mi + 1;
