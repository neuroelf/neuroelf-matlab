function [t, b, ix, p, r, w, tpreds, d] = mdm_ComputeVOIGLM(xo, voi, opts)
% MDM::ComputeGLM  - compute a VOI-based GLM from an MDM and VOI file
%
% FORMAT:       [t, b, ix, p, r, w, pn, d] = mdm.ComputeVOIGLM(voi, [options])
%
% Input fields:
%
%       options     optional 1x1 struct with fields
%        .combsub   flag, combine each subject's data (default: false)
%        .confounds flag, also return confound betas (default: false)
%        .motpars   motion parameters (Sx1 cell array with sdm/txt files)
%        .motparsd  also add diff of motion parameters (default: false)
%        .motparsq  also add squared motion parameters (default: false)
%        .ndcreg    if set > 0, perform deconvolution (only with PRTs!)
%        .orthconf  orthogonalize confounds (and motion parameters, true)
%        .ppicond   list of regressors (or differences) to interact
%        .ppirob    perform robust regression on VOI timecourse and remove
%                   outliers from timecourse/model (threshold, default: 0)
%        .ppivoi    VOI used to extract time-course from
%        .prtpnorm  normalize parameters of PRT.Conds (true)
%        .regdiff   flag, regress first discreet derivatives (diff) instead
%        .restcond  remove rest condition (rest cond. name, default: '')
%        .robust    perform robust instead of OLS regression
%        .savesdms  token, if not empty, save on-the-fly SDMs (e.g. '.sdm')
%        .showsdms  token, passed to SDM::ShowDesign (if valid)
%        .subsel    cell array with subject IDs to work on
%        .subvois   subject specific VOIs, either 'sub_', '_sub', {'voi'}
%        .tfilter   add filter regressors to SDMs (cut-off in secs)
%        .tfilttype temporal filter type (either 'dct' or {'fourier'})
%        .unique    flag, only extract unique functional voxels (false)
%        .voisel    cell array with sub-VOI selection to use
%        .vweight   combine runs/studies variance-weighted (default: false)
%
% Output fields:
%
%       t           time courses (copy of MDM::VOITimeCourses)
%       b           SxP beta weights, either S = studies or subjects
%       ix          inverse design matrices (can be used for weighting)
%       p           predicted time-courses
%       r           residual time-courses
%       w           time-courses weights
%       pn          predictor names
%       d           design matrices
%
% Using: emptystruct, findfirst, makelabel, multimatch.

% Version:  v1.1
% Build:    16020116
% Date:     Feb-01 2016, 4:53 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2012, 2014, 2015, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'mdm') || ...
    numel(voi) ~= 1 || ~xffisobject(voi, true, 'voi')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end

% try to extract time courses
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
try
    [t, tf, tvi, tr] = mdm_VOITimeCourses(xo, voi, opts);
catch xfferror
    rethrow(xfferror);
end
bc = xo.C;
rtv = bc.RunTimeVars;

% check options
if ~isfield(opts, 'ar1') || ~islogical(opts.ar1) || isempty(opts.ar1)
    opts.ar1 = false;
else
    opts.ar1 = opts.ar1(1);
end
if ~isfield(opts, 'combsub') || ~islogical(opts.combsub) || numel(opts.combsub) ~= 1
    opts.combsub = false;
end
if ~isfield(opts, 'confounds') || ~islogical(opts.confounds) || numel(opts.confounds) ~= 1
    opts.confounds = false;
end
if isfield(opts, 'motpars') && islogical(opts.motpars) && numel(opts.motpars) == 1 && ...
    opts.motpars && isfield(rtv, 'MotionParameters') && iscell(rtv.MotionParameters) && ...
    numel(rtv.MotionParameters) == size(bc.XTC_RTC, 1)
    opts.motpars = rtv.MotionParameters;
end
if ~isfield(opts, 'motpars') || ~iscell(opts.motpars)
    opts.motpars = {};
end
if ~isfield(opts, 'motparsd') || ~islogical(opts.motparsd) || numel(opts.motparsd) ~= 1
    opts.motparsd = false;
end
if ~isfield(opts, 'motparsq') || ~islogical(opts.motparsq) || numel(opts.motparsq) ~= 1
    opts.motparsq = false;
end
if ~isfield(opts, 'ndcreg') || ~isa(opts.ndcreg, 'double') || numel(opts.ndcreg) ~= 1 || ...
    isinf(opts.ndcreg) || isnan(opts.ndcreg) || opts.ndcreg < 0 || opts.ndcreg > 20
    opts.ndcreg = 0;
else
    opts.ndcreg = floor(opts.ndcreg);
end
if ~isfield(opts, 'orthconf') || ~islogical(opts.orthconf) || numel(opts.orthconf) ~= 1
    opts.orthconf = true;
end
if ~isfield(opts, 'ppicond') || (~iscell(opts.ppicond)) && (~ischar(opts.ppicond) || isempty(opts.ppicond))
    opts.ppicond = {};
    opts.ppirob = false;
    opts.ppivoi = [];
elseif ischar(opts.ppicond)
    opts.ppicond = {opts.ppicond(:)'};
else
    for pcc = numel(opts.ppicond):-1:1
        if ~ischar(opts.ppicond{pcc}) || isempty(opts.ppicond{pcc})
            error('neuroelf:xff:badArgument', 'Invalid ppicond list.');
        end
    end
end
if ~isfield(opts, 'ppirob') || ~islogical(opts.ppirob) || numel(opts.ppirob) ~= 1
    opts.ppirob = false;
end
if ~isfield(opts, 'ppivoi') || numel(opts.ppivoi) ~= 1 || ...
   (~xffisobject(opts.ppivoi, true, 'poi') && ~xffisobject(opts.ppivoi, true, 'voi'))
    opts.ppivoi = [];
end
if ~isfield(opts, 'prtpnorm') || ~islogical(opts.prtpnorm) || numel(opts.prtpnorm) ~= 1
    opts.prtpnorm = true;
end
if ~isfield(opts, 'regdiff') || ~islogical(opts.regdiff) || numel(opts.regdiff) ~= 1
    opts.regdiff = false;
end
if ~isfield(opts, 'restcond') || (~ischar(opts.restcond) && ~iscell(opts.restcond)) || ...
    isempty(opts.restcond)
    opts.restcond = {};
elseif ischar(opts.restcond)
    opts.restcond = {lower(opts.restcond(:)')};
else
    for rcc = numel(opts.restcond):-1:1
        if ~ischar(opts.restcond{rcc}) || isempty(opts.restcond{rcc})
            opts.restcond(rcc) = [];
        else
            opts.restcond{rcc} = lower(opts.restcond{rcc}(:)');
        end
    end
end
if ~isfield(opts, 'robust') || numel(opts.robust) ~= 1 || ~islogical(opts.robust)
    opts.robust = false;
end
if ~isfield(opts, 'savesdms') || ~ischar(opts.savesdms) || numel(opts.savesdms) < 4
    opts.savesdms = '';
else
    opts.savesdms = opts.savesdms(:)';
    if ~strcmpi(opts.savesdms(end-3:end), '.sdm')
        error('neuroelf:xff:badArgument', '.savesdms token must end in ''.sdm''.');
    end
end
if ~isfield(opts, 'seppred') || ~isa(opts.seppred, 'double') || numel(opts.seppred) ~= 1 || ...
   ~any(0:2 == opts.seppred)
    if numel(bc.SeparatePredictors) == 1 && isa(bc.SeparatePredictors, 'double') && any(0:2 == bc.SeparatePredictors)
        opts.seppred = bc.SeparatePredictors;
    else
        opts.seppred = 0;
    end
end
if ~isfield(opts, 'showsdms') || ~ischar(opts.showsdms) || isempty(opts.showsdms) || ...
   ~any(lower(opts.showsdms(1)) == 'iop')
    opts.showsdms = '';
else
    opts.showsdms = lower(opts.showsdms(1));
end
if ~isfield(opts, 'subsel') || ~iscell(opts.subsel) || isempty(opts.subsel)
    opts.subsel = mdm_Subjects(xo);
else
    opts.subsel = opts.subsel(:);
    try
        ssm = ne_methods.multimatch(opts.subsel, mdm_Subjects(xo));
    catch xfferror
        rethrow(xfferror);
    end
    if any(ssm < 1)
        error('neuroelf:xff:badArgument', 'Invalid subject ID in selection.');
    end
end
if ~isfield(opts, 'tfilter') || ~isa(opts.tfilter, 'double') || numel(opts.tfilter) ~= 1 || ...
    isnan(opts.tfilter) || opts.tfilter < 60
    opts.tfilter = Inf;
end
if ~isfield(opts, 'tfilttype') || ~ischar(opts.tfilttype) || isempty(opts.tfilttype) || ...
    lower(opts.tfilttype(1)) ~= 'd'
    opts.tfilttype = 'fourier';
else
    opts.tfilttype = 'dct';
end
if ~isfield(opts, 'trans') || ~ischar(opts.trans) || isempty(opts.trans) || ...
   ~any('npz' == lower(opts.trans(1)))
    opts.trans = 'n';
    if bc.PSCTransformation > 0
        opts.trans = 'p';
    end
    if bc.zTransformation > 0
        opts.trans = 'z';
    end
end
if ~isfield(opts, 'vweight') || ~islogical(opts.vweight) || numel(opts.vweight) ~= 1
    opts.vweight = false;
end
if ~isfield(opts, 'xconfound') || ~iscell(opts.xconfound)
    opts.xconfound = {};
end
subjids = opts.subsel;
if isempty(subjids)
    error('neuroelf:xff:badArgument', 'Invalid or missing subject IDs supplied.');
end
if any(strcmp(subjids, ''))
    error('neuroelf:xff:invalidObject', 'Invalid subject IDs for some subjects.');
end

% let MDM::SDMs do most of the work (for the modelling part)
rfiles = bc.XTC_RTC;
rfobjs = cell(size(rfiles));
mfiles = rfiles(:, end);
for stc = 1:numel(mfiles)
    if isempty(mfiles{stc}) || (~ischar(mfiles{stc}) && ~xffisobject(mfiles{stc}, true, {'prt', 'sdm'}))
        error('neuroelf:xff:invalidModels', 'Invalid model references.');
    end
    if xffisobject(mfiles{stc}, true)
        mfiles{stc} = mfiles{stc}.F;
    end
end
try
    mdmsubst = mdm_Subjects(xo, true);
    [rfobjs(:, end), stlist, sdmtr, rfobjs(:, end-1)] = mdm_SDMs(xo, opts);
catch xfferror
    rethrow(xfferror);
end

% get required PRT/SDM objects
rlines = ne_methods.multimatch(tf, bc.XTC_RTC(:, 1));
rfiles = bc.XTC_RTC(rlines, :);
mdmsubst = mdmsubst(rlines);
rfobjo = rfobjs;
rfobjs = rfobjs(rlines, :);
sdmtr = sdmtr(rlines);
if any(tr(:) ~= sdmtr(:))
    clearxffobjects(rfobjo(:));
    error('neuroelf:xff:processingError', 'Invalid TR combination between methods.');
end
numstudy = size(rfiles, 1);

% figure handle for showsdms
showsdmf = [];

% check PRT/SDMs
numtp = 0;
makelabel = ne_methods.makelabel;
try
    predcol = struct;
    preds = cell(numstudy, 1);
    predr = cell(numstudy, 1);
    ststr = ne_methods.emptystruct({ ...
        'NrOfTimePoints', 'NameOfAnalyzedFile', 'NameOfSSMFile', 'NameOfSDMFile', ...
        'NrOfConfounds', 'SubjectID'}, [1, numstudy]);

    % check studies that make it into the dataset
    for stc = 1:numstudy

        % check objects
        xtcsc = t{stc};
        sdmsc = rfobjs{stc, end};
        stnumtp = size(xtcsc, 1);
        if stnumtp ~= sdmsc.C.NrOfDataPoints
            error('neuroelf:xff:layoutMismatch', ...
                'Number of time points between SDM and VOI-timecourse must match.');
        end
        sdmm = sdmsc.C.SDMMatrix;
        ststr(stc).NameOfAnalyzedFile = tf{stc};
        ststr(stc).NameOfSDMFile = sdmsc.F;
        ststr(stc).NrOfConfounds = max(1, size(sdmm, 2) - (sdmsc.C.FirstConfoundPredictor - 1));
        ststr(stc).NrOfTimePoints = stnumtp;
        ststr(stc).SubjectID = mdmsubst{stc};
        numtp = numtp + stnumtp;

        % and keep track of names (that are not confounds)
        preds{stc} = sdmsc.C.PredictorNames(:);
        predr{stc} = preds{stc}(1:max(1, min(numel(preds{stc}), ...
            sdmsc.C.FirstConfoundPredictor - 1)));

        % then re-copy sub-portion
        predr{stc} = preds{stc}(1:max(1, min(numel(preds{stc}), ...
            sdmsc.C.FirstConfoundPredictor - 1)));

        % fill in colors
        for pcc = 1:numel(preds{stc})
            predlab = makelabel(preds{stc}{pcc});
            if ~isfield(predcol, predlab)
                predcol.(predlab) = [sdmsc.C.PredictorColors(pcc, :); zeros(3, 3)];
            end
        end
    end
catch xfferror
    clearxffobjects(rfobjo(:));
    rethrow(xfferror);
end

% reorder ststr to comply with selected subjects
pstc = 1;
ordstud = zeros(1, numstudy);
stsubid = {ststr.SubjectID};
stsubid = stsubid(:);
for ssc = 1:numel(subjids)

    % find studies of that subject
    addstud = find(strcmpi(stsubid, subjids{ssc}));
    if ~isempty(addstud)
        ordstud(pstc:pstc+numel(addstud)-1) = addstud;
        pstc = pstc + numel(addstud);
    end
end
if pstc <= numstudy
    clearxffobjects(rfobjo(:));
    error('neuroelf:xff:invalidOperation', 'Error matching subject IDs of studies.');
end
predr = predr(ordstud);
preds = preds(ordstud);
% ststr = ststr(ordstud);
t = t(ordstud);
% opreds = preds;

% get predictors right
tpreds = cat(1, predr{:});
tpredk = true(size(tpreds));
for fc = numel(tpreds):-1:2
    if any(strcmpi(tpreds{fc}, tpreds(1:fc-1)))
        tpredk(fc) = false;
    end
end
tpreds = tpreds(tpredk);

% take care of baseline
tpredc = find(strcmpi('constant', tpreds));
if isempty(tpredc)
    tpreds{end+1} = 'Constant';
else
    tpreds = [tpreds(1:tpredc-1); tpreds(tpredc+1:end); {'Constant'}];
end
tpredn = numel(tpreds);

% prepare outputs
b = zeros(numel(t), tpredn, size(t{1}, 2));
ix = cell(size(t));
p = cell(size(t));
r = cell(size(t));
w = cell(size(t));
if nargout > 7
    d = cell(size(t));
end

% options
cbopt = struct('regdiff', opts.regdiff, 'robust', opts.robust, 'tdim', 1, 'trans', 'none');

% iterate over time courses
findfirst = ne_methods.findfirst;
for tc = 1:numel(t)

    % compute single-study GLM
    sdmc = rfobjs{tc, 2}.C;

    % show design
    if ~isempty(opts.showsdms)
        if ~isempty(showsdmf) && ishandle(showsdmf)
            delete(showsdmf);
        end
        showsdmf = sdm_ShowDesign(rfobjs{tc, 2}, struct('type', opts.showsdms));
    end

    % computation done by SDM::CalcBetas
    [betas, iXX, ptc, se, ws] = sdm_CalcBetas(rfobjs{tc, 2}, t{tc}, cbopt);

    % iterate over predictor names we need
    tpredi = zeros(1, tpredn);
    for bc = 1:tpredn

        % find matching predictor for this study
        predi = findfirst(strcmpi(preds{tc}, tpreds{bc}));

        % if not found, assume confound for last entry
        if isempty(predi)
            if bc < tpredn
                continue;
            end
            predi = size(betas, 1);
        end
        tpredi(bc) = predi;
    end

    % insert into outputs
    tpredu = find(tpredi > 0);
    b(tc, tpredu, :) = betas(:, tpredi(tpredu))';
    ix{tc} = iXX(tpredi(tpredu), tpredi(tpredu));
    p{tc} = sdmc.SDMMatrix(:, tpredi(tpredu)) * betas(:, tpredi(tpredu))';
    r{tc} = t{tc} - sdmc.SDMMatrix * betas';
    w{tc} = ws;
    if nargout > 7
        d{tc} = sdmc.SDMMatrix(:, tpredi(tpredu));
    end
end

% clean up
clearxffobjects(rfobjo(:));
if ~isempty(showsdmf) && ishandle(showsdmf)
    delete(showsdmf);
end

% combine subjects?
if opts.combsub

    % make backup of outputs
    st = t;
    sb = b;
    si = ix;
    sp = p;
    sr = r;
    sw = w;

    % generate new betas and outputs
    b = zeros(numel(opts.subsel), tpredn, size(st{1}, 2));
    t = cell(size(b, 1), 1);
    ix = cell(size(t));
    p = cell(size(t));
    r = cell(size(t));
    w = cell(size(t));
    if nargout > 7
        sd = d;
        d = cell(size(t));
    end

    % get subject IDs
    sids = mdm_Subjects(xo, true);
    sids = sids(rlines);

    % find matches
    for sc = 1:numel(opts.subsel)
        sid = opts.subsel{sc};
        rid = find(strcmpi(sids, sid));
        tps = zeros(numel(rid), 1);
        for ic = 1:numel(rid)
            tps(ic) = size(st{rid(ic)}, 1);
        end

        % weighting
        if opts.vweight
            wv = zeros(numel(rid), tpredn, size(b, 3));
        else
            wv = tps(:, ones(1, tpredn), ones(1, size(b, 3)));
        end

        % re-set time-series
        t{sc} = cat(1, st{rid});
        nix = zeros(numel(rid) * tpredn);
        if nargout > 7
            nd = zeros(sum(tps), numel(rid) * tpredn);
        end
        for ic = 1:numel(rid)
            nix((ic-1)*tpredn+1:ic*tpredn, (ic-1)*tpredn+1:ic*tpredn) = si{rid(ic)};
            if opts.vweight
                dsi = (1 ./ diag(si{rid(ic)}));
                dsi(isinf(dsi) | isnan(dsi)) = 0;
                wv(ic, :, :) = dsi * (sum(sw{rid(ic)}) ./ tps(ic));
            end
            if nargout > 7
                nd(1+sum(tps(1:ic-1)):sum(tps(1:ic)), (ic-1)*tpredn+1:ic*tpredn) = sd{rid(ic)};
            end
            b(sc, :, :) = b(sc, :, :) + sb(rid(ic), :, :) .* wv(ic, :, :);
        end
        b(sc, :, :) = b(sc, :, :) ./ sum(wv);
        ix{sc} = nix;
        p{sc} = cat(1, sp{rid});
        r{sc} = cat(1, sr{rid});
        w{sc} = cat(1, sw{rid});
        if nargout > 7
            d{sc} = nd;
        end
    end
end
