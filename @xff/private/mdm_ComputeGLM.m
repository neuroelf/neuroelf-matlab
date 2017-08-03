function glm = mdm_ComputeGLM(xo, opts)
% MDM::ComputeGLM  - compute a GLM from an MDM file
%
% FORMAT:       glm = mdm.ComputeGLM([options])
%
% Input fields:
%
%       options     optional 1x1 struct with fields
%        .bbox      bounding box for GLM (default: full MNI)
%        .bcomp     HRF-boost computation type (default: 'none')
%        .bftype    HRF-boost basis-function type, either {'bfs'} or 'sdm'
%        .bwtype    HRF-boost weighting type (default: 'varfract')
%        .discard   Sx1 cell array with volumes to discard (in addition)
%        .discmpfd  discard-MPFD threshold (default: Inf)
%        .discmpfdv volumes to discard for MPFD threshold (default: 0,
%                   where, e.g., 2 means detected plus two following!)
%        .ffxwread  read weights from .RunTimeVars.FFXWeights.(NAME)
%                   if non-existant, use column of 1's (equal weighting)
%        .globsigd  also add diff of global signals as nuisance regressors
%        .globsigs  add global signals as confound, one of
%                   0 - none
%                   1 - entire dataset (above threshold/within mask)
%                   2 - two (one per hemisphere, split at BV Z=128)
%                   3 or more, perform PCA of time courses and first N
%                   xff object(s), extract average time course from masks
%        .imeth     interpolation method (where necessary, default: cubic)
%        .ithresh   intensity threshold, default: 100
%        .loadglm   boolean flag, load GLM file named in .outfile
%        .mask      optional masking, default: no mask (for now only VTC)
%        .maxiter   maximum number of robust iterations (default: 30)
%        .motpars   motion parameters (Sx1 cell array with sdm/txt files)
%        .motparsd  also add diff of motion parameters (default: false)
%        .motparsq  also add squared motion parameters (default: false)
%        .ndcreg    if set > 0, perform deconvolution (only with PRTs!)
%        .orthconf  orthogonalize confounds (and motion parameters, true)
%        .outfile   output filename of GLM file, default: no saving
%        .ppicond   list of regressors (or differences) to interact
%        .ppirob    perform robust regression on VOI timecourse and remove
%                   outliers from timecourse/model (threshold, default: 0)
%        .ppitfilt  temporally filter PPI VOI timecourse (default: true)
%        .ppivoi    VOI object used to extract time-course from
%        .ppivoiidx intra-VOI-object index (default: 1)
%        .progress  either {true} or a 1x1 xfigure::progress or xprogress
%        .prtpnorm  normalize parameters of PRT.Conds (true)
%        .redo      selected subjects will be overwritten (default: false)
%        .regdiff   flag, regress first discreet derivatives (diff) instead
%        .res       resolution for GLM (default: 3)
%        .restcond  remove rest condition (rest cond. name, default: '')
%        .robust    perform robust instead of OLS regression
%        .savesdms  token, if not empty, save on-the-fly SDMs (e.g. '.sdm')
%        .showsdms  token, passed to SDM::ShowDesign (if valid)
%        .shuflab   PRT labels (conditions names) to shuffle
%        .shuflabm  minimum number of onsets per label (1x1 or 1xL)
%        .sngtpool  pool all but single trial to one (lower model DF, false)
%                   see http://ncbi.nlm.nih.gov/pmc/articles/PMC3251697
%        .sngtrial  single-trial GLM (only with PRTs, default: false)
%        .sngtskip  condition list to skip during single-trial conversion
%        .sortregs  sort regressors (default: true for PRTs, false for SDMs)
%        .subsel    cell array with subject IDs to work on
%        .tfilter   add filter regressors to SDMs (cut-off in secs)
%        .tfilttype temporal filter type (one of {'dct'}, 'fourier', 'poly')
%        .tmaps     store t-maps rather than betas (e.g. for RSA, false)
%        .transio   boolean flag, if true, save GLM and use transio
%        .vweight   combine runs/studies variance-weighted (default: true)
%        .wdvarsfd  weigh volumes by DVARS and FD (if available, false)
%        .writeres  write residual VTCs, default: '' (off, e.g. '%_res.vtc')
%        .xconfound just as motpars, but without restriction on number
%
% Output fields:
%
%       glm         GLM object
%
% Note: if RFX flag in MDM is set to true, predictor separation will be
%       set to "Subjects".
%       if VTC files are used (at least one VTC) or the GLM is loaded
%       the .bbox and .res fields are not taken into consideration
%       if .outfile is given, GLM is saved after each subject for robust
%       regression models (to allow later continuation of crashed/broken
%       computation using the .subsel field)
%       for .ppi to work, the model filenames must be PRTs!
%       all additional fields for the call to PRT::CreateSDM are supported!
%
% Using: bvcoordconv, ddeblank, emptystruct, findfirst, fitrobustbisquare_img,
%        glmtstat, hrfboost, hsvconv, importvtcfromanalyze, lsqueeze,
%        makelabel, multimatch, poolnonsingletrial, psctrans, resampleaa,
%        robustmean, varc, ztrans.

% TODO:
%        .ar1       boolean flag, correction for AR(1), default: no

% Version:  v1.1
% Build:    17061217
% Date:     Jun-12 2017, 5:39 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2015, 2016, 2017, Jochen Weber
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

% neuroelf library (set in this workspace)
using(neuroelf, {'bvcoordconv', 'ddeblank', 'emptystruct', 'findfirst', ...
    'fitrobustbisquare_img', 'glmtstat', 'hrfboost', 'hsvconv', ...
    'importvtcfromanalyze', 'lsqueeze', 'makelabel', 'multimatch', ...
    'poolnonsingletrial', 'psctrans', 'resampleaa', 'robustmean', 'varc', ...
    'ztrans'});

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'mdm')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
rtv = bc.RunTimeVars;

% check options
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'ar1') || ~islogical(opts.ar1) || isempty(opts.ar1)
    opts.ar1 = false;
else
    opts.ar1 = opts.ar1(1);
end
if ~isfield(opts, 'bbox') || ~isa(opts.bbox, 'double') || ~isequal(size(opts.bbox), [2, 3]) || ...
    any(isinf(opts.bbox(:)) | isnan(opts.bbox(:)) | opts.bbox(:) < 0 | opts.bbox(:) > 256 | ...
    opts.bbox(:) ~= fix(opts.bbox(:))) || any(diff(opts.bbox) < 0)
    opts.bbox = [44, 38, 44; 242, 194, 212];
end
if ~isfield(opts, 'bcomp') || ~ischar(opts.bcomp) || isempty(opts.bcomp) || ...
   ~any('abhmnp' == lower(opts.bcomp(1)))
    opts.bcomp = 'n';
else
    opts.bcomp = lower(opts.bcomp(1));
    if opts.bcomp == 'h'
        opts.bcomp = 'b';
    end
end
if ~isfield(opts, 'bftype') || ~ischar(opts.bftype) || isempty(opts.bftype) || ...
   ~any('bs' == lower(opts.bftype(1)))
    opts.bftype = 'b';
else
    opts.bftype = lower(opts.bftype(1));
end
if ~isfield(opts, 'bwtype') || ~ischar(opts.bwtype) || isempty(opts.bwtype) || ...
   ~any('rv' == lower(opts.bwtype(1)))
    opts.bwtype = 'v';
else
    opts.bwtype = lower(opts.bwtype(1));
end
if isfield(opts, 'discard') && islogical(opts.discard) && numel(opts.discard) == 1 && ...
    opts.discard && isfield(rtv, 'Discard')
    opts.discard = rtv.Discard;
elseif isfield(opts, 'discard') && isa(opts.discard, 'double') && ...
   ~any(isinf(opts.discard(:)) | isnan(opts.discard(:)) | opts.discard(:) < 1)
    opts.discard = repmat({unique(round(opts.discard(:)))}, size(bc.XTC_RTC, 1), 1);
end
if ~isfield(opts, 'discard') || ~iscell(opts.discard) || numel(opts.discard) ~= size(bc.XTC_RTC, 1)
    opts.discard = {};
else
    opts.discard = opts.discard(:);
end
if ~isfield(opts, 'discmpfd') || ~isa(opts.discmpfd, 'double') || numel(opts.discmpfd) ~= 1 || ...
    isnan(opts.discmpfd) || opts.discmpfd <= 0
    opts.discmpfd = Inf;
end
if ~isfield(opts, 'discmpfdv') || ~isa(opts.discmpfdv, 'double') || numel(opts.discmpfdv) ~= 1 || ...
    isinf(opts.discmpfdv) || isnan(opts.discmpfdv) || opts.discmpfdv < 0
    opts.discmpfdv = 0;
else
    opts.discmpfdv = floor(opts.discmpfdv);
end
if ~isfield(opts, 'ffxwread') || ~ischar(opts.ffxwread) || isempty(opts.ffxwread) || ...
   ~isvarname(opts.ffxwread(:)')
    opts.ffxwread = '';
else
    opts.ffxwread = opts.ffxwread(:)';
end
cleargs = false;
if isfield(opts, 'globsigs') && ischar(opts.globsigs) && ...
    any(strcmpi(opts.globsigs(:)', {'cb2', 'cb23', 'cb3', 'cb33'}))
    mp = neuroelf_path('masks');
    switch (lower(opts.globsigs(:)'))
        case {'cb2'}
            if exist([mp '/colin_brain_ICBMnorm_brain2mm.msk'], 'file') == 2
                try
                    opts.globsigs = {xff([mp '/colin_brain_ICBMnorm_brain2mm.msk'])};
                    cleargs = true;
                catch xfferror
                    fprintf('Error loading 2mm Colin brain mask: %s.\n', xfferror.message);
                end
            end
        case {'cb23'}
            if exist([mp '/colin_brain_ICBMnorm_gray2mm.msk'], 'file') == 2 && ...
                exist([mp '/colin_brain_ICBMnorm_white2mm.msk'], 'file') == 2 && ...
                exist([mp '/colin_brain_ICBMnorm_csfx2mm.msk'], 'file') == 2
                try
                    opts.globsigs = {xff([mp '/colin_brain_ICBMnorm_gray2mm.msk']); ...
                        xff([mp '/colin_brain_ICBMnorm_white2mm.msk']); ...
                        xff([mp '/colin_brain_ICBMnorm_csfx2mm.msk'])};
                    cleargs = true;
                catch xfferror
                    fprintf('Error loading 2mm Colin brain masks: %s.\n', xfferror.message);
                end
            end
        case {'cb3'}
            if exist([mp '/colin_brain_ICBMnorm_brain3mm.msk'], 'file') == 2
                try
                    opts.globsigs = {xff([mp '/colin_brain_ICBMnorm_brain3mm.msk'])};
                    cleargs = true;
                catch xfferror
                    fprintf('Error loading 3mm Colin brain mask: %s.\n', xfferror.message);
                end
            end
        case {'cb33'}
            if exist([mp '/colin_brain_ICBMnorm_gray3mm.msk'], 'file') == 2 && ...
                exist([mp '/colin_brain_ICBMnorm_white3mm.msk'], 'file') == 2 && ...
                exist([mp '/colin_brain_ICBMnorm_csfx3mm.msk'], 'file') == 2
                try
                    opts.globsigs = {xff([mp '/colin_brain_ICBMnorm_gray3mm.msk']); ...
                        xff([mp '/colin_brain_ICBMnorm_white3mm.msk']); ...
                        xff([mp '/colin_brain_ICBMnorm_csfx3mm.msk'])};
                    cleargs = true;
                catch xfferror
                    fprintf('Error loading 3mm Colin brain masks: %s.\n', xfferror.message);
                end
            end
    end
    if ischar(opts.globsigs)
        warning('neuroelf:xff:fileNorFound', 'Required mask file(s) not found.');
        opts.globsigs = [];
    end
end
if ~isfield(opts, 'imeth') || ~ischar(opts.imeth) || ...
   ~any(strcmpi(opts.imeth(:)', {'cubic', 'lanczos3', 'linear', 'nearest'}))
    opts.imeth = 'cubic';
else
    opts.imeth = lower(opts.imeth(:)');
end
if ~isfield(opts, 'ithresh') || ~isnumeric(opts.ithresh) || numel(opts.ithresh) ~= 1
    if isfield(rtv, 'IntensityThreshold') && numel(rtv.IntensityThreshold) == 1
        opts.ithresh = rtv.IntensityThreshold;
    else
        opts.ithresh = 100;
    end
else
    opts.ithresh = double(real(opts.ithresh(1)));
end
opts.mskf = '';
if ~isfield(opts, 'loadglm') || ~islogical(opts.loadglm) || isempty(opts.loadglm)
    opts.loadglm = false;
else
    opts.loadglm = opts.loadglm(1);
end
if ~isfield(opts, 'mask') || isempty(opts.mask) || (~ischar(opts.mask) && ...
    (numel(opts.mask) ~= 1 || ~xffisobject(opts.mask, true, 'msk')))
    opts.mask = [];
elseif ischar(opts.mask)
    try
        opts.mskf = opts.mask(:)';
        opts.mask = xff(opts.mskf);
        if ~xffisobject(opts.mask, true, 'msk')
            error('INVALID_MASK');
        end
    catch xfferror;
        warning('neuroelf:xff:badArgument', 'Invalid mask file specified: %s.', opts.mskf);
        neuroelf_lasterr(xfferror);
        opts.mask = [];
        opts.mskf = '';
    end
end
if ~isempty(opts.mask)
    mbc = opts.mask.C;
    opts.mask = (mbc.Mask(:, :, :) ~= 0);
end
if ~isfield(opts, 'maxiter') || ~isa(opts.maxiter, 'double') || numel(opts.maxiter) ~= 1 || ...
    isinf(opts.maxiter) || isnan(opts.maxiter) || opts.maxiter <= 1
    opts.maxiter = 30;
else
    opts.maxiter = min(100, max(2, round(opts.maxiter)));
end
if isfield(opts, 'motpars') && islogical(opts.motpars) && numel(opts.motpars) == 1 && ...
    opts.motpars && isfield(rtv, 'MotionParameters') && iscell(rtv.MotionParameters) && ...
    numel(rtv.MotionParameters) == size(bc.XTC_RTC, 1)
    opts.motpars = rtv.MotionParameters;
end
if ~isfield(opts, 'motparsd') || ~islogical(opts.motparsd) || numel(opts.motparsd) ~= 1
    opts.motparsd = false;
end
if ~isfield(opts, 'motparsq') || ~islogical(opts.motparsq) || numel(opts.motparsq) ~= 1
    opts.motparsq = false;
end
if ~isfield(opts, 'nderiv') || ~isa(opts.nderiv, 'double') || ...
    any(isinf(opts.nderiv(:)) | isnan(opts.nderiv(:))) || ...
   ~all(opts.nderiv(:) == 1 | opts.nderiv(:) == 2)
    opts.nderiv = [];
else
    opts.nderiv = unique(opts.nderiv(:))';
end
if ~isfield(opts, 'outfile') || ~ischar(opts.outfile) || numel(opts.outfile) < 5 || ...
   ~strcmpi(opts.outfile(end-3:end), '.glm')
    opts.outfile = '';
    opts.loadglm = false;
    opts.transio = false;
else
    opts.outfile = opts.outfile(:)';
end
if ~isfield(opts, 'progress') || numel(opts.progress) ~= 1 || ...
   (~islogical(opts.progress) && ~isa(opts.progress, 'xprogress') && ~isxfigure(opts.progress, true))
    opts.progress = true;
end
if ~isfield(opts, 'redo') || ~islogical(opts.redo) || numel(opts.redo) ~= 1
    opts.redo = false;
end
if ~isfield(opts, 'regdiff') || ~islogical(opts.regdiff) || numel(opts.regdiff) ~= 1
    opts.regdiff = false;
end
if ~isfield(opts, 'res') || ~isa(opts.res, 'double') || numel(opts.res) ~= 1 || ...
    isinf(opts.res) || isnan(opts.res) || opts.res ~= fix(opts.res) || opts.res < 1 || opts.res > 12
    opts.res = 3;
end
if ~isfield(opts, 'rfx') || numel(opts.rfx) ~= 1 || ~islogical(opts.rfx)
    opts.rfx = (bc.RFX_GLM > 0);
end
if ~opts.rfx && opts.loadglm
    warning('neuroelf:xff:notSupported', 'Only RFX GLMs can be computed incrementally.');
    opts.loadglm = false;
end
if ~isfield(opts, 'robust') || numel(opts.robust) ~= 1 || ~islogical(opts.robust)
    opts.robust = false;
end
if opts.robust && ~isempty(opts.outfile)
    opts.transio = true;
end
if ~isfield(opts, 'seppred') || ~isa(opts.seppred, 'double') || numel(opts.seppred) ~= 1 || ...
   ~any(0:2 == opts.seppred)
    if numel(bc.SeparatePredictors) == 1 && isa(bc.SeparatePredictors, 'double') && ...
        any(0:2 == bc.SeparatePredictors)
        opts.seppred = bc.SeparatePredictors;
    else
        opts.seppred = 0;
    end
end
if opts.rfx
    opts.seppred = 2;
end
if ~isfield(opts, 'showsdms') || ~ischar(opts.showsdms) || ...
    isempty(opts.showsdms) || ~any(lower(opts.showsdms(1)) == 'iop')
    opts.showsdms = '';
else
    opts.showsdms = lower(opts.showsdms(1));
end
if ~isfield(opts, 'sngtpool') || ~islogical(opts.sngtpool) || numel(opts.sngtpool) ~= 1
    opts.sngtpool = false;
end
if ~isfield(opts, 'sortregs') || ~islogical(opts.sortregs) || numel(opts.sortregs) ~= 1
    opts.sortregs = any(~cellfun('isempty', regexpi(bc.XTC_RTC(:, end), '\.prt$')));
end
mdmsubs = mdm_Subjects(xo);
if ~isfield(opts, 'subsel') || ~iscell(opts.subsel)
    opts.subsel = mdmsubs;
else
    opts.subsel = opts.subsel(:);
    for ssc = numel(opts.subsel):-1:1
        if ~ischar(opts.subsel{ssc}) || isempty(opts.subsel{ssc}) || ...
           ~any(strcmpi(opts.subsel{ssc}(:)', mdmsubs))
            opts.subsel(ssc) = [];
        else
            opts.subsel{ssc} = opts.subsel{ssc}(:)';
        end
    end
end
if ~isfield(opts, 'tfilter') || ~isa(opts.tfilter, 'double') || numel(opts.tfilter) ~= 1 || ...
    isnan(opts.tfilter) || opts.tfilter < 60
    opts.tfilter = Inf;
end
if ~isfield(opts, 'tfilttype') || ~ischar(opts.tfilttype) || isempty(opts.tfilttype) || ...
   ~any(strcmpi(opts.tfilttype(:)', {'dct', 'fourier', 'poly'}))
    opts.tfilttype = 'dct';
else
    opts.tfilttype = lower(opts.tfilttype(:)');
end
if ~isfield(opts, 'tmaps') || ~islogical(opts.tmaps) || numel(opts.tmaps) ~= 1
    opts.tmaps = false;
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
if ~isfield(opts, 'transio') || ~islogical(opts.transio) || numel(opts.transio) ~= 1
    opts.transio = false;
end
if ~isfield(opts, 'vweight') || ~islogical(opts.vweight) || numel(opts.vweight) ~= 1
    opts.vweight = true;
end
if ~isfield(opts, 'wdvarsfd') || ~islogical(opts.wdvarsfd) || numel(opts.wdvarsfd) ~= 1
    opts.wdvarsfd = false;
end
if ~isfield(opts, 'writeres') || ~ischar(opts.writeres) || isempty(opts.writeres) || ...
   ~any(opts.writeres(:) == '%')
    opts.writeres = '';
end

% get and check subject IDs
subjids = opts.subsel;
if isempty(subjids)
    error('neuroelf:xff:badArgument', 'Invalid or missing subject IDs supplied.');
end
if any(strcmp(subjids, ''))
    error('neuroelf:xff:invalidObject', 'Invalid subject IDs for some subjects.');
end
csubjids = subjids;

% let MDM::SDMs do most of the work (for the modelling part)
rfiles = bc.XTC_RTC;
rfobjs = cell(size(rfiles));
mfiles = rfiles(:, end);
for stc = 1:numel(mfiles)
    if isempty(mfiles{stc}) || (~ischar(mfiles{stc}) && ...
        ~xffisobject(mfiles{stc}, true, {'prt', 'sdm'}))
        error('neuroelf:xff:invalidModels', 'Invalid model references.');
    end
    if xffisobject(mfiles{stc}, true)
        mfiles{stc} = mfiles{stc}.F;
    end
end

% determine progress bar capabilities
try
    closepbar = false;
    if islogical(opts.progress)
        if opts.progress
            pbar = xprogress;
            xprogress(pbar, 'setposition', [80, 200, 640, 36]);
            xprogress(pbar, 'settitle', 'Computing multi-study GLM...');
            xprogress(pbar, 0, 'Checking referenced files...', 'visible', 0, 1);
            closepbar = true;
        else
            pbar = [];
        end
    else
        pbar = opts.progress;
        pbarvis = pbar.Visible;
        pbar.Progress(0, 'Checking referenced files...');
        pbar.Visible = 'on';
    end
catch xfferror
    pbar = [];
    neuroelf_lasterr(xfferror);
end
opts.pbar = pbar;

try
    [rfobjs(:, end), stlist, sdmtr, rfobjs(:, end-1), bfs, prtcs] = ...
        mdm_SDMs(xo, opts);
catch xfferror
    if ~isempty(pbar)
        if closepbar
            closebar(pbar);
        else
            pbar.Visible = pbarvis;
        end
    end
    xfferr2 = xfferror;
    if isfield(opts, 'globsigs') && iscell(opts.globsigs) && ~isempty(opts.globsigs) && ...
        xffisobject(opts.globsigs{1}, true)
        try
            if cleargs
                clearxffobjects(opts.globsigs);
            end
        catch xfferror
            neuroelf_lasterr(xfferror);
        end
    end
    rethrow(xfferr2);
end
if isfield(opts, 'globsigs') && iscell(opts.globsigs) && ~isempty(opts.globsigs) && ...
    xffisobject(opts.globsigs{1}, true)
    try
        if cleargs
            clearxffobjects(opts.globsigs);
        end
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end
if size(rfobjs, 2) > 2
    try
        for stc = 1:size(rfobjs, 1)
            rfobjs{stc, 1} = xff(rfiles{stc, 1});
            if ~xffisobject(rfobjs{stc, 1}, true, {'ssm', 'tsm'})
                error('neuroelf:xff:badArgument', ...
                    'MTC-based GLM requires SSM/TSM objects.');
            end
        end
    catch xfferror
        clearxffobjects(rfobjs(:));
        if ~isempty(pbar)
            if closepbar
                closebar(pbar);
            else
                pbar.Visible = pbarvis;
            end
        end
        rethrow(xfferror);
    end
end

% discard based on MPFD alone?
if ~isinf(opts.discmpfd) && isempty(opts.discard)
    opts.discard = cell(size(rfobjs, 1), 1);
    for stc = 1:size(rfobjs, 1)
        xtcc = rfobjs{stc, end - 1}.C;
        if isfield(xtcc.RunTimeVars, 'MPFD')
            discmpfd = 1 + find(xtcc.RunTimeVars.MPFD > opts.discmpfd);
            for dscc = 1:opts.discmpfdv
                discmpfd = union(discmpfd, 1 + dscc + find(xtcc.RunTimeVars.MPFD > opts.discmpfd));
            end
            opts.discard{stc} = union(xtcc.RunTimeVars.Discard, discmpfd(:));
        end
    end
end

% smart-weighing DVARS and FD information
dvars = {};
if opts.wdvarsfd

    % get "general picture" of DVARS and FD
    dvars = cell(size(rfobjs, 1), 1);
    fds = dvars;
    singlews = dvars;
    for stc = 1:numel(fds)
        xtcc = rfobjs{stc, end - 1}.C;
        if isfield(xtcc.RunTimeVars, 'DVARS')
            dvars{stc} = xtcc.RunTimeVars.DVARS{1, 2};
        end
        if isfield(xtcc.RunTimeVars, 'MPFD')
            fds{stc} = xtcc.RunTimeVars.MPFD;
        end
    end

    % overall
    stc = cat(1, dvars{:});
    mdvars = robustmean(stc);
    sdvars = std(stc);
    stc = cat(1, fds{:});
    mfds = robustmean(stc);
    sfds = std(stc);

    % some volumes will be discarded completely!
    if isempty(opts.discard)
        opts.discard = cell(numel(fds), 1);
    end

    % now parse weights
    for stc = 1:numel(fds)
        if ~isa(opts.discard{stc}, 'double') || any(isinf(opts.discard{stc}(:)) | ...
            isnan(opts.discard{stc}(:)) | opts.discard{stc}(:) < 1)
            opts.discard{stc} = [];
        end

        % generate weights
        singlew = 1;
        if ~isempty(dvars{stc})
            xw = max(0, 1 - (abs((dvars{stc} - mdvars) ./ (1.65 * sdvars)) .^ 2));
            singlew = min(min([xw; 1], [1; xw]), [1; 1; xw(2:end)]);
        end
        if ~isempty(fds{stc})
            xw = max(0, 1 - (abs((fds{stc} - mfds) ./ (1.65 * sfds)) .^ 2));
            xw = min(min([xw; 1], [1; xw]), [1; 1; xw(2:end)]);
            singlew = min(singlew, xw);
        end

        % update discard and store?
        if numel(singlew) > 1
            opts.discard{stc} = union(opts.discard{stc}(:), find(singlew(:) == 0));
            singlews{stc} = singlew;
        end
    end
end

% add discarded volumes to objects
if ~isempty(opts.discard)
    for stc = 1:size(rfobjs, 1)
        if ~isempty(opts.discard{stc}) && isa(opts.discard{stc}, 'double') && ...
           ~any(isinf(opts.discard{stc}(:)) | isnan(opts.discard{stc}(:)) | opts.discard{stc}(:) < 1)
            xtcc = rfobjs{stc, end - 1}.C;
            if ~isfield(xtcc.RunTimeVars, 'Discard')
                xtcc.RunTimeVars.Discard = unique(round(opts.discard{stc}(:)));
            else
                xtcc.RunTimeVars.Discard = union(xtcc.RunTimeVars.Discard(:), ...
                    unique(round(opts.discard{stc}(:))));
            end
            rfobjs{stc, end - 1}.C = xtcc;
        end
    end
end

% load a pre-defined GLM
xffroot = xff();
tioosz = root_TransIOSize(xffroot, 1e5);
oglm = {};
olay = [];
otrf = [];
ogns = 0;
if opts.loadglm && exist(opts.outfile, 'file') == 2
    try
        if ~isempty(pbar)
            pbar.Progress(0, ['Reading in existing GLM: ' opts.outfile]);
        end
        oglm = {xff(opts.outfile)};
        if ~xffisobject(oglm{1}, true, 'glm')
            error('neuroelf:xff:BadArgument', ...
                'The GLM named in .outfile is not a GLM file.');
        end

        % get GLM layout, subjects, and load/update mask
        glmc = oglm{1}.C;
        if glmc.ProjectType == 1
            olay = aft_Layout(oglm{1});
            otrf = glmc.RunTimeVars.TrfPlus;
        elseif glmc.ProjectType == 2
            olay = glmc.NrOfVertices;
        end
        glmsubs = glm_Subjects(oglm{1});
        ogns = numel(glmsubs);
        glmspred = glm_SubjectPredictors(oglm{1});
        if isempty(opts.mskf)
            if ~isempty(glmc.CortexBasedStatisticsMaskFile)
                try
                    opts.mskf = glmc.CortexBasedStatisticsMaskFile;
                    opts.mask = xff(opts.mskf);
                    if ~xffisobject(opts.mask, true, 'msk')
                        error('INVALID_MASK');
                    end
                catch xfferror;
                    warning('neuroelf:xff:badArgument', ...
                        'Invalid mask file in GLM specified: %s.', opts.mskf);
                    neuroelf_lasterr(xfferror);
                    opts.mskf = '';
                    if glmc.ProjectTypeRFX > 0
                        opts.mask = (glmc.GLMData.RFXGlobalMap(:, :, :) > 0);
                    else
                        opts.mask = any(glmc.GLMData.BetaMaps(:, :, :, :) ~= 0, 4);
                    end
                end
            end
        else
            glmc.CortexBasedStatisticsMaskFile = opts.mskf;
        end

        % unless redo
        if ~opts.redo

            % remove double subjects from list!
            subjids(multimatch(subjids, glmsubs) > 0) = [];
        end
    catch xfferror
        clearxffobjects(oglm);
        clearxffobjects(rfobjs(:));
        root_TransIOSize(xffroot, tioosz);
        if ~isempty(pbar)
            if closepbar
                closebar(pbar);
            else
                pbar.Visible = pbarvis;
            end
        end
        rethrow(xfferror);
    end

    % nothing left to do
    if isempty(subjids)
        disp('All subjects done.');
        clearxffobjects(rfobjs(:));
        glm = oglm{1};
        if ~isempty(pbar)
            if closepbar
                closebar(pbar);
            else
                pbar.Visible = pbarvis;
            end
        end
        return;
    end
end

% apply study selection
if ~isfield(rtv, 'Subjects') || ~iscell(rtv.Subjects) || ...
    numel(rtv.Subjects) ~= size(rfiles, 1) || ~all(cellfun(@ischar, rtv.Subjects(:))) || ...
    any(cellfun('isempty', rtv.Subjects(:))) || ...
    any(~cellfun('isempty', regexpi(rtv.Subjects(:), '[_ ]')))
    mdmsubst = mdm_Subjects(xo, true);
else
    mdmsubst = rtv.Subjects(:);
end
mdmsubls = mdmsubst;
for stc = 1:numel(mdmsubls)
    mdmsubls{stc} = makelabel(ddeblank(mdmsubls{stc}));
end
mdmsnmat = struct;
mdmtrfpl = struct;
studies = find(multimatch(mdmsubst, subjids) > 0);
numstudy = numel(studies);

% check ?TCs
numtp = 0;
vtclay = [];
hvsz = zeros(numstudy, 3);
vtrf = cell(numstudy, 1);
reqi = false(numstudy, 1);
hashdr = false;
try
    if size(rfobjs, 2) == 2
        fltypes = {'fmr', 'hdr', 'vtc', 'mtc'};
        layout = aft_Layout(rfobjs{studies(1), end - 1});
        if layout(end) == 875706624
            layout = olay;
        end
    else
        fltypes = {'mtc'};
        layout = rfobjs{studies(1), 1}.C.NrOfTargetVertices;
    end
    li = [1:3, 5:13, 15];
    if (numel(olay) > 1 && any(olay(li(1:end-1)) ~= layout(li(1:end-1)))) || ...
       (numel(olay) == 1 && layout ~= olay)
        error('neuroelf:xff:layoutMismatch', ...
            'Loaded GLM and XTCs must match in layout!');
    end
    predcol = struct;
    preds = cell(numstudy, 1);
    predr = cell(numstudy, 1);
    ststr = emptystruct({'NrOfTimePoints', 'NameOfAnalyzedFile', 'NameOfSSMFile', ...
        'NameOfSDMFile', 'NrOfConfounds', 'RunTimeVars', 'SubjectID'}, [1, numstudy]);

    % tally FFX data
    if ~opts.rfx
        ffxpreds = 0;
        ffxconfs = 0;
        ffxstntp = zeros(numstudy, 1);
    end

    % check studies that make it into the dataset
    for stc = 1:numstudy

        % subject label
        stslab = mdmsubls{stc};

        % update progress bar
        fc = studies(stc);
        if ~isempty(pbar)
            pbar.Progress(0, sprintf('Checking dims of study %d/%d...', stc, numstudy));
        end

        % check objects
        xtcsc = rfobjs{fc, end - 1};
        xtcc = xtcsc.C;
        sdmsc = rfobjs{fc, end};
        sdmc = sdmsc.C;
        if isfield(xtcc, 'NrOfVolumes')
            stnumtp = xtcc.NrOfVolumes;
        elseif isfield(xtcc, 'NrOfTimePoints')
            stnumtp = xtcc.NrOfTimePoints;
        else
            stnumtp = xtcc.ImgDim.Dim(5);
        end
        if ~isfield(xtcc, 'NrOfTimePoints')
            vtrf{stc} = xtcc.RunTimeVars.TrfPlus;
        else
            vtrf{stc} = [];
        end
        if size(rfobjs, 2) > 2
            ssmc = rfobjs{fc, 1}.C;
            if ssmc.NrOfTargetVertices ~= layout
                error('neuroelf:xff:layoutMismatch', ...
                    'NrOfTargetVertices must match across studies.');
            end
            ststr(stc).NameOfSSMFile = rfobjs{fc, 1}.F;
        else
            tlayout = aft_Layout(xtcsc);
            if tlayout(end) ~= 875706624 && ~isequal(layout(li), tlayout(li))
                error('neuroelf:xff:layoutMismatch', ...
                    'Data layout must be consistent across studies.');
            elseif tlayout(end) == 875706624
                hdrcf = hdr_CoordinateFrame(xtcsc);
                hdrsz = size(xtcc.VoxelData);
                if numel(hdrsz) < 3
                    hdrsz(3) = 1;
                end
                vtrf{stc} = hdrcf.Trf * [0, -1, 0, 129; 0, 0, -1, 129; -1, 0, 0, 129; 0, 0, 0, 1];
                hvsz(stc, :) = hdrsz(1:3);
                reqi(stc) = true;
            end
            if tlayout(end) == 993669504 && isempty(vtclay)
                vtclay = tlayout;
            end
            ststr(stc).NameOfSSMFile = '';
        end
        if ~any(strcmpi(xtcsc.S.Extensions{1}, fltypes))
            error('neuroelf:xff:notSupported', ...
                'Filetype %s not supported as time-course file.', ...
                upper(xtcsc.S.Extensions{1}) ...
            );
        end
        if lower(xtcsc.S.Extensions{1}(1)) == 'h'
            hashdr = true;

        % for VTCs, also check and store SPMsn and TrfPlus
        elseif lower(xtcsc.S.Extensions{1}(1)) == 'v'
            xtcrtv = xtcc.RunTimeVars;
            if isfield(xtcrtv, 'SPMsn') && isstruct(xtcrtv.SPMsn) && ...
                numel(xtcrtv.SPMsn) == 1 && isfield(xtcrtv.SPMsn, 'VG') && ...
                isstruct(xtcrtv.SPMsn.VG) && ~isempty(xtcrtv.SPMsn.VG) && ...
                isfield(xtcrtv.SPMsn.VG, 'dim') && isfield(xtcrtv.SPMsn.VG, 'mat') && ...
                isfield(xtcrtv.SPMsn, 'VF') && isstruct(xtcrtv.SPMsn.VF) && ...
                numel(xtcrtv.SPMsn.VF) == 1 && isfield(xtcrtv.SPMsn.VF, 'dim') && ...
                isfield(xtcrtv.SPMsn.VF, 'mat') && isfield(xtcrtv.SPMsn, 'Tr') && ...
                isfield(xtcrtv.SPMsn, 'Affine')
                spmsn = xtcrtv.SPMsn;
                if isfield(mdmsnmat, stslab)
                    if ~isequal(mdmsnmat.(stslab).VG(1).dim, spmsn.VG(1).dim) || ...
                       ~isequal(mdmsnmat.(stslab).VG(1).mat, spmsn.VG(1).mat) || ...
                       ~isequal(mdmsnmat.(stslab).VF.dim, spmsn.VF.dim) || ...
                       ~isequal(mdmsnmat.(stslab).VF.mat, spmsn.VF.mat) || ...
                       ~isequal(mdmsnmat.(stslab).Tr, spmsn.Tr) || ...
                       ~isequal(mdmsnmat.(stslab).Affine, spmsn.Affine)
                        error('neuroelf:xff:SPMsnMismatch', ...
                            'SPM-based SN-mat content must match within subject (%s).', stslab);
                    end
                else
                    mdmsnmat.(stslab) = spmsn;
                end
            end
            if isfield(mdmtrfpl, stslab)
                if ~isequal(mdmtrfpl.(stslab), vtrf{stc})
                    error('neuroelf:xff:trfPlusMismatch', ...
                        'TrfPlus must match within subject (%s).', stslab);
                end
            elseif ~isequal(vtrf{stc}, eye(4))
                mdmtrfpl.(stslab) = vtrf{stc};
            end
        end
        if stnumtp ~= sdmc.NrOfDataPoints
            error('neuroelf:xff:layoutMismatch', ...
                'Number of time points between SDM and VTC must match.');
        end
        sdmm = sdmc.SDMMatrix;
        ststr(stc).NameOfAnalyzedFile = xtcsc.F;
        ststr(stc).NameOfSDMFile = sdmsc.F;
        if isempty(sdmsc.F) && isfield(sdmc.RunTimeVars, 'CreatedFromPRT')
            ststr(stc).NameOfSDMFile = ...
                sprintf('<created from %s>', sdmc.RunTimeVars.CreatedFromPRT);
        end
        ststr(stc).NrOfConfounds = ...
            max(1, size(sdmm, 2) - (sdmc.FirstConfoundPredictor - 1));
        ststr(stc).NrOfTimePoints = stnumtp;
        if size(sdmm, 1) > size(sdmm, 2)
            if all(sum(abs(sdmm)) > 1e-4)
                sdmmi = inv(sdmm' * sdmm);
            else
                sdmmi = pinv(sdmm' * sdmm);
            end
        else
            sdmnzi = find(any(sdmm ~= 0, 1));
            sdmmi = zeros(size(sdmm, 1));
            if numel(sdmnzi) >= size(sdmm, 1)
                fprintf('Warning: too many columns in SDMMatrix (%s).\n', ...
                    ststr(stc).NameOfSDMFile);
            end
            sdmmi(sdmnzi, sdmnzi) = pinv(sdmm(:, sdmnzi)' * sdmm(:, sdmnzi));
            sdmmi = sparse(sdmmi);
        end
        ststr(stc).RunTimeVars = struct( ...
            'FWHMResEst',    [Inf, Inf, Inf], ...
            'FromFile',      mfiles{stc}, ...
            'NrOfConfounds', ststr(stc).NrOfConfounds, ...
            'PRTContent',    prtcs{stc}, ...
            'Predictors',    {sdmc.PredictorNames}, ...
            'SDMMatrix',     sdmm, ...
            'SDMMatrixInv',  sdmmi);
        if isfield(xtcc.RunTimeVars, 'Discard')
            ststr(stc).RunTimeVars.TCDiscard = xtcc.RunTimeVars.Discard;
        end
        if ~isempty(dvars)
            ststr(stc).RunTimeVars.DVARS = dvars{stc};
            ststr(stc).RunTimeVars.MPFD = fds{stc};
            ststr(stc).RunTimeVars.SingleWs = singlews{stc};
        end
        ststr(stc).SubjectID = mdmsubst{fc};
        numtp = numtp + stnumtp;

        % and keep track of names (that are not confounds)
        preds{stc} = sdmc.PredictorNames(:);
        predr{stc} = preds{stc}(1:max(1, min(numel(preds{stc}), ...
            sdmc.FirstConfoundPredictor - 1)));

        % tally FFX and alter names
        if ~opts.rfx

            % keep track of number of time points
            ffxstntp(stc) = stnumtp;

            % depending on separation status
            switch (opts.seppred)
                case {0}
                    ffxpreds = max(ffxpreds, numel(predr{stc}));
                case {1}
                    ffxpreds = ffxpreds + numel(predr{stc});
                    for pcc = 1:numel(predr{stc})
                        preds{stc}{pcc} = sprintf( ...
                            'Study %d: %s', stc, preds{stc}{pcc});
                    end
                case {2}
                    for pcc = 1:numel(predr{stc})
                        preds{stc}{pcc} = sprintf( ...
                            'Subject %s: %s', mdmsubst{fc}, preds{stc}{pcc});
                    end
            end
            ffxconfs = ffxconfs + ststr(stc).NrOfConfounds;

            % update
            for pcc = (numel(predr{stc}) + 1):numel(preds{stc})
                preds{stc}{pcc} = sprintf('Study %d: %s', stc, preds{stc}{pcc});
            end
            sdmc.PredictorNames = preds{stc};
            sdmsc.C = sdmc;
        end

        % then re-copy sub-portion
        predr{stc} = preds{stc}(1:max(1, min(numel(preds{stc}), ...
            sdmc.FirstConfoundPredictor - 1)));

        % fill in colors
        for pcc = 1:numel(preds{stc})
            predlab = makelabel(preds{stc}{pcc});
            if ~isfield(predcol, predlab)
                predcol.(predlab) = [sdmc.PredictorColors(pcc, :); zeros(3, 3)];
            end
        end
    end

    % check that TrfPlus are also in line!
    mixhdr = true;
    if any(hvsz(:) > 0)
        if any(any(diff(hvsz(hvsz(:, 1) > 0, :), 1, 1) ~= 0))
            error('neuroelf:xff:layoutMismatch', 'Loaded HDRs must match in size!');
        end
        if all(hvsz(:) > 0) && (isempty(oglm) || isequal(otrf, vtrf{1}))
            mixhdr = false;
            opts.bbox = [0, 0, 0; hvsz(1, :)];
        end
    end

    % multiple subjects normalization NOT FOR FFX!
    if ~opts.rfx && (numel(fieldnames(mdmsnmat)) > 1 || numel(fieldnames(mdmtrfpl)) > 1)
        error('neuroelf:xff:notSupported', ...
            'Multiple normalizations not supported for FFX-GLMs.');
    end
catch xfferror
    root_TransIOSize(xffroot, tioosz);
    if ~isempty(pbar);
        if closepbar
            closebar(pbar);
        else
            pbar.Visible = pbarvis;
        end
    end
    clearxffobjects(oglm);
    clearxffobjects(rfobjs(:));
    rethrow(xfferror);
end
root_TransIOSize(xffroot, tioosz);

% reorder ststr to comply with selected subjects
pstc = 1;
ordstud = zeros(1, numstudy);
stfname = {ststr.NameOfAnalyzedFile};
stsubid = {ststr.SubjectID};
stsubid = stsubid(:);
for ssc = 1:numel(subjids)

    % find studies of that subject
    addstud = find(strcmpi(stsubid, subjids{ssc}));
    ordstud(pstc:pstc+numel(addstud)-1) = addstud;
    pstc = pstc + numel(addstud);
end
if pstc <= numstudy || numel(unique(stfname)) ~= numel(stfname)
    clearxffobjects(rfobjs(:));
    clearxffobjects(oglm);
    if ~isempty(pbar)
        if closepbar
            closebar(pbar);
        else
            pbar.Visible = pbarvis;
        end
    end
    error('neuroelf:xff:invalidOperation', 'Error matching subject IDs of studies.');
end
predr = predr(ordstud);
preds = preds(ordstud);
ststr = ststr(ordstud);
if opts.wdvarsfd
    singlews = singlews(ordstud);
end
opreds = preds;
ndsti = {};

% add HRF-boost combination if requested
if ~isempty(opts.nderiv) && opts.rfx && opts.bcomp ~= 'n'

    % resample bfs (if used)
    if opts.bftype == 'b'
        bfs = resampleaa(bfs, size(bfs, 1) / 100);
        if all(bfs(end, :) == 0)
            bfs(end, :) = [];
        end
    end

    % hrfboost options
    hrfbopt = struct('bf', bfs, 'comp', opts.bcomp, 'wcutoff', 0.01, 'wtype', opts.bwtype);

    % matched pattern
    ndpat = sprintf(' - %d.deriv', opts.nderiv(end));

    % source/target indices for derivatives
    ndsti = cell(numel(predr), 1);

    % parse all predictor structs
    for ssc = 1:numel(predr)

        % create required space
        ndsti{ssc} = cell(numel(preds{ssc}), 1);

        % go through list of predictors
        for spc = numel(predr{ssc}):-1:1

            % last derivative?
            if ~isempty(strfind(predr{ssc}{spc}, ndpat))

                % add boosted parameter after
                predr{ssc} = [predr{ssc}(1:spc); ...
                    {strrep(predr{ssc}{spc}, ndpat, ' - HRFboost')}; predr{ssc}(spc+1:end)];
            end
        end

        % do the same for the list of regressors
        for spc = numel(preds{ssc}):-1:1

            % last derivative?
            if ~isempty(strfind(preds{ssc}{spc}, ndpat))

                % add boosted parameter after
                preds{ssc} = [preds{ssc}(1:spc); ...
                    {strrep(preds{ssc}{spc}, ndpat, ' - HRFboost')}; preds{ssc}(spc+1:end)];
            end
        end

        % and then determine the source/target indices for HRFboost
        hrfbt = find(~cellfun('isempty', regexp(preds{ssc}, ' - HRFboost$')));
        for spc = hrfbt(:)'

            % find predictors that match
            hrfbs = find(~cellfun('isempty', regexp(opreds{ssc}, ...
                strrep(preds{ssc}{spc}, ' - HRFboost', ''))));

            % store source indices in target
            if numel(hrfbs) == size(bfs, 2)
                ndsti{ssc}{spc} = hrfbs(1:end);

                % requires predictor color
                if ~isfield(predcol, makelabel(preds{ssc}{spc}))

                    % get HSV values
                    predhsv = hsvconv(uint8(predcol.(makelabel(opreds{ssc}{hrfbs(1)}))), 2);

                    % set saturation and value closer to 1
                    predhsv(1, 2:3) = 1 - 0.5 * (1 - predhsv(1, 2:3));

                    % store boosted color
                    predcol.(makelabel(preds{ssc}{spc})) = double(hsvconv(predhsv));
                end
            end
        end
    end
end

% get predictors right
tpreds = cat(1, predr{:});
if ~opts.rfx
    ffxstntp = ffxstntp(ordstud);
    tpreds = [tpreds; cat(1, preds{:})];
    ffxstud = rfobjs(studies(ordstud), :);
end
tpredk = true(size(tpreds));
for fc = numel(tpreds):-1:2
    if any(strcmpi(tpreds{fc}, tpreds(1:fc-1)))
        tpredk(fc) = false;
    end
end
tpreds = tpreds(tpredk);

% attempt to sort according to stlist
if opts.sortregs
    tpredi = 1;
    for ssc = 1:numel(stlist)
        tpredk = ~cellfun('isempty', regexp(tpreds(tpredi:end), ['^' stlist{ssc}]));
        if any(tpredk)
            tpredc = find(tpredk) + (tpredi - 1);
            tprednc = find(~tpredk) + (tpredi - 1);
            tpreds = [tpreds(1:tpredi-1); tpreds(tpredc); tpreds(tprednc)];
            tpredi = tpredi + numel(tpredc);
        end
        if tpredi > numel(tpreds)
            break;
        end
    end
end

% take care of baseline
if opts.rfx
    tpredc = find(strcmpi('constant', tpreds));
    if isempty(tpredc)
        tpreds{end+1} = 'Constant';
    else
        tpreds = [tpreds(1:tpredc-1); tpreds(tpredc+1:end); {'Constant'}];
    end
end
tpredn = numel(tpreds);

% GLM loaded?
if ~isempty(olay)

    % get handle
    glm = oglm{1};

    % compare subject predictors
    if numel(tpreds) ~= numel(glmspred) || ~all(strcmpi(tpreds, glmspred))

        % check that at least no new predictors are present!
        for spc = 1:tpredn
            if ~any(strcmpi(tpreds{spc}, glmspred))

                % clear objects and bail out
                clearxffobjects(rfobjs(:));
                clearxffobjects(oglm);
                if ~isempty(pbar)
                    if closepbar
                        closebar(pbar);
                    else
                        pbar.Visible = pbarvis;
                    end
                end
                error('neuroelf:xff:invalidOperation', ...
                    'Adding new predictors to existing GLM prohibited.');
            end
        end
    end

    % do NOT update colors
    nsp = glmc.NrOfSubjectPredictors;
    predcols = fieldnames(predcol);
    for pcc = 1:nsp
        predname = makelabel(regexprep(glmc.Predictor(pcc).Name2, '^.*:\s*', ''));
        predcolm = findfirst(strcmpi(predcols, predname));
        if ~isempty(predcolm)
            predcol.(predcols{predcolm}) = glmc.Predictor(pcc).RGB;
        else
            predcol.(predname) = glmc.Predictor(pcc).RGB;
        end
    end
    predcols = fieldnames(predcol);

    % output size
    outsize = size(glmc.GLMData.RFXGlobalMap);
    boutsz = [outsize, tpredn];
    outsz3 = ones(1, 3);
    outsz3(1:numel(outsize)) = outsize;

    % iterate over subjects in MDM selection
    for ssc = 1:numel(subjids)

        % find studies of that subject
        addstud = find(strcmpi(subjids{ssc}, {ststr.SubjectID}));
        newststr = ststr(addstud);

        % already in GLM?
        if any(strcmpi(subjids{ssc}, glmsubs))

            % re-do
            if opts.redo

                % remove!
                glm_RemoveSubject(glm, subjids{ssc});
                glmc = glm.C;
                glmsubs = glm_Subjects(glm);

            % leave untouched
            else

                % clear from list
                subjids{ssc} = '';

                % and don't go on
                continue;
            end
        end

        % extend GLM!
        glmsubs{end+1} = subjids{ssc};
        nsp = glmc.NrOfSubjectPredictors;
        fnsp = (nsp - 1) * glmc.NrOfSubjects + 1;
        glmc.Predictor = glmc.Predictor( ...
            [1:(fnsp-1), (fnsp+1-nsp):(fnsp-1), fnsp:end, end]);
        for stc = fnsp:numel(glmc.Predictor)
            glmc.Predictor(stc).Name1 = sprintf('Predictor: %d', stc);
        end
        for stc = fnsp:(fnsp + nsp - 2)
            glmc.Predictor(stc).Name2 = ...
                sprintf('Subject %s: %s', subjids{ssc}, glmspred{stc+1-fnsp});
            predcolm = findfirst(strcmpi(predcols, makelabel(glmspred{stc+1-fnsp})));
            if ~isempty(predcolm)
                glmc.Predictor(stc).RGB = predcol.(predcols{predcolm});
            else
                glmc.Predictor(stc).RGB = [floor(255.999 .* rand(3, 1)); zeros(3, 3)];
            end
        end
        glmc.Predictor(end).Name2 = sprintf('Subject %s: Constant', subjids{ssc});
        glmc.NrOfSubjects = glmc.NrOfSubjects + 1;
        glmc.NrOfTimePoints = glmc.NrOfTimePoints + sum(cat(1, newststr.NrOfTimePoints));
        glmc.NrOfPredictors = numel(glmc.Predictor);
        glmc.NrOfStudies = glmc.NrOfStudies + numel(addstud);
        glmc.NrOfStudiesWithConfounds = glmc.NrOfStudiesWithConfounds + numel(addstud);
        glmc.NrOfConfoundsPerStudy = ...
            [glmc.NrOfConfoundsPerStudy, cat(2, newststr.NrOfConfounds)];
        for stc = 1:numel(addstud)
            glmc.Study(end+1).NrOfTimePoints = newststr(stc).NrOfTimePoints;
            glmc.Study(end).NameOfAnalyzedFile = newststr(stc).NameOfAnalyzedFile;
            glmc.Study(end).NameOfSSMFile = newststr(stc).NameOfSSMFile;
            glmc.Study(end).NameOfSDMFile = newststr(stc).NameOfSDMFile;
        end
        glmc.GLMData.Subject(end+1).BetaMaps = single(zeros([outsize, tpredn]));
        glm.C = glmc;
    end

    % remove non-redos from list
    subjids(cellfun('isempty', subjids)) = [];

    % replace with new value
    numsubs = numel(subjids);

    % global map
    if glmc.ProjectTypeRFX > 0
        grfx = (0.75 * glmc.NrOfSubjects) * double(glmc.GLMData.RFXGlobalMap(:, :, :) > 0);
    end

% create GLM object and empty contents
else

    glm = xff('new:glm');
    glmc = glm.C;
    numsubs = numel(subjids);
    glmsubs = subjids;

    % some initial settings
    switch lower(xtcsc.S.Extensions{1}(1))
        case {'f'}
            ptype = 0;
            outsize = [xtcc.ResolutionX, xtcc.ResolutionY, xtcc.NrOfSlices];
        case {'h', 'v'}
            ptype = 1;
            if ~isempty(vtclay)
                outres = vtclay(11);
                outsize = vtclay(1:3);
                opts.bbox = [vtclay(5:7); vtclay(8:10)];
            elseif ~mixhdr
                outres = 1;
                outsize = opts.bbox(2, :);
            else
                outres = opts.res;
                opts.bbox = [min(opts.bbox); min(opts.bbox) + ...
                    outres .* ceil(abs(diff(opts.bbox) - 0.5) ./ outres)];
                outsize = diff(opts.bbox) ./ outres;
            end
        case {'m'}
            ptype = 2;
            outsize = size(xtcc.MTCData, 2); % doesn't take SSM/TSM into account!
        otherwise
            clearxffobjects(oglm);
            clearxffobjects(rfobjs(:));
            if ~isempty(pbar)
                if closepbar
                    closebar(pbar);
                else
                    pbar.Visible = pbarvis;
                end
            end
            error('neuroelf:xff:badObject', 'ProjectType unsupported.');
    end
    glmc.ProjectType = ptype;
    boutsz = [outsize, tpredn];
    outsz3 = ones(1, 3);
    outsz3(1:numel(outsize)) = outsize;
    glmc.NrOfTimePoints = numtp;
    snumpred = tpredn * numsubs;
    glmc.Predictor = cell2struct(cell(1, snumpred, 3), {'Name1', 'Name2', 'RGB'}, 3);
    if opts.rfx
        glmc.ProjectTypeRFX = 1;
        glmc.NrOfSubjects = numsubs;
        glmc.NrOfSubjectPredictors = tpredn;
        prdc = 1;
        for sc = 1:numsubs
            for pc = 1:(tpredn-1)
                glmc.Predictor(prdc).Name1 = sprintf('Predictor: %d', prdc);
                glmc.Predictor(prdc).Name2 =  ...
                    sprintf('Subject %s: %s', subjids{sc}, tpreds{pc});
                glmc.Predictor(prdc).RGB = predcol.(makelabel(tpreds{pc}));
                prdc = prdc + 1;
            end
        end
        for sc = 1:numsubs
            glmc.Predictor(prdc).Name1 = sprintf('Predictor: %d', prdc);
            glmc.Predictor(prdc).Name2 = ...
                sprintf('Subject %s: Constant', subjids{sc});
            glmc.Predictor(prdc).RGB = [255; 0; 0; 0] * ones(1, 3);
            prdc = prdc + 1;
        end
        glmc.GLMData.RFXGlobalMap = single(zeros(outsz3));
        glmc.GLMData.Subject = repmat( ...
            struct('BetaMaps', single(zeros([outsize, tpredn]))), 1, numsubs);
    else
        glmc.ProjectTypeRFX = 0;
    end
    glmc.NrOfPredictors = snumpred;
    glmc.NrOfStudies = numstudy;
    glmc.NrOfStudiesWithConfounds = numstudy;
    glmc.NrOfConfoundsPerStudy = cat(2, ststr.NrOfConfounds);
    glmc.NrOfConfounds = sum(glmc.NrOfConfoundsPerStudy);
    glmc.SeparatePredictors = opts.seppred;
    switch (opts.trans)
        case {'n'}
            glmc.TransformationType = 0;
        case {'p'}
            glmc.TransformationType = 3;
        case {'z'}
            glmc.TransformationType = 1;
    end
    glmc.SerialCorrelation = opts.ar1;
    glmc.MeanARPre = 0;
    glmc.MeanARPost = 0;
    if ptype == 2
        glmc.Resolution = 3;
        glmc.NrOfVertices = outsize;
    elseif ptype == 1
        glmc.Resolution = outres;
        glmc.XStart = opts.bbox(1, 1);
        glmc.XEnd   = opts.bbox(2, 1);
        glmc.YStart = opts.bbox(1, 2);
        glmc.YEnd   = opts.bbox(2, 2);
        glmc.ZStart = opts.bbox(1, 3);
        glmc.ZEnd   = opts.bbox(2, 3);
    else
        if ~isempty(pbar)
            if closepbar
                closebar(pbar);
            else
                pbar.Visible = pbarvis;
            end
        end
        error('neuroelf:xff:notYetImplemented', ...
            'GLM computation of FMR-based data not supported.');
    end
    if isempty(opts.mask)
        glmc.CortexBasedStatistics = 0;
        glmc.CortexBasedStatisticsMaskFile = '';
    else
        glmc.CortexBasedStatistics = 1;
        glmc.CortexBasedStatisticsMaskFile = opts.mskf;
    end
    ststr = rmfield(ststr, 'NrOfConfounds');
    glmc.Study = rmfield(ststr, 'SubjectID');
end

% generating FFX design matrix
if ~opts.rfx
    if ~isempty(pbar)
        pbar.Progress(0, 'Creating FFX design matrix...', 'visible', 0, 1);
    end
    ffxXn = sum(ffxstntp);
    ffxX = zeros(ffxXn, tpredn);
    ffxstntpc = 1 + [0; cumsum(ffxstntp)];
    prstr = emptystruct({'Name1', 'Name2', 'RGB'}, [1, tpredn]);
    discarded = [];
    discardls = cell(1, numstudy);
    ffxw = [];
    if ~isempty(opts.ffxwread)
        ffxw = zeros(ffxXn, 1);
    end
    for stc = 1:numstudy
        ffxi1 = ffxstntpc(stc):(ffxstntpc(stc+1)-1);
        sdmsc = ffxstud{stc, end};
        xtcsc = ffxstud{stc, end-1};
        if isfield(xtcsc.C.RunTimeVars, 'Discard') && ...
           ~isempty(xtcsc.C.RunTimeVars.Discard) && ...
            isa(xtcsc.C.RunTimeVars.Discard, 'double') && ...
           ~any(isinf(xtcsc.C.RunTimeVars.Discard(:)) & isnan(xtcsc.C.RunTimeVars.Discard(:)))
            discardls{stc} = ...
                unique(min(ffxstntp(stc), max(1, round(xtcsc.C.RunTimeVars.Discard(:)))));
            discarded = [discarded(:); discardls{stc} + ffxstntpc(stc) - 1];
        end
        for pcc = 1:numel(preds{stc})
            ffxX(ffxi1, findfirst(strcmp(tpreds, preds{stc}{pcc}))) = sdmsc.C.SDMMatrix(:, pcc);
        end
        if ~isempty(ffxw)
            if ~isfield(xtcsc.C.RunTimeVars, 'FFXWeights') || ...
               ~isstruct(xtcsc.C.RunTimeVars.FFXWeights) || ...
               ~isfield(xtcsc.C.RunTimeVars.FFXWeights, opts.ffxwread)
                ffxw(ffxi1) = ones(numel(ffxi1), 1);
            else
                ffxw(ffxi1) = xtcsc.C.RunTimeVars.FFXWeights.(opts.ffxwread);
            end
        end
    end
    for pcc = 1:tpredn
        prstr(pcc).Name1 = sprintf('Predictor: %d', pcc);
        prstr(pcc).Name2 = tpreds{pcc};
        prstr(pcc).RGB = predcol.(makelabel(tpreds{pcc}));
    end

    % remove discarded time points
    undiscarded = true(size(ffxX, 1), 1);
    if ~isempty(discarded)
        ffxX(discarded, :) = [];
        glmc.NrOfTimePoints = size(ffxX, 1);
        undiscarded(discarded) = false;
        if ~isempty(ffxw)
            ffxw(discarded) =[];
        end
    end

    % removing empty regressors
    unusebets = all(ffxX == 0, 1);
    anusebets = any(unusebets);
    if ~isempty(pbar)
        pbar.Progress(0.25, 'Inverting FFX design matrix...');
    end
    glmc.DesignMatrix = ffxX;
    if sum(ffxX(:) ~= 0) < (0.05 * numel(ffxX)) || opts.sngtpool
        ffxX = sparse(ffxX);
    end
    if anusebets
        dousebets = ~unusebets;
        susebets = sum(dousebets);
        glmc.iXX = zeros(size(ffxX, 2));
        douseiXX = inv(ffxX(:, dousebets)' * ffxX(:, dousebets));
        glmc.iXX(dousebets, dousebets) = douseiXX;
    else
        glmc.iXX = inv(ffxX' * ffxX);
    end
    ixxt = (glmc.iXX ~= 0 & abs(glmc.iXX) < sqrt(eps));
    glmc.iXX(ixxt) = 0;
    if sum(glmc.iXX(:) ~= 0) < (0.1 * numel(glmc.iXX))
        glmc.iXX = sparse(glmc.iXX);
    end
    glmc.Predictor = prstr;
    glmc.NrOfPredictors = numel(prstr);

    % apply weights (for single w-OLS)
    if ~isempty(ffxw) && ~opts.robust
        glmc.RunTimeVars.FFXWeights = ffxw;
        sffxw = sqrt(ffxw);
        wffxX = ffxX .* sffxw(:, ones(1, size(ffxX, 2)));
        if anusebets
            wiXX = zeros(size(ffxX, 2));
            wiXX(dousebets, dousebets) = inv(wffxX(:, dousebets)' * wffxX(:, dousebets));
        else
            wiXX = inv(wffxX' * wffxX);
        end
        ixxt = (wiXX ~= 0 & abs(wiXX) < sqrt(eps));
        wiXX(ixxt) = 0;
        glmc.RunTimeVars.FFXiXX = wiXX;
    else
        wffxX = ffxX;
        wiXX = glmc.iXX;
    end

    % part 1 single trial code HERE (create iXX for trials)
    if opts.sngtpool

        % if less than 3 _00X condition
        stpn = {prstr.Name2};
        stpm = (~cellfun('isempty', regexp(stpn, '_T\d+$')));
        if sum(stpm) < 3

            % disable, as those two need to be separated anyway
            opts.sngtpool = false;
        else

            % find conditions that need to be "re-pooled"
            stpc = find(stpm);
            nstpc = numel(stpc);

            % and then, for each unique condition
            stpx = cell(numel(stpc), 3);
            if ~isempty(pbar)
                pbar.Progress(0.35, 'Inverting pooled single-trial FFX design matrices...');
            end
            for pc = 1:nstpc

                % combine SDM accordingly
                [stpx{pc, :}] = poolnonsingletrial(wffxX, stpn, stpc(pc));
                stpx{pc, 2} = stpx{pc, 3} * stpx{pc, 2};
            end
        end
    end
end

% auto-store RunTimeVars
glmc.RunTimeVars.AutoSave = true;

% transio access (only for RFX for now)
if opts.transio && opts.rfx

    % empty filename?
    glm.C = glmc;
    aft_SaveAs(glm, opts.outfile);

    % reload with transio
    delete(glm);
    glmtio = root_TransIOSize(xffroot, 'glm', 1e5);
    try
        glm = xff(opts.outfile);
        oglm = {glm};
        glmc = glm.C;
    catch xfferror
        root_TransIOSize(xffroot, 'glm', glmtio);
        clearxffobjects(oglm);
        clearxffobjects(rfobjs(:));
        rethrow(xfferror);
    end

    % reset transiosize
    root_TransIOSize(xffroot, 'glm', glmtio);
end

% new progress bar and variable for ShowDesign handle
if ~isempty(pbar)
    if closepbar
        closebar(pbar);
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 200, 640, 36]);
        xprogress(pbar, 'settitle', 'Computing multi-study GLM...');
        xprogress(pbar, 0, 'Regressing timecourses...', 'visible', 0, 1);
        opts.pbar = pbar;
    else
        pbar.Progress(0, 'Regressing timecourses...');
    end
end
showsdmf = [];

% calc depends of Separate predictors / RFX -> currently always 12!
switch (opts.seppred + 10 * opts.rfx)

    % FFX
    case {0, 1, 2}

        % force transio for FFX if exceeds 2GB
        if (8 * prod([outsize, tpredn])) > 2e9
            if isempty(opts.outfile)
                clearxffobjects(oglm);
                clearxffobjects(rfobjs(:));
                if ~isempty(pbar)
                    if closepbar
                        closebar(pbar);
                    else
                        pbar.Visible = pbarvis;
                    end
                end
                error('neuroelf:xff:badArgument', ...
                    'Large FFX computations require outfile to be set.');
            end
            opts.transio = true;
        end

        % compute TAL coords ?
        if hashdr
            [hcx, hcy, hcz] = ndgrid(1:outsize(1), 1:outsize(2), 1:outsize(3));
            if mixhdr
                glm.C = glmc;
                hcx = bvcoordconv([hcx(:), hcy(:), hcz(:)], 'bvc2tal', aft_BoundingBox(glm));
            else
                hcx = [hcx(:), hcy(:), hcz(:), ones(numel(hcx), 1)];
                hcx = hcx * hdrcf.Trf';
                hcx(:, 4) = [];
            end
        end

        % common data
        glmc.GLMData.MultipleRegressionR = single(zeros([outsize, 1]));
        glmc.GLMData.MCorrSS = single(zeros([outsize, 1]));
        ffxbetaf = '';
        ffxxytcf = '';
        if opts.transio
            try
                ffxbetaf = [tempname '.tio'];
                ffxbetas = transio(ffxbetaf, 'ieee-le', ...
                    'single', 0, [prod(outsize), tpredn], true);
                ffxxytcf = [tempname '.tio'];
                ffxxytcs = transio(ffxxytcf, 'ieee-le', ...
                    'single', 0, [prod(outsize), tpredn], true);
                glmc.GLMData.BetaMaps = ffxbetas;
                glmc.GLMData.XY = ffxxytcs;
            catch xfferror
                xfferr2 = xfferror;
                if ~isempty(ffxbetaf) && exist(ffxbetaf, 'file') > 0
                    try
                        delete(ffxbetaf);
                    catch xfferror
                        neuroelf_lasterr(xfferror);
                    end
                end
                if ~isempty(ffxxytcf) && exist(ffxxytcf, 'file') > 0
                    try
                        delete(ffxxytcf);
                    catch xfferror
                        neuroelf_lasterr(xfferror);
                    end
                end
                if ~isempty(pbar)
                    if closepbar
                        closebar(pbar);
                    else
                        pbar.Visible = pbarvis;
                    end
                end
                clearxffobjects(oglm);
                clearxffobjects(rfobjs(:));
                rethrow(xfferr2);
            end
        else
            glmc.GLMData.BetaMaps = repmat(single(0), [prod(outsize), tpredn]);
            glmc.GLMData.XY = repmat(single(0), [prod(outsize), tpredn]);
        end
        glmc.GLMData.TimeCourseMean = single(zeros([outsize, 1]));

        % compute how many timecourses to regress at any given time
        numtcs = prod(outsize);
        stptcs = min(4096, max(64, floor(1e9 / numel(ffxX))));
        idxtcs = [1:stptcs:numtcs, numtcs + 1];
        if isempty(opts.mask)
            tottcs = numtcs;
        else
            tottcs = sum(opts.mask(:));
        end
        stptcs = 0;

        % for XY computation, keep transpose of X
        tffxX = ffxX';
        if ~opts.robust && ~isempty(ffxw)
            wtffxX = wffxX';
        else
            wtffxX = tffxX;
        end

        % for robust
        if opts.robust

            % use lower maxiter flag
            robopts = struct('maxiter', opts.maxiter, 'pbar', pbar, 'prange', [0, 1]);
            stptct = 1 / tottcs;

            % and precompute h and perm
            if anusebets
                [stepc, robh, robperm] = qr(ffxX(:, dousebets), 0);
            else
                [stepc, robh, robperm] = qr(ffxX, 0);
            end
            robh = {robh, robperm};
        end

        % for t-maps
        if opts.tmaps
            ffxnval = glmc.NrOfTimePoints - glmc.NrOfPredictors;
            tmapc = eye(size(ffxX, 2));
        end

        % initialize progress
        if ~isempty(pbar)
            pbar.Progress(stptcs / tottcs, sprintf( ...
                'Regressing timecourses on %d-by-%d matrix (%d/%d done)...', ...
                size(ffxX, 1), size(ffxX, 2), stptcs, tottcs));
        end

        % iterate over steps
        for stepc = 1:(numel(idxtcs)-1)

            % generate indices to address
            tcidxs = idxtcs(stepc):(idxtcs(stepc+1)-1);
            tcarray = zeros(ffxXn, numel(tcidxs));
            maski = true(numel(tcidxs), 1);

            % pre-masking
            stptco = stptcs;
            if ~isempty(opts.mask)
                maski = maski & lsqueeze(opts.mask(tcidxs));
                stptcs = stptcs + sum(maski);
                if ~any(maski)
                    continue;
                end
            else
                stptcs = stptcs + sum(maski);
            end

            % depending on functional datatype
            switch (ptype)

                % VTC
                case {1}

                    % requires TAL coords
                    if hashdr
                        hcy = hcx(tcidxs(maski), :);
                    end

                    % iterare over studies
                    for stc = 1:numstudy

                        % get objects target time indices
                        ttidx = ffxstntpc(stc):(ffxstntpc(stc+1)-1);

                        % get object
                        xtcsc = ffxstud{stc, end - 1};

                        % depending on type
                        switch (lower(xtcsc.S.Extensions{1}(1)))
                            case {'h'}
                                for hvolc = 1:ffxstntp
                                    tcarray(ttidx(hvolc), maski) = ...
                                        aft_SampleData3D(xtcsc, hcy, struct('mapvol', hvolc))';
                                end
                            case {'v'}
                                tcarray(ttidx, maski) = xtcsc.C.VTCData(:, tcidxs(maski));
                        end
                    end

                % MTC
                case {2}

                    % iterare over studies
                    for stc = 1:numstudy

                        % get objects target time indices
                        ttidx = ffxstntpc(stc):(ffxstntpc(stc+1)-1);

                        % get object
                        xtcsc = ffxstud{stc, end - 1};

                        % without SSM/TSMs
                        if size(ffxstud, 2) == 2

                            % copy the data
                            tcarray(ttidx, maski) = xtcsc.C.MTCData(:, tcidxs(maski));

                        % with SSM/TSM objects
                        else

                            % get xSM content
                            xsmsc = ffxstud{stc, 1};

                            % SSM
                            if lower(xsmsc.Extensions{1}(1)) == 's'

                                % we require the unique indices first
                                [ssmui, ssmua, ssmub] = unique( ...
                                    xsmsc.C.SourceOfTarget(tcidxs(maski)));

                                % get unique required time courses
                                ssmut = xtcsc.C.MTCData(:, ssmui);

                                % apply SSM
                                tcarray(ttidx, maski) = ssmut(:, ssmub);

                            % TSM
                            else

                                % NOT YET IMPLEMENTED!!

                            end
                        end
                    end

            end

            % implicit masking
            if opts.ithresh > 0
                maski = maski & (mean(tcarray, 1) > opts.ithresh)';
                if ~any(maski)
                    if ~isempty(pbar)
                        pbar.Progress(stptcs / tottcs, sprintf( ...
                            'Regressing timecourses on %d-by-%d matrix (%d/%d done)...', ...
                            size(ffxX, 1), size(ffxX, 2), stptcs, tottcs));
                    end
                    continue;
                end
            end

            % transform
            if opts.trans == 'p'
                for stc = 1:numstudy
                    tcarray(ffxstntpc(stc):(ffxstntpc(stc+1)-1), maski) = ...
                        psctrans(tcarray(ffxstntpc(stc):(ffxstntpc(stc+1)-1), maski), 1);
                end
            elseif opts.trans == 'z'
                for stc = 1:numstudy
                    tcarray(ffxstntpc(stc):(ffxstntpc(stc+1)-1), maski) = ...
                        ztrans(tcarray(ffxstntpc(stc):(ffxstntpc(stc+1)-1), maski), 1);
                end
            end

            % compute regression stuff
            if opts.robust
                robopts.prange = stptct .* [stptco, stptcs];
                if anusebets
                    ffxbeta = zeros(sum(maski), size(ffxX, 2));
                    if isempty(discarded)
                        [ffxbeta(:, dousebets), robw] = fitrobustbisquare_img(ffxX(:, dousebets), ...
                            tcarray(:, maski)', [], [], robopts, robh, douseiXX, ffxw);
                        ffxbeta = ffxbeta';
                        % SINGLE-TRIAL CODE, PART 2
                        ptc = (ffxX * ffxbeta) .* robw' + (tcarray(:, maski) .* (1 - robw'));
                    else
                        if isempty(ffxw)
                            [ffxbeta(:, dousebets), robw] = fitrobustbisquare_img(ffxX(:, dousebets), ...
                                tcarray(undiscarded, maski)', [], [], robopts, robh, douseiXX);
                        else
                            [ffxbeta(:, dousebets), robw] = fitrobustbisquare_img(ffxX(:, dousebets), ...
                                tcarray(undiscarded, maski)', [], [], robopts, robh, douseiXX, ffxw);
                        end
                        ffxbeta = ffxbeta';
                        % SINGLE-TRIAL CODE, PART 2
                        ptc = (ffxX * ffxbeta) .* robw' + ...
                            (tcarray(undiscarded, maski) .* (1 - robw'));
                    end
                else
                    if isempty(discarded)
                        [ffxbeta, robw] = fitrobustbisquare_img(ffxX, ...
                            tcarray(:, maski)', [], [], robopts, robh, glmc.iXX, ffxw);
                        ffxbeta = ffxbeta';
                        % SINGLE-TRIAL CODE, PART 2
                        ptc = (ffxX * ffxbeta) .* robw' + (tcarray(:, maski) .* (1 - robw'));
                    else
                        if isempty(ffxw)
                            [ffxbeta, robw] = fitrobustbisquare_img(ffxX, ...
                                tcarray(undiscarded, maski)', [], [], robopts, robh, glmc.iXX);
                        else
                            [ffxbeta, robw] = fitrobustbisquare_img(ffxX, ...
                                tcarray(undiscarded, maski)', [], [], robopts, robh, glmc.iXX, ffxw);
                        end
                        ffxbeta = ffxbeta';
                        % SINGLE-TRIAL CODE, PART 2
                        ptc = (ffxX * ffxbeta) .* robw' + ...
                            (tcarray(undiscarded, maski) .* (1 - robw'));
                    end
                end
            else
                % weighting
                if isempty(ffxw)
                    wtcarray = tcarray;
                else
                    wtcarray = tcarray .* sffxw(:, ones(1, size(tcarray, 2)));
                end
                if isempty(discarded)
                    ffxbeta = wiXX * (wtffxX * wtcarray(:, maski));
                    if opts.sngtpool
                        for pc = 1:nstpc
                            sbetas = stpx{pc, 2} * wtcarray(:, maski);
                            ffxbeta(stpc(pc), :) = sbetas(1, :);
                        end
                    end
                else
                    ffxbeta = wiXX * (wtffxX * wtcarray(undiscarded, maski));
                    if opts.sngtpool
                        for pc = 1:nstpc
                            sbetas = stpx{pc, 2} * wtcarray(undiscarded, maski);
                            ffxbeta(stpc(pc), :) = sbetas(1, :);
                        end
                    end
                end
                ptc = ffxX * ffxbeta;
            end

            % additional computations
            if isempty(discarded)
                vartc = varc(tcarray(:, maski));
                ptcxy = tffxX * tcarray(:, maski);
            else
                vartc = varc(tcarray(undiscarded, maski));
                ptcxy = tffxX * tcarray(undiscarded, maski);
            end

            % then fill in elements
            ffxmrr = std(ptc) ./ sqrt(vartc);
            glmc.GLMData.MultipleRegressionR(tcidxs(maski)) = ffxmrr;
            ffxmcss = (glmc.NrOfTimePoints - 1) .* vartc;
            glmc.GLMData.MCorrSS(tcidxs(maski)) = ffxmcss;
            if opts.tmaps
                ffxse = sqrt((1 - (double(ffxmrr) .^ 2)) .* double(ffxmcss) / ffxnval);
                ffxse(ffxse == 0) = Inf;
                for pc = 1:size(ffxX, 2)
                    ffxbeta(pc, :) = glmtstat(tmapc(pc, :), ffxbeta', glmc.iXX, ffxse)';
                end
            end
            glmc.GLMData.BetaMaps(tcidxs(maski), :) = ffxbeta';
            glmc.GLMData.XY(tcidxs(maski), :) = ptcxy';
            if isempty(discarded)
                glmc.GLMData.TimeCourseMean(tcidxs(maski)) = ...
                    (1 / glmc.NrOfTimePoints) .* sum(tcarray(:, maski), 1)';
            else
                glmc.GLMData.TimeCourseMean(tcidxs(maski)) = ...
                    (1 / glmc.NrOfTimePoints) .* sum(tcarray(undiscarded, maski), 1);
            end

            % progress
            if ~isempty(pbar)
                pbar.Progress(stptcs / tottcs, sprintf( ...
                    'Regressing timecourses on %d-by-%d matrix (%d/%d done)...', ...
                    size(ffxX, 1), size(ffxX, 2), stptcs, tottcs));
            end
        end

        % reshaping ND elements
        glmc.GLMData.BetaMaps = reshape(glmc.GLMData.BetaMaps, [outsize, tpredn]);
        glmc.GLMData.XY = reshape(glmc.GLMData.XY, [outsize, tpredn]);

        % ensure iXX is not sparse
        if issparse(glmc.iXX)
            glmc.iXX = full(glmc.iXX);
        end

        % save output
        if ~isempty(opts.outfile)
            try
                glm.C = glmc;
                aft_SaveAs(glm, opts.outfile);
                aft_SaveRunTimeVars(glm);
                glmc = glm.C;
                if ~isempty(ffxbetaf)
                    delete(ffxbetaf);
                end
                if ~isempty(ffxxytcf)
                    delete(ffxxytcf);
                end
            catch xfferror
                warning('neuroelf:xff:saveFailed', ...
                    'Error saving GLM: %s.', xfferror.message);
                neuroelf_lasterr(xfferror);
            end
        end

    % RFX
    case {12}

        % if not loaded
        if isempty(olay)

            % global rfx map
            grfx = zeros(outsz3);
        end

        % options
        cbopt = struct( ...
            'mask',     opts.mask, ...
            'maxiter',  opts.maxiter, ...
            'pbar',     pbar, ...
            'regdiff',  opts.regdiff, ...
            'robust',   opts.robust, ...
            'singlew',  [], ...
            'sngtpool', opts.sngtpool, ...
            'tdim',     [], ...
            'thresh',   opts.ithresh, ...
            'tmaps',    opts.tmaps, ...
            'trans',    opts.trans, ...
            'tssm',     [], ...
            'writeres', opts.writeres);

        % iterate over subjects
        pstc = 1;
        for sc = 1:numsubs

            % target in GLMData
            glmtarg = findfirst(strcmpi(subjids{sc}, glmsubs));

            % get studies for subjects
            substud = find(strcmp(mdmsubst, subjids{sc}));

            % create beta maps and weighting with number of runs
            bm = zeros([outsz3, tpredn]);
            iXw = zeros(tpredn, 1);

            % weighting array for HRFboost
            if ~isempty(ndsti)
                hbwst = cell(tpredn, 1);
            end

            % iterate over studies
            for stc = 1:numel(substud)

                % study number
                sti = substud(stc);

                % pass in correct SSM/TSM
                if size(rfobjs, 2) > 2
                    cbopt.tssm = rfobjs{sti, 1};
                end

                % and weights
                if opts.wdvarsfd
                    cbopt.singlew = singlews{sti};
                end

                % progress
                if ~isempty(pbar)
                    pbar.Progress((pstc - 1) / numstudy, sprintf( ...
                        'Working on study %d/%d...', pstc, numstudy));
                    cbopt.prange = [pstc - 1, pstc] ./ numstudy;
                end

                % show design
                if ~isempty(opts.showsdms)
                    if ~isempty(showsdmf) && ishandle(showsdmf)
                        delete(showsdmf);
                    end
                    showsdmf = sdm_ShowDesign(rfobjs{sti, end}, struct('type', opts.showsdms));
                end

                % temporary interpolation is required
                if reqi(sti)
                    rfobjo = rfobjs(sti, end - 1);
                    try
                        if mixhdr
                            rfobjs{sti, end - 1} = importvtcfromanalyze( ...
                                rfobjo, opts.bbox, opts.res, opts.imeth);
                        else
                            rfobjs{sti, end - 1} = hdr_Dyn3DToVTC(rfobjo{1});
                        end
                    catch xfferror
                        clearxffobjects(oglm);
                        clearxffobjects(rfobjs(:));
                        if ~isempty(pbar)
                            if closepbar
                                closebar(pbar);
                            else
                                pbar.Visible = pbarvis;
                            end
                        end
                        error('neuroelf:xff:badObject', ...
                            'Error creating temporary VTC of HDR/NII object: %s.', ...
                            xfferror.message);
                    end
                end

                % computation done by SDM::CalcBetas
                xtcsc = rfobjs{sti, end - 1};
                if opts.writeres
                    [betas, iXX, bc, cbse, cbw, rsm] = ...
                        sdm_CalcBetas(rfobjs{sti, end}, xtcsc, cbopt);
                    xtcsci = find(strcmpi(stfname, xtcsc.F));
                    if numel(xtcsci) == 1
                        glmc.Study(xtcsci).RunTimeVars.FWHMResEst = rsm{1};
                        glmc.Study(xtcsci).RunTimeVars.FWHMResImg = single(rsm{2});
                    end
                else
                    [betas, iXX] = sdm_CalcBetas(rfobjs{sti, end}, xtcsc, cbopt);
                end

                % reshape betas to 4-D matrix, regardless
                betas = reshape(betas, [outsz3, size(betas, ndims(betas))]);

                % free temporary object and set original object back
                if reqi(sti)
                    delete(rfobjs{sti, end - 1});
                    rfobjs(sti, end - 1) = rfobjo;
                end

                % get weights correct
                iXX = 1 ./ diag(iXX);
                iXX(isinf(iXX) | isnan(iXX)) = 0;
                if ~opts.vweight
                    iXX = ststr(pstc).NrOfTimePoints .* double(iXX ~= 0);
                end

                % iterate over predictor names we need
                for bc = 1:tpredn

                    % find matching predictor for this study
                    predi = findfirst(strcmpi(opreds{pstc}, tpreds{bc}));

                    % if not found
                    if isempty(predi)

                        % not yet confound
                        if bc < tpredn

                            % HRFboost?
                            if ~isempty(ndsti) && ~isempty(ndsti{sti}{bc})

                                % replace basis functions
                                if opts.bftype == 's'
                                    sdmsc = rfobjs{sti, end};
                                    hrfbopt.bf = sdmsc.C.SDMMatrix(:, ndsti{sti}{bc});
                                end

                                % perform HRFboost
                                [hrfbm, hrfbw] = hrfboost( ...
                                    betas(:, :, :, ndsti{sti}{bc}), hrfbopt);

                                % make sure bad voxels are 0!
                                hrfbw(isinf(hrfbw) | isnan(hrfbw) | ...
                                    isinf(hrfbm) | isnan(hrfbm) | hrfbm == 0) = 0;
                                hrfbm(hrfbw == 0) = 0;

                                % add to betas
                                hrfbw = iXX(ndsti{sti}{bc}(1)) .* hrfbw;
                                bm(:, :, :, bc) = bm(:, :, :, bc) + hrfbw .* hrfbm;

                                % set or add to weighting average
                                if isempty(hbwst{bc})
                                    hbwst{bc} = hrfbw;
                                else
                                    hbwst{bc} = hbwst{bc} + hrfbw;
                                end
                            end

                            % don't do anything else
                            continue;
                        end

                        % for last, use confound
                        predi = size(betas, 4);
                    end

                    % add to betas and weighting
                    if iXX(predi) ~= 0
                        bm(:, :, :, bc) = bm(:, :, :, bc) + ...
                            iXX(predi) .* betas(:, :, :, predi);
                        iXw(bc) = iXw(bc) + iXX(predi);
                    end
                end

                % increase counter
                pstc = pstc + 1;
            end

            % HRFboost maps in use...
            if ~isempty(ndsti)
                for bc = 1:tpredn
                    if ~isempty(hbwst{bc})

                        % reweight
                        bm(:, :, :, bc) = bm(:, :, :, bc) ./ hbwst{bc};

                        % and no further reweighting!
                        iXw(bc) = 1;
                    end
                end
            end

            % mask betas
            bm(isinf(bm) | isnan(bm)) = 0;

            % sum and weight and put into output as single
            glmc.GLMData.Subject(glmtarg).BetaMaps(:, :, :, :) = ...
                reshape(single(repmat(reshape(1 ./ iXw, [1, 1, 1, tpredn]), ...
            	outsz3) .* bm), boutsz);

            % keep track of well behaved data
            grfx = grfx + any(bm ~= 0, 4);

            % update global RFX field
            glmc.GLMData.RFXGlobalMap = single(grfx >= 0.5 * sc + ogns);

            % for robust regression, try saving
            if opts.robust && ...
               ~istransio(glmc.GLMData.Subject(glmtarg).BetaMaps) && ~isempty(opts.outfile)
                try
                    glm.C = glmc;
                    aft_SaveAs(glm, opts.outfile);
                catch xfferror
                    fprintf('Couldn''t interim-save GLM file: %s.\n', xfferror.message);
                    neuroelf_lasterr(xfferror);
                end
            end
        end
end

% clean up
clearxffobjects(rfobjs(:));
if ~isempty(showsdmf) && ishandle(showsdmf)
    delete(showsdmf);
end

% update object
glm.C = glmc;
glmsubjids = glm_Subjects(glm);
gsmatch = multimatch(mdmsubs, glmsubjids);
msmatch = multimatch(glmsubjids, mdmsubs);

% Groups defined in mdm.RunTimeVars?
if isfield(rtv, 'Groups') && iscell(rtv.Groups) && ~isempty(rtv.Groups) && size(rtv.Groups, 2) == 2

    % make a matchkey between subject lists
    try
        for gc = 1:size(rtv.Groups, 1)
            rtv.Groups{gc, 2} = gsmatch(rtv.Groups{gc, 2});
            rtv.Groups{gc, 2}(rtv.Groups{gc, 2} < 1) = [];
        end
        glmc.RunTimeVars.Groups = rtv.Groups;
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end

% also keep other fields
bc = xo.C;
if isfield(rtv, 'Contrasts') && isfield(rtv, 'ContrastColors')
    glmc.RunTimeVars.Contrasts = rtv.Contrasts;
    glmc.RunTimeVars.ContrastColors = rtv.ContrastColors;
end
if isfield(rtv, 'CovariatesData') && isfield(rtv, 'CovariatesNames') && ...
    isa(rtv.CovariatesData, 'double') && iscell(rtv.CovariatesNames) && ...
    size(rtv.CovariatesData, 2) == numel(rtv.CovariatesNames) && ...
    all(cellfun(@ischar, rtv.CovariatesNames(:))) && ...
   ~any(cellfun('isempty', rtv.CovariatesNames(:))) && ...
    numel(mdmsubs) == size(rtv.CovariatesData, 1) && ...
    isfield(glmc.RunTimeVars, 'CovariatesData') && ...
    isempty(glmc.RunTimeVars.CovariatesData) && all(msmatch > 0)
    glmc.RunTimeVars.CovariatesData = rtv.CovariatesData(msmatch, :);
    glmc.RunTimeVars.CovariatesNames = rtv.CovariatesNames(:);
end
if isfield(opts, 'motpars') && iscell(opts.motpars) && numel(opts.motpars) == size(bc.XTC_RTC, 1)

    % match study filenames
    mstf = bc.XTC_RTC(:, end-1);
    gstf = {glmc.Study.NameOfAnalyzedFile};
    gstf = gstf(:);
    for bc = 1:numel(mstf)
        [nullp, mstf{bc}] = fileparts(mstf{bc});
    end
    for bc = 1:numel(gstf)
        [nullp, gstf{bc}] = fileparts(gstf{bc});
    end
    stmatch = multimatch(gstf, mstf);
    if all(stmatch > 0)
        glmc.RunTimeVars.MotionParameters = opts.motpars(stmatch);
    else
        glmc.RunTimeVars.MotionParameters = opts.motpars(:);
    end
    glmc.RunTimeVars.MotionParamsDiff = opts.motparsd;
    glmc.RunTimeVars.MotionParamsSquared = opts.motparsq;
end
if opts.robust
    glmc.RunTimeVars.RobustRegression = true;
else
    glmc.RunTimeVars.RobustRegression = false;
end
if ~isinf(opts.tfilter)
    glmc.RunTimeVars.TempFilterCutoff = opts.tfilter;
    glmc.RunTimeVars.TempFilterType = opts.tfilttype;
end

% add subject-specific SPMsn and TrfPlus fields
if ~isfield(glmc.RunTimeVars, 'SubjectSPMsn') || ...
   ~isstruct(glmc.RunTimeVars.SubjectSPMsn) || ...
    numel(glmc.RunTimeVars.SubjectSPMsn) ~= 1
    glmc.RunTimeVars.SubjectSPMsn = struct;
end
stslab = fieldnames(mdmsnmat);
for stc = 1:numel(stslab)
    glmc.RunTimeVars.SubjectSPMsn.(stslab{stc}) = mdmsnmat.(stslab{stc});
end
if ~isfield(glmc.RunTimeVars, 'SubjectTrfPlus') || ...
   ~isstruct(glmc.RunTimeVars.SubjectTrfPlus) || ...
    numel(glmc.RunTimeVars.SubjectTrfPlus) ~= 1
    glmc.RunTimeVars.SubjectTrfPlus = struct;
end
stslab = fieldnames(mdmtrfpl);
for stc = 1:numel(stslab)
    if ~isequal(mdmtrfpl.(stslab{stc}), eye(4))
        glmc.RunTimeVars.SubjectTrfPlus.(stslab{stc}) = mdmtrfpl.(stslab{stc});
    end
end

% update again
glm.C = glmc;

% close pbar
if ~isempty(pbar)
    if closepbar
        closebar(pbar);
    else
        pbar.Visible = pbarvis;
    end
end

% try to save?
if ~isempty(opts.outfile)
    try
        aft_SaveAs(glm, opts.outfile);
        aft_SaveRunTimeVars(glm);
    catch xfferror
        warning('neuroelf:xff:saveFailed', 'Error saving GLM: %s.', xfferror.message);
        neuroelf_lasterr(xfferror);
    end
end
