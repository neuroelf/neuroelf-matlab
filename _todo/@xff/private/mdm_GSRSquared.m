function vmp = mdm_GSRSquared(xo, opts)
% MDM::WMRSquared  - compute one RSquared map for each run showing GS R2
%
% FORMAT:       glm = mdm.GSRSquared([options])
%
% Input fields:
%
%       options     optional 1x1 struct with fields
%        .discard   Sx1 cell array with volumes to discard (in addition)
%        .globsigd  also add diff of global signals as nuisance regressors
%        .globsigs  add global signals as confound, one of
%                   1 - entire dataset (above threshold/within mask)
%                   2 - two (one per hemisphere, split at BV Z=128)
%                   3 or more, perform PCA of time courses and first N
%                   xff object(s), extract average time course from masks
%                   default is to use the white matter mask matching
%                   the VTC resolution
%        .motpars   motion parameters (Sx1 cell array with sdm/txt files)
%        .motparsd  also add diff of motion parameters (default: false)
%        .motparsq  also add squared motion parameters (default: false)
%        .mpmaps    also create (partial) R2 maps for motion params (false)
%        .ndcreg    if set > 0, perform deconvolution (only with PRTs!)
%        .ppicond   list of regressors (or differences) to interact
%        .ppitfilt  temporally filter PPI VOI timecourse (default: true)
%        .ppivoi    VOI object used to extract time-course from
%        .ppivoiidx intra-VOI-object index (default: 1)
%        .restcond  remove rest condition (rest cond. name, default: '')
%        .sngtrial  single-trial GLM (only with PRTs, default: false)
%        .sngtskip  condition list to skip during single-trial conversion
%        .subsel    cell array with subject IDs to work on
%        .tfilter   add filter regressors to SDMs (cut-off in secs)
%        .tfilttype temporal filter type (one of {'dct'}, 'fourier', 'poly')
%        .tfmaps    also create (partial) R2 maps for temp filters (false)
%        .xconfound just as motpars, but without restriction on number
%
% Output fields:
%
%       vmp         VMP object with one map per VTC
%
% Note: only works with VTC as time course files
%       for .ppi to work, the model filenames must be PRTs!
%       all additional fields for the call to PRT::CreateSDM are supported!
%
% Using: modelcomp, newnatresvmp.

% Version:  v1.1
% Build:    16020917
% Date:     Feb-09 2016, 5:21 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'mdm')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
hc = xo.H;
bc = xo.C;
rtv = bc.RunTimeVars;

% only valid for VTC
if ~strcmpi(bc.TypeOfFunctionalData, 'vtc') || ...
    any(cellfun('isempty', regexpi(bc.XTC_RTC(:, 1), '\.vtc$')))
    error('neuroelf:xff:badArgument', 'Only valid for VTC-based MDM.');
end

% check list of files
cfs = struct('autofind', true, 'silent', true);
try
    if ~isfield(hc, 'FilesChecked') || ~islogical(hc.FilesChecked) || ~hc.FilesChecked
        mdm_CheckFiles(xo, cfs);
        bc = xo.C;
    end
catch xfferror
    rethrow(xfferror);
end

% load first VTC header
try
    vtchead = xff(bc.XTC_RTC{1, 1}, 'h');
    vtcres = vtchead.Resolution;
catch xfferror
    rethrow(xfferror);
end

% check options
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
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
mp = neuroelf_path('masks');
if isfield(opts, 'globsigs') && ischar(opts.globsigs) && ...
    any(strcmpi(opts.globsigs(:)', {'cb2', 'cb23', 'cb3', 'cb33'}))
    switch (lower(opts.globsigs(:)'))
        case 'cb2'
            if exist([mp '/colin_brain_ICBMnorm_brain2mm.msk'], 'file') == 2
                try
                    opts.globsigs = {xff([mp '/colin_brain_ICBMnorm_brain2mm.msk'])};
                catch xfferror
                    neuroelf_lasterr(xfferror);
                end
            end
        case 'cb23'
            if exist([mp '/colin_brain_ICBMnorm_gray2mm.msk'], 'file') == 2 && ...
                exist([mp '/colin_brain_ICBMnorm_white2mm.msk'], 'file') == 2 && ...
                exist([mp '/colin_brain_ICBMnorm_csfx2mm.msk'], 'file') == 2
                try
                    opts.globsigs = {xff([mp '/colin_brain_ICBMnorm_gray2mm.msk']); ...
                        xff([mp '/colin_brain_ICBMnorm_white2mm.msk']); ...
                        xff([mp '/colin_brain_ICBMnorm_csfx2mm.msk'])};
                catch xfferror
                    neuroelf_lasterr(xfferror);
                end
            end
        case 'cb3'
            if exist([mp '/colin_brain_ICBMnorm_brain3mm.msk'], 'file') == 2
                try
                    opts.globsigs = {xff([mp '/colin_brain_ICBMnorm_brain3mm.msk'])};
                catch xfferror
                    neuroelf_lasterr(xfferror);
                end
            end
        case 'cb33'
            if exist([mp '/colin_brain_ICBMnorm_gray3mm.msk'], 'file') == 2 && ...
                exist([mp '/colin_brain_ICBMnorm_white3mm.msk'], 'file') == 2 && ...
                exist([mp '/colin_brain_ICBMnorm_csfx3mm.msk'], 'file') == 2
                try
                    opts.globsigs = {xff([mp '/colin_brain_ICBMnorm_gray3mm.msk']); ...
                        xff([mp '/colin_brain_ICBMnorm_white3mm.msk']); ...
                        xff([mp '/colin_brain_ICBMnorm_csfx3mm.msk'])};
                catch xfferror
                    neuroelf_lasterr(xfferror);
                end
            end
    end
    if ischar(opts.globsigs)
        warning('neuroelf:xff:fileNorFound', 'Required mask file(s) not found.');
        opts.globsigs = [];
    end
else
    opts.globsigs = [];
end
if isfield(opts, 'motpars') && islogical(opts.motpars) && numel(opts.motpars) == 1 && ...
    opts.motpars && isfield(rtv, 'MotionParameters')
    opts.motpars = rtv.MotionParameters;
end
if ~isfield(opts, 'motparsd') || ~islogical(opts.motparsd) || numel(opts.motparsd) ~= 1
    opts.motparsd = false;
end
if ~isfield(opts, 'motparsq') || ~islogical(opts.motparsq) || numel(opts.motparsq) ~= 1
    opts.motparsq = false;
end
if ~isfield(opts, 'motpars') || ~iscell(opts.motpars) || isempty(opts.motpars)
    opts.mpmaps = false;
end
if ~isfield(opts, 'mpmaps') || ~islogical(opts.mpmaps) || numel(opts.mpmaps) ~= 1
    opts.mpmaps = false;
end
if ~isfield(opts, 'nderiv') || ~isa(opts.nderiv, 'double') || ...
    any(isinf(opts.nderiv(:)) | isnan(opts.nderiv(:))) || ...
   ~all(opts.nderiv(:) == 1 | opts.nderiv(:) == 2)
    opts.nderiv = [];
else
    opts.nderiv = unique(opts.nderiv(:))';
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
    opts.tfmaps = false;
end
if ~isfield(opts, 'tfilttype') || ~ischar(opts.tfilttype) || isempty(opts.tfilttype) || ...
   ~any(strcmpi(opts.tfilttype(:)', {'dct', 'fourier', 'poly'}))
    opts.tfilttype = 'dct';
else
    opts.tfilttype = lower(opts.tfilttype(:)');
end
if ~isfield(opts, 'tfmaps') || ~islogical(opts.tfmaps) || numel(opts.tfmaps) ~= 1
    opts.tfmaps = false;
end
opts.trans = 'n';

% get and check subject IDs
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

% determine progress bar capabilities
try
    pbar = xprogress;
    xprogress(pbar, 'setposition', [80, 200, 640, 36]);
    xprogress(pbar, 'settitle', 'Computing global signal R2 maps...');
    xprogress(pbar, 0, 'Checking referenced files...', 'visible', 0, 1);
catch xfferror
    pbar = [];
    neuroelf_lasterr(xfferror);
end
opts.pbar = pbar;

% determine global signals (late)
mskfile = {[]};
if isempty(opts.globsigs) || (iscell(opts.globsigs) && isempty(opts.globsigs{1}))
    try
        mskfile{1} = xff(sprintf('%s/colin_brain_ICBMnorm_white%dmm.msk', mp, vtcres));
        opts.globsigs = mskfile;
    catch xfferror
        if ~isempty(pbar)
            closebar(pbar);
        end
        clearxffobjects(mskfile);
        rethrow(xfferror);
    end
end

% then load/create SDMs
try
    [rfobjs(:, end), stlist, sdmtr, rfobjs(:, end-1), bfs] = mdm_SDMs(xo, opts);
catch xfferror
    if ~isempty(pbar)
        closebar(pbar);
    end
    clearxffobjects(mskfile);
    rethrow(xfferror);
end
clearxffobjects(mskfile);

% create VMP
vmp = ne_methods.newnatresvmp();
vmpc = vmp.C;
vmpc.Resolution = vtcres;
vmpc.XStart = vtchead.XStart;
vmpc.XEnd = vtchead.XEnd;
vmpc.YStart = vtchead.YStart;
vmpc.YEnd = vtchead.YEnd;
vmpc.ZStart = vtchead.ZStart;
vmpc.ZEnd = vtchead.ZEnd;
nummaps = size(rfobjs, 1) * (1 + double(opts.mpmaps) + double(opts.tfmaps));
vmpc.Map.Type = 2;
vmpc.Map.LowerThreshold = 0.25;
vmpc.Map.UpperThreshold = 1;
vmpc.Map = vmpc.Map(ones(nummaps, 1));
mti = 1;

% iterate over time course files
modelcomp = ne_methods.modelcomp;
for sc = 1:size(rfobjs, 1)

    % progress
    if ~isempty(pbar)

    end

    % load VTC data
    vtcc = rfobjs{sc, 1}.C;
    vtcd = vtcc.VTCData(:, :, :, :);

    % get SDM
    sdmc = rfobjs{sc, 2}.C;

    % find global signal regressor(s)
    gsi = find(~cellfun('isempty', regexpi(sdmc.PredictorNames, '^gs\dd?$')));

    % compute F maps
    rmi = 1:size(sdmc.SDMMatrix, 2);
    rmi(gsi) = [];
    [f, df1, df2] = modelcomp(sdmc.SDMMatrix, sdmc.SDMMatrix(:, rmi), 1);
    f(isinf(f) | isnan(f)) = 0;

    % convert to partial Eta2 (R2)

    % store in VMP
    vmpc.Map(mti).VMPData = single(f);

    % increase counter
    mti = mti + 1;

    % also compute for motion parameters/temporal filters
end

% clear object
clearxffobjects(rfobjs(:));

% set in object
vmp.C = vmpc;
