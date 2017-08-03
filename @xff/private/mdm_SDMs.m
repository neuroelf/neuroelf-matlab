function [sdms, stlist, sdmtr, tfiles, bfs, prtcs] = mdm_SDMs(xo, opts)
% MDM::SDMs  - create SDMs from an MDM file
%
% FORMAT:       [sdms, preds, sdmtr, tfiles, bfs, prtcs] = mdm.SDMs([options])
%
% Input fields:
%
%       options     optional 1x1 struct with fields
%        .asstruct  return structs instead of objects (default: false)
%        .globsigd  also add diff of global signals as nuisance regressors
%        .globsigf  filter global signals with confounds (default: true)
%        .globsigo  orthogonalize global signals (default: true)
%        .globsigs  add global signals as confound, one of
%                   0 - none
%                   1 - entire dataset (above threshold/within mask)
%                   2 - two (one per hemisphere, split at BV Z=128)
%                   3 or more, perform PCA of time courses and first N
%                   xff object(s), extract average time course from masks
%        .globsigup force-update of global signal estimates (default: false)
%        .motpars   motion parameters (Sx1 cell array with sdm/txt files)
%        .motparsd  also add diff of motion parameters (default: false)
%        .motparsq  also add squared motion parameters (default: false)
%        .ndcreg    if set > 0, perform deconvolution (only with PRTs!)
%        .orthconf  orthogonalize confounds (and motion parameters, true)
%        .partfft   partial FFT (cell array with list of conditions, {})
%        .partfftn  number of frequencies (default: 1)
%        .pbar      progress bar object (to show progress, default: [])
%        .pbrange   progress range (default 0 .. 1)
%        .pnames    provide names for parametric regressors (PRT only)
%        .ppicond   list of regressors (or differences) to interact
%        .ppirob    perform robust regression on VOI timecourse and remove
%                   outliers from timecourse/model (threshold, default: 0)
%        .ppitfilt  temporally filter PPI VOI timecourse (default: true)
%        .ppivoi    VOI object used to extract time-course from
%        .ppivoiidx intra-VOI-object index (default: 1)
%        .prtr      1x1 or Sx1 TR (in ms) for PRT::CreateSDM
%        .prtpnorm  normalize parameters of PRT.Conds (true)
%        .remodisis remodel ISIs using PRT::RemodelISIs function ({})
%        .restcond  remove rest condition (rest cond. name, default: '')
%        .savesdms  token, if not empty, save on-the-fly SDMs (e.g. '.sdm')
%        .shuflab   PRT labels (conditions names) to shuffle
%        .shuflabm  minimum number of onsets per label (1x1 or 1xL)
%        .sngtrial  single-trial SDMs
%        .sngtskip  condition list to skip during single-trial conversion
%        .tfilter   add filter regressors to SDMs (cut-off in secs)
%        .tfilttype temporal filter type (one of {'dct'}, 'fourier', 'poly')
%        .xconfound just as motpars, but without restriction on number
%
% Output fields:
%
%       sdms        1xS cell array with (new) SDM objects (or structs)
%       preds       list of predictor names
%       sdmtr       TR assumed for use with each SDM
%       tfiles      time-course file objects (transio access)
%       bfs         basis function set (if PRTs are used, empty otherwise)
%       prtcs       PRT contents (prior to PRT::CreateSDM calls)
%
% Note: all additional fields for the call to PRT::CreateSDM are supported!
%
% Using: meannoinfnan, multimatch, ne_fastica, orthvecs, psctrans,
%        singletrialprts, tffio, varc, ztrans.

% Version:  v1.1
% Build:    16082613
% Date:     Aug-26 2016, 1:15 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% global access to xffcont for RTC/SDM format description
global xffsngl;

% import functions from neuroelf library into workspace
using(neuroelf, {'bvcoordconv', 'findfirst', 'meannoinfnan', 'multimatch', 'ne_fastica', ...
    'orthvecs', 'psctrans', 'singletrialprts', 'tffio', 'varc', 'ztrans'});

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'mdm')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'asstruct') || ~islogical(opts.asstruct) || numel(opts.asstruct) ~= 1
    opts.asstruct = false;
end
if ~isfield(opts, 'pbar') || numel(opts.pbar) ~= 1 || ...
   (~isa(opts.pbar, 'xfigure') && ~isa(opts.pbar, 'xprogress'))
    pb = [];
else
    pb = opts.pbar;
end
if ~isfield(opts, 'pbrange') || ~isa(opts.pbrange, 'double') || numel(opts.pbrange) ~= 2 || ...
    any(isinf(opts.pbrange) | isnan(opts.pbrange) | opts.pbrange < 0)
    pbr = [0, 1];
else
    pbr = sort(opts.pbrange(:))';
end
pbm = pbr(1);
pbd = eps + pbr(2) - pbm;

% get object content (including format and handles)
bc = xo.C;
hc = xo.H;
bfs = [];

% check list of files
try
    if ~isfield(hc, 'FilesChecked') || ~islogical(hc.FilesChecked) || ~hc.FilesChecked
        if ~isempty(pb)
            pb.Progress(pbm, 'Checking file locations...');
        end
        mdm_CheckFiles(xo, struct('autofind', true, 'silent', true));
        bc = xo.C;
        hc = xo.H;
    end
catch xfferror
    rethrow(xfferror);
end

% get list of design files and of time course files and prepare output
rfobjs = bc.XTC_RTC;
numstudy = size(rfobjs, 1);
rfiles = rfobjs(:, end);
tfiles = rfobjs(:, end-1);
sdms = cell(size(rfiles));
xffroot = xff();

% check options
if ~isfield(opts, 'collapse') || ~iscell(opts.collapse) || ...
   ~any([2, 3] == size(opts.collapse, 2)) || ndims(opts.collapse) ~= 2 || isempty(opts.collapse)
    collapse = cell(0, 2);
else
    collapse = opts.collapse;
end
if ~isfield(opts, 'globsigd') || ~islogical(opts.globsigd) || numel(opts.globsigd) ~= 1
    opts.globsigd = false;
end
if ~isfield(opts, 'globsigf') || ~islogical(opts.globsigf) || numel(opts.globsigf) ~= 1
    opts.globsigf = true;
end
if ~isfield(opts, 'globsigo') || ~islogical(opts.globsigo) || numel(opts.globsigo) ~= 1
    opts.globsigo = true;
end
if ~isfield(opts, 'globsigs') || ((~isa(opts.globsigs, 'double') || ...
     numel(opts.globsigs) ~= 1 || isinf(opts.globsigs) || isnan(opts.globsigs) || opts.globsigs < 0) && ...
    (~iscell(opts.globsigs) || isempty(opts.globsigs)) && ...
    (~xffisobject(opts.globsigs, true, {'hdr', 'msk', 'vmr', 'voi'})))
    opts.globsigs = 0;
elseif isa(opts.globsigs, 'double')
    opts.globsigs = floor(opts.globsigs);
elseif numel(opts.globsigs) == 1 && xffisobject(opts.globsigs, true)
    opts.globsigs = {opts.globsigs};
else
    opts.globsigs = opts.globsigs(:)';
    if iscell(opts.globsigs)
        opts.globsigs(cellfun('isempty', opts.globsigs)) = [];
    end
end
if ~isfield(opts, 'globsigup') || ~islogical(opts.globsigup) || numel(opts.globsigup) ~= 1
    opts.globsigup = false;
end
if isfield(opts, 'motpars') && islogical(opts.motpars) && numel(opts.motpars) == 1 && opts.motpars
    if isfield(bc.RunTimeVars, 'MotionParameters') && ...
        iscell(bc.RunTimeVars.MotionParameters) && ...
        numel(bc.RunTimeVars.MotionParameters) == size(bc.XTC_RTC, 1)
        opts.motpars = bc.RunTimeVars.MotionParameters;
    else
        opts.motpars = repmat({'RTV'}, size(bc.XTC_RTC, 1), 1);
    end
end
if ~isfield(opts, 'motpars') || ~iscell(opts.motpars) || numel(opts.motpars) ~= size(bc.XTC_RTC, 1)
    opts.motpars = cell(size(bc.XTC_RTC, 1), 1);
else
    opts.motpars = opts.motpars(:);
end
if ~isfield(opts, 'motparsd') || ~islogical(opts.motparsd) || numel(opts.motparsd) ~= 1
    opts.motparsd = false;
end
if ~isfield(opts, 'motparsq') || ~islogical(opts.motparsq) || numel(opts.motparsq) ~= 1
    opts.motparsq = false;
end
if ~isfield(opts, 'ndcreg') || ~isa(opts.ndcreg, 'double') || numel(opts.ndcreg) ~= 1 || ...
    isinf(opts.ndcreg) || isnan(opts.ndcreg) || opts.ndcreg < 0 || opts.ndcreg > 120
    opts.ndcreg = 0;
else
    opts.ndcreg = floor(opts.ndcreg);
end
if opts.ndcreg > 0
    sdmtype = 'fir';
else
    sdmtype = 'hrf';
end
if isfield(hc, 'NrOfVolumes') && numel(hc.NrOfVolumes) == numstudy
    opts.nvol = hc.NrOfVolumes(:);
else
    opts.nvol = [];
end
if ~isfield(opts, 'orthconf') || ~islogical(opts.orthconf) || numel(opts.orthconf) ~= 1
    opts.orthconf = true;
end
if ~isfield(opts, 'partfft') || ~iscell(opts.partfft) || isempty(opts.partfft)
    opts.partfft = {};
else
    opts.partfft = opts.partfft(:);
    for ppc = numel(opts.partfft):-1:1
        if ~ischar(opts.partfft{ppc}) || isempty(opts.partfft{ppc})
            opts.partfft(ppc) = [];
        else
            opts.partfft{ppc} = opts.partfft{ppc}(:)';
        end
    end
    if ~isempty(opts.partfft)
        opts.partfft = unique(opts.partfft);
    end
end
if ~isfield(opts, 'partfftn') || ~isa(opts.partfftn, 'double') || ...
    numel(opts.partfftn) ~= 1 || isinf(opts.partfftn) || isnan(opts.partfftn) || opts.partfftn < 1
    opts.partfftn = 1;
else
    opts.partfftn = min(8, round(opts.partfftn));
end
if ~isfield(opts, 'pnames') || ~iscell(opts.pnames) || isempty(opts.pnames) || ...
   ~all(cellfun(@ischar, opts.pnames(:))) || any(cellfun('isempty', opts.pnames(:)))
    opts.pnames = {};
else
    opts.pnames = opts.pnames(:)';
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
if ~isfield(opts, 'ppitfilt') || ~islogical(opts.ppitfilt) || numel(opts.ppitfilt) ~= 1
    opts.ppitfilt = true;
end
if ~isfield(opts, 'ppivoi') || numel(opts.ppivoi) ~= 1 || ...
   (~xffisobject(opts.ppivoi, true, 'poi') && ~xffisobject(opts.ppivoi, true, 'voi'))
    opts.ppivoi = [];
end
if ~isfield(opts, 'ppivoiidx') || numel(opts.ppivoiidx) ~= 1 || ~isa(opts.ppivoiidx, 'double') || ...
    isinf(opts.ppivoiidx) || isnan(opts.ppivoiidx) || opts.ppivoiidx < 1
    opts.ppivoiidx = 1;
else
    opts.ppivoiidx = floor(opts.ppivoiidx);
end
if ~isa(opts.globsigs, 'double') || opts.globsigs ~= 0 || ...
   (~isempty(opts.ppicond) && ~isempty(opts.ppivoi)) || nargout > 3
    try
        tiosz = root_TransIOSize(xffroot, 1e5);
        for stc = 1:numstudy
            if ~isempty(pb)
                pb.Progress(pbm, sprintf('Accessing XTC %d...', stc));
            end
            tfiles{stc} = xff(tfiles{stc});
            if ~xffisobject(tfiles{stc}, true, {'hdr', 'mtc', 'vtc'})
                error('neuroelf:xff:badArgument', ...
                    'PPI/global signals currently only supported for VTCs.');
            end
            tfilec = tfiles{stc}.C;
            if ischar(opts.motpars{stc}) && strcmpi(opts.motpars{stc}, 'rtv')
                opts.motpars{stc} = [];
                if isfield(tfilec.RunTimeVars, 'MotionParameters') && ...
                   ~isempty(tfilec.RunTimeVars.MotionParameters) && ...
                    size(tfilec.RunTimeVars.MotionParameters, 2) == 6
                    opts.motpars{stc} = tfilec.RunTimeVars.MotionParameters;
                end
            end
        end
    catch xfferror
        root_TransIOSize(xffroot, tiosz);
        clearxffobjects(tfiles);
        rethrow(xfferror);
    end
    root_TransIOSize(xffroot, tiosz);
end
if ~isfield(opts, 'prtr') || ~isa(opts.prtr, 'double') || ...
   ~any([1, numstudy] == numel(opts.prtr)) || ...
    any(isinf(opts.prtr(:)) | isnan(opts.prtr(:)) | opts.prtr(:) <= 0)
    if isfield(hc, 'TR') && numel(hc.TR) == numstudy
        opts.prtr = hc.TR(:);
    else
        opts.prtr = [];
    end
else
    opts.prtr = opts.prtr(:);
    if numel(opts.prtr) == 1
        opts.prtr = opts.prtr .* ones(numstudy, 1);
    end
end
if ~isfield(opts, 'prtpnorm') || ~islogical(opts.prtpnorm) || numel(opts.prtpnorm) ~= 1
    opts.prtpnorm = true;
end
if ~isfield(opts, 'remodisis') || ~iscell(opts.remodisis)
    opts.remodisis = {};
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
if ~isfield(opts, 'savesdms') || ~ischar(opts.savesdms) || isempty(regexpi(opts.savesdms(:)', '\.sdm$'))
    opts.savesdms = '';
end
if ~isfield(opts, 'shuflab') || ~iscell(opts.shuflab)
    opts.shuflab = {};
else
    opts.shuflab = opts.shuflab(:)';
    for lc = numel(opts.shuflab):-1:1
        if ~ischar(opts.shuflab{lc}) || isempty(opts.shuflab{lc}) || (lc > 1 && ...
            any(strcmpi(opts.shuflab{lc}, opts.shuflab(1:lc-1))))
            opts.shuflab(lc) = [];
        else
            opts.shuftlab{lc} = opts.shuflab{lc}(:)';
        end
    end
end
if ~isfield(opts, 'shuflabm') || ~isa(opts.shuflabm, 'double') || ...
    any(isinf(opts.shuflabm(:)) | isnan(opts.shuflabm(:))) || ...
   ~any([1, numel(opts.shuflab)] == numel(opts.shuflabm))
    opts.shuflabm = 1;
end
if ~isfield(opts, 'sngtrial') || ~islogical(opts.sngtrial) || numel(opts.sngtrial) ~= 1
    opts.sngtrial = false;
elseif opts.sngtrial
    if opts.ndcreg > 0
        warning('neuroelf:xff:invalidOption', ...
            'Single-Trial SDMs cannot be combined with FIR modeling.');
        opts.ndcreg = 0;
        sdmtype = 'hrf';
    end
end
if ~isfield(opts, 'sngtskip') || ~iscell(opts.sngtskip)
    opts.sngtskip = {};
else
    opts.sngtskip = opts.sngtskip(:)';
end
if ~isfield(opts, 'tfilter') || ~isa(opts.tfilter, 'double') || ...
    numel(opts.tfilter) ~= 1 || isnan(opts.tfilter) || opts.tfilter < 60
    opts.tfilter = Inf;
end
if ~isfield(opts, 'tfilttype') || ~ischar(opts.tfilttype) || isempty(opts.tfilttype) || ...
   ~any(strcmpi(opts.tfilttype(:)', {'dct', 'fourier', 'poly'}))
    opts.tfilttype = 'dct';
else
    opts.tfilttype = lower(opts.tfilttype(:)');
end
if ~isfield(opts, 'xconfound') || ~iscell(opts.xconfound)
    opts.xconfound = {};
else
    opts.xconfound = opts.xconfound(:);
end

% no valid NrOfVolumes/TR given -> read info from files
if isempty(opts.nvol) || isempty(opts.prtr)
    try
        opts.nvol = zeros(numstudy, 1);
        if isempty(opts.prtr)
            opts.prtr = zeros(numstudy, 1);
        end
        for stc = 1:numel(opts.prtr)
            if ischar(tfiles{stc})
                if ~isempty(pb)
                    pb.Progress(pbm, sprintf('Accessing XTC %d...', stc));
                end
                trstr = xff(tfiles{stc}, 'h');
                trstr = trstr.C;
            else
                trstr = tfiles{stc}.C;
            end
            if isfield(trstr, 'NrOfVolumes')
                opts.nvol(stc) = trstr.NrOfVolumes;
            elseif isfield(trstr, 'NrOfTimePoints')
                opts.nvol(stc) = trstr.NrOfTimePoints;
            else
                opts.nvol(stc) = trstr.ImgDim.Dim(5);
            end
            if opts.prtr(stc) == 0
                if isfield(trstr, 'TR')
                    opts.prtr(stc) = trstr.TR;
                end
            end
        end
    catch xfferror
        rethrow(xfferror);
    end
end

% make sure motpars and xconfound are correctly sized
if numel(opts.motpars) ~= numstudy
    opts.motpars = cell(numstudy, 1);
end
if numel(opts.xconfound) ~= numstudy
    opts.xconfound = cell(numstudy, 1);
end
mpnb = {'trX',  'trY',  'trZ',  'rotX',  'rotY',  'rotZ' };
mpnd = {'trXd', 'trYd', 'trZd', 'rotXd', 'rotYd', 'rotZd'};
mpn2 = {'trX2', 'trY2', 'trZ2', 'rotX2', 'rotY2', 'rotZ2'};

% load objects (PRTs and/or SDMs)
try
    for stc = 1:numstudy
        if ~isempty(pb)
            pb.Progress(pbm, sprintf('Reading design for study %d...', stc));
        end
        sdms{stc} = xff(rfiles{stc});
        if ~xffisobject(sdms{stc}, true, {'prt', 'sdm'})
            error('neuroelf:xff:badArgument', ...
                'Invalid design file ''%s'' for study %d.', rfiles{stc}, stc);
        end
    end
catch xfferror
    clearxffobjects(tfiles);
    clearxffobjects(sdms);
    rethrow(xfferror);
end

% for single-trial GLMs
if opts.sngtrial

    % check that all files ARE PRTs, otherwise, error and end
    for stc = 1:numstudy
        if ~xffisobject(sdms{stc}, true, 'prt')
            clearxffobjects(tfiles);
            clearxffobjects(sdms);
            error('neuroelf:xff:badArgument', ...
                'Single-trial SDMs only possible with all-PRTs in XTC_RTC.');
        end
        if ~isempty(opts.remodisis)
            prt_RemodelISIs(sdms{stc}, opts.remodisis);
        end
        if ~isempty(collapse)
            for coc = 1:size(collapse, 1)
                prt_Collapse(sdms{stc}, collapse{coc, :});
            end
        end
    end

    % convert and get list
    try
        stlist = singletrialprts(sdms, mdm_Subjects(xo, true), opts.sngtskip);
    catch xfferror
        clearxffobjects(tfiles);
        clearxffobjects(sdms);
        rethrow(xfferror);
    end

% regular mode
else

    % if required, create list
    stlist = cell(0, 1);
    for stc = 1:numstudy
        if strcmpi(sdms{stc}.S.Extensions{1}, 'prt')
            if ~isempty(opts.remodisis)
                prt_RemodelISIs(sdms{stc}, opts.remodisis);
            end
            if ~isempty(collapse)
                for coc = 1:size(collapse, 1)
                    prt_Collapse(sdms{stc}, collapse{coc, :});
                end
            end
            condlist = cat(1, sdms{stc}.C.Cond.ConditionName);
        else
            if sdms{stc}.C.FirstConfoundPredictor > 0
                condlist = sdms{stc}.C.PredictorNames(1:sdms{stc}.C.FirstConfoundPredictor-1);
                condlist = condlist(:);
            else
                condlist = sdms{stc}.C.PredictorNames(:);
            end
        end
        condlist(multimatch(condlist, opts.restcond) > 0) = [];
        if isempty(stlist)
            stlist = condlist;
        else
            for clc = numel(condlist):-1:1
                if any(strcmp(condlist{clc}, stlist))
                    condlist(clc) = [];
                end
            end
            if ~isempty(condlist)
                stlist = cat(1, stlist, condlist);
            end
        end
    end

    % make sure constant is not in the list!
    stlist(strcmpi('constant', stlist)) = [];
end

% lookup coordinates for global time course masks
if iscell(opts.globsigs)

    % number of components
    if numel(opts.globsigs) > 1 && isa(opts.globsigs{end}, 'double') && ...
        numel(opts.globsigs{end}) == 1 && ~isinf(opts.globsigs{end}) && ...
       ~isnan(opts.globsigs{end}) && opts.globsigs{end} >= 1
        ngs = round(opts.globsigs{end});
        opts.globsigs(end) = [];
    else
        ngs = 1;
    end

    % only valid for VTC files
    if ~xffisobject(tfiles{1}, true, 'vtc')
        clearxffobjects(tfiles);
        clearxffobjects(sdms);
        error('neuroelf:xff:badArgument', 'Mask-based global signals only valid for VTCs.');
    end

    % for each object, find coordinates
    for oc = 1:numel(opts.globsigs)

        % get object's content
        gso = opts.globsigs{oc}.C;
        xtcc = tfiles{1}.C;
        if ~isempty(pb)
            if numel(gso) == 1 && xffisobject(opts.globsigs{oc}, true, {'hdr', 'head', 'msk', 'vmr'})
                gsf = opts.globsigs{oc}.F;
                if isempty(gsf)
                    gsf = sprintf('global signal object %d', oc);
                end
            elseif ischar(gso)
                for fc = 1:numel(tfiles)
                    if ~any(strcmpi(tfiles{fc}.C.RunTimeVars.GlobSigs(:, 3), [gso(:)' ',1']))
                        clearxffobjects(tfiles);
                        clearxffobjects(sdms(:));
                        error('neuroelf:xff:maskNotFound', ...
                            'Requested mask %s not found for %s.', gso(:)', bc.XTC_RTC{fc, 1});
                    end
                end
            else
                clearxffobjects(tfiles);
                clearxffobjects(sdms(:));
                error('neuroelf:xff:invalidOption', 'Invalid global signal option.');
            end
            pb.Progress(pbm, sprintf('Getting voxel indices for %s...', gsf));
        end

        % for VMRs
        if xffisobject(opts.globsigs{oc}, true, 'vmr')

            % find voxels > 0
            [gv1, gv2, gv3] = ind2sub(size(gso.VMRData), find(gso.VMRData(:, :, :) > 0));
            if gso.OffsetX ~= 0
                gv1 = gv1 + gso.OffsetX;
            end
            if gso.OffsetY ~= 0
                gv2 = gv2 + gso.OffsetY;
            end
            if gso.OffsetZ ~= 0
                gv3 = gv3 + gso.OffsetZ;
            end

        % for VOIs
        elseif xffisobject(opts.globsigs{oc}, true, 'voi')

            % allow SPMsn-content
            gv1 = unique(cat(1, gso.VOI.Voxels), 'rows');
            if isfield(xtcc.RunTimeVars, 'SPMsn') && isstruct(xtcc.RunTimeVars.SPMsn) && ...
                isfield(xtcc.RunTimeVars.SPMsn, 'VG') && isfield(xtcc.RunTimeVars.SPMsn, 'VF') && ...
                isfield(xtcc.RunTimeVars.SPMsn, 'Tr') && isfield(xtcc.RunTimeVars.SPMsn, 'Affine') && ...
               ~isempty(xtcc.RunTimeVars.SPMsn.VG) && ~isempty(xtcc.RunTimeVars.SPMsn.VF) && ...
               ~isempty(xtcc.RunTimeVars.SPMsn.Tr) && ~isempty(xtcc.RunTimeVars.SPMsn.Affine)
                opts.globsigs{oc} = gv1;
                continue;
            end
            
            % recompute into VTC voxels
            gv3 = 128 - gv1(:, 1);
            gv2 = 128 - gv1(:, 3);
            gv1 = 128 - gv1(:, 2);

        % for MSKs
        elseif xffisobject(opts.globsigs{oc}, true, 'msk')
            
            % allow SPMsn-content
            if isfield(xtcc.RunTimeVars, 'SPMsn') && isstruct(xtcc.RunTimeVars.SPMsn) && ...
                isfield(xtcc.RunTimeVars.SPMsn, 'VG') && isfield(xtcc.RunTimeVars.SPMsn, 'VF') && ...
                isfield(xtcc.RunTimeVars.SPMsn, 'Tr') && isfield(xtcc.RunTimeVars.SPMsn, 'Affine') && ...
               ~isempty(xtcc.RunTimeVars.SPMsn.VG) && ~isempty(xtcc.RunTimeVars.SPMsn.VF) && ...
               ~isempty(xtcc.RunTimeVars.SPMsn.Tr) && ~isempty(xtcc.RunTimeVars.SPMsn.Affine)
                opts.globsigs{oc} = bvcoordconv(find(gso.Mask(:)), 'bvx2tal', ...
                    aft_BoundingBox(opts.globsigs{oc}));
                continue;
            end

            % ensure compatibility
            if ~isequal(gso.Resolution, xtcc.Resolution) || ...
               ~isequal(gso.XStart, xtcc.XStart) || ~isequal(gso.YStart, xtcc.YStart) || ...
               ~isequal(gso.ZStart, xtcc.ZStart) || ~isequal(gso.XEnd, xtcc.XEnd) || ...
               ~isequal(gso.YEnd, xtcc.YEnd) || ~isequal(gso.ZEnd, xtcc.ZEnd)
                clearxffobjects(tfiles);
                clearxffobjects(sdms);
                error('neuroelf:xff:badArgument', 'Global signal MSK and VTC spatial mismatch.');
            end

            % find voxels
            opts.globsigs{oc} = find(gso.Mask(:));
            continue;

        % for HDR/NIIs
        else

            % find voxels > 0
            [gv1, gv2, gv3] = ind2sub(size(gso.VoxelData), find(gso.VoxelData(:, :, :, 1) >= 0.5));

            % get transformation matrix
            vtrf = hdr_CoordinateFrame(opts.globsigs{oc});
            vtrf = vtrf.Trf;

            % allow SPMsn-content
            gv1 = [gv1, gv2, gv3, ones(numel(gv1), 1)] * vtrf';
            gv1(:, 4) = [];
            if isfield(xtcc.RunTimeVars, 'SPMsn') && ...
                isstruct(xtcc.RunTimeVars.SPMsn) && ...
                isfield(xtcc.RunTimeVars.SPMsn, 'VG') && ...
                isfield(xtcc.RunTimeVars.SPMsn, 'VF') && ...
                isfield(xtcc.RunTimeVars.SPMsn, 'Tr') && ...
                isfield(xtcc.RunTimeVars.SPMsn, 'Affine') && ...
               ~isempty(xtcc.RunTimeVars.SPMsn.VG) && ...
               ~isempty(xtcc.RunTimeVars.SPMsn.VF) && ...
               ~isempty(xtcc.RunTimeVars.SPMsn.Tr) && ...
               ~isempty(xtcc.RunTimeVars.SPMsn.Affine)
                opts.globsigs{oc} = gv1;
                continue;
            end

            % compute voxel coordinates in BV space
            gv3 = 128 - gv1(:, 1);
            gv2 = 128 - gv1(:, 3);
            gv1 = 128 - gv1(:, 2);
        end

        % get offset, size, and resolution
        voff = [xtcc.XStart, xtcc.YStart, xtcc.ZStart];
        vsiz = size(xtcc.VTCData);
        vsiz(1) = [];
        vres = xtcc.Resolution;

        % compute coordinates to extract data from
        gv1 = unique(round(1 + (1 / vres) .* ([gv1, gv2, gv3] - ...
            repmat(voff, [numel(gv1), 1]))), 'rows');
        gv1(any(gv1 < 1, 2) | gv1(:, 1) > vsiz(1) | ...
            gv1(:, 2) > vsiz(2) | gv1(:, 3) > vsiz(3), :) = [];
        opts.globsigs{oc} = sub2ind(vsiz, gv1(:, 1), gv1(:, 2), gv1(:, 3));
    end
end

% get RTC/SDM format for access of motion parameter/confound files
sdmtff = xffsngl.FF.tff(xffsngl.EXT.sdm{2});

% copy of options for SDM::CreatePRT
psdmopts = opts;

% onsets, parameters
if nargout > 5
    prtcs = cell(numstudy, 1);
end

% work on each study
ppitc = [];
ppitf = {};
try
    for fc = 1:numstudy

        % progress bar
        if ~isempty(pb)
            pb.Progress(pbm + pbd * (fc / numstudy), ...
                sprintf('Processing design for study %d...', fc));
        end

        % set TR
        sttr = opts.prtr(fc);
        if sttr == 0 || opts.nvol(fc) == 0
            if ischar(tfiles{fc})
                tch = xff(tfiles{fc}, 'h');
                tch = tch.C;
            else
                tch = tfiles{fc}.C;
            end
        end
        if sttr == 0
            sttr = tch.TR;
        end
        if opts.nvol(fc) == 0
            if isfield(tch, 'NrOfVolumes')
                opts.nvol(fc) = tch.NrOfVolumes;
            elseif isfield(tch, 'NrOfTimePoints')
                opts.nvol(fc) = tch.NrOfTimePoints;
            elseif isfield(tch, 'NIIFileType')
                opts.nvol(fc) = tch.ImgDim.Dim(5);
            elseif isfield(tch, 'Brick')
                opts.nvol(fc) = numel(tch.Brick);
            end
        end

        % extract PPI VOI timecourse
        if ~isempty(opts.ppivoi)
            if size(rfobjs, 2) > 2
                ppitc = mtc_POITimeCourse(tfiles{fc}, opts.ppivoi);
            else
                ppitc = aft_VOITimeCourse(tfiles{fc}, opts.ppivoi);
            end
            stnumtp = size(ppitc, 1);

            % select time course
            if opts.ppivoiidx > size(ppitc, 2)
                ppitc = ppitc(:, 1);
            else
                ppitc = ppitc(:, opts.ppivoiidx);
            end

            % normalize (PSC!)
            ppitc = psctrans(ppitc);
            ppitc = ppitc - mean(ppitc);

            % invalid TC?
            if any(isinf(ppitc) | isnan(ppitc))

                % don't use for this study!
                ppitc = [];
                ppitf = {};

            % otherwise
            else

                % apply correct filter settings
                if ~isinf(opts.tfilter) && opts.ppitfilt
                    if opts.tfilttype(1) == 'd'
                        ppitf = {'tempdct', 1000 * opts.tfilter};
                    elseif opts.tfilttype(1) == 'f'
                        ppitf = {'tempsc', floor(0.001 * sttr * stnumtp / opts.tfilter)};
                    else
                        ppitf = {'temppoly', 2 + floor(0.002 * sttr * stnumtp / opts.tfilter)};
                    end
                end
            end
        end

        % create SDM from PRT if necessary
        prtsc = '';
        if strcmpi(sdms{fc}.S.Extensions{1}, 'prt')
            try
                prt = sdms{fc};
                prtsc = sdms{fc}.F;
                prtc = prt_ConditionNames(prt);
                if ~isempty(opts.restcond)
                    rcond = false(1, numel(prtc));
                    for rcc = 1:numel(prtc)
                        rcond(rcc) = any(strcmpi(prtc{rcc}, opts.restcond));
                    end
                    rcond = find(rcond);
                else
                    rcond = [];
                end
                psdmopts.nvol = opts.nvol(fc);
                psdmopts.pnorm = opts.prtpnorm;
                psdmopts.ppitc = ppitc;
                psdmopts.ppitf = ppitf;
                psdmopts.prtr = sttr;
                psdmopts.rcond = rcond;
                psdmopts.type = sdmtype;

                % shuffle condition labels
                if ~isempty(opts.shuflab)

                    % create copied PRT with shuffled labels
                    sprt = prt_ShuffleLabels(prt, opts.shuflab, opts.shuflabm);

                    % copy content to original PRT
                    prt.C = sprt.C;

                    % remove new PRT
                    delete(sprt);
                end

                % add parameter names to PRT
                if prt.C.ParametricWeights > 0 && numel(opts.pnames) >= prt.C.ParametricWeights
                    prt.C.RunTimeVars.ParameterNames = opts.pnames(1:prt.C.ParametricWeights);
                end

                % store PRT content
                fftprt = prt.C;
                if nargout > 5
                    prtcs{fc} = fftprt;
                end

                % create SDM from PRT
                [sdms{fc}, bfs] = prt_CreateSDM(prt, psdmopts);
                aft_ClearObject(prt);
                sdmc = sdms{fc}.C;
                
                % replace onsets with FFT-style list
                if ~isempty(opts.partfft)
                    
                    % replace which columns
                    fftcols = multimatch(opts.partfft, sdmc.PredictorNames(:));
                    fftcols = sort(fftcols(fftcols > 0));
                    
                    % get PRT condition names
                    fftprtcs = cat(1, fftprt.Cond.ConditionName);
                    
                    % iterate over replacement conditions
                    for rcc = numel(fftcols):-1:1
                        
                        % find PRT condition
                        fftprtci = find(strcmpi(fftprtcs, sdmc.PredictorNames{fftcols(rcc)}));
                        if isempty(fftprtci)
                            error('neuroelf:xff:conditionNotFound', ...
                                'Condition ''%s'' not found for FFT-replacement', ...
                                sdmc.PredictorNames{fftcols(rcc)});
                        end
                        
                        % get first onset and compute average difference
                        fftfo = fftprt.Cond(fftprtci(1)).OnOffsets;
                        fftcn = fftprt.Cond(fftprtci(1)).ConditionName{1};
                        fftfod = mean(diff(fftfo(:, 1))) / sttr;
                        fftfo = floor(fftfo(1) / sttr);
                        fftfl = size(sdmc.SDMMatrix, 1) - fftfo;
                        
                        % compute sin/cos of correct length/periods
                        fftscs = zeros(fftfl, 2 * opts.partfftn);
                        fftscn = cell(1, 2 * opts.partfftn);
                        for rccs = 1:opts.partfftn
                            fftscs(:, 2*rccs-1) = sin(0:(rccs*2*pi/fftfod):(fftfl-0.5)*(rccs*2*pi/fftfod))';
                            fftscs(:, 2*rccs) = cos(0:(rccs*2*pi/fftfod):(fftfl-0.5)*(rccs*2*pi/fftfod))';
                            fftscn{2*rccs-1} = sprintf('FFT_%s_sin%d', fftcn, rccs);
                            fftscn{2*rccs} = sprintf('FFT_%s_cos%d', fftcn, rccs);
                        end
                        
                        % replace
                        sdmc.NrOfPredictors = sdmc.NrOfPredictors + numel(fftscn) - 1;
                        sdmc.FirstConfoundPredictor = ...
                            sdmc.FirstConfoundPredictor + numel(fftscn) - 1;
                        sdmc.PredictorNames = [sdmc.PredictorNames(1:fftcols(rcc)-1), ...
                            fftscn, sdmc.PredictorNames(fftcols(rcc)+1:end)];
                        sdmc.PredictorColors = [ ...
                            sdmc.PredictorColors(1:fftcols(rcc)-1, :); ...
                            sdmc.PredictorColors(fftcols(rcc) .* ones(1, 2 * opts.partfftn), :); ...
                            sdmc.PredictorColors(fftcols(rcc)+1:end, :)];
                        sdmc.SDMMatrix = [sdmc.SDMMatrix(:, 1:fftcols(rcc)-1), ...
                            [zeros(fftfo, 2 * opts.partfftn); fftscs], ...
                            sdmc.SDMMatrix(:, fftcols(rcc)+1:end)];
                    end
                    sdms{fc}.C = sdmc;
                end
            catch xfferror
                error('neuroelf:xff:internalError', ...
                    'Error converting PRT to SDM: %s.', xfferror.message);
            end
        elseif ~strcmpi(sdms{fc}.S.Extensions{1}, 'sdm')
            error('neuroelf:xff:badArgument', ...
                'PRT or SDM needed to run study %d.', fc);
        end
        stnumtp = size(sdms{fc}.C.SDMMatrix, 1);
        ffreg = size(sdms{fc}.C.SDMMatrix, 2) + 1;

        % add filters
        fltmx = [];
        if ~isinf(opts.tfilter)
            if opts.tfilttype(1) == 'd'
                stfilt = floor(0.002 * sttr * stnumtp / opts.tfilter);
            elseif opts.tfilttype(1) == 'f'
                stfilt = floor(0.001 * sttr * stnumtp / opts.tfilter);
            else
                stfilt = 2 + floor(0.002 * sttr * stnumtp / opts.tfilter);
            end
            if stfilt > 0
                sdm_AddFilters(sdms{fc}, struct('ftype',  opts.tfilttype, 'number', stfilt));
                if opts.orthconf
                    fltmx = sdms{fc}.C.SDMMatrix(:, ffreg:end);
                    fltiv = pinv(fltmx' * fltmx) * fltmx';
                end
            end
        end

        % add motion parameters
        if ~isempty(opts.motpars{fc})
            mp = [];
            if ischar(opts.motpars{fc}) && numel(opts.motpars{fc}) > 5 && ...
                any(strcmpi(opts.motpars{fc}(end-3:end), {'.rtc', '.sdm'}))
                try
                    mpsdm = tffio(opts.motpars{fc}(:)', sdmtff);
                    mp = mpsdm.SDMMatrix;
                    if ~isequal(size(mp), [stnumtp, 6]) || any(isinf(mp(:)) | isnan(mp(:)))
                        mp = [];
                    end
                catch xfferror
                    fprintf('Motion parameters file %s not found/readable.\n', opts.motpars{fc}(:)');
                    neuroelf_lasterr(xfferror);
                    mp = [];
                end
            elseif ischar(opts.motpars{fc}) && numel(opts.motpars{fc}) > 5 && ...
                strcmpi(opts.motpars{fc}(end-3:end), '.txt')
                try
                    mp = load(opts.motpars{fc}(:)');
                    if ~isequal(size(mp), [stnumtp, 6]) || any(isinf(mp(:)) | isnan(mp(:)))
                        fprintf('Motion parameters file %s doesn''t match.\n', opts.motpars{fc}(:)');
                        mp = [];
                    end
                catch xfferror
                    fprintf('Motion parameters file %s not found/readable.\n', opts.motpars{fc}(:)');
                    neuroelf_lasterr(xfferror);
                    mp = [];
                end
            elseif isnumeric(opts.motpars{fc}) && isequal(size(opts.motpars{fc}), [stnumtp, 6])
                mp = double(opts.motpars{fc});
            end
            if ~isempty(mp)
                if opts.motparsd
                    dmp = [zeros(1, size(mp, 2)); diff(mp)];
                end
                if opts.motparsd && opts.motparsq
                    mp = [ztrans(mp), ztrans(dmp), ztrans(mp .* mp)];
                    mpn = [mpnb, mpnd, mpn2];
                    mpc = floor(255.999 * rand(18, 3));
                elseif opts.motparsd
                    mp = [ztrans(mp), ztrans(dmp)];
                    mpn = [mpnb, mpnd];
                    mpc = floor(255.999 * rand(12, 3));
                elseif opts.motparsq
                    mp = [ztrans(mp), ztrans(mp .* mp)];
                    mpn = [mpnb, mpn2];
                    mpc = floor(255.999 * rand(12, 3));
                else
                    mp = ztrans(mp);
                    mpc = floor(255.999 * rand(6, 3));
                    mpn = mpnb;
                end
                if opts.orthconf
                    if ~isempty(fltmx)
                        mp = mp - fltmx * (fltiv * mp);
                    end
                    mp = orthvecs(mp);
                end
                sdms{fc}.C.NrOfPredictors = sdms{fc}.C.NrOfPredictors + size(mp, 2);
                sdms{fc}.C.PredictorNames = [sdms{fc}.C.PredictorNames(:)', mpn];
                sdms{fc}.C.PredictorColors = [sdms{fc}.C.PredictorColors; mpc];
                sdms{fc}.C.SDMMatrix = [sdms{fc}.C.SDMMatrix, mp];
            end
        end

        % add additional confounds
        if ~isempty(opts.xconfound{fc})
            xc = [];
            if ischar(opts.xconfound{fc}) && numel(opts.xconfound{fc}) > 5 && ...
                any(strcmpi(opts.xconfound{fc}(end-3:end), {'.rtc', '.sdm'}))
                try
                    xcsdm = tffio(opts.xconfound{fc}(:)', sdmtff);
                    xc = xcsdm.SDMMatrix;
                    if ndims(xc) > 2 || size(xc, 1) ~= stnumtp
                        xc = [];
                    else
                        xc(:, any(isinf(xc) | isnan(xc))) = [];
                        xc(:, sum(abs(diff(xc))) == 0) = [];
                    end
                catch xfferror
                    neuroelf_lasterr(xfferror);
                    xc = [];
                end
            elseif ischar(opts.xconfound{fc}) && numel(opts.xconfound{fc}) > 5 && ...
                strcmpi(opts.xconfound{fc}(end-3:end), '.txt')
                try
                    xc = load(opts.xconfound{fc}(:)');
                    if ndims(xc) > 2 || size(xc, 1) ~= stnumtp || any(isinf(xc(:)) | isnan(xc(:)))
                        xc = [];
                    end
                catch xfferror
                    neuroelf_lasterr(xfferror);
                    xc = [];
                end
            end
            if ~isempty(xc)
                xc = ztrans(xc);
                xcc = floor(255.999 * rand(size(xc, 2), 3));
                xcn = cell(1, size(xc, 2));
                for xcnc = 1:numel(xcn)
                    xcn{xcnc} = sprintf('xc%03d', xcnc);
                end
                if opts.orthconf
                    if ~isempty(fltmx)
                        xc = xc - fltmx * (fltiv * xc);
                    end
                    if size(xc, 2) > 1
                        xc = orthvecs(xc);
                    end
                end
                sdms{fc}.C.NrOfPredictors = sdms{fc}.C.NrOfPredictors + size(xc, 2);
                sdms{fc}.C.PredictorNames = [sdms{fc}.C.PredictorNames(:)', xcn];
                sdms{fc}.C.PredictorColors = [sdms{fc}.C.PredictorColors; xcc];
                sdms{fc}.C.SDMMatrix = [sdms{fc}.C.SDMMatrix, xc];
            end
        end

        % add global signals
        glsigs = [];
        if iscell(opts.globsigs)
            glsigup = opts.globsigup(1, ones(1, numel(opts.globsigs)));
            xtcc = tfiles{fc}.C;
            glsigs = zeros(stnumtp, numel(opts.globsigs) * ngs);
            for oc = 1:numel(opts.globsigs)
                switch lower(tfiles{fc}.S.Extensions{1})
                    case {'fmr', 'mtc'}
                        error('neuroelf:xff:badArgument', ...
                            'Mask-based extraction of FMR not supported.');
                    case {'vtc'}
                        glsiguse = 0;
                        if isfield(xtcc.RunTimeVars, 'GlobSigs') && ...
                           ~isempty(xtcc.RunTimeVars.GlobSigs)
                            for occ = 1:size(xtcc.RunTimeVars.GlobSigs, 1)
                                if size(xtcc.RunTimeVars.GlobSigs, 2) < 3 || ...
                                    isempty(xtcc.RunTimeVars.GlobSigs{occ, 3})
                                    xtcc.RunTimeVars.GlobSigs{occ, 3} = ...
                                        sprintf('pca-%dvox,%d', ...
                                        numel(xtcc.RunTimeVars.GlobSigs{occ, 1}), ...
                                        size(xtcc.RunTimeVars.GlobSigs{occ, 2}, 2));
                                end
                                if ~opts.globsigup
                                    if ischar(opts.globsigs{oc})
                                        glsigposs = findfirst(strcmpi(xtcc.RunTimeVars.GlobSigs(:, 3), ...
                                            sprintf('%s,%d', opts.globsigs{oc}(:)', ngs)));
                                        if ~isempty(glsigposs)
                                            glsiguse = occ;
                                        end
                                    elseif isequal(xtcc.RunTimeVars.GlobSigs{occ, 1}, opts.globsigs{oc}) && ...
                                        size(xtcc.RunTimeVars.GlobSigs{occ, 2}, 2) == ngs
                                        glsiguse = occ;
                                        break;
                                    end
                                end
                            end
                        else
                            xtcc.RunTimeVars.GlobSigs = cell(0, 3);
                        end
                        if glsiguse > 0
                            glsiglist = xtcc.RunTimeVars.GlobSigs{glsiguse, 2};
                        else
                            glsigup(oc) = true;
                            if ischar(opts.globsigs{oc})
                                % % % to implement
                            elseif size(opts.globsigs{oc}, 2) > 1
                                vtcdsz = size(xtcc.VTCData);
                                [firstvol, glsigcrd] = aft_SampleData3D( ...
                                    tfiles{fc}, opts.globsigs{oc}, struct('snmat', true));
                                glsigcrd(any(glsigcrd <= 0.5, 2) | ...
                                    glsigcrd(:, 1) > vtcdsz(2) | ...
                                    glsigcrd(:, 2) > vtcdsz(3) | ...
                                    glsigcrd(:, 3) > vtcdsz(4), :) = [];
                                glsigcrd = 1 + (unique(round(glsigcrd), 'rows') - 1) * ...
                                    [1; vtcdsz(2); vtcdsz(2) * vtcdsz(3)];
                                if numel(glsigcrd) > 1000
                                    glsiglist = xtcc.VTCData(:, :, :, :);
                                    glsiglist = glsiglist(:, glsigcrd);
                                else
                                    glsiglist = xtcc.VTCData(:, glsigcrd);
                                end
                            elseif numel(opts.globsigs{oc}) > 1000
                                glsiglist = xtcc.VTCData(:, :, :, :);
                                glsiglist = glsiglist(:, opts.globsigs{oc});
                            else
                                glsiglist = xtcc.VTCData(:, opts.globsigs{oc});
                            end
                        end
                end
                if opts.globsigf
                    fltmx = [ztrans(sdms{fc}.C.SDMMatrix(:, ffreg:end)), ones(stnumtp, 1)];
                    fltmx(:, all(fltmx == 0)) = [];
                    glsiglist = glsiglist - fltmx * ((fltmx' * fltmx) \ fltmx' * glsiglist);
                end
                if ngs == 1
                    glsigs(:, oc) = meannoinfnan(glsiglist, 2);
                    if glsigup(oc)
                        xtcc.RunTimeVars.GlobSigs(end+1, :) = ...
                            {opts.globsigs{oc}, glsigs(:, oc), sprintf('mean-%dvox,1', numel(opts.globsigs{oc}))};
                    end
                elseif glsigup(oc)
                    glsiglist(:, any(isinf(glsiglist) | isnan(glsiglist), 1) | ...
                        all(diff(glsiglist, 1, 1) == 0)) = [];
                    glsigs(:, (oc*ngs):-1:(1+(oc-1)*ngs)) = ztrans(ne_fastica(double( ...
                        glsiglist - ones(size(glsiglist, 1), 1) * mean(glsiglist, 1)), ...
                        struct('step', 'pca', 'eign', ngs)));
                    if glsigup(oc)
                        xtcc.RunTimeVars.GlobSigs(end+1, :) = ...
                            {opts.globsigs{oc}, glsigs(:, (oc*ngs):-1:(1+(oc-1)*ngs)), ...
                            sprint('pca-%dvox,%d', numel(opts.globsigs{oc}), ngs)};
                    end
                else
                    glsigs(:, (oc*ngs):-1:(1+(oc-1)*ngs)) = glsiglist;
                end
            end
            if any(glsigup)
                tfiles{fc}.C = xtcc;
                try
                    aft_SaveRunTimeVars(tfiles{fc});
                catch xfferror
                    neuroelf_lasterr(xfferror);
                end
            end
            glsigs(:, varc(glsigs) < sqrt(eps)) = [];
        elseif opts.globsigs > 0
            xtcc = tfiles{fc}.C;
            switch lower(tfiles{fc}.S.Extensions{1})
                case {'fmr'}
                    if xtcc.DataStorageFormat == 1
                        globd = xtcc.Slice(1).STCData(:, :, :, :);
                        globd(1, 1, 1, numel(xtcc.Slice)) = globd(1, 1, 1, 1);
                        for slc = 2:size(globd, 4)
                            globd(:, :, :, slc) = xtcc.Slice(slc).STCData(:, :, :);
                        end
                    else
                        globd = xtcc.Slice(1).STCData(:, :, :, :);
                    end
                    globd = permute(globd, [3, 2, 4, 1]);
                    globd = globd(:, :);
                case {'mtc'}
                    globd = xtcc.MTCData(:, :);
                case {'vtc'}
                    globd = xtcc.VTCData(:, :, :, :);
                    globd = globd(:, :);
            end
            switch opts.globsigs

                % one global signal: mean
                case {1}

                    glsigs = double(mean(globd, 2));

                % two global signals: left/right (for MTC: just one)
                case {2}
                    if strcmpi(tfiles{fc}.S.Extensions{1}, 'mtc')
                        glsigs = double(mean(globd, 2));
                    else
                        glsigs = double([mean(globd(:, 1:floor(0.5 * size(globd, 2))), 2), ...
                            mean(globd(:, ceil(1 + 0.5 * size(globd, 2)):end), 2)]);
                    end

                % more than 2
                otherwise

                    % perform PCA
                    glsigs = ne_fastica(double(globd), struct('step', 'pca'));

                    % take last N
                    glsigs = glsigs(:, end:-1:max(1,size(glsigs,2)+1-opts.globsigs));
            end
        end

        % global signals exist
        if ~isempty(glsigs)

            % also add derivatives
            if opts.globsigd
                dglsigs = [2 .* glsigs(1, :) - glsigs(2, :); glsigs; 2 .* glsigs(end - 1, :) - glsigs(end, :)];
                glsigs = [glsigs, 0.5 .* (diff(dglsigs(1:end-1, :)) + diff(dglsigs(2:end, :)))];
            end

            % transform
            glsigs = glsigs - (ones(size(glsigs, 1), 1) * mean(glsigs));

            % orthogonalize
            if opts.globsigo && size(glsigs, 2) > 1
                glsigs = orthvecs(glsigs);
            end

            % melt-down to correlations of < 0.7
            if size(glsigs, 2) > 1
                glcorr = corrcoef(glsigs);
                glcorr(1:(size(glcorr, 1)+1):end) = 0;
                [glchi1, glchi2] = ind2sub(size(glcorr), find(abs(glcorr(:)) > sqrt(0.5)));
                while ~isempty(glchi2)
                    if glcorr(glchi1(1), glchi2(1)) > 0
                        glsigs(:, glchi1(1)) = 0.5 .* ...
                            (glsigs(:, glchi1(1)) + glsigs(:, glchi2(1)));
                    else
                        glsigs(:, glchi1(1)) = 0.5 .* ...
                            (glsigs(:, glchi1(1)) - glsigs(:, glchi2(1)));
                    end
                    glsigs(:, glchi1(2)) = [];
                    glcorr = corrcoef(glsigs);
                    glcorr(1:(size(glcorr, 1)+1):end) = 0;
                    [glchi1, glchi2] = ind2sub(size(glcorr), find(abs(glcorr(:)) > sqrt(0.5)));
                end
            end

            % add to SDMs
            xcn = cell(1, size(glsigs, 2));
            for xcnc = 1:numel(xcn)
                xcn{xcnc} = sprintf('gs%d', xcnc);
            end
            sdms{fc}.C.NrOfPredictors = sdms{fc}.C.NrOfPredictors + numel(xcn);
            sdms{fc}.C.PredictorColors = ...
                [sdms{fc}.C.PredictorColors; floor(255.999 .* rand(numel(xcn), 3))];
            sdms{fc}.C.PredictorNames = [sdms{fc}.C.PredictorNames, xcn];
            sdms{fc}.C.SDMMatrix = [sdms{fc}.C.SDMMatrix, glsigs];
        end

        % make sure the matrix has exactly one mean confound
        confounds = find(all(sdms{fc}.C.SDMMatrix == ...
            (ones(size(sdms{fc}.C.SDMMatrix, 1), 1) * sdms{fc}.C.SDMMatrix(1, :))) & ...
            sdms{fc}.C.SDMMatrix(1, :) ~= 0);

        % confound(s) found
        if ~isempty(confounds)

            % make sure all of them are confounds
            if any(confounds < sdms{fc}.C.FirstConfoundPredictor)
                error('neuroelf:xff:badArgument', ...
                    'SDM %d has mean/constant predictor before confounds.', fc);
            end

            % remove from design
            sdms{fc}.C.NrOfPredictors = sdms{fc}.C.NrOfPredictors - numel(confounds);
            sdms{fc}.C.PredictorColors(confounds, :) = [];
            sdms{fc}.C.PredictorNames(confounds) = [];
            sdms{fc}.C.SDMMatrix(:, confounds) = [];
        end

        % then add *at the end* of the matrix/names
        sdms{fc}.C.IncludesConstant = 1;
        sdms{fc}.C.PredictorColors(end+1, :) = [255, 255, 255];
        sdms{fc}.C.PredictorNames{end+1} = 'Constant';
        sdms{fc}.C.SDMMatrix(:, end+1) = 1;
        sdms{fc}.C.NrOfPredictors = size(sdms{fc}.C.SDMMatrix, 2);

        % remove all empty confound regressors
        emptyregs = find(all(sdms{fc}.C.SDMMatrix == 0));
        emptyregs(emptyregs < sdms{fc}.C.FirstConfoundPredictor) = [];
        if ~isempty(emptyregs)
            sdms{fc}.C.NrOfPredictors = sdms{fc}.C.NrOfPredictors - numel(emptyregs);
            sdms{fc}.C.PredictorColors(emptyregs, :) = [];
            sdms{fc}.C.PredictorNames(emptyregs) = [];
            sdms{fc}.C.SDMMatrix(:, emptyregs) = [];
        end

        % save SDM?
        if ~isempty(opts.savesdms) && ~isempty(prtsc) && ischar(prtsc)
            try
                aft_SaveAs(sdms{fc}, [prtsc(1:end-4) opts.savesdms]);
            catch xfferror
                neuroelf_lasterr(xfferror);
            end
        end
    end

    % get structures?
    if opts.asstruct
        sdmstr = sdms;
        for fc = 1:numstudy
            sdmstr{fc} = sdms{fc}.C;
        end
    end

    % return as structs
    if opts.asstruct
        clearxffobjects(sdms(:));
        sdms = sdmstr;
    end

    % return tr?
    if nargout > 2
        sdmtr = opts.prtr;
    end

% deal with errors
catch xfferror
    clearxffobjects(tfiles);
    clearxffobjects(sdms(:));
    rethrow(xfferror);
end
