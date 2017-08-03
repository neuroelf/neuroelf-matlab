function [jobs, jobhelps, funcfiles, rff, rfs, rf] = spmx_preprojobs(sf, fp, ap, opts)
% spmx_preprojobs  - create jobs and jobhelps from options
%
% FORMAT:       [jobs, jobhelps, ff] = spmx_preprojobs(sf, fp, ap, opts)
% FORMAT:       spmx_preprojobs;
%
% Input fields:
%
%       sf          subject (base) folder
%       fp          functional folder pattern (e.g. f*/r*)
%       ap          anatomical folder pattern (e.g. a*/*)
%       opts        1x1 struct with options
%        .afp       anatomical file pattern (default: '*.img' / '*.nii')
%        .extbrain  extract brain from anatomical images (default: true)
%        .ffp       anatomical file pattern (default: '*.img' / '*.nii')
%        .fsingle   force images to single precision (default: true)
%        .fun2str   boolean flag, coregister func2struct (default: true)
%        .jobfile   save output into a mat file (full or partial filename)
%        .jobrun    flag, run jobs?
%        .normsess  normalize sessions separately (only epi, default: false)
%        .normtype  normalization, one of 'anat', 'epi', 'dartel', {'seg'}
%        .reganat   register anatomical scan to T1 template (default: true)
%        .reslice   re-slice images are realignment (default: true)
%        .rinterpe  realignment estimation spline degree (default: 2)
%        .rinterpr  realignment reslicing spline degree (default: 4)
%        .rqual     realignment quality (accuracy trade-off default: 0.9)
%        .rrtm      realignment register-to-mean flag (default: true)
%        .rsep      realignment sampling grid resolution (default: 3mm)
%        .skip      skip N volumes at beginning of each run (default: 0)
%        .smk       smoothing kernel (either 1x1 or 1x3 double, default: 6)
%        .spmv      which SPM version (if empty, auto detect)
%        .st        flag, perform slice time correction (default: true)
%        .sto       slice acquisition order, either of
%                   'a', 'aie', {'aio'}, 'd', 'die', 'dio' or numeric list
%        .str       reference slice (default: 1)
%        .stta      TA (default: auto)
%        .sttr      TR (default: 2 sec)
%        .t1temp    T1 template filename (if empty, use from SPM templates)
%        .wbb       normalize writing bounding box
%                   (default: [-84, -114, -66; 84, 84, 90])
%        .winterp   normalize writing interpolation spline (default: 1)
%        .wsna      write spatially-normalize anatomical (default: true)
%        .wvox      normalize writing voxel size (default: [3, 3, 3])
%
% Output fields:
%
%       jobs        cell array with jobs suitable for spm_jobman
%       jobhelps    cell array with helps (required when saved as jobs.mat)
%       ff          cell array with functional files (after preprocessing)
%
% Note: without inputs, the UI-based mode will be used.
%
% This function loads a pre-compiled list of jost from
%
% NEUROELF_PATH/_files/spm/spmX_prepro.mat
%
% which contains the following jobs
%
% - temporal.st (slice-timing; which is dropped if not used)
% - spatial.realign (estwrite sub-type)
% - util.imcalc (to manually compute a mean image)
% - spatial.coreg (estimate coregistration between mean-func and anat)
% - spatial.preproc (segmentation)
% - spatial.normalise (write, functional, using segmentation sn_mat)
% - spatial.smooth
%
% alternatively, the last two will be replace by set of DARTEL jobs
%
% and optionally, the anatomical is spatially normalised with an extra job

% Version:  v1.1
% Build:    16060816
% Date:     Jun-08 2016, 4:30 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% global UI variable
global ne_ui;

% check options
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'afp') || (~ischar(opts.afp) && ~iscell(opts.afp)) || isempty(opts.afp)
    opts.afp = {'*.img', '*.nii'};
elseif iscell(opts.afp)
    opts.afp = opts.afp(:)';
else
    opts.afp = {opts.afp(:)'};
end
for fpc = numel(opts.afp):-1:1
    if ~ischar(opts.afp{fpc}) || isempty(opts.afp{fpc}) || ~any(opts.afp{fpc}(:) == '*')
        opts.afp(fpc) = [];
    else
        opts.afp{fpc} = opts.afp{fpc}(:)';
    end
end
if isempty(opts.afp)
    opts.afp = {'*.img', '*.nii'};
end
if ~isfield(opts, 'extbrain') || ~islogical(opts.extbrain) || numel(opts.extbrain) ~= 1
    opts.extbrain = true;
end
if ~isfield(opts, 'ffp') || (~ischar(opts.ffp) && ~iscell(opts.ffp)) || isempty(opts.ffp)
    opts.ffp = {'*.img', '*.nii'};
elseif iscell(opts.ffp)
    opts.ffp = opts.ffp(:)';
else
    opts.ffp = {opts.ffp(:)'};
end
for fpc = numel(opts.ffp):-1:1
    if ~ischar(opts.ffp{fpc}) || isempty(opts.ffp{fpc}) || ~any(opts.ffp{fpc}(:) == '*')
        opts.ffp(fpc) = [];
    else
        opts.ffp{fpc} = opts.ffp{fpc}(:)';
    end
end
if isempty(opts.ffp)
    opts.ffp = {'*.img', '*.nii'};
end
if ~isfield(opts, 'fsingle') || ~islogical(opts.fsingle) || numel(opts.fsingle) ~= 1
    opts.fsingle = true;
end
if ~isfield(opts, 'fun2str') || numel(opts.fun2str) ~= 1 || ~islogical(opts.fun2str)
    opts.fun2str = true;
end
if ~isfield(opts, 'jobfile') || isempty(opts.jobfile) || ~ischar(opts.jobfile)
    opts.jobfile = '';
else
    opts.jobfile = opts.jobfile(:)';
end
if ~isfield(opts, 'jobrun') || numel(opts.jobrun) ~= 1 || ...
   (~isa(opts.jobrun, 'double') && ~islogical(opts.jobrun))
    opts.jobrun = false;
end
if ~isfield(opts, 'normsess') || ~islogical(opts.normsess) || numel(opts.normsess) ~= 1
    opts.normsess = false;
end
if ~isfield(opts, 'normtype') || ~ischar(opts.normtype) || isempty(opts.normtype) || ...
   ~any(lower(opts.normtype(1)) == 'ades')
    opts.normtype = 's';
else
    opts.normtype = lower(opts.normtype(1));
end
if opts.normtype ~= 'e'
    opts.normsess = false;
else
    opts.reganat = false;
end
if ~isfield(opts, 'reganat') || ~islogical(opts.reganat) || numel(opts.reganat) ~= 1
    opts.reganat = true;
end
if ~isfield(opts, 'reslice') || ~islogical(opts.reslice) || numel(opts.reslice) ~= 1
    opts.reslice = true;
end
if ~isfield(opts, 'rinterpe') || ~isa(opts.rinterpe, 'double') || ...
    numel(opts.rinterpe) ~= 1 || isinf(opts.rinterpe) || isnan(opts.rinterpe) || ...
   ~any((1:7) == opts.rinterpe)
    opts.rinterpe = 2;
end
if ~isfield(opts, 'rinterpr') || ~isa(opts.rinterpr, 'double') || ...
    numel(opts.rinterpr) ~= 1 || isinf(opts.rinterpr) || isnan(opts.rinterpr) || ...
   ~any((1:7) == opts.rinterpr)
    opts.rinterpr = 4;
end
if ~isfield(opts, 'rqual') || ~isa(opts.rqual, 'double') || numel(opts.rqual) ~= 1 || ...
    isinf(opts.rqual) || isnan(opts.rqual) || opts.rqual <= 0 || opts.rqual > 1
    opts.rqual = 0.9;
end
if ~isfield(opts, 'rrtm') || numel(opts.rrtm) ~= 1 || ...
   (~isa(opts.rrtm, 'double') && ~islogical(opts.rrtm))
    opts.rrtm = true;
elseif isa(opts.rrtm, 'double')
    opts.rrtm = (opts.rrtm ~= 0);
end
opts.rrtm = 0 + opts.rrtm;
if ~isfield(opts, 'rsep') || ~isa(opts.rsep, 'double') || numel(opts.rsep) ~= 1 || ...
    isinf(opts.rsep) || isnan(opts.rsep) || opts.rsep < 1 || opts.rsep > 8
    opts.rsep = 3;
end
if ~isfield(opts, 'skip') || ~isa(opts.skip, 'double') || numel(opts.skip) ~= 1 || ...
    isinf(opts.skip) || isnan(opts.skip) || opts.skip < 0
    opts.skip = 0;
else
    opts.skip = round(opts.skip);
end
if ~isfield(opts, 'smk') || ~isa(opts.smk, 'double') || ~any([1, 3] == numel(opts.smk)) || ...
    any(isinf(opts.smk) | isnan(opts.smk) | opts.smk < 0 | opts.smk > 16)
    opts.smk = 6;
end
if numel(opts.smk) == 1
    opts.smk = opts.smk([1, 1, 1]);
end
if ~isfield(opts, 'spmv') || ~isa(opts.spmv, 'double') || numel(opts.spmv) ~= 1 || ...
    isinf(opts.spmv) || isnan(opts.spmv) || ~any([5, 8, 12] == opts.spmv)
    opts.spmv = [];
end
if ~isfield(opts, 'st') || numel(opts.st) ~= 1 || ...
   (~isa(opts.st, 'double') && ~islogical(opts.st))
    opts.st = true;
end
if isa(opts.st, 'double')
    opts.st = (opts.st ~= 0);
end
if ~isfield(opts, 'sto') || (~ischar(opts.sto) && ~isa(opts.sto, 'double'))
    opts.sto = 'aio';
end
if ischar(opts.sto)
    if ~any(strcmpi(opts.sto(:)', {'a', 'aie', 'aio', 'd', 'die', 'dio'}))
        opts.sto = 'aio';
    else
        opts.sto = lower(opts.sto(:)');
    end
else
    opts.sto = opts.sto(:)';
    if any(isinf(opts.sto) | isnan(opts.sto) | opts.sto < 1 | opts.sto ~= fix(opts.sto)) || ...
        numel(opts.sto) ~= numel(unique(opts.sto))
        warning('neuroelf:general:badArgument', ...
            'Invalid slice order argument given.');
        opts.sto = 'aio';
    end
end
if ~isfield(opts, 'str') || ~isa(opts.str, 'double') || numel(opts.str) ~= 1 || ...
    isinf(opts.str) || isnan(opts.str) || opts.str < 1 || opts.str ~= fix(opts.str) || ...
   (isa(opts.sto, 'double') && opts.str > max(opts.sto))
    opts.str = 1;
end
if ~isfield(opts, 'stta') || ~isa(opts.stta, 'double') || numel(opts.stta) ~= 1 || ...
    isinf(opts.stta) || isnan(opts.stta) || opts.stta <= 0
    opts.stta = [];
end
if ~isfield(opts, 'sttr') || ~isa(opts.sttr, 'double') || numel(opts.sttr) ~= 1 || ...
    isinf(opts.sttr) || isnan(opts.sttr) || opts.sttr <= 0
    opts.sttr = 2;
end
if ~isempty(opts.stta) && opts.stta >= opts.sttr
    opts.stta = [];
end
if ~isfield(opts, 't1temp') || ~ischar(opts.t1temp) || exist(opts.t1temp(:)', 'file') ~= 2
    opts.t1temp = '';
end
if ~isfield(opts, 'wbb') || ~isa(opts.wbb, 'double') || ~isequal(size(opts.wbb), [2, 3]) || ...
    any(isinf(opts.wbb(:)) | isnan(opts.wbb(:))) || any(diff(opts.wbb) < 0)
    opts.wbb = [-84, -114, -66; 84, 84, 90];
end
if ~isfield(opts, 'winterp') || ~isa(opts.winterp, 'double') || ...
    numel(opts.winterp) ~= 1 || ~any((1:4) == opts.winterp)
    opts.winterp = 1;
end
if ~isfield(opts, 'wsna') || ~islogical(opts.wsna) || numel(opts.wsna) ~= 1
    opts.wsna = true;
end
if ~isfield(opts, 'wvox') || ~isa(opts.wvox, 'double') || ~any([1, 3] == numel(opts.wvox)) || ...
    any(isinf(opts.wvox) | isnan(opts.wvox) | opts.wvox < 0.1 | opts.wvox > 16)
    opts.wvox = 3;
end
if numel(opts.wvox) == 1
    opts.wvox = opts.wvox([1, 1, 1]);
end

% preset runflag
rf = false;

% check SPM is available
try
    spmver = lower(spm('ver'));
    if ~any(strcmpi(spmver, {'spm5', 'spm8', 'spm12b', 'spm12'}))
        error('WRONG_SPM_VERSION');
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error('neuroelf:SPMError:versionUnsupported', ...
        'No supported SPM version available.');
end
spmv = str2double(regexprep(spmver, '^spm(\d+).*$', '$1'));
if isempty(opts.spmv)
    opts.spmv = spmv;
end

% UI call
if nargin < 1 || ~ischar(sf) || isempty(sf)

    % make sure TPM files are available
    try
        spmp = fileparts(which('spm'));
        if opts.spmv <= 8
            tpm = { ...
                [spmp filesep 'tpm' filesep 'grey.nii']; ...
                [spmp filesep 'tpm' filesep 'white.nii']; ...
                [spmp filesep 'tpm' filesep 'csf.nii']};
            if exist(tpm{1}, 'file') ~= 2 || ...
                exist(tpm{2}, 'file') ~= 2 || ...
                exist(tpm{3}, 'file') ~= 2
                error('FILE_NOT_FOUND');
            end
        else
            tpm = { ...
                [spmp filesep 'tpm' filesep 'TPM.nii,1']; ...
                [spmp filesep 'tpm' filesep 'TPM.nii,2']; ...
                [spmp filesep 'tpm' filesep 'TPM.nii,3']};
            if exist(tpm{1}(1:end-2), 'file') ~= 2
                error('FILE_NOT_FOUND');
            end
        end
        if isempty(opts.t1temp)
            if opts.spmv <= 8
                t1temp = [spmp filesep 'templates' filesep 'T1.nii,1'];
            else
                t1temp = [spmp filesep 'canonical' filesep 'avg305T1.nii,1'];
            end
        else
            t1temp = opts.t1temp(:)';
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error('neuroelf:SPMError:fileNotFound', ...
            'Some file required for SPM was not found.');
    end

    % load figure
    try
        hFig = xfigure([neuroelf_path('tfg') '/spmx_preprojobs.tfg']);
    catch ne_eo;
        error('neuroelf:xfigure:openError', ...
            'Error creating UI for spmx_preprojobs: %s.', ne_eo.message);
    end

    % get tags
    hTag = hFig.TagStruct;

    % set into global variable
    ne_ui.spmx_prepro = struct( ...
        'af',   '', ...         % anatomical folder
        'afp',  {{}}, ...       % anatomical file pattern
        'agsc', 0, ...          % advanced global signal components
        'agsr', false, ...      % advanced global signal regression
        'arhp', 0, ...          % advanced raw high-pass filter
        'arpn', 0, ...          % advanced replace noise setting
        'ashp', 0, ...          % advanced subject-space high-pass filter
        'aslp', 0, ...          % advanced subject-space low-pass FWHM
        'astt', 'none', ...     % advanced subject-space trans (none/psc/z)
        'asvp', '_filt.vtc', ...% advanced subject-space filter VTC pattern
        'atms', 'none', ...     % advanced apply spatial mask (none/g/gw)
        'avsd', Inf, ...        % advanced vessel voxels SD threshold
        'avsl', 'none', ...     % advanced vessel voxels (all/atlas/none)
        'asdn', Inf, ...        % advanced SD noise threshold
        'cupa', false, ...      % clean-up slice-timing (a...) files
        'cupk', false, ...      % clean-up skip-volumes (k...) files
        'cupr', false, ...      % clean-up realigned+resliced (r...) files
        'cups', false, ...      % clean-up smoothed (s...) files
        'cupw', false, ...      % clean-up normalized (w...) files
        'cupx', false, ...      % clean-up brain-extracted (x...) files
        'exbr', true, ...       %
        'ff',   '', ...         % functional folder
        'ffp',  {{}}, ...       % functional file pattern
        'func', {{}}, ...       %
        'hFig', hFig, ...       % handle of xfigure to UI
        'hTag', hTag, ...       % struct with xfigure handles to controls
        'jobs', {{}}, ...       % jobs (SPM)
        'jobh', {{}}, ...       % job helps (SPM)
        'nslc', 28, ...         % number of slices
        'qual', false, ...      % run fmriquality on unprocessed data?
        'quav', false, ...      % run fmriquality on VTCs
        'rff',  {{}}, ...       % 
        'rfs',  {{}}, ...       %
        'runj', false, ...      % run jobs (instead of just creating)?
        'rvtc', false, ...      % 
        'sfld', {{}}, ...       % source folders
        'skip', 0, ...          % number of skipped volumes
        'spmv', opts.spmv, ...  % SPM version
        't1t',  t1temp, ...     % T1 template (for T1-based normalization)
        'vids', {{}}, ...       % 
        'vtci', 1, ...          % VTC-naming start index
        'vtcp', '%s/%s_RUN%02d_MNI.vtc', ... % VTC filename pattern
        'vtcs', false, ...      % generate VTCs
        'wsna', true);          % 

    % initialize controls
    hTag.DD_spmxprepro_skipvols.String = ...
        {'0'; '1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'};
    hTag.DD_spmxprepro_skipvols.Value = 1;
    hTag.LB_spmxprepro_sf_found.Max = 3;
    hTag.LB_spmxprepro_sf_found.Value = [];
    hTag.LB_spmxprepro_sf_found.String = {};
    hTag.DD_spmxprepro_sto.String = { ...
        '1,2,3,4,5,6, ...'; ...
        '30,29,28,27, ...'; ...
        '1,3,5, ... 2,4,6, ...'; ...
        '2,4,6, ... 1,3,5, ...'; ...
        '29,27, ... 30,28, ...'; ...
        '30,28, ... 29,27, ...'; ...
        '1,6, ... 26,2,7, ...'};
    hTag.DD_spmxprepro_sto.Value = 3;
    hTag.DD_spmxprepro_str.String = ...
        splittocellc(deblank(sprintf('%d ', 1:30)), ' ');
    hTag.DD_spmxprepro_rinterpe.String = {'1'; '2'; '3'; '4'; '5'; '6'; '7'};
    hTag.DD_spmxprepro_rinterpe.Value = 2;
    hTag.DD_spmxprepro_rinterpr.String = {'1'; '2'; '3'; '4'; '5'; '6'; '7'};
    hTag.DD_spmxprepro_rinterpr.Value = 4;
    ntypes = hTag.DD_spmxprepro_nrmtype.String;
    if ~iscell(ntypes)
        ntypes = cellstr(ntypes);
    end
    if exist([neuroelf_path('spm') filesep 'dartel' filesep 'Template1mm'], 'dir') ~= 7 || ...
        exist([neuroelf_path('spm') filesep 'dartel' filesep 'Template1mm' filesep 'Template1mm_7.nii'], 'file') ~= 2 || ...
        isequal(opts.spmv, 5)
        ntypes(2) = [];
        hTag.DD_spmxprepro_nrmtype.String = ntypes;
    end
    hTag.ED_spmxprepro_t1temp.String = t1temp;
    hTag.DD_spmxprepro_wvox.String = ...
        {'0.2'; '0.3'; '0.5'; '0.75'; '1.0'; '1.5'; '2.0'; '3.0'; '4.0'; '5.0'; '6.0'; '8.0'};
    hTag.DD_spmxprepro_wvox.Value = 8;
    hTag.PB_spmxprepro_progress.Visible = 'off';

    % make sure OK folder/file groups are disabled
    hFig.SetGroupEnabled('FoldOK', 'off');
    hFig.SetGroupEnabled('FilesOK', 'off');

    % set callbacks
    hTag.ED_spmxprepro_sf.Callback = @spp_sf_edit;
    hTag.BT_spmxprepro_sf_browse.Callback = @spp_sf_browse;
    hTag.BT_spmxprepro_sf_search.Callback = @spp_sf_search;
    hTag.CB_spmxprepro_advanced.Callback = @spp_advanced;
    hTag.BT_spmxprepro_delsubj.Callback = @spp_sf_delete;
    hTag.DD_spmxprepro_sto.Callback = @spp_sto;
    hTag.CB_spmxprepro_stom.Callback = @spp_stom;
    hTag.ED_spmxprepro_stol.Callback = @spp_stol;
    hTag.ED_spmxprepro_sttr.Callback = @spp_sttr;
    hTag.CB_spmxprepro_sttaa.Callback = @spp_sttaa;
    hTag.CB_spmxprepro_reslice.Callback = @spp_reslice;
    hTag.DD_spmxprepro_nrmtype.Callback = @spp_nrmtype;
    hTag.BT_spmxprepro_t1_browse.Callback = @spp_t1browse;
    hTag.CB_spmxprepro_autovtc.Callback = @spp_dovtcs;
    hTag.ED_spmxprepro_vtcfpat.Callback = @spp_vtcfpat;
    hTag.ED_spmxprepro_vtcfidx.Callback = @spp_vtcfidx;
    hTag.BT_spmxprepro_cancel.Callback = @spp_closeui;
    hTag.BT_spmxprepro_loadcfg.Callback = @spp_loadcfg;
    hTag.BT_spmxprepro_savecfg.Callback = @spp_savecfg;
    hTag.BT_spmxprepro_createjob.Callback = {@spp_create, 0};
    hTag.BT_spmxprepro_candrjob.Callback = {@spp_create, 1};
    hFig.CloseRequestFcn = @spp_closeui;

    % no "create jobs" button
    if nargin > 0 && isstruct(sf) && numel(sf) == 1 && ...
        isfield(sf, 'nocreate') && islogical(sf.nocreate) && ...
        numel(sf.nocreate) == 1 && sf.nocreate
        hTag.BT_spmxprepro_createjob.Visible = 'off';
    end

    % set visible and modal
    hFig.HandleVisibility = 'callback';
    hFig.Visible = 'on';

    % wait for dialog to be "uiresumed"
    uiwait(hFig.MLHandle);

    % copy over things
    rf = ne_ui.spmx_prepro;
    jobs = rf.jobs;
    jobhelps = rf.jobh;
    funcfiles = rf.func;
    rff = rf.rff;
    rfs = rf.rfs;

    % then return if runflag is requested
    if nargout > 5
        if ~rf.qual
            rf.sfls = [];
        end
        hFig.Delete;
        ne_ui.spmx_prepro = [];
        return;
    end

    % otherwise, if runflag is set, RUN
    if rf.runj

        % make buttons invisible
        hTag.BT_spmxprepro_loadcfg.Visible = 'off';
        hTag.BT_spmxprepro_savecfg.Visible = 'off';
        hTag.BT_spmxprepro_createjob.Visible = 'off';
        hTag.BT_spmxprepro_candrjob.Visible = 'off';
        hTag.BT_spmxprepro_cancel.Visible = 'off';
        pbar = hTag.PB_spmxprepro_progress;
        pbar.Visible = 'on';
        hFig.SetGroupEnabled('AllUIC',  'off');
        hFig.SetGroupEnabled('FoldOK',  'off');
        hFig.SetGroupEnabled('FilesOK', 'off');
        hFig.SetGroupEnabled('NrT1',    'off');
        hFig.SetGroupEnabled('NrEPI',   'off');
        hFig.SetGroupEnabled('NrANA',   'off');
        hFig.SetGroupEnabled('DoVTCs',  'off');
        drawnow;

        % no motion correction or sheet for now
        qo = struct('fftsd', true, 'motcor', false, 'qasheet', false, ...
            'savefq', true, 'skip', rf.skip, 'storefilt', true);
        
        % compute overall progress counter
        [jobcost, jobnames] = spp_jobcost(jobs);
        totalsessions = sum(cellfun('prodofsize', rf.sfls));
        ptotal = ...
            0.5 * double(rf.qual) * totalsessions + ...
            sum(jobcost) + ...
            0.125 * double(rf.rvtc) * totalsessions + ...
            0.5 * double(rf.rvtc) * double(rf.quav) * totalsessions + ...
            0.125 * double(rf.vtcs) * totalsessions + ...
            0.5 * double(rf.vtcs) * double(rf.quav) * totalsessions;
        pcount = 0;

        % quality check
        if rf.qual

            % error handling
            try

                % iterate over subjects
                for sc = 1:numel(rf.sfls)

                    % iterate over sessions
                    for sec = 1:numel(rf.sfls{sc})

                        % first file -> quality file
                        fmrifirst = rf.sfls{sc}{sec}{1};
                        fqfile = regexprep(fmrifirst, '\.(hdr|img|nii)(\,\d+)?$', ...
                            '.fq', 'ignorecase');
                        if numel(fmrifirst) ~= numel(fqfile) && ...
                            exist(fqfile, 'file') > 0
                            pcount = pcount + 0.5;
                            continue;
                        end

                        % run fmriquality and save info
                        pbar.Progress(pcount / ptotal, sprintf( ...
                            'fMRI quality checking subject %d, session %d...', sc, sec));
                        drawnow;
                        q = fmriquality(rf.sfls{sc}{sec}, qo);
                        pcount = pcount + 0.5;
                        fprintf('fmriqasheet output: %s\n', q.Filename);

                        % print some general information numbers
                        qfh = fmriqasheet(q);
                        set(qfh, 'CloseRequestFcn', '', 'DeleteFcn', '');
                        udt = get(qfh, 'UserData');
                        if isstruct(udt) && isfield(udt, 'object') && ...
                            numel(udt.object) == 1 && isxff(udt.object, true)
                            udt.object.ClearObject;
                        end
                        delete(qfh);
                    end
                end

            % error happened
            catch ne_eo;
                hFig.Delete;
                rethrow(ne_eo);
            end
        end
        qo.fftsd = false;
        qo.skip = 0;
        qo.storefilt = false;

        % run preprocessing
        try
            if rf.spmv < 8
                spm_defaults;
            else
                spm('defaults', 'FMRI');
                spm_jobman('initcfg');
            end
            for sc = 1:numel(jobs)
                pbar.Progress(pcount / ptotal, sprintf( ...
                    'Running %s job... (%d/%d jobs done...)', ...
                    jobnames{sc}, (sc - 1), numel(jobs)));
                spm_jobman('run', jobs(sc));
                pcount = pcount + jobcost(sc);
            end
        catch ne_eo;
            hFig.Delete;
            rethrow(ne_eo);
        end

        % backup the func files for cleanup
        cupfiles = rf.func;

        % convert to VTCs (after realignment)
        if rf.rvtc && ~isempty(rf.rff)

            % iterate over number of target sessions
            for sc = 1:numel(rf.rff)

                % import VTC
                try
                    vtc = cell(1, 1);
                    if iscell(rf.rff{sc})
                        rfffile = rf.rff{sc}{1};
                    else
                        rfffile = rf.rff{sc};
                    end
                    rfffile = regexprep(rfffile, '\.(hdr|head|img|nii)\,?\d*$', '', 'ignorecase');
                    rfffile = [fnpp(rfffile, [rf.vids{sc} '_']) '.vtc'];
                    pbar.Progress(pcount / ptotal, sprintf('Importing %s...', rfffile));
                    vtc{1} = ...
                        importvtcfromanalyze(rf.rff{sc}, [], 2, 'nearest', [], ...
                        struct('raw', true, 'snmat', rf.rfs{sc}));
                    vtc{1}.TR = floor(1000 * rf.sttr);
                    vtc{1}.RunTimeVars.AutoSave = true;
                    vtc{1}.SaveAs(rfffile);
                    pcount = pcount + 0.1;
                    rf.rff{sc} = vtc{1}.FilenameOnDisk;
                    if rf.quav
                        pbar.Progress(pcount / ptotal, sprintf( ...
                            'Running fMRI quality check on %s...', rf.rff{sc}));
                        q = fmriquality(rf.rff{sc}, qo);
                        fprintf(' -> fmriqasheet output for %s:\n', q.Filename);
                        qfh = fmriqasheet(q);
                        set(qfh, 'CloseRequestFcn', '', 'DeleteFcn', '');
                        udt = get(qfh, 'UserData');
                        if isstruct(udt) && isfield(udt, 'object') && ...
                            numel(udt.object) == 1 && isxff(udt.object, true)
                            udt.object.ClearObject;
                        end
                        delete(qfh);
                        pcount = pcount + 0.5;
                    end
                catch ne_eo;
                    hFig.Delete;
                    rethrow(ne_eo);
                end

                % remove object from memory
                clearxffobjects(vtc);
            end
        end

        % convert to VTCs (final product)
        if rf.vtcs

            % iterate over number of target sessions
            lastid = '::INVALID::';
            lastidc = rf.vtci;
            for sc = 1:numel(rf.func)

                % import VTC
                try
                    vtc = cell(1, 1);
                    if ~strcmp(lastid, rf.vids{sc})
                        lastid = rf.vids{sc};
                        lastidc = rf.vtci;
                    end
                    vtcfile = sprintf(rf.vtcp, rf.vidf{sc}, rf.vids{sc}, lastidc);
                    pbar.Progress(pcount / ptotal, sprintf('Importing %s...', vtcfile));
                    vtc{1} = ...
                        importvtcfromanalyze(rf.func{sc}, rf.wbb, rf.wvox, 'cubic');
                    vtc{1}.TR = floor(1000 * rf.sttr);
                    vtc{1}.RunTimeVars.AutoSave = true;
                    vtc{1}.SaveAs(vtcfile);
                    pcount = pcount + 0.1;
                    rf.func{sc} = vtc{1}.FilenameOnDisk;

                    % fMRIquality on VTCs
                    if rf.quav
                        pbar.Progress(pcount / ptotal, sprintf( ...
                            'Running fMRI quality check on %s...', rf.func{sc}));
                        q = fmriquality(rf.func{sc}, qo);
                        fprintf(' -> fmriqasheet output for %s:\n', q.Filename);
                        qfh = fmriqasheet(q);
                        set(qfh, 'CloseRequestFcn', '', 'DeleteFcn', '');
                        udt = get(qfh, 'UserData');
                        if isstruct(udt) && isfield(udt, 'object') && ...
                            numel(udt.object) == 1 && isxff(udt.object, true)
                            udt.object.ClearObject;
                        end
                        delete(qfh);
                        pcount = pcount + 0.5;
                    end
                catch ne_eo;
                    hFig.Delete;
                    rethrow(ne_eo);
                end

                % remove object from memory
                clearxffobjects(vtc);

                % increase session counter
                lastidc = lastidc + 1;
            end
        end

        % clean-up (functional files)
        if rf.cupa || rf.cupk || rf.cupr || rf.cupw || rf.cups
            [cupfold, cupfile1] = fileparts(cupfiles{1}{1});
            delimgs = ~isempty(regexpi(cupfiles{1}{1}, '.img(\,\d+)?$'));
            if cupfile1(1) == 's'
                if rf.cups
                    for sc = 1:numel(cupfiles)
                        delfiles = unique(regexprep(cupfiles{sc}, '\,\d+$', ''));
                        if exist(delfiles{1}, 'file') ~= 2
                            break;
                        end
                        mdelete(delfiles);
                        if delimgs
                            mdelete(regexprep(delfiles, '\.img$', '.hdr', 'preservecase'));
                        end
                    end
                end
                cupfiles = fnxp(cupfiles);
                [cupfold, cupfile1] = fileparts(cupfiles{1}{1});
            end
            if cupfile1(1) == 'w'
                if rf.cupw
                    for sc = 1:numel(cupfiles)
                        delfiles = unique(regexprep(cupfiles{sc}, '\,\d+$', ''));
                        if exist(delfiles{1}, 'file') ~= 2
                            break;
                        end
                        mdelete(delfiles);
                        if delimgs
                            mdelete(regexprep(delfiles, '\.img$', '.hdr', 'preservecase'));
                        end
                    end
                end
                cupfiles = fnxp(cupfiles);
                [cupfold, cupfile1] = fileparts(cupfiles{1}{1});
            end
            if cupfile1(1) == 'r'
                if rf.cupr
                    for sc = 1:numel(cupfiles)
                        delfiles = unique(regexprep(cupfiles{sc}, '\,\d+$', ''));
                        mdelete(delfiles);
                        if delimgs
                            mdelete(regexprep(delfiles, '\.img$', '.hdr', 'preservecase'));
                        end
                    end
                end
                cupfiles = fnxp(cupfiles);
                [cupfold, cupfile1] = fileparts(cupfiles{1}{1});
            end
            if cupfile1(1) == 'a'
                if rf.cupa
                    for sc = 1:numel(cupfiles)
                        delfiles = unique(regexprep(cupfiles{sc}, '\,\d+$', ''));
                        mdelete(delfiles);
                        if delimgs
                            mdelete(regexprep(delfiles, '\.img$', '.hdr', 'preservecase'));
                        end
                    end
                end
                cupfiles = fnxp(cupfiles);
                [cupfold, cupfile1] = fileparts(cupfiles{1}{1});
            end
            if cupfile1(1) == 'k' && rf.cupk
                for sc = 1:numel(cupfiles)
                    delfiles = unique(regexprep(cupfiles{sc}, '\,\d+$', ''));
                    mdelete(delfiles);
                    if delimgs
                        mdelete(regexprep(delfiles, '\.img$', '.hdr', 'preservecase'));
                    end
                end
            end
        end
    end

    % and then return
    hFig.Delete;
    ne_ui.spmx_prepro = [];
    return;
end

% for non-UI mode, argument check
if nargin < 3 || ...
   ~ischar(sf) || ...
    isempty(sf) || ...
    exist(sf(:)', 'dir') ~= 7 || ...
   ~ischar(fp) || ...
    isempty(fp) || ...
   ~ischar(ap) || ...
    isempty(ap)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument supplied.' ...
    );
end
sf = sf(:)';
fp = fp(:)';
ap = ap(:)';

% load job specification file
try
    if opts.spmv < 8
        load([neuroelf_path('spm') '/spm5_prepro.mat']);
    else
        load([neuroelf_path('spm') '/spm8_prepro.mat']);
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:BadFile', ...
        'SPM5/SPM8 jobs template file not found or corrupt.' ...
    );
end

% set required SPM files
try
    spmp = fileparts(which('spm'));
    if opts.spmv <= 8
        tpm = { ...
            [spmp filesep 'tpm' filesep 'grey.nii']; ...
            [spmp filesep 'tpm' filesep 'white.nii']; ...
            [spmp filesep 'tpm' filesep 'csf.nii']};
        if exist(tpm{1}, 'file') ~= 2 || ...
            exist(tpm{2}, 'file') ~= 2 || ...
            exist(tpm{3}, 'file') ~= 2
            error('FILE_NOT_FOUND');
        end
    else
        tpm = { ...
            [spmp filesep 'tpm' filesep 'TPM.nii,1']; ...
            [spmp filesep 'tpm' filesep 'TPM.nii,2']; ...
            [spmp filesep 'tpm' filesep 'TPM.nii,3']};
        if exist(tpm{1}(1:end-2), 'file') ~= 2
            error('FILE_NOT_FOUND');
        end
    end
    if isempty(opts.t1temp)
        if opts.spmv <= 8
            t1temp = [spmp filesep 'templates' filesep 'T1.nii,1'];
        else
            t1temp = [spmp filesep 'canonical' filesep 'avg305T1.nii,1'];
        end
    else
        t1temp = opts.t1temp(:)';
    end
    if opts.spmv < 8
        jobs{4}.spatial{2}.preproc.opts.tpm = tpm;
    else
        jobs{5}.spm.spatial.preproc.opts.tpm = tpm;
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:FileNotFound', ...
        'SPM5/SPM8 not installed or required file/s not found.' ...
    );
end

% normalization check 1
if opts.normtype == 'd' && ...
   (~isequal(spmv, 8) || ...
    exist([neuroelf_path('spm') filesep 'dartel' filesep 'Template1mm'], 'dir') ~= 7 || ...
    exist([neuroelf_path('spm') filesep 'dartel' filesep 'Template1mm' filesep 'Template1mm_7.nii'], 'file') ~= 2)
    error( ...
        'neuroelf:BadOption', ...
        'DARTEL normalization only available with SPM8 and precomputed template.' ...
    );
end

% DARTEL usage
dartel = false;
try
    if opts.normtype == 'd'

        % load additional job files
        dimport = load([neuroelf_path('spm') '/spm8_dartel_import.mat']);
        dimport = dimport.jobs;
        dstwarp = load([neuroelf_path('spm') '/spm8_dartel_warp_to_template.mat']);
        dstwarp = dstwarp.jobs;
        dwtomni = load([neuroelf_path('spm') '/spm8_dartel_mni_normalize_anat.mat']);
        dwtomni = dwtomni.jobs;
        dtpimgs = findfiles([neuroelf_path('spm') '/dartel/Template1mm/'], '*.nii', 'depth=1');
        if numel(dtpimgs) == 8
            dartel = true;
            opts.reslice = true;
            for tc = 1:8
                dstwarp{1}.spm.tools.dartel.warp1.settings.param(tc).template = dtpimgs(tc);
            end
            dwtomni{1}.spm.tools.dartel.mni_norm.template = dtpimgs(end);
        else
            error('INVALID_DARTEL_FILES');
        end
        dftomni = dwtomni;
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:FileNotFound', ...
        'DARTEL-related file not found or invalid.' ...
    );
end

% make settings (without further checks)
if opts.spmv < 8
    jobs{1}.temporal{1}.st.tr = opts.sttr;
    jobs{1}.temporal{1}.st.refslice = opts.str;
    jobs{2}.spatial{1}.realign{1}.estwrite.eoptions.interp = opts.rinterpe;
    jobs{2}.spatial{1}.realign{1}.estwrite.eoptions.quality = opts.rqual;
    jobs{2}.spatial{1}.realign{1}.estwrite.eoptions.rtm = opts.rrtm;
    jobs{2}.spatial{1}.realign{1}.estwrite.eoptions.sep = opts.rsep;
    jobs{2}.spatial{1}.realign{1}.estwrite.roptions.interp = opts.rinterpr;
    jobs{3}.util{1}.imcalc.options.interp = opts.rinterpr;
    jobs{4}.spatial{3}.normalise{1}.write.roptions.bb = opts.wbb;
    jobs{4}.spatial{3}.normalise{1}.write.roptions.interp = opts.winterp;
    jobs{4}.spatial{3}.normalise{1}.write.roptions.vox = opts.wvox;
    jobs{5}.spatial{1}.smooth.fwhm = opts.smk;
    if ~opts.reslice
        jobs{2}.spatial{1}.realign{1}.estimate = ...
            jobs{2}.spatial{1}.realign{1}.estwrite;
        jobs{2}.spatial{1}.realign{1} = ...
            rmfield(jobs{2}.spatial{1}.realign{1}, 'estwrite');
    end
else
    jobs{1}.spm.temporal.st.tr = opts.sttr;
    jobs{1}.spm.temporal.st.refslice = opts.str;
    jobs{2}.spm.spatial.realign.estwrite.eoptions.interp = opts.rinterpe;
    jobs{2}.spm.spatial.realign.estwrite.eoptions.quality = opts.rqual;
    jobs{2}.spm.spatial.realign.estwrite.eoptions.rtm = opts.rrtm;
    jobs{2}.spm.spatial.realign.estwrite.eoptions.sep = opts.rsep;
    jobs{2}.spm.spatial.realign.estwrite.roptions.interp = opts.rinterpr;
    jobs{3}.spm.util.imcalc.options.interp = opts.rinterpr;
    jobs{6}.spm.spatial.normalise.write.roptions.bb = opts.wbb;
    jobs{6}.spm.spatial.normalise.write.roptions.interp = opts.winterp;
    jobs{6}.spm.spatial.normalise.write.roptions.vox = opts.wvox;
    if ~opts.reslice
        jobs{2}.spm.spatial.realign.estimate = ...
            jobs{2}.spm.spatial.realign.estwrite;
        jobs{2}.spm.spatial.realign = ...
            rmfield(jobs{2}.spm.spatial.realign, 'estwrite');
    end
end

% register anatomical to T1 first?
if opts.reganat

    % load additional job file
    try
        if opts.spmv < 8
            crjob = load([neuroelf_path('spm') '/spm5_coreg.mat']);
        else
            crjob = load([neuroelf_path('spm') '/spm8_coreg.mat']);
        end
        crjob = crjob.jobs;
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error( ...
            'neuroelf:BadFile', ...
            'SPM5/SPM8 coregistration job template file not found or corrupt.' ...
        );
    end
end

% normalization check 2
if opts.normtype ~= 's' && ...
    opts.normtype ~= 'd'

    % load additional job file
    try
        if opts.spmv < 8
            nrjob = load([neuroelf_path('spm') '/spm5_normest.mat']);
        else
            nrjob = load([neuroelf_path('spm') '/spm8_normest.mat']);
        end
        nrjob = nrjob.jobs;
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error( ...
            'neuroelf:BadFile', ...
            'SPM5/SPM8 normalization job template file not found or corrupt.' ...
        );
    end
end

% find folders and files
anatfolders = spp_findfolders(sf, ap);
if isempty(anatfolders)
    error( ...
        'neuroelf:FolderNotFound', ...
        'No anatomical folders found...' ...
    );
end
anatfiles = spp_findfiles(anatfolders, opts.afp);
if isempty(anatfiles)
    error( ...
        'neuroelf:NoFilesFound', ...
        'No anatomical file(s) found.' ...
    );
end

% extract brain?
if opts.extbrain
    try
        spmx_extract_brain(anatfiles{1});
        anatfiles = fnpp(anatfiles, 'x');
    catch ne_eo;
        rethrow(ne_eo);
    end
end

% patch files
anatfiles{1} = [anatfiles{1} ',1'];
funcfolders = spp_findfolders(sf, fp);
[funcfiles, nslc, fourd] = spp_findfiles(funcfolders, opts.ffp, opts.skip);
if isempty(funcfiles)
    error( ...
        'neuroelf:BadArgument', ...
        'No functional files/folders found.' ...
    );
end
if any(diff(nslc))
    error( ...
        'neuroelf:BadArgument', ...
        'Functional files must all have the same number of slices.' ...
    );
end
nslc = nslc(1);

% disable norm-per session if not needed
nff = numel(funcfiles);
if nff < 2
    opts.normsess = false;
end

% get number of slices and possibly auto settings for ST
if isempty(opts.stta)
    opts.stta = opts.sttr * (nslc - 1) / nslc;
end
if ischar(opts.sto)
    switch (opts.sto)
        case {'a'}
            opts.sto = 1:nslc;
        case {'aie'}
            opts.sto = [2:2:nslc, 1:2:nslc];
        case {'aio'}
            opts.sto = [1:2:nslc, 2:2:nslc];
        case {'d'}
            opts.sto = nslc:-1:1;
        case {'die'}
            fslc = 2 * floor(0.5 * nslc);
            sslc = 2 * floor(0.5 * (nslc - 1)) + 1;
            opts.sto = [fslc:-2:1, sslc:-2:1];
        case {'dio'}
            fslc = 2 * floor(0.5 * (nslc - 1)) + 1;
            sslc = 2 * floor(0.5 * nslc);
            opts.sto = [fslc:-2:1, sslc:-2:1];
    end
end
if opts.spmv < 8
    jobs{1}.temporal{1}.st.nslices = nslc;
    jobs{1}.temporal{1}.st.scans = funcfiles;
    jobs{1}.temporal{1}.st.so = opts.sto;
    jobs{1}.temporal{1}.st.ta = opts.stta;
else
    jobs{1}.spm.temporal.st.nslices = nslc;
    jobs{1}.spm.temporal.st.scans = funcfiles;
    jobs{1}.spm.temporal.st.so = opts.sto;
    jobs{1}.spm.temporal.st.ta = opts.stta;
end

% rewrite in single precision
if opts.fsingle

    % iterate over all sessions/files
    for sc = 1:nff
        for fc = 1:numel(funcfiles{sc})

            % load object
            h = [];
            fd = false;
            try
                h = xff(regexprep(regexprep(funcfiles{sc}{fc}, ...
                    ',\d+$', ''), '\.img$', '.hdr', 'preservecase'));
                if ~isa(h.VoxelData(1), 'single')
                    if h.ImgDim.ScalingSlope ~= 1 || ...
                        h.ImgDim.ScalingIntercept ~= 0
                        if h.ImgDim.ScalingSlope == 0 || ...
                            isinf(h.ImgDim.ScalingSlope) || ...
                            isnan(h.ImgDim.ScalingSlope)
                            h.ImgDim.ScalingSlope = 1;
                        end
                        if isinf(h.ImgDim.ScalingIntercept) || ...
                            isinf(h.ImgDim.ScalingIntercept)
                            h.ImgDim.ScalingIntercept = 0;
                        end
                        h.VoxelData = h.ImgDim.ScalingSlope .* ...
                            (double(h.VoxelData) + h.ImgDim.ScalingIntercept);
                    end
                    h.VoxelData = single(h.VoxelData);
                    h.ImgDim.ScalingSlope = 1;
                    h.ImgDim.ScalingIntercept = 0;
                    h.ImgDim.DataType = 16;
                    h.ImgDim.BitsPerPixel = 32;
                    if size(h.VoxelData, 4) > 1
                        fd = true;
                    elseif h.NIIFileType < 2
                        vdname = h.FilenameOnDisk;
                        vdfile = sprintf('%simg_%06x', vdname(1:end-3), ceil(16777215 * rand(1)));
                        vdfid = fopen(vdfile, 'wb', h.Endian);
                        if vdfid < 1
                            error( ...
                                'neuroelf:FileOpenError', ...
                                'Error opening file for altering datatype.' ...
                            );
                        end
                        fseek(vdfid, 0, 'bof');
                        if h.ImgDim.VoxOffset > 0
                            fwrite(vdfid, ...
                                char(zeros(1, floor(h.ImgDim.VoxOffset))), ...
                                'char');
                        end
                        if fwrite(vdfid, h.VoxelData(:), 'single') ~= numel(h.VoxelData)
                            fclose(vdfid);
                            error( ...
                                'neuroelf:FileWriteError', ...
                                'Error writing single-precision data to file.' ...
                            );
                        end
                        fclose(vdfid);
                        if ~movefile(vdfile, vdfile(1:end-7))
                            error( ...
                                'neuroelf:FileMoveError', ...
                                'Error moving file to target location.' ...
                            );
                        end
                    end
                    h.Save;
                end
                h.ClearObject;
            catch ne_eo;
                clearxffobjects({h});
                rethrow(ne_eo);
            end
            if fd
                break;
            end
        end
    end
end

% begin to build scans/data into jobs -> slice timing
if opts.st
    funcfiles = fnpp(funcfiles, 'a');
end

% realignment
if opts.reslice
    if opts.spmv < 8
        jobs{2}.spatial{1}.realign{1}.estwrite.data = funcfiles;
    else
        jobs{2}.spm.spatial.realign.estwrite.data = funcfiles;
    end
    funcfiles = fnpp(funcfiles, 'r');
    rff = funcfiles;
else
    if opts.spmv < 8
        jobs{2}.spatial{1}.realign{1}.estimate.data = funcfiles;
    else
        jobs{2}.spm.spatial.realign.estimate.data = funcfiles;
    end
    rff = {};
end

% mean image (imcalc)
if opts.spmv < 8
    jobs{3}.util{1}.imcalc.outdir{1} = anatfolders{1};
    if opts.normsess
        jobs{3}.util = jobs{3}.util(ones(1, nff));
        for sc = 1:nff
            jobs{3}.util{sc}.imcalc.input = funcfiles{sc};
            jobs{3}.util{sc}.imcalc.output =  ...
                sprintf('mean_functional%02d.img', sc);
        end
    else
        jobs{3}.util{1}.imcalc.input = cat(1, funcfiles{:});
    end
else
    jobs{3}.spm.util.imcalc.outdir{1} = anatfolders{1};
    if opts.normsess
        jobs{3}.spm = jobs{3}.spm(ones(1, nff));
        for sc = 1:nff
            jobs{3}.spm(sc).util.imcalc.input = funcfiles{sc};
            jobs{3}.spm(sc).util.imcalc.output =  ...
                sprintf('mean_functional%02d.img', sc);
        end
    else
        jobs{3}.spm.util.imcalc.input = cat(1, funcfiles{:});
    end
end

% coregistration
if opts.fun2str
    if opts.spmv < 8
        jobs{4}.spatial{1}.coreg{1}.estimate.other = jobs{3}.util{1}.imcalc.input;
        jobs{4}.spatial{1}.coreg{1}.estimate.ref{1} = anatfiles{1};
        jobs{4}.spatial{1}.coreg{1}.estimate.source{1} = ...
            [anatfolders{1} filesep 'mean_functional.img,1'];
    else
        jobs{4}.spm.spatial.coreg.estimate.other = jobs{3}.spm(1).util.imcalc.input;
        jobs{4}.spm.spatial.coreg.estimate.ref{1} = anatfiles{1};
        jobs{4}.spm.spatial.coreg.estimate.source{1} = ...
            [anatfolders{1} filesep 'mean_functional.img,1'];
    end
else
    if opts.spmv < 8
        jobs{4}.spatial{1}.coreg{1}.estimate.ref{1} = ...
            [anatfolders{1} filesep 'mean_functional.img,1'];
        jobs{4}.spatial{1}.coreg{1}.estimate.source{1} = anatfiles{1};
    else
        jobs{4}.spm.spatial.coreg.estimate.ref{1} = ...
            [anatfolders{1} filesep 'mean_functional.img,1'];
        jobs{4}.spm.spatial.coreg.estimate.source{1} = anatfiles{1};
    end
end

% segmentation (normalization estimation)
if opts.spmv < 8
    jobs{4}.spatial{2}.preproc.data{1} = anatfiles{1};
else
    jobs{5}.spm.spatial.preproc.data{1} = anatfiles{1};
end

% normalization reslicing
rfs = [anatfiles{1}(1:end-6) '_seg_sn.mat'];
if opts.spmv < 8
    jobs{4}.spatial{3}.normalise{1}.write.subj.matname{1} = rfs;
    jobs{4}.spatial{3}.normalise{1}.write.subj.resample = ...
        cat(1, funcfiles{:}, {[anatfolders{1} filesep 'mean_functional.img,1']});
else
    jobs{6}.spm.spatial.normalise.write.subj.matname{1} = rfs;
    jobs{6}.spm.spatial.normalise.write.subj.resample = ...
        cat(1, funcfiles{:}, {[anatfolders{1} filesep 'mean_functional.img,1']});
end
if opts.normsess
    normfuncfiles = funcfiles;
end
funcfiles = fnpp(funcfiles, 'w');

% smoothing
if opts.smk > 0
    if opts.spmv < 8
        jobs{5}.spatial{1}.smooth.data = cat(1, funcfiles{:});
    else
        jobs{7}.spm.spatial.smooth.data = cat(1, funcfiles{:});
    end
    funcfiles = fnpp(funcfiles, 's');
else
    jobs(end) = [];
    jobhelps(end) = [];
end

% type of normalization (default: segmentation, nothing to do)
if opts.normtype ~= 's' && ...
    opts.normtype ~= 'd'

    % use anatomical with T1.nii -> normalization
    if opts.normtype == 'a'

        % parameter file
        snpfile = [anatfiles{1}(1:end-6) '_sn.mat'];
        if opts.spmv < 8
            nrjob{1}.spatial{1}.normalise{1}.est.subj.source{1} = anatfiles{1};
            nrjob{1}.spatial{1}.normalise{1}.est.eoptions.template{1} = t1temp;
            jobs{4}.spatial{2} = nrjob{1}.spatial{1};
            jobs{4}.spatial{3}.normalise{1}.write.subj.matname{1} = snpfile;
        else
            nrjob{1}.spm.spatial.normalise.est.subj.source{1} = anatfiles{1};
            if opts.spmv > 8
                nrjob{1}.spm.spatial.normalise.est.subj.vol = ...
                    nrjob{1}.spm.spatial.normalise.est.subj.source;
                nrjob{1}.spm.spatial.normalise.est.subj = rmfield( ...
                    nrjob{1}.spm.spatial.normalise.est.subj, 'source');
                nrjob{1}.spm.spatial.normalise.est.subj = rmfield( ...
                    nrjob{1}.spm.spatial.normalise.est.subj, 'wtsrc');
            end
            nrjob{1}.spm.spatial.normalise.est.eoptions.template{1} = t1temp;
            jobs{5}.spm.spatial = nrjob{1}.spm.spatial;
            jobs{6}.spm.spatial.normalise.write.subj.matname{1} = snpfile;
        end

    % use mean functional with EPI.nii -> normalization
    else

        % parameter file
        snpfile = [anatfolders{1} filesep 'mean_functional_sn.mat'];
        if opts.spmv < 8
            nrjob{1}.spatial{1}.normalise{1}.est.subj.source{1} = ...
                [anatfolders{1} filesep 'mean_functional.img,1'];
            nrjob{1}.spatial{1}.normalise{1}.est.eoptions.template{1} = ...
                [spmp filesep 'templates' filesep 'EPI.nii,1'];
            jobs{4}.spatial{2} = nrjob{1}.spatial{1};
            jobs{4}.spatial{3}.normalise{1}.write.subj.matname{1} = snpfile;
            jobs{4}.spatial(1) = [];
            if opts.normsess
                snpfile = '';
                jobs{4}.spatial = ...
                    jobs{4}.spatial([ones(1, nff), 2 .* ones(1, nff)]);
                for sc = 1:nff
                    jobs{4}.spatial{sc}.normalise{1}.est.subj.source{1} = ...
                        sprintf('%s%smean_functional%02d.img', ...
                        anatfolders{1}, filesep, sc);
                    jobs{4}.spatial{sc+nff}.normalise{1}.write.subj.matname{1} = ...
                        sprintf('%s%smean_functional%02d_sn.mat', ...
                        anatfolders{1}, filesep, sc);
                    jobs{4}.spatial{sc+nff}.normalise{1}.write.subj.resample = ...
                        normfuncfiles{sc};
                end
            end
        else
            nrjob{1}.spm.spatial.normalise.est.subj.source{1} = ...
                [anatfolders{1} filesep 'mean_functional.img,1'];
            if opts.spmv > 8
                nrjob{1}.spm.spatial.normalise.est.subj.vol = ...
                    nrjob{1}.spm.spatial.normalise.est.subj.source;
                nrjob{1}.spm.spatial.normalise.est.subj = rmfield( ...
                    nrjob{1}.spm.spatial.normalise.est.subj, 'source');
                nrjob{1}.spm.spatial.normalise.est.subj = rmfield( ...
                    nrjob{1}.spm.spatial.normalise.est.subj, 'wtsrc');
            end
            nrjob{1}.spm.spatial.normalise.est.eoptions.template{1} = ...
                [spmp filesep 'templates' filesep 'EPI.nii,1'];
            jobs{5}.spm.spatial = nrjob{1}.spm.spatial;
            jobs{6}.spm.spatial.normalise.write.subj.matname{1} = snpfile;
            jobs(4) = [];
            if opts.normsess
                snpfile = '';
                jobs = [jobs(1:3), jobs(4 * ones(1, nff)), jobs(5 * ones(1, nff)), jobs(6:end)];
                for sc = 1:nff
                    jobs{3+sc}.spm.spatial.normalise.est.subj.source{1} = ...
                        sprintf('%s%smean_functional%02d.img', ...
                        anatfolders{1}, filesep, sc);
                    if opts.spmv > 8
                        jobs{3+sc}.spm.spatial.normalise.est.subj.vol = ...
                            jobs{3+sc}.spm.spatial.normalise.est.subj.source;
                        jobs{3+sc}.spm.spatial.normalise.est.subj = rmfield( ...
                            jobs{3+sc}.spm.spatial.normalise.est.subj, 'source');
                        jobs{3+sc}.spm.spatial.normalise.est.subj = rmfield( ...
                            jobs{3+sc}.spm.spatial.normalise.est.subj, 'wtsrc');
                    end
                    jobs{3+sc+nff}.spm.spatial.normalise.write.subj.matname{1} = ...
                        sprintf('%s%smean_functional%02d_sn.mat', ...
                        anatfolders{1}, filesep, sc);
                    jobs{3+sc+nff}.spm.spatial.normalise.write.subj.resample = ...
                        normfuncfiles{sc};
                end
            end
        end
    end

% default: normalization via segmentation
else

    % parameter file
    snpfile = [anatfiles{1}(1:end-6) '_seg_sn.mat'];
end
if isempty(rff)
    rfs = {};
else
    rfs = repmat({snpfile}, nff, 1);
end

% DARTEL
if dartel
    flowfield = fnpp([anatfiles{1}(1:end-6) '.nii'], 'u_rc1');
    dimport{1}.spm.tools.dartel.initial.matnames = {snpfile};
    dimport{1}.spm.tools.dartel.initial.odir{1} = fileparts(snpfile);
    dstwarp{1}.spm.tools.dartel.warp1.images{1} = ...
        {fnpp([anatfiles{1}(1:end-6) '.nii'], 'rc1')};
    dstwarp{1}.spm.tools.dartel.warp1.images{2} = ...
        {fnpp([anatfiles{1}(1:end-6) '.nii'], 'rc2')};
    dwtomni{1}.spm.tools.dartel.mni_norm.data.subj.flowfield = {flowfield};
    dwtomni{1}.spm.tools.dartel.mni_norm.data.subj.images = ...
        {fnpp([anatfiles{1}(1:end-6) '.nii'], 'r')};
    dftomni{1}.spm.tools.dartel.mni_norm.data.subj.flowfield = {flowfield};
    urff = rff;
    for sc = 1:numel(urff)
        urff{sc} = unique(regexprep(rff{sc}, ',\d+$', ''));
    end
    dftomni{1}.spm.tools.dartel.mni_norm.data.subj.images = ...
        cat(1, urff{:});
    dftomni{1}.spm.tools.dartel.mni_norm.bb = opts.wbb;
    dftomni{1}.spm.tools.dartel.mni_norm.vox = opts.wvox;
    if opts.smk > 0
        dftomni{1}.spm.tools.dartel.mni_norm.fwhm = opts.smk;
    else
        dftomni{1}.spm.tools.dartel.mni_norm.fwhm = opts.wvox;
    end
    jobs = [jobs(1:5), dimport, dstwarp, dwtomni, dftomni];
    jobhelps = [jobhelps(1:5), {'', '', '', ''}];
end

% unfold mean image computations if necessary
if opts.normsess && ...
    opts.spmv >= 8
    jobs = [jobs(1:2), jobs(3 * ones(1, nff)), jobs(4:end)];
    for sc = 1:nff
        jobs{2 + sc}.spm = jobs{2 + sc}.spm(sc);
    end
end

% slice scantiming
if ~opts.st
    jobs(1) = [];
    jobhelps(1) = [];
end

% write normalized anatomical image
if opts.wsna && ...
   ~isempty(snpfile) && ...
   ~isempty(anatfiles{1})

    % add job
    if ~dartel
        jobs(end+1) = writesn(1);
    end
    if opts.spmv < 8
        jobs{end}.spatial{1}.normalise{1}.write.subj.matname = {snpfile};
        jobs{end}.spatial{1}.normalise{1}.write.subj.resample = anatfiles(1);
    elseif ~dartel
        jobs{end}.spm.spatial.normalise.write.subj.matname = {snpfile};
        jobs{end}.spm.spatial.normalise.write.subj.resample = anatfiles(1);
    end
end

% coregister anatomical to T1 first
if opts.reganat
    if opts.spmv < 8
        crjob{1}.spatial{1}.coreg{1}.estimate.ref{1} = ...
            [spmp filesep 'templates' filesep 'T1.nii,1'];
        crjob{1}.spatial{1}.coreg{1}.estimate.source{1} = anatfiles{1};
    else
        crjob{1}.spm.spatial.coreg.estimate.ref{1} = ...
            [spmp filesep 'templates' filesep 'T1.nii,1'];
        crjob{1}.spm.spatial.coreg.estimate.source{1} = anatfiles{1};
    end
    jobs = [crjob, jobs];
    jobhelps = [{''}, jobhelps];
end

% save job?
if ~isempty(opts.jobfile)
    try
        save(opts.jobfile, 'jobs', 'jobhelps', '-v6');
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        try
            save(opts.jobfile, 'jobs', 'jobhelps');
        catch ne_eo;
            warning( ...
                'neuroelf:MatlabError', ...
                'Error saving jobfile: %s', ...
                ne_eo.message ...
            );
        end
    end
end

% run job?
if opts.jobrun
    try
        if opts.spmv < 8
            spm_defaults;
        else
            spm('defaults', 'FMRI');
            spm_jobman('initcfg');
        end
        spm_jobman('run', jobs);
    catch ne_eo;
        rethrow(ne_eo);
    end
end

% if funcfiles is output, remove ",VOL" parts
if nargout > 2
    for sc = 1:nff

        % four-D file, only one file anyway
        if fourd(sc)

            % take first filename without ",1"
            funcfiles{sc} = {funcfiles{sc}{1}(1:end-2)};
            if ~isempty(rff{sc})
                rff{sc} = {rff{sc}{1}(1:end-2)};
            end

        % regular 3D files
        else

            % depending on rff availability
            if isempty(rff)

                % iterate over files
                for fc = 1:numel(funcfiles{sc})

                    % and remove ",VOL" part
                    funcfiles{sc}{fc} = regexprep(funcfiles{sc}{fc}, ',\d+$', '');
                end
            else
                for fc = 1:numel(funcfiles{sc})
                    funcfiles{sc}{fc} = regexprep(funcfiles{sc}{fc}, ',\d+$', '');
                    rff{sc}{fc} = regexprep(rff{sc}{fc}, ',\d+$', '');
                end
            end
        end
    end
end



% UI functions



function spp_closeui(varargin)
global ne_ui;
uiresume(ne_ui.spmx_prepro.hFig.MLHandle);


function spp_sf_edit(varargin)
global ne_ui;
cf = ne_ui.spmx_prepro.hFig;
ch = ne_ui.spmx_prepro.hTag;

% current configuration
fs = ch.ED_spmxprepro_sf.String;

% if empty, set certain controls disabled and return
if isempty(fs)
    cf.SetGroupEnabled('FoldOK', 'off');
    return;
end

% removing trailing fileseps
while ~isempty(fs) && ...
    any('\/' == fs(end))
    fs(end) = [];
end

% split into path and name
[fp, fn, fe] = fileparts(fs);
fn = [fn, fe];

% replace empty path
if isempty(fp)
    fp = pwd;
end
if isempty(fn)
    fn = '*';
end

% then look for subject folders
sf = findfiles(fp, fn, 'dirs', 'depth=1');
sf = sf(:);
ne_ui.spmx_prepro.sfld = sf;
vids = findfiles(fp, fn, 'dirs', 'depth=1', 'relative=');
vids = vids(:);
if numel(unique(vids)) ~= numel(vids)
    vidp = sf;
    for sc = 1:numel(sf)
        vidp{sc} = fileparts(vidp{sc});
    end
    splitfld = true;
    while splitfld
        for sc = 1:numel(sf)
            olds = vids{sc};
            [vidp{sc}, vids{sc}] = fileparts(vidp{sc});
            vids{sc} = [vids{sc} '_' olds];
        end
        if any(cellfun('isempty', vidp)) || ...
            numel(unique(vids)) == numel(vids)
            splitfld = false;
        end
    end
end
ne_ui.spmx_prepro.vids = [sf(:), vids(:)];

% plural form
if numel(sf) ~= 1
    pls = 's';
else
    pls = '';
end
ch.TX_spmxprepro_subffound.String = ...
    sprintf('(%d subject folder%s found)', numel(sf), pls);

% if some folders found
if ~isempty(sf)

    % enabled folder OK but not files
    cf.SetGroupEnabled('FoldOK', 'on');
    cf.SetGroupEnabled('FilesOK', 'off');
end


function spp_sf_browse(varargin)
global ne_ui;
cf = ne_ui.spmx_prepro.hFig;
ch = ne_ui.spmx_prepro.hTag;

% request dir name
nd = uigetdir(pwd, 'Please select a subject folder...');

% cancelled or empty
if isequal(nd, 0) || ...
    isempty(nd)

    % simply return
    return;
end

% otherwise set edit field
ch.ED_spmxprepro_sf.String = nd;

% and check whether exists
if exist(nd, 'dir') > 0
    ff = 'on';
    nd = {nd};
    nf = '(1 subject folder found)';
    [sfsf, sf] = fileparts(nd{1});
    sf = {sf};
else
    ff = 'off';
    nd = {};
    nf = '(0 subject folders found)';
    sf = {};
end

% then set flags
cf.SetGroupEnabled('FoldOK', ff);
cf.SetGroupEnabled('FilesOK', 'off');
ch.TX_spmxprepro_subffound.String = nf;
ne_ui.spmx_prepro.sfld = nd;
ne_ui.spmx_prepro.vids = [nd(:), sf(:)];


function spp_sf_search(varargin)
global ne_ui;
sx = ne_ui.spmx_prepro;
cf = sx.hFig;
ch = sx.hTag;

% set cursor to hourglass
cf.Pointer = 'watch';
drawnow;

% default: disable Create buttons
cf.SetGroupEnabled('FilesOK', 'off');

% get anat/func file patterns
af = ch.ED_spmxprepro_ap.String;
ff = ch.ED_spmxprepro_fp.String;
afp = ch.ED_spmxprepro_afp.String;
ffp = ch.ED_spmxprepro_ffp.String;
ddbl = (ch.CB_spmxprepro_dedouble.Value > 0);

% number of skipped volumes (for computation)
sx.skip = ch.DD_spmxprepro_skipvols.Value - 1;

% make sure they are valid
af(af == ' ') = [];
ff(ff == ' ') = [];
afp(afp == ' ') = [];
ffp(ffp == ' ') = [];

% and expand to lists
afp = splittocellc(afp, ',;', true, true);
ffp = splittocellc(ffp, ',;', true, true);

% iterate over folders
sf = sx.sfld;
vi = sx.vids;
sff = sf;
sfr = sf;
nsub = 0;
nses = 0;
anslc = 0;
snslc = zeros(numel(sf), 1);
for fc = numel(sf):-1:1

    % de-double?
    if ddbl && ...
       ~isempty(spp_findfiles(sf{fc}, {'*.vtc'}))
        afld = {};

    % need to do the work
    else

        % find anatomical folders
        afld = spp_findfolders(sf{fc}, af);
    end

    % if no folders
    if isempty(afld)

        % remove from list and continue
        sf(fc) = [];
        sff(fc) = [];
        sfr(fc) = [];
        vi(fc, :) = [];
        snslc(fc) = [];
        continue;
    end

    % find anatomical files
    afls = spp_findfiles(afld, afp);

    % if no file
    if isempty(afls)

        % equally remove from list and go on
        sf(fc) = [];
        sff(fc) = [];
        sfr(fc) = [];
        vi(fc, :) = [];
        snslc(fc) = [];
        continue;
    end

    % find functional files and test number of slices
    ffld = spp_findfolders(sf{fc}, ff);
    if isempty(ffld)
        sf(fc) = [];
        sff(fc) = [];
        sfr(fc) = [];
        vi(fc, :) = [];
        snslc(fc) = [];
        continue;
    end
    [ffld, nslc] = spp_findfiles(ffld, ffp, 0);
    if isempty(ffld) || ...
        any(diff(nslc))
        sf(fc) = [];
        sff(fc) = [];
        sfr(fc) = [];
        vi(fc, :) = [];
        snslc(fc) = [];
        continue;
    end
    sfr{fc} = ffld;
    snslc(fc) = nslc(1);

    % add to number of runs
    nsub = nsub + 1;
    nses = nses + numel(ffld);

    % and make a nice filename
    filec = 0;
    for sc = 1:numel(ffld)
        filec = filec + numel(ffld{sc}) - sx.skip;
    end
    sff{fc} = sprintf('%s (%d volumes)', ffld{1}{1}(1:end-2), filec);

    % update VTC IDs
    if ~iscell(vi{fc, 1})
        vi{fc, 1} = repmat(vi(fc, 1), numel(ffld), 1);
        vi{fc, 2} = repmat(vi(fc, 2), numel(ffld), 1);
    else
        vi{fc, 1} = repmat(vi{fc, 1}(1), numel(ffld), 1);
        vi{fc, 2} = repmat(vi{fc, 2}(1), numel(ffld), 1);
    end

    % check number of slices
    if anslc == 0
        anslc = nslc(1);
    else
        if anslc ~= nslc(1)
            anslc = -1;
        end
    end
end

% update configuration
sx.af = af;
sx.ff = ff;
sx.afp = afp;
sx.ffp = ffp;
sx.nslc = anslc;
sx.sfld = sf;
sx.sfls = sfr;
sx.vids = vi;
ne_ui.spmx_prepro = sx;

% update UI
ch.TX_spmxprepro_subffound.String = ...
    sprintf('(%d valid subjects/%d runs found)', numel(sf), nses);
ch.LB_spmxprepro_sf_found.String = sff(:);

% update slice options
if anslc > 0
    ch.DD_spmxprepro_str.Value = max(1, min(ch.DD_spmxprepro_str.Value, anslc));
    dstr = '%d,%d, ... %d,%d, ...';
    if mod(anslc, 2) == 0
        dio = sprintf(dstr, anslc - 1, anslc - 3, anslc, anslc - 2);
        die = sprintf(dstr, anslc, anslc - 2, anslc - 1, anslc - 3);
    else
        dio = sprintf(dstr, anslc, anslc - 2, anslc - 1, anslc - 3);
        die = sprintf(dstr, anslc - 1, anslc - 3, anslc, anslc - 2);
    end
    slg = round(sqrt(anslc));
    if slg > 2
        slg1 = 1:slg:anslc;
        slg2 = 2:slg:anslc;
        slge = slg:slg:anslc;
        dsqr = sprintf('%d,%d, ... %d,%d,%d ... %d', ...
            slg1(1), slg1(2), slg1(end), slg2(1), slg2(2), slge(end));
    else
        slg1 = 1:slg:anslc;
        slg2 = 2:slg:anslc;
        dsqr = sprintf('%d,%d, ... %d,%d,%d ... %d', ...
            slg1(1), slg1(2), slg1(end), slg2(1), slg2(2), slg2(end));
    end
    ch.DD_spmxprepro_sto.String = {
        '1,2,3,4,5,6, ...'; ...
        sprintf('%d,%d,%d,%d, ...', anslc - [0, 1, 2, 3]); ...
        '1,3,5, ... 2,4,6, ...'; ...
        '2,4,6, ... 1,3,5, ...'; ...
        dio; die; dsqr};
    % ch.DD_spmxprepro_sto.Value = 3; (don't set to default!)
    ch.DD_spmxprepro_str.String = ...
        splittocellc(deblank(sprintf('%d ', 1:anslc)), ' ');

    cf.SetGroupEnabled('FilesOK', 'on');

% warn of things
elseif anslc == 0
    uiwait(warndlg('No valid runs found!', 'spmx_preprojobs - Info', 'modal'));
else
    uslc = unique(snslc);
    usln = zeros(numel(uslc), 1);
    eslc = cell(numel(uslc), 1);
    for sc = 1:numel(uslc)
        usln(sc) = sum(snslc == uslc(sc));
        idslc = sprintf('  %s\n', sf{snslc == uslc(sc)});
        eslc{sc} = sprintf('%d subjects with %d slices:\n%s', usln(sc), uslc(sc), idslc);
    end
    mslc = uslc(maxpos(usln));
    eslc(maxpos(usln)) = [];
    slmerror = sprintf(['Number of slices mismatch!\n\n', ...
        '%d subjects with %d slices and\n%s\nPlease preprocess separately!'], ...
        max(usln), mslc, gluetostring(eslc, ''));
    uiwait(warndlg(slmerror, 'spmx_preprojobs - Info', 'modal'));
end

% set pointer back to arrow
cf.Pointer = 'arrow';


function spp_advanced(varargin)
global ne_ui;
sx = ne_ui.spmx_prepro;
ch = sx.hTag;

if ch.CB_spmxprepro_advanced.Value > 0
    sx.hFig.Position(3) = 960;
else
    sx.hFig.Position(3) = 634;
end

function spp_sf_delete(varargin)
global ne_ui;
sx = ne_ui.spmx_prepro;
ch = sx.hTag;

todel = ch.LB_spmxprepro_sf_found.Value;
tdstr = ch.LB_spmxprepro_sf_found.String;
if ~iscell(tdstr)
    tdstr = cellstr(tdstr);
end
tdstr(todel) = [];
sx.sfld(todel) = [];
sx.sfls(todel) = [];
sx.vids(todel, :) = [];
ne_ui.spmx_prepro = sx;
ch.LB_spmxprepro_sf_found.ListboxTop = max(1, min(min(todel) - 9, numel(tdstr)));
ch.LB_spmxprepro_sf_found.String = tdstr;
ch.LB_spmxprepro_sf_found.Value = [];


function spp_sto(varargin)
global ne_ui;
ch = ne_ui.spmx_prepro.hTag;

% get selected slice order
stolist = ch.DD_spmxprepro_sto.String;
if ~iscell(stolist)
    stolist = cellstr(stolist);
end
stolist = stolist{ch.DD_spmxprepro_sto.Value};

% get first element
stolist = regexprep(stolist, ',.*$', '');

% set as reference slice
ch.DD_spmxprepro_str.Value = str2double(stolist);


function spp_stom(varargin)
global ne_ui;
ch = ne_ui.spmx_prepro.hTag;

% set enabled state
if ch.CB_spmxprepro_stom.Value > 0
    ef = 'on';
    es = deblank(sprintf('%d ', ...
        spp_sliceorder(ne_ui.spmx_prepro.nslc, ch.DD_spmxprepro_sto.Value)));
else
    ef = 'off';
    es = ' <pre-specified>';
end
ch.ED_spmxprepro_stol.Enable = ef;
ch.ED_spmxprepro_stol.String = es;


function spp_stol(varargin)
global ne_ui;
ch = ne_ui.spmx_prepro.hTag;

% check string
es = ch.ED_spmxprepro_stol.String;
ns = ne_ui.spmx_prepro.nslc;
try
    esn = round(min(ns, max(1, eval(['[' es ']']))));
    if numel(unique(esn)) ~= ns
        error('Wrong number of slices.');
    end
    ch.ED_spmxprepro_stol.String = deblank(sprintf('%d ', esn));
catch ne_eo;
    rethrow(ne_eo);
end


function spp_sttaa(varargin)
global ne_ui;
ch = ne_ui.spmx_prepro.hTag;

% depending on setting
if ch.CB_spmxprepro_sttaa.Value > 0
    ef = 'off';
    es = '<auto>';
else
    ef = 'on';
    try
        tr = str2double(ch.ED_spmxprepro_sttr.String);
        ta = tr * (1 - 1 / ne_ui.spmx_prepro.nslc);
        es = sprintf('%.3f', ta);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        es = '<bad-TR>';
    end
end
ch.ED_spmxprepro_stta.String = es;
ch.ED_spmxprepro_stta.Enable = ef;


function spp_sttr(varargin)
global ne_ui;
ch = ne_ui.spmx_prepro.hTag;

% check string
es = ch.ED_spmxprepro_sttr.String;
try
    esn = str2double(es);
    if numel(esn) ~= 1 || ...
        isinf(esn) || ...
        isnan(esn) || ...
        esn <= 0 || ...
        esn > 120
        error('Bad TR value.');
    end
    ch.ED_spmxprepro_sttr.String = sprintf('%.3f', esn);
catch ne_eo;
    rethrow(ne_eo);
end


function spp_reslice(varargin)
global ne_ui;
ch = ne_ui.spmx_prepro.hTag;

% make settings
if ch.CB_spmxprepro_reslice.Value > 0
    ch.DD_spmxprepro_rinterpr.Enable = 'on';
    ch.CB_spmxprepro_cupr.Enable = 'on';
else
    ch.DD_spmxprepro_rinterpr.Enable = 'off';
    ch.CB_spmxprepro_cupr.Enable = 'off';
end


function spp_nrmtype(varargin)
global ne_ui;
sx = ne_ui.spmx_prepro;
cf = sx.hFig;
ch = sx.hTag;

% set groups status
ntypes = ch.DD_spmxprepro_nrmtype.String;
if ~iscell(ntypes)
    ntypes = cellstr(ntypes);
end
if ch.DD_spmxprepro_nrmtype.Value < numel(ntypes)
    cf.SetGroupEnabled('NrANA', 'on');
    cf.SetGroupEnabled('NrEPI', 'off');
    ch.CB_spmxprepro_nrmsess.Value = 0;
    if ch.DD_spmxprepro_nrmtype.Value == (numel(ntypes) - 1)
        cf.SetGroupEnabled('NrT1', 'on');
    else
        cf.SetGroupEnabled('NrT1', 'off');
    end
    if numel(ntypes) == 4
        if ch.DD_spmxprepro_nrmtype.Value == 2
            ch.ED_spmxprepro_vtcfpat.String = strrep( ...
                ch.ED_spmxprepro_vtcfpat.String, '_MNI.', '_DMNI.');
            ch.CB_spmxprepro_reavtc.Enable = 'off';
            ch.CB_spmxprepro_reavtc.Value = 0;
        else
            ch.ED_spmxprepro_vtcfpat.String = strrep( ...
                ch.ED_spmxprepro_vtcfpat.String, '_DMNI.', '_MNI.');
            ch.CB_spmxprepro_reavtc.Enable = 'on';
        end
    end
else
    cf.SetGroupEnabled('NrANA', 'off');
    cf.SetGroupEnabled('NrEPI', 'on');
    ch.CB_spmxprepro_writesn.Value = 0;
    ch.CB_spmxprepro_reavtc.Enable = 'on';
    if numel(ntypes) == 4
        ch.ED_spmxprepro_vtcfpat.String = strrep( ...
            ch.ED_spmxprepro_vtcfpat.String, '_DMNI.', '_MNI.');
    end
end


function spp_t1browse(varargin)
global ne_ui;
sx = ne_ui.spmx_prepro;
ch = sx.hTag;

% request dir name
[t1file, t1path] = uigetfile( ...
    {'*.img;*.nii', 'Analyze/NIftI files (*.img, *.nii)'}, ...
    'Please select a T1 template file...');

% cancelled or empty
if isempty(t1file) || ...
    isequal(t1file, 0)

    % simply return
    return;
end
if isempty(t1path)
    t1path = pwd;
end
t1file = strrep(strrep([t1path '/' t1file], '\', '/'), '//', '/');

% exists
if exist(t1file, 'file') ~= 2
    return;
end

% otherwise set edit field
ch.ED_spmxprepro_t1temp.String = t1file;


function spp_loadcfg(varargin)
global ne_ui;
sx = ne_ui.spmx_prepro;
hFig = sx.hFig;
hTag = sx.hTag;

% file dialog
[cfgfile, cfgpath] = uigetfile( ...
    {'*.ini', 'SPM preprocessing configuration (*.ini)'}, ...
    'Please select the preprocessing configuration...');
if isequal(cfgfile, 0) || ...
    isequal(cfgpath, 0) || ...
    isempty(cfgfile)
    return;
end
if isempty(cfgpath)
    cfgpath = pwd;
end
cfgfile = [cfgpath '/' cfgfile];
if exist(cfgfile, 'file') ~= 2
    return;
end

% try to load file
cfgini = [];
try
    cfgini = xini(cfgfile, 'convert');
    cc = cfgini.GetComplete;
    cfgini.Release;
    cfgini = [];
    if ~isfield(cc, 'f') || ...
       ~isstruct(cc.f) || ...
        numel(cc.f) ~= 1 || ...
       ~isfield(cc, 'x') || ...
       ~isstruct(cc.x) || ...
        numel(cc.x) ~= 1 || ...
       ~isfield(cc, 'k') || ...
       ~isstruct(cc.k) || ...
        numel(cc.k) ~= 1 || ...
       ~isfield(cc, 'a') || ...
       ~isstruct(cc.a) || ...
        numel(cc.a) ~= 1 || ...
       ~isfield(cc, 'r') || ...
       ~isstruct(cc.r) || ...
        numel(cc.r) ~= 1 || ...
       ~isfield(cc, 'w') || ...
       ~isstruct(cc.w) || ...
        numel(cc.w) ~= 1 || ...
       ~isfield(cc, 's') || ...
       ~isstruct(cc.s) || ...
        numel(cc.s) ~= 1 || ...
       ~isfield(cc, 'o') || ...
       ~isstruct(cc.o) || ...
        numel(cc.o) ~= 1
        error( ...
            'neuroelf:BadIniContent', ...
            'Invalid SPM preprocessing configuration.' ...
        );
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    if ~isempty(cfgini)
        cfgini.Release;
    end
    return;
end

% make settings
try
    ne_ui.spmx_prepro = struct( ...
        'af',   cc.f.ap, ...
        'afp',  {cc.f.afp}, ...
        'cupa', cc.a.cupa, ...
        'cupk', cc.k.cupk, ...
        'cupr', cc.r.cupr, ...
        'cups', cc.s.cups, ...
        'cupw', cc.w.cupw, ...
        'cupx', cc.x.cupx, ...
        'exbr', cc.x.bex, ...
        'ff',   cc.f.fp, ...
        'ffp',  {cc.f.ffp}, ...
        'func', {{}}, ...
        'hFig', hFig, ...
        'hTag', hTag, ...
        'jobs', {{}}, ...
        'jobh', {{}}, ...
        'nslc', cc.a.nslc, ...
        'qual', cc.o.qual, ...
        'quav', false, ...
        'rff',  {{}}, ...
        'rfs',  {{}}, ...
        'runj', false, ...
        'rvtc', cc.o.rvtc, ...
        'sfld', {{}}, ...
        'skip', cc.k.n, ...
        'spmv', sx.spmv, ...
        'vids', {{}}, ...
        'vtci', cc.o.vtci, ...
        'vtcp', cc.o.vtcf, ...
        'vtcs', cc.o.vtc, ...
        'wsna', cc.w.wsna);

    % override T1 template
    if ~isempty(cc.w.t1t)
        ne_ui.spmx_prepro.t1t = cc.w.t1t;
    end

    % initialize controls; folders
    hTag.ED_spmxprepro_sf.String = ddeblank(cc.f.sf);
    hTag.ED_spmxprepro_ap.String = ['  ' cc.f.ap];
    hTag.ED_spmxprepro_afp.String = ['  ' gluetostringc(cc.f.afp, ', ')];
    hTag.ED_spmxprepro_fp.String = ['  ' cc.f.fp];
    hTag.ED_spmxprepro_ffp.String = ['  ' gluetostringc(cc.f.ffp, ', ')];
    hTag.ED_spmxprepro_dp.String = ['  ' cc.f.dp];
    hTag.ED_spmxprepro_dfp.String = ['  ' gluetostringc(cc.f.dfp, ', ')];

    % data handling
    if isfield(cc.x, 'ddbl')
        hTag.CB_spmxprepro_dedouble.Value = double(cc.x.ddbl);
    else
        hTag.CB_spmxprepro_dedouble.Value = 1;
    end
    hTag.CB_spmxprepro_single.Value = double(cc.x.sng);
    hTag.CB_spmxprepro_bextract.Value = double(cc.x.bex);
    hTag.CB_spmxprepro_reganat.Value = double(cc.x.rt1);
    hTag.CB_spmxprepro_cupx.Value = double(cc.x.cupx);

    % skipping volumes
    hTag.DD_spmxprepro_skipvols.Value = max(1, min(11, 1 + cc.k.n));
    hTag.CB_spmxprepro_cupk.Value = double(cc.k.cupk);

    % set found folders to empty
    hTag.TX_spmxprepro_subffound.String = '(0 subjects found)';
    hTag.LB_spmxprepro_sf_found.Max = 3;
    hTag.LB_spmxprepro_sf_found.Value = [];
    hTag.LB_spmxprepro_sf_found.String = {};

    % slice-timing
    hTag.CB_spmxprepro_st.Value = double(cc.a.do);
    hTag.CB_spmxprepro_cupa.Value = double(cc.a.cupa);
    nslc = cc.a.nslc;
    oslc = 1:2:nslc;
    eslc = 2:2:nslc;
    nslcsq = round(sqrt(nslc));
    sslc = 1:nslcsq:nslc;
    hTag.DD_spmxprepro_sto.String = { ...
        '1,2,3,4,5,6, ...'; ...
        sprintf('%d,%d,%d,%d, ...', nslc - [0, 1, 2, 3]); ...
        '1,3,5, ... 2,4,6, ...'; ...
        '2,4,6, ... 1,3,5, ...'; ...
        sprintf('%d,%d, ... %d,%d, ...', oslc(end), oslc(end-1), eslc(end), eslc(end-1)); ...
        sprintf('%d,%d, ... %d,%d, ...', eslc(end), eslc(end-1), oslc(end), oslc(end-1)); ...
        sprintf('%d,%d, ... %d,%d,%d, ...', 1, 1 + nslcsq, sslc(end), 2, 2 + nslcsq)};
    switch lower(cc.a.o)
        case {'a'}
            hTag.DD_spmxprepro_sto.Value = 1;
        case {'aie'}
            hTag.DD_spmxprepro_sto.Value = 4;
        case {'aio'}
            hTag.DD_spmxprepro_sto.Value = 3;
        case {'d'}
            hTag.DD_spmxprepro_sto.Value = 2;
        case {'die'}
            hTag.DD_spmxprepro_sto.Value = 6;
        case {'dio'}
            hTag.DD_spmxprepro_sto.Value = 5;
        case {'sqrt'}
            hTag.DD_spmxprepro_sto.Value = 7;
        otherwise
            hTag.DD_spmxprepro_sto.Value = 3;
    end
    hTag.CB_spmxprepro_stom.Value = double(cc.a.om);
    feval(hTag.CB_spmxprepro_stom.Callback);
    if cc.a.om
        sllist = sprintf('%d, ', cc.a.ol);
        hTag.ED_spmxprepro_stol.String = ['  ' sllist(1:end-1)];
    else
        hTag.ED_spmxprepro_stol.String = ' <pre-specified>';
    end
    hTag.ED_spmxprepro_sttr.String = sprintf('  %.3f', cc.a.tr);
    hTag.CB_spmxprepro_sttaa.Value = double(cc.a.taa);
    feval(hTag.CB_spmxprepro_sttaa.Callback);
    if cc.a.taa
        hTag.ED_spmxprepro_stta.String = '<auto>';
    else
        hTag.ED_spmxprepro_stta.String = sprintf('%.3f', cc.a.ta);
    end
    hTag.DD_spmxprepro_str.String = ...
        splittocellc(deblank(sprintf('%d ', 1:nslc)), ' ');
    hTag.DD_spmxprepro_str.Value = cc.a.r;
    hTag.DD_spmxprepro_mux.Value = cc.a.mux;

    % realignment
    if cc.r.dir(end) == 'c'
        hTag.DD_spmxprepro_regtype.Value = 2;
    else
        hTag.DD_spmxprepro_regtype.Value = 1;
    end
    hTag.ED_spmxprepro_rqual.String = sprintf('%g', cc.r.qual);
    hTag.ED_spmxprepro_rsep.String = sprintf('%g', cc.r.sep);
    hTag.CB_spmxprepro_rrtm.Value = double(cc.r.rtm);
    hTag.DD_spmxprepro_rinterpe.Value = cc.r.ie;
    hTag.DD_spmxprepro_rinterpr.Value = cc.r.ir;
    hTag.CB_spmxprepro_reslice.Value = double(cc.r.res);
    hTag.CB_spmxprepro_cupr.Value = double(cc.r.cupr);

    % normalization
    switch lower(cc.w.type)
        case {'dartel'}
            if numel(hTag.DD_spmxprepro_nrmtype.String) > 3
                hTag.DD_spmxprepro_nrmtype.Value = 2;
            else
                hTag.DD_spmxprepro_nrmtype.Value = 1;
            end
        case {'normana'}
            if numel(hTag.DD_spmxprepro_nrmtype.String) > 3
                hTag.DD_spmxprepro_nrmtype.Value = 3;
            else
                hTag.DD_spmxprepro_nrmtype.Value = 2;
            end
        case {'normepi'}
            if numel(hTag.DD_spmxprepro_nrmtype.String) > 3
                hTag.DD_spmxprepro_nrmtype.Value = 4;
            else
                hTag.DD_spmxprepro_nrmtype.Value = 3;
            end
        case {'seg'}
            hTag.DD_spmxprepro_nrmtype.Value = 1;
    end
    feval(hTag.DD_spmxprepro_nrmtype.Callback);
    if ~isempty(cc.w.t1t)
        hTag.ED_spmxprepro_t1temp.String = ['  ' cc.w.t1t];
    end
    hTag.ED_spmxprepro_wbb.String = ...
        sprintf('  [%d, %d, %d; %d, %d, %d]', lsqueeze(cc.w.bb'));
    wvsz = hTag.DD_spmxprepro_wvox.String;
    if ~iscell(wvsz)
        wvsz = cellstr(wvsz);
    end
    for sc = 1:numel(wvsz)
        wvsz{sc} = str2double(ddeblank(wvsz{sc}));
    end
    wvsz = cat(1, wvsz{:});
    hTag.DD_spmxprepro_wvox.Value = minpos(abs(wvsz - cc.w.vox));
    hTag.CB_spmxprepro_writesn.Value = double(cc.w.wsna);
    hTag.CB_spmxprepro_nrmsess.Value = double(cc.w.sep);
    hTag.CB_spmxprepro_cupw.Value = double(cc.w.cupw);

    % smoothing
    hTag.ED_spmxprepro_smk.String = sprintf('%.1f', cc.s.k);
    hTag.CB_spmxprepro_sm.Value = double(cc.s.do);
    hTag.CB_spmxprepro_cups.Value = double(cc.s.cups);

    % quality and import
    hTag.CB_spmxprepro_fquality.Value = double(cc.o.qual);
    hTag.CB_spmxprepro_fqualvtc.Value = double(cc.o.quav);
    hTag.CB_spmxprepro_reavtc.Value = double(cc.o.rvtc);
    hTag.CB_spmxprepro_autovtc.Value = double(cc.o.vtc);
    hTag.ED_spmxprepro_vtcfpat.String = cc.o.vtcf;
    hTag.ED_spmxprepro_vtcfidx.String = sprintf('%d', round(max(1, cc.o.vtci)));

    % make sure OK folder/file groups are disabled
    hFig.SetGroupEnabled('FoldOK', 'off');
    hFig.SetGroupEnabled('FilesOK', 'off');

    % evaluate callback for search
    feval(hTag.ED_spmxprepro_sf.Callback);
catch ne_eo;
    uiwait(warndlg(ne_eo.message, 'NeuroElf - error', 'modal'));
    return;
end


function spp_savecfg(varargin)
global ne_ui;
sx = ne_ui.spmx_prepro;
hTag = sx.hTag;

% file dialog
[cfgfile, cfgpath] = uiputfile( ...
    {'*.ini', 'SPM preprocessing configuration (*.ini)'}, ...
    'Please select the preprocessing configuration...', 'spmx_prepro.ini');
if isequal(cfgfile, 0) || ...
    isequal(cfgpath, 0) || ...
    isempty(cfgfile)
    return;
end
if isempty(cfgpath)
    cfgpath = pwd;
end
cfgfile = [cfgpath '/' cfgfile];

% read original file
try
    cfgini = xini( ...
        [neuroelf_path filesep '_core' filesep 'config' filesep 'spmx_prepro.ini'], ...
        'convert');

    % save under new name
    cfgini.SaveAs(cfgfile);

    % and then make settings
    cfgini.f.sf = ddeblank(hTag.ED_spmxprepro_sf.String);
    cfgini.f.ap = ddeblank(hTag.ED_spmxprepro_ap.String);
    cfgini.f.afp = ddeblank(splittocellc(hTag.ED_spmxprepro_afp.String, ',;', true, true));
    cfgini.f.fp = ddeblank(hTag.ED_spmxprepro_fp.String);
    cfgini.f.ffp = ddeblank(splittocellc(hTag.ED_spmxprepro_ffp.String, ',;', true, true));
    cfgini.f.dp = ddeblank(hTag.ED_spmxprepro_dp.String);
    cfgini.f.dfp = ddeblank(splittocellc(hTag.ED_spmxprepro_dfp.String, ',;', true, true));
    cfgini.x.ddbl = (hTag.CB_spmxprepro_dedouble.Value > 0);
    cfgini.x.sng = (hTag.CB_spmxprepro_single.Value > 0);
    cfgini.x.bex = (hTag.CB_spmxprepro_bextract.Value > 0);
    cfgini.x.rt1 = (hTag.CB_spmxprepro_reganat.Value > 0);
    cfgini.x.cupx = (hTag.CB_spmxprepro_cupx.Value > 0);
    cfgini.k.n = hTag.DD_spmxprepro_skipvols.Value - 1;
    cfgini.k.cupk = (hTag.CB_spmxprepro_cupk.Value > 0);
    cfgini.a.do = (hTag.CB_spmxprepro_st.Value > 0);
    cfgini.a.nslc = sx.nslc;
    slorders = {'a', 'd', 'aio', 'aie', 'dio', 'die', 'sqrt'};
    cfgini.a.o = slorders{hTag.DD_spmxprepro_sto.Value};
    cfgini.a.om = (hTag.CB_spmxprepro_stom.Value > 0);
    if cfgini.a.om
        cfgini.a.ol = u8str2double(hTag.ED_spmxprepro_stol.String);
    else
        cfgini.a.ol = [];
    end
    cfgini.a.tr = str2double(ddeblank(hTag.ED_spmxprepro_sttr.String));
    if hTag.CB_spmxprepro_sttaa.Value > 0
        cfgini.a.ta = [];
        cfgini.a.taa = true;
    else
        cfgini.a.ta = str2double(ddeblank(hTag.ED_spmxprepro_stta.String));
        cfgini.a.taa = false;
    end
    cfgini.a.r = hTag.DD_spmxprepro_str.Value;
    cfgini.a.mux = hTag.DD_spmxprepro_mux.Value;
    cfgini.a.cupa = (hTag.CB_spmxprepro_cupa.Value > 0);
    rdirs = {'toana', 'tofunc'};
    cfgini.r.dir = rdirs{hTag.DD_spmxprepro_regtype.Value};
    cfgini.r.qual = str2double(ddeblank(hTag.ED_spmxprepro_rqual.String));
    cfgini.r.sep = str2double(ddeblank(hTag.ED_spmxprepro_rsep.String));
    cfgini.r.rtm = (hTag.CB_spmxprepro_rrtm.Value > 0);
    cfgini.r.ie = hTag.DD_spmxprepro_rinterpe.Value;
    cfgini.r.ir = hTag.DD_spmxprepro_rinterpr.Value;
    cfgini.r.res = (hTag.CB_spmxprepro_reslice.Value > 0);
    cfgini.r.cupr = (hTag.CB_spmxprepro_cupr.Value > 0);
    if numel(hTag.DD_spmxprepro_nrmtype.String) > 3
        nrmtypes = {'seg', 'dartel', 'normana', 'normepi'};
    else
        nrmtypes = {'seg', 'normana', 'normepi'};
    end
    cfgini.w.type = nrmtypes{hTag.DD_spmxprepro_nrmtype.Value};
    if strcmpi(hTag.ED_spmxprepro_t1temp.Enable, 'on')
        cfgini.w.t1t = hTag.ED_spmxprepro_t1temp.String;
    else
        cfgini.w.t1t = '';
    end
    cfgini.w.bb = u8str2double(hTag.ED_spmxprepro_wbb.String);
    wvsz = hTag.DD_spmxprepro_wvox.String;
    if ~iscell(wvsz)
        wvsz = cellstr(wvsz);
    end
    cfgini.w.vox = str2double(ddeblank(wvsz{hTag.DD_spmxprepro_wvox.Value}));
    cfgini.w.wsna = (hTag.CB_spmxprepro_writesn.Value > 0);
    cfgini.w.sep = (hTag.CB_spmxprepro_nrmsess.Value > 0);
    cfgini.w.cupw = (hTag.CB_spmxprepro_cupw.Value > 0);
    cfgini.s.do = (hTag.CB_spmxprepro_sm.Value > 0);
    cfgini.s.k = str2double(ddeblank(hTag.ED_spmxprepro_smk.String));
    cfgini.s.cups = (hTag.CB_spmxprepro_cups.Value > 0);
    cfgini.o.qual = (hTag.CB_spmxprepro_fquality.Value > 0);
    cfgini.o.quav = (hTag.CB_spmxprepro_fqualvtc.Value > 0);
    cfgini.o.rvtc = (hTag.CB_spmxprepro_reavtc.Value > 0);
    cfgini.o.vtc = (hTag.CB_spmxprepro_autovtc.Value > 0);
    cfgini.o.vtcf = ddeblank(hTag.ED_spmxprepro_vtcfpat.String);
    cfgini.o.vtci = max(1, round(str2double(ddeblank(hTag.ED_spmxprepro_vtcfidx.String))));

    % then save and release
    cfgini.Save;
    cfgini.Release;
catch ne_eo;
    uiwait(warndlg(['Error saving config file: ' ne_eo.message], ...
        'NeuroElf - error', 'modal'));
    return;
end


function spp_create(varargin)
global ne_ui;
sx = ne_ui.spmx_prepro;
ch = sx.hTag;

% nothing to do after all
if isempty(sx.sfld)
    return;
end

% plausibility check
try
    rq = str2double(ch.ED_spmxprepro_rqual.String);
    if numel(rq) ~= 1 || ...
        isinf(rq) || ...
        isnan(rq) || ...
        rq < 0.1 || ...
        rq > 1
        error('Bad rqual value.');
    end
    rs = str2double(ch.ED_spmxprepro_rsep.String);
    if numel(rs) ~= 1 || ...
        isinf(rs) || ...
        isnan(rs) || ...
        rs < 1 || ...
        rs > 16
        error('Bad rsep value.');
    end
    sk = str2double(ch.ED_spmxprepro_smk.String);
    if ch.CB_spmxprepro_sm.Value <= 0
        sk = 0;
    end
    if numel(sk) ~= 1 || ...
        isinf(sk) || ...
        isnan(sk) || ...
        sk < 0 || ...
        sk > 15
        error('Bad smk value.');
    end
    tr = str2double(ch.ED_spmxprepro_sttr.String);
    if numel(tr) ~= 1 || ...
        isinf(tr) || ...
        isnan(tr) || ...
        tr < 0.025 || ...
        tr > 60
        error('Bad TR value.');
    end
    st = (ch.CB_spmxprepro_st.Value > 0);
    if st
        if ch.CB_spmxprepro_stom.Value <= 0
            so = spp_sliceorder(sx.nslc, ch.DD_spmxprepro_sto.Value);
        else
            try
                so = round(max(1, min(sx.nslc, ...
                    eval(['[' ch.ED_spmxprepro_stol.String ']']))));
            catch ne_eo;
                rethrow(ne_eo);
            end
            if numel(unique(so)) ~= sx.nslc || ...
                numel(so) ~= sx.nslc
                error('Non-unique or missing slices in slice order.');
            end
        end
        if ch.CB_spmxprepro_sttaa.Value > 0
            ta = tr * (1 - 1 / sx.nslc);
        else
            ta = str2double(ch.ED_spmxprepro_stta.String);
            if numel(ta) ~= 1 || ...
                isinf(ta) || ...
                isnan(ta) || ...
                tr < 0.1 || ...
                ta >= tr
                error('Bad TA value.');
            end
        end
    else
        so = 1:sx.nslc;
        ta = [];
    end
    try
        wb = eval(['[' ch.ED_spmxprepro_wbb.String ']']);
    catch ne_eo;
        rethrow(ne_eo);
    end
    if ~isequal(size(wb), [2, 3]) || ...
        any(isinf(wb(:)) | isnan(wb(:)) | wb(:) < -128 | wb(:) > 128) || ...
        any(diff(wb) <= 0)
        error('Invalid wbb definition.');
    end
catch ne_eo;
    uiwait(warndlg(ne_eo.message, 'spmx_preprojobs - Info', 'modal'));
    return;
end

% set run flag
sx.runj = (varargin{3} > 0);

% options
ntypes = ch.DD_spmxprepro_nrmtype.String;
if ~iscell(ntypes)
    ntypes = cellstr(ntypes);
end
if numel(ntypes) == 4
    ntypes = {'seg', 'dartel', 'anat', 'epi'};
else
    ntypes = {'seg', 'anat', 'epi'};
end
opts = struct( ...
    'afp', {sx.afp}, ...
    'extbrain', (ch.CB_spmxprepro_bextract.Value > 0), ...
    'ffp', {sx.ffp}, ...
    'fsingle', (ch.CB_spmxprepro_single.Value > 0), ...
    'fun2str', (ch.DD_spmxprepro_regtype.Value == 1), ...
    'jobfile', '', ...
    'jobrun', false, ...
    'normsess', (ch.CB_spmxprepro_nrmsess.Value > 0), ...
    'normtype', ntypes{ch.DD_spmxprepro_nrmtype.Value}, ...
    'reganat', (ch.CB_spmxprepro_reganat.Value > 0), ...
    'reslice', (ch.CB_spmxprepro_reslice.Value > 0), ...
    'rinterpe', ch.DD_spmxprepro_rinterpe.Value, ...
    'rinterpr', ch.DD_spmxprepro_rinterpr.Value, ...
    'rqual', rq, ...
    'rrtm', (ch.CB_spmxprepro_rrtm.Value > 0), ...
    'rsep', rs, ...
    'skip', ch.DD_spmxprepro_skipvols.Value - 1, ...
    'smk', sk, ...
    'st', st, ...
    'sto', so, ...
    'str', ch.DD_spmxprepro_str.Value, ...
    'stta', ta, ...
    'sttr', tr, ...
    't1temp', ch.ED_spmxprepro_t1temp.String, ...
    'wbb', wb, ...
    'winterp', 1, ...
    'wsna', (ch.CB_spmxprepro_writesn.Value > 0), ...
    'wvox', str2double(ch.DD_spmxprepro_wvox.String{ch.DD_spmxprepro_wvox.Value}));

% initiate outputs
jobs = {};
jobhelps = {};
funcfiles = {};
rff = {};
rfs = {};

% set arrow to hour glass
sx.hFig.Pointer = 'watch';

% make buttons invisible
ch.BT_spmxprepro_loadcfg.Visible = 'off';
ch.BT_spmxprepro_savecfg.Visible = 'off';
ch.BT_spmxprepro_createjob.Visible = 'off';
ch.BT_spmxprepro_candrjob.Visible = 'off';
ch.BT_spmxprepro_cancel.Visible = 'off';
pbar = ch.PB_spmxprepro_progress;
pbar.Visible = 'on';
sx.hFig.SetGroupEnabled('AllUIC',  'off');
sx.hFig.SetGroupEnabled('FoldOK',  'off');
sx.hFig.SetGroupEnabled('FilesOK', 'off');
sx.hFig.SetGroupEnabled('NrT1',    'off');
sx.hFig.SetGroupEnabled('NrEPI',   'off');
sx.hFig.SetGroupEnabled('NrANA',   'off');
sx.hFig.SetGroupEnabled('DoVTCs',  'off');

% iterate over folders
try
    for fc = 1:numel(sx.sfld)
        fldname = sx.sfld{fc};
        if numel(fldname) > 80
            fldname = [fldname(1:20) '...' fldname(end-55:end)];
        end
        ch.LB_spmxprepro_sf_found.Value = fc;
        ch.LB_spmxprepro_sf_found.ListboxTop = max(1, fc - 8);
        pbar.Progress((fc - 1) / numel(sx.sfld), sprintf('preparing %s...', fldname));
        drawnow;
        [p1, p2, p3, p4, p5] = spmx_preprojobs(sx.sfld{fc}, sx.ff, sx.af, opts);
        jobs = [jobs(:)', p1(:)'];
        jobhelps = [jobhelps(:)', p2(:)'];
        funcfiles = [funcfiles(:)', p3(:)'];
        rff = [rff(:)', p4(:)'];
        rfs = [rfs(:)', p5(:)'];
    end
catch ne_eo;
    sx.hFig.Pointer = 'arrow';
    uiwait(warndlg(ne_eo.message, 'spmx_preprojobs - Error', 'modal'));
    spp_closeui;
    return;
end

% store back in global struct
sx.cupa = (ch.CB_spmxprepro_cupa.Value > 0);
sx.cupk = (ch.CB_spmxprepro_cupk.Value > 0);
sx.cupr = (ch.CB_spmxprepro_cupr.Value > 0);
sx.cups = (ch.CB_spmxprepro_cups.Value > 0);
sx.cupw = (ch.CB_spmxprepro_cupw.Value > 0);
sx.cupx = (ch.CB_spmxprepro_cupx.Value > 0);
sx.func = funcfiles;
sx.jobh = jobhelps;
sx.jobs = jobs;
sx.qual = (ch.CB_spmxprepro_fquality.Value > 0);
sx.quav = (ch.CB_spmxprepro_fqualvtc.Value > 0);
sx.rff  = rff;
sx.rfs  = rfs;
sx.rvtc = (ch.CB_spmxprepro_reavtc.Value > 0);
sx.skip = opts.skip;
sx.sttr = opts.sttr;
vidf = sx.vids(:, 1);
sx.vidf = cat(1, vidf{:});
vids = sx.vids(:, 2);
sx.vids = cat(1, vids{:});
sx.vtci = str2double(ch.ED_spmxprepro_vtcfidx.String);
sx.vtcp = ch.ED_spmxprepro_vtcfpat.String;
sx.vtcs = (ch.CB_spmxprepro_autovtc.Value > 0);
if ~sx.vtcs
    if sk > 0
        sx.cups = false;
    else
        sx.cupw = false;
    end
end
sx.wbb = wb;
sx.wvox = opts.wvox;
ne_ui.spmx_prepro = sx;

% resume
ch.LB_spmxprepro_sf_found.ListboxTop = 1;
ch.LB_spmxprepro_sf_found.Value = [];
sx.hFig.Pointer = 'arrow';
uiresume(sx.hFig.MLHandle);


% change DoVTCs group enabled status
function spp_dovtcs(varargin)
global ne_ui;
sx = ne_ui.spmx_prepro;
cf = sx.hFig;
ch = sx.hTag;
if ch.CB_spmxprepro_autovtc.Value > 0
    cf.SetGroupEnabled('DoVTCs', 'on');
else
    cf.SetGroupEnabled('DoVTCs', 'off');
    if ch.CB_spmxprepro_sm.Value > 0
        ch.CB_spmxprepro_cups.Value = 0;
    else
        ch.CB_spmxprepro_cupw.Value = 0;
    end
end


% check VTC filename index
function spp_vtcfidx(varargin)
global ne_ui;
sx = ne_ui.spmx_prepro;
ch = sx.hTag;
vtcfi = ddeblank(ch.ED_spmxprepro_vtcfidx.String);
if isempty(vtcfi) || ...
    any(vtcfi < '0' | vtcfi > '9')
    vtcfi = '1';
end
ch.ED_spmxprepro_vtcfidx.String = vtcfi;



% check VTC filename pattern
function spp_vtcfpat(varargin)
global ne_ui;
sx = ne_ui.spmx_prepro;
ch = sx.hTag;
vtcfp = ddeblank(ch.ED_spmxprepro_vtcfpat.String);
if sum(vtcfp == '%') ~= 3 || ...
    isempty(regexpi(vtcfp, '.*\%s.+\%s.+\%\d*d.*\.vtc$'))
    vtcfp = '%s/%s_RUN%02d_MNI.vtc';
end
ch.ED_spmxprepro_vtcfpat.String = vtcfp;



% internal functions



% prepend filenames with letter (string)
function f = fnpp(f, p)

% for list of files
if iscell(f)

    % apply to each and return
    for cc = 1:numel(f)
        f{cc} = fnpp(f{cc}, p);
    end
    return;
end

% for single file get parts and prepend
[fp{1:3}] = fileparts(f);
if isempty(fp{1})
    fp{1} = pwd;
end
f = [fp{1}, filesep, p, fp{2}, fp{3}];


% remove first letter of filename
function f = fnxp(f)

% for list of files
if iscell(f)

    % apply to each (filename or cell of filenames) and return
    for cc = 1:numel(f)
        f{cc} = fnxp(f{cc});
    end
    return;
end

% for single file
[fp{1:3}] = fileparts(f);
if isempty(fp{1})
    f = [fp{2}(2:end), fp{3}];
else
    f = [fp{1}, '/', fp{2}(2:end), fp{3}];
end


% find folders according to spec
function ff = spp_findfolders(sf, fp)

% find folders and files
if any(fp == filesep)
    [ffp, fp, fpe] = fileparts(fp);
    fp = [fp, fpe];
    ffp = [filesep ffp];
else
    ffp = '';
end
ff = findfiles([sf ffp], fp, 'dirs', 'depth=1');


% find files according to spec
function [ff, ns, fd] = spp_findfiles(ff, fp, skip)

% make sure folders are cell
if ~iscell(ff)
    ff = {ff};
end
nf = numel(ff);

% for anatomical files
if nargout == 1

    % iterate over folders
    for sc = 1:nf
        files = findfiles(ff{sc}, fp, 'depth=1');
        if ~isempty(files)
            ff = files(1);
            return;
        end
    end

    % obviously nothing found
    ff = {};

% find all files
else

    % default skip
    if nargin < 3
        skip = 0;
    end

    % iterater over folders
    for sc = 1:nf

        % and replace folder name with found files
        ff{sc} = spp_short(findfiles(ff{sc}, fp, 'depth=1'));
    end

    % remove empty file lists
    ff(cellfun('isempty', ff)) = [];

    % default is 0 slices, and not 4D
    ns = zeros(size(ff));
    fd = false(size(ff));

    % for each file list
    for sc = numel(fd):-1:1

        % number of to-be skipped volumes
        toskip = skip;

        % try to read first file
        try
            first = {};
            fname = ff{sc}{1};

            % and make sure to not use ".img" with xff
            if strcmpi(fname(end-2:end), 'img')
                fname = regexprep(fname, 'img$', 'hdr', 'preservecase');
            end
            first = {xff(fname)};

            % get number of slices and 4-D flag
            ns(sc) = size(first{1}.VoxelData, 3);
            fd(sc) = (size(first{1}.VoxelData, 4) > 1);

            % 4D-volume/s -> expand name/s
            if fd(sc)

                % skip requested
                if toskip > 0

                    % check if file with skipped volumes exists
                    ff{sc}{1} = fnpp(ff{sc}{1}, 'k');
                    skipok = false;

                    % try reading and
                    skipped = {};
                    try
                        skipped = {xff(ff{sc}{1})};
                        if size(skipped{1}.VoxelData, 4) == (size(first{1}.VoxelData, 4) - toskip)
                            skipok = true;
                        end
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                    end
                    clearxffobjects(skipped);

                    % otherwise rewrite file and change filename
                    if ~skipok
                        try
                            first{1}.LoadVoxelData;
                            first{1}.VoxelData(:, :, :, 1:toskip) = [];
                            first{1}.ImgDim.Dim(5) = size(first{1}.VoxelData, 4);
                            first{1}.SaveAs(ff{sc}{1});
                            if isfield(first{1}.RunTimeVars, 'Discard') && ...
                               ~isempty(first{1}.RunTimeVars.Discard)
                                first{1}.RunTimeVars.Discard = ...
                                    first{1}.RunTimeVars.Discard - toskip;
                                first{1}.SaveRunTimeVars;
                            end
                            toskip = 0;
                        catch ne_eo;
                            rethrow(ne_eo);
                        end
                    end
                end

                % repeat first name as often as volumes
                ff{sc} = ff{sc}(ones(size(first{1}.VoxelData, 4) - toskip, 1));

                % and append ",VOL"
                for fc = 1:numel(ff{sc})
                    ff{sc}{fc} = sprintf('%s,%d', ff{sc}{fc}, fc);
                end

            % otherwise, add ',1' to each file
            else
                for fc = 1:numel(ff{sc})
                    ff{sc}{fc} = sprintf('%s,1', ff{sc}{fc});
                end
            end

        % on error
        catch ne_eo;
            neuroelf_lasterr(ne_eo);

            % remove session from list
            ff(sc) = [];
            ns(sc) = [];
            fd(sc) = [];
        end

        % skip files?
        if toskip > 0 && ...
           ~fd(sc)
            ff{sc}(1:min(skip, numel(ff{sc}))) = [];
        end

        % remove object (if necessary)
        clearxffobjects(first);
    end
end


function f = spp_short(f)

% get lengths of filenames
nf = numel(f);
lf = zeros(nf, 1);
for fc = 1:nf
    lf(fc) = numel(f{fc});
end

% return if all the same
if ~any(diff(lf))
    return;
end

% get shortest name
[sp, sn] = fileparts(f{minpos(lf)});
sc = sn(1:4);

% and remove files for which the first letter is different
for fc = nf:-1:1
    [sp, sn] = fileparts(f{fc});
    if any(sn(1:4) ~= sc)
        f(fc) = [];
    end
end


function so = spp_sliceorder(ns, so)
sos = {'a', 'd', 'aio', 'aie', 'dio' , 'die', 'dsq'};
if ~ischar(so)
    so = sos{so};
end

% what slice order
switch (lower(so(:)'))

    % ascending
    case {'a'}
        so = 1:ns;

    % descending
    case {'d'}
        so = ns:-1:1;

    % ascending, interleaved, odd-numbered slices first
    case {'aio'}
        so = [1:2:ns, 2:2:ns];

    % ascending, interleaved, even-numbered slices first
    case {'aie'}
        so = [2:2:ns, 1:2:ns];

    % descending, interleaved, odd-numbered slices first
    case {'dio'}
        if mod(ns, 2) == 1
            so = [ns:-2:1, (ns-1):-2:1];
        else
            so = [(ns-1):-2:1, ns:-2:1];
        end

    % descending, interleaved, even-numbered slices first
    case {'die'}
        if mod(ns, 2) == 1
            so = [(ns-1):-2:1, ns:-2:1];
        else
            so = [ns:-2:1, (ns-1):-2:1];
        end

    % Philips interleaved
    case {'dsq'}
        sns = round(sqrt(ns));
        so = lsqueeze(reshape(1:(sns*(sns+2)), sns, sns + 2)');
        so(so > ns) = [];
        so = so(:)';
end


function [jc, jn, jf] = spp_jobcost(jobs)
jc = ones(numel(jobs), 1);
jn = repmat({'SPM job'}, size(jc));
jf = repmat({'functional images'}, size(jc));
for c = 1:numel(jc)
    sj = jobs{c};
    if ~isstruct(sj) || ...
        numel(sj) ~= 1
        continue;
    end
    if isfield(sj, 'spm')
        sj = sj.spm;
        if ~isstruct(sj) || ...
            numel(sj) ~= 1
            continue;
        end
    end
    fn = fieldnames(sj);
    if numel(fn) ~= 1
        continue;
    end
    fn = fn{1};
    if iscell(sj.(fn)) && ...
        numel(sj.(fn)) == 1
        sj.(fn) = sj.(fn){1};
    end
    if ~isstruct(sj.(fn)) || ...
        numel(sj.(fn)) ~= 1 || ...
        numel(fieldnames(sj.(fn))) < 1
        continue;
    end
    sfn = fieldnames(sj.(fn));
    sfn = lower(sfn{1});
    jn{c} = sprintf('SPM %s', fn);
    switch lower(fn)
        case {'spatial'}
            switch sfn
                case {'coreg'}
                    if numel(fieldnames(sj.(fn))) > 1
                        jc(c) = 2.5;
                        jn{c} = 'coregistration + segmentation';
                    else
                        jc(c) = 0.5;
                        jn{c} = 'coregistration';
                    end
                case {'normalise'}
                    jc(c) = 0.75;
                    jn{c} = 'normalization';
                case {'realign'}
                    jc(c) = 3;
                    jn{c} = 'realignment';
                case {'preproc'}
                    jc(c) = 1.5;
                    jn{c} = 'segmentation';
                case {'smooth'}
                    jc(c) = 0.5;
                    jn{c} = 'smoothing';
            end
        case {'temporal'}
            jc(c) = 2;
            jn{c} = 'slice-timing';
        case {'tools'}
            switch sfn
                case {'dartel'}
                    dtype = fieldnames(sj.(fn).dartel);
                    switch lower(dtype{1})
                        case {'initial'}
                            jc(c) = 0.5 * numel(sj.(fn).dartel.initial.matnames);
                            jn{c} = 'DARTEL import';
                        case {'warp'}
                            jc(c) = 0.5 * numel(sj.(fn).dartel.warp.images{1}) * ...
                                sum(cat(1, sj.(fn).dartel.warp.settings.param.its) .* ...
                                (2 .^ cat(1, sj.(fn).dartel.warp.settings.param.K)));
                            jn{c} = 'DARTEL create-template';
                        case {'warp1'}
                            jc(c) = 0.25 * numel(sj.(fn).dartel.warp1.images{1}) * ...
                                sum(cat(1, sj.(fn).dartel.warp1.settings.param.its) .* ...
                                (0.5 + log(1 + cat(1, sj.(fn).dartel.warp1.settings.param.K))));
                            jc(c) = 0.25 * round(4 * jc(c));
                            jn{c} = 'DARTEL apply-template';
                    end
            end
        case {'util'}
            switch sfn
                case {'imcalc'}
                    jc(c) = 0.5;
                    jn{c} = 'SPM imcalc';
            end
    end
end
