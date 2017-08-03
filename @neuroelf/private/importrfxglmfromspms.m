function glm = importrfxglmfromspms(spms, opts)
% importrfxglmfromspms  - import SPM beta maps to a BrainVoyager GLM file
%
% FORMAT:       glm = importrfxglmfromspms(spms, opts)
%
% Input fields:
%
%       spms        list of SPM.mat filenames to use for import
%       opts        options for the import
%        .bbox      bounding box to use (in BVS notation)
%        .cond      1xC structure with fields (default: auto detect)
%         .bvname   condition name for the BV file
%         .color    1x3 RGB color for condition
%         .spmname  SPM based name (or pattern)
%        .filename  output GLM filename (otherwise unsaved)
%        .imeth     interpolation 'cubic', 'lanczos3', {'linear'}, 'nearest'
%        .impvtc    pattern for VTCs (default: '', e.g. '#_run%02d.vtc')
%        .keepbf    keep basis-function (*bf(X)) part when detecting names
%        .mmaskfac  mean-masking factor, default 0.5
%                   results in masking out voxels where constant < 0.25
%                   of its mean overall value (per subject)
%        .pbar      either xprogress or xfigure:XProgress object
%        .res       resolution (default: floor of SPM resolution)
%        .subjids   subject identifiers (default: auto detect)
%        .trans     either {'none'}, 'psc'
%                   PSC-transform uses estimated mean (constant) as 100!
%        .vweight   variance-based weighting (false, uses SPM.xX.Bcov)
%
% Output fields:
%
%       glm         GLM object (saved if .filename is given)
%
% Note: if VTC import is desired, the VTC filename pattern (.impvtc)
%       must contain one or two occurrences of the hash mark (which will be
%       replaced by the subject ID) and one occurrence %d or %XXd for the
%       Sess number; other arguments (bbox, imeth, res) will be passed on

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:21 AM EST
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

% variable for UI stuff
global ne_ui;

% allow UI call
if nargin > 0 && ...
    ischar(spms) && ...
    strcmpi(spms(:)', 'ui')

    % open figure
    try
        hFig = xfigure([neuroelf_path('tfg') '/importrfxglmfromspms.tfg']);
    catch ne_eo;
        error( ...
            'neuroelf:xfigureError', ...
            'Error creating UI for GLM import: %s.', ...
            ne_eo.message ...
        );
    end

    % initialize global variable
    hTag = hFig.TagStruct;
    ne_ui.importrfxglmfromspms = struct( ...
        'ccol', zeros(0, 3), ...
        'glm',  [], ...
        'hFig', hFig, ...
        'hTag', hTag, ...
        'pbar', []);

    % prepare listboxes
    hTag.LB_importrfxglm_subids.Value = [];
    hTag.LB_importrfxglm_subids.String = {};
    hTag.LB_importrfxglm_files.Value = [];
    hTag.LB_importrfxglm_files.String = {};
    hTag.LB_importrfxglm_spmcond.Value = [];
    hTag.LB_importrfxglm_spmcond.String = {};
    hTag.LB_importrfxglm_glmcond.Value = [];
    hTag.LB_importrfxglm_glmcond.String = {};

    % set callbacks
    hTag.BT_importrfxglm_folder.Callback = @ig_btbrowse;
    hTag.BT_importrfxglm_search.Callback = @ig_btsearch;
    hTag.LB_importrfxglm_subids.Callback = @ig_subidselect;
    hTag.LB_importrfxglm_files.Callback = @ig_spmmatselect;
    hTag.BT_importrfxglm_fdown.Callback = @ig_btfiledown;
    hTag.BT_importrfxglm_fup.Callback = @ig_btfileup;
    hTag.BT_importrfxglm_subjprop.Callback = @ig_btsubjprop;
    hTag.BT_importrfxglm_fplus.Callback = @ig_btfileadd;
    hTag.BT_importrfxglm_fminus.Callback = @ig_btfileremove;
    hTag.LB_importrfxglm_spmcond.Callback = @ig_lbspmcselect;
    hTag.LB_importrfxglm_glmcond.Callback = @ig_lbglmcselect;
    hTag.BT_importrfxglm_cdown.Callback = @ig_btconddown;
    hTag.BT_importrfxglm_cup.Callback = @ig_btcondup;
    hTag.BT_importrfxglm_glmcprop.Callback = @ig_btcondprop;
    hTag.BT_importrfxglm_glmccol.Callback = @ig_btcondcol;
    hTag.BT_importrfxglm_cplus.Callback = @ig_btcondadd;
    hTag.BT_importrfxglm_cminus.Callback = @ig_btcondremove;
    hTag.CB_importrfxglm_autcond.Callback = @ig_toggleautocond;
    hTag.CB_importrfxglm_autores.Callback = @ig_toggleautores;
    hTag.BT_importrfxglm_saveas.Callback = @ig_btsaveas;
    hTag.CB_importrfxglm_impvtc.Callback = @ig_toggleimpvtc;
    hTag.BT_importrfxglm_cancel.Callback = @ig_closeui;
    hTag.BT_importrfxglm_import.Callback = @ig_btimport;

    % progress bar option given after all?
    if nargin > 1 && ...
        isstruct(opts) && ...
        numel(opts) == 1
        if isfield(opts, 'pbar')
            ne_ui.importrfxglmfromspms.pbar = opts.pbar;
        end
    end

    % make figure modal
    hFig.HandleVisibility = 'callback';
    hFig.Visible = 'on';
    hFig.WindowStyle = 'modal';

    % wait for figure to close
    uiwait(hFig.mlhandle);

    % assign GLM output
    glm = ne_ui.importrfxglmfromspms.glm;
    ne_ui.importrfxglmfromspms = [];

    % return
    return;
end

% argument check
if nargin < 1 || ...
   ~iscell(spms) || ...
    isempty(spms)
    error( ...
        'neuroelf:BadArgument', ...
        'At least three subjects needed for RFX GLM.' ...
    );
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'bbox') || ...
   ~isa(opts.bbox, 'double') || ...
   ~isequal(size(opts.bbox), [2, 3]) || ...
    any(isnan(opts.bbox(:)) | opts.bbox(:) < 0 | opts.bbox(:) > 255)
    opts.bbox = [];
else
    opts.bbox = round(opts.bbox);
end
if ~isfield(opts, 'cond') || ...
   ~isstruct(opts.cond) || ...
    isempty(opts.cond) || ...
   ~isfield(opts.cond, 'bvname') || ...
   ~isfield(opts.cond, 'color') || ...
   ~isfield(opts.cond, 'spmname')
    opts.cond = emptystruct({'bvname', 'color', 'spmname'});
else
    spmnames = {opts.cond.spmname};
    glmnames = {opts.cond.bvname};
    if any(cellfun('isempty', spmnames)) || ...
        numel(unique(lower(spmnames))) ~= numel(spmnames) || ...
        any(cellfun('isempty', glmnames)) || ...
        numel(unique(lower(glmnames))) ~= numel(glmnames)
        error( ...
            'neuroelf:BadArgument', ...
            'Bad condition names request.' ...
        );
    end
end
if ~isfield(opts, 'filename') || ...
   ~ischar(opts.filename) || ...
    isempty(opts.filename) || ...
    numel(opts.filename) > 255
    opts.filename = '';
else
    opts.filename = opts.filename(:)';
end
if ~isfield(opts, 'imeth') || ...
   ~ischar(opts.imeth) || ...
   ~any(strcmpi(opts.imeth(:)', {'cubic', 'lanczos3', 'linear', 'nearest'}))
    opts.imeth = 'linear';
else
    opts.imeth = lower(opts.imeth(:)');
end
if ~isfield(opts, 'impvtc') || ...
   ~ischar(opts.impvtc) || ...
   ~any([1, 2] == sum(opts.impvtc(:) == '#')) || ...
    isempty(regexpi(opts.impvtc(:)', '\%\d*d.*\.vtc$'))
    opts.impvtc = '';
else
    [ivp, ivf] = fileparts(opts.impvtc(:)');
    if isempty(ivp)
        if ~isempty(opts.filename)
            ivp = fileparts(opts.filename);
            if isempty(ivp)
                ivp = strrep(pwd, '\', '/');
            end
        else
            ivp = strrep(pwd, '\', '/');
        end
        ivp = strrep(ivp, '\', '/');
    else
        ivp = strrep(ivp, '\', '/');
        if ivp(1) ~= '/'
            if ~isempty(opts.filename)
                ivpp = strrep(fileparts(opts.filename), '\', '/');
                if isempty(ivpp)
                    ivpp = strrep(pwd, '\', '/');
                end
                ivp = strrep([ivpp '/' ivp], '//', '/');
            else
                ivp = strrep([strrep(pwd, '\', '/') '/' ivp], '//', '/');
            end
        end
    end
    opts.impvtc = strrep([ivp '/' ivf], '//', '/');
end
if ~isfield(opts, 'keepbf') || ...
   ~islogical(opts.keepbf) || ...
    numel(opts.keepbf) ~= 1
    opts.keepbf = false;
end
if ~isfield(opts, 'mmaskfac') || ...
   ~isa(opts.mmaskfac, 'double') || ...
    numel(opts.mmaskfac) ~= 1 || ...
    isinf(opts.mmaskfac) || ...
    isnan(opts.mmaskfac) || ...
    opts.mmaskfac < 0 || ...
    opts.mmaskfac > 2
    opts.mmaskfac = 0.25;
end
if ~isfield(opts, 'pbar') || ...
    numel(opts.pbar) ~= 1 || ...
   ~any(strcmpi(class(opts.pbar), {'xfigure', 'xprogress'}))
    opts.pbar = [];
end
if ~isfield(opts, 'res') || ...
   ~isa(opts.res, 'double') || ...
    numel(opts.res) ~= 1 || ...
   ~any((1:4) == opts.res)
    opts.res = [];
end
if ~isfield(opts, 'subjids') || ...
   ~iscell(opts.subjids) || ...
    numel(opts.subjids) ~= numel(spms)
    opts.subjids = {};
end
if numel(spms) == 1 && ...
    isempty(opts.subjids)
    error( ...
        'neuroelf:BadArgument', ...
        'For single subject, a subject ID must be given.' ...
    );
end
for sc = 1:numel(opts.subjids)
    if ~ischar(opts.subjids{sc}) || ...
        isempty(opts.subjids{sc}) || ...
        numel(opts.subjids{sc}) > 63
        opts.subjids = {};
        break;
    else
        opts.subjids{sc} = opts.subjids{sc}(:)';
    end
end
if ~isempty(opts.subjids)
    if numel(opts.subjids) ~= numel(spms) || ...
        numel(unique(lower(opts.subjids))) ~= numel(spms)
        opts.subjids = {};
    end
end
if ~isfield(opts, 'trans') || ...
   ~ischar(opts.trans) || ...
   ~any(strcmpi(opts.trans, {'none', 'psc'}))
    opts.trans = 'none';
else
    opts.trans = lower(opts.trans(:)');
end
if ~isfield(opts, 'vweight') || ...
   ~islogical(opts.vweight) || ...
    numel(opts.vweight) ~= 1
    opts.vweight = false;
end

% try load SPMs and store essential information
spmc = emptystruct({ ...
    'filename', 'filepath', 'Vbeta', 'nSess', 'iXXdiag', 'swd', 'TR', ...
    'xXsize', 'xXname', 'xXrgbf', 'xXsess', 'xXuname', 'xY', 'xYnvol'}, ...
    [1, numel(spms)]);

for sc = 1:numel(spms)
    try
        if ~ischar(spms{sc}) || ...
            isempty(spms{sc})
            spms{sc} = sprintf('No/empty string at position %d.', sc);
            error('BADSTRING');
        end
        spms{sc} = spms{sc}(:)';
        s = load(spms{sc});
        if ~isfield(s, 'SPM')
            s = struct('SPM', s);
        end
        if ~isfield(s.SPM, 'Sess') || ...
           ~isstruct(s.SPM.Sess) || ...
           ~isfield(s.SPM, 'Vbeta') || ...
           ~isfield(s.SPM, 'xX') || ...
           ~isstruct(s.SPM.xX)
            error('INVALID_SPM');
        end
        if ~isstruct(s.SPM.Vbeta)
            if iscell(s.SPM.Vbeta)
                for bc = 1:numel(s.SPM.Vbeta)
                    if ~isstruct(s.SPM.Vbeta{bc})
                        if ischar(s.SPM.Vbeta{bc})
                            s.SPM.Vbeta{bc} = ...
                                struct('fname', s.SPM.Vbeta{bc}(:)');
                        else
                            error('INVALID_SPM');
                        end
                    end
                end
                s.SPM.Vbeta = cat(1, s.SPM.Vbeta{:});
            else
                error('INVALID_SPM');
            end
        end
        if ~isfield(s.SPM.xX, 'name')
            if isfield(s.SPM.xX, 'Xnames')
                s.SPM.xX.name = s.SPM.xX.Xnames;
            else
                error('NO_xX.name');
            end
        end
        spmc(sc).filename = spms{sc};
        spmc(sc).filepath = [fileparts(spms{sc}), filesep];
        spmc(sc).iXXdiag = 1 ./ sqrt(diag(s.SPM.xX.Bcov));
        try
            spmc(sc).swd = s.SPM.swd;
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            spmc(sc).swd = strrep(fileparts(spms{sc}), '\', '/');
        end
        spmc(sc).TR = round(1000 * s.SPM.xBF.T * s.SPM.xBF.dt);
        spmc(sc).Vbeta = strrep({s.SPM.Vbeta(:).fname}, 'beta', ...
            [spmc(sc).filepath, filesep, 'beta']);
        spmc(sc).xXsize = size(s.SPM.xX.X, 1);
        spmc(sc).xXname = s.SPM.xX.name;
        if ischar(s.SPM.xY.P)
            spmc(sc).xY = cellstr(s.SPM.xY.P);
        else
            spmc(sc).xY = s.SPM.xY.P;
        end
        spmc(sc).xYnvol = s.SPM.nscan;
        if numel(spmc(sc).Vbeta) ~= numel(spmc(sc).xXname)
            error('BAD_NROFxX.name');
        end
        spmc(sc).Vbeta = regexprep(spmc(sc).Vbeta, '\.img', '.hdr', 'preservecase');
        spmc(sc).xXrgbf = regexprep(spmc(sc).xXname, '^.*\*bf\((\d+)\)', '$1');
        spmc(sc).xXsess = regexprep(spmc(sc).xXname, '^.*Sn\((\d+)\).*$', '$1');
        if opts.keepbf
            spmc(sc).xXname = regexprep(spmc(sc).xXname, '^.*Sn\(\d+\)\s*', '');
        else
            spmc(sc).xXname = regexprep(regexprep(spmc(sc).xXname, ...
                '^.*Sn\(\d+\)\s*', ''), '\s*\*bf\(\d+\).*$', '');
        end
        [spmc(sc).xXuname, foi, spmc(sc).xXuidx] = unique(spmc(sc).xXname);
        spmc(sc).nSess = str2double(spmc(sc).xXsess{end});
    catch ne_eo;
        error( ...
            'neuroelf:SPMError', ...
            'Invalid SPM.mat file: %s (%s).', ...
            spms{sc}, ne_eo.message ...
        );
    end
end

% subject ID's ?
if isempty(opts.subjids)
    subjids = char(spms);
    subjids = cellstr(subjids(:, findfirst(any(diff(subjids))):end));
    subjidl = 0;
    for sc = 1:numel(subjids)
        subjids{sc} = fileparts(subjids{sc});
        subjidl = max(subjidl, numel(subjids{sc}));
    end
    if numel(opts.subjids) > 1
        subjidl = reshape(sprintf(sprintf('%%%ds', subjidl), subjids{:}), ...
            subjidl, numel(subjids))';
        subjidl = subjidl(:, 1:findfirst(any(diff(subjidl)), -1));
        for sc = 1:numel(subjids)
            subjids{sc} = strrep(subjidl(sc, :), ' ', '');
        end
    end
    opts.subjids = subjids;
end

% time points
timepoints = {spmc(:).xXsize};
timepoints = cat(1, timepoints{:});
nroftotaltimepoints = sum(timepoints);

% check spatial layout/dimension
Vbetaobj = cell(1, numel(spmc));
for sc = 1:numel(spmc)
    try
        Vbetaobj{sc} = xff(spmc(sc).Vbeta{1});
    catch ne_eo;
        clearxffobjects(Vbetaobj);
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid Analyze/NIftI beta_0001 file for %s: %s', ...
            spmc(sc).filename, ne_eo.message ...
        );
    end
end
VResDims = zeros(numel(spmc), 19);
for sc = 1:numel(spmc)
    VResDims(sc, :) = [lsqueeze(Vbetaobj{sc}.CoordinateFrame.Trf(:))', ...
        Vbetaobj{sc}.CoordinateFrame.Dimensions(1:3)];
end
if any(any(diff(VResDims, 1, 1)))
    warning( ...
        'neuroelf:SPMWarning', ...
        'Incompatible map sizes. Output must be taken carefully!!' ...
    );
end

% create list of regressors from first subject (in order of appearance)
uname = spmc(1).xXuname;
sname = uname;
uidx = spmc(1).xXuidx;
for cc = 1:numel(spmc(1).xXuname)
    sname{cc} = uname{uidx(1)};
    uidx(uidx == uidx(1)) = [];
end

% and ensure regressors are available in all other subjects!
for sc = 2:numel(spmc)
    for cc = numel(sname):-1:1
        if ~any(strcmpi(spmc(sc).xXuname, sname{cc}))
            sname(cc) = [];
        end
    end
end

% no common regressors?
if isempty(sname)
    clearxffobjects(Vbetaobj);
    error( ...
        'neuroelf:BadArgument', ...
        'Condition names don''t match between subject.' ...
    );
end

% no conditions named
if isempty(opts.cond)

    % put into cond list
    opts.cond(numel(sname)).bvname = '';
    for cc = 1:numel(sname)
        opts.cond(cc).bvname = sname{cc};
        opts.cond(cc).color = [floor(255.999 * rand(1, 3)); zeros(3, 3)];
        opts.cond(cc).spmname = sname{cc};
    end
    opts.cond(end).color = 255 .* [ones(1, 3); zeros(3, 3)];

% for named conditions
else

    % check they exist in all subjects
    for cc = 1:numel(opts.cond)
        if ~any(strcmpi(opts.cond(cc).spmname, sname))
            clearxffobjects(Vbetaobj);
            error( ...
                'neuroelf:BadArgument', ...
                'Requested condition (%s) not found in SPM.mat files.', ...
                opts.cond(cc).spmname ...
            );
        end
    end
end

% get spm condition names
spmcn = {opts.cond(:).spmname};
consti = find(~cellfun('isempty', regexpi(spmcn, 'constant')));
if ~isempty(consti)
    opts.cond(consti) = [];
end
spmcn = {opts.cond(:).spmname};
glmcn = {opts.cond(:).bvname};
glmcc = cat(1, opts.cond(:).color);
opts.cond(end + 1).bvname = 'Constant';
opts.cond(end).color = [255, 255, 255];
opts.cond(end).spmname = 'constant';

% make sure SPM names are patterns alright
for cc = 1:numel(opts.cond)
    if opts.cond(cc).spmname(1) ~= '^'
        opts.cond(cc).spmname = ['^' opts.cond(cc).spmname];
    end
    if opts.cond(cc).spmname(end) ~= '$'
        opts.cond(cc).spmname(end+1) = '$';
    end
end

% check resolution
if isempty(opts.res)
    opts.res = floor(0.05 + mean(sqrt(sum( ...
        Vbetaobj{1}.CoordinateFrame.Trf(1:3, 1:3) .^ 2))));
end

% temporary clean-up
clearxffobjects(Vbetaobj);

% create GLM structure
glm = xff('new:glm');

% make some initial settings
glm.ProjectType = 1;
glm.ProjectTypeRFX = 1;
glm.NrOfSubjects = numel(spmc);
glm.NrOfSubjectPredictors = numel(opts.cond);
glm.NrOfTimePoints = nroftotaltimepoints;
glm.NrOfPredictors = numel(spmc) * glm.NrOfSubjectPredictors;
glm.NrOfConfounds = numel(spmc);
glm.NrOfStudies = numel(spmc);
glm.NrOfStudiesWithConfounds = numel(spmc);
glm.NrOfConfoundsPerStudy = ones(1, numel(spmc));
glm.SeparatePredictors = 2;
glm.Resolution = opts.res;
glm.SerialCorrelation = 0;
glm.Study(numel(spmc)).NrOfTimePoints = 0;
glm.Predictor(glm.NrOfPredictors).Name1 = '';
prc = 1;
studyxtcs = {};
studyprts = {};
studysdms = {};
studysdmc = [];
studysdmn = [];

% initialize progress bar
if isempty(opts.pbar)
    try
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 264, 640, 36]);
        xprogress(pbar, 'settitle', ...
            sprintf('Importing %d subjects'' SPM.mat to RFX-GLM...', numel(spmc)));
        xprogress(pbar, 0, 'Importing subject ...', 'visible', 0, 1);
        pbarn = '';
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        pbar = [];
    end
else
    if isxfigure(opts.pbar, true)
        opts.pbar.Visible = 'on';
    end
    pbar = opts.pbar;
    pbar.Progress(0, ...
        sprintf('importrfxglmfromspms: Importing subject 1/%d...', numel(spmc)));
    pbarn = 'importrfxglmfromspms: ';
end

% iterate over subjects
for sc = 1:numel(spmc)

    % update progress bar
    if ~isempty(pbar)
        pbar.Progress((sc - 1) / numel(spmc), ...
            sprintf('%sImporting subject ''%s'' (%d/%d)...', pbarn, ...
            opts.subjids{sc}, sc, numel(spmc)));
    end

    % also import data
    if ~isempty(opts.impvtc)

        % generate name for this subject
        impvtc = strrep(opts.impvtc, '#', opts.subjids{sc});

        % models (with error handling)
        try

            % start with protocols
            prts = {};
            sdms = {};
            s = load(spmc(sc).filename);
            [prt, prts] = spmmat2prt(s.SPM, [tempname '.prt']);

            % for each protocol
            for pc = 1:numel(prts)

                % delete file
                delete(prts{pc}.FilenameOnDisk);

                % match to requested conditions
                cmatch = multimatch(lower(spmcn), lower(prts{pc}.ConditionNames));
                if any(cmatch < 1)
                    error('MISMATCHED_CONDITION');
                end

                % update names and colors
                for cc = 1:numel(cmatch)
                    prts{pc}.Cond(cmatch(cc)).ConditionName = glmcn(cc);
                    prts{pc}.Cond(cmatch(cc)).Color = glmcc(cc, :);
                end

                % remove other conditions
                prts{pc}.Cond(setdiff(1:numel(prts{pc}.Cond), cmatch(:)')) = [];
                prts{pc}.NrOfConditions = numel(prts{pc}.Cond);

                % and store as desired file
                tfile = sprintf([impvtc '.prt'], pc);
                if exist(fileparts(tfile), 'dir') ~= 7
                    mkadir(fileparts(tfile));
                    if exist(fileparts(tfile), 'dir') ~= 7
                        error('DIR_NOT_CREATED');
                    end
                end
                prts{pc}.Experiment = ...
                    sprintf('%s session %d', opts.subjids{sc}, pc);
                prts{pc}.SaveAs(tfile);

                % then keep filename, but clear object
                pfile = prts{pc}.FilenameOnDisk;
                prts{pc}.ClearObject;
                prts{pc} = pfile;
            end

            % then do SDMs
            [sdm, sdms] = spmmat2sdm(s.SPM, [tempname '.sdm'], spmcn);
            sdmc = ones(numel(sdms), 1);
            sdmn = spmc(sc).xYnvol(:);

            % for each SDM
            for pc = 1:numel(sdms)

                % delete original file
                delete(sdms{pc}.FilenameOnDisk);

                % match to requested conditions
                cmatch = multimatch(lower(spmcn), lower(sdms{pc}.PredictorNames));
                if any(cmatch(:)' ~= (1:numel(cmatch)))
                    error('MISMATCHED_CONDITION');
                end

                % re-color
                sdms{pc}.PredictorColors(1:numel(spmcn), :) = glmcc;

                % and store as desired file
                tfile = sprintf([impvtc '.sdm'], pc);
                if exist(fileparts(tfile), 'dir') ~= 7
                    mkadir(fileparts(tfile));
                    if exist(fileparts(tfile), 'dir') ~= 7
                        error('DIR_NOT_CREATED');
                    end
                end
                sdms{pc}.SaveAs(tfile);
                sdmc(pc) = size(sdms{pc}.SDMMatrix, 2) + 1 - sdms{pc}.FirstConfoundPredictor;

                % then keep filename, but clear object
                pfile = sdms{pc}.FilenameOnDisk;
                sdms{pc}.ClearObject;
                sdms{pc} = pfile;
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            clearxffobjects(prts);
            clearxffobjects(sdms);
            prts = {};
            sdms = {};
            sdmc = [];
            sdmn(:) = 0;
            warning('neuroelf:warning', ne_eo.message);
        end

        % import VTCs
        volc = 1;
        vtcs = cell(spmc(sc).nSess, 1);
        kprt = false(1, numel(vtcs));
        for pc = 1:spmc(sc).nSess

            % filenames to import
            inames = regexprep( ...
                spmc(sc).xY(volc:(volc + spmc(sc).xYnvol(pc) - 1)), ...
                '\.img', '.hdr', 'preservecase');
            volc = volc + spmc(sc).xYnvol(pc);

            % files don't exist
            iname1 = regexprep(inames{1}, '\,\d+$', '');
            if exist(iname1, 'file') ~= 2

                % location of SPM mat file
                spmndir = strrep(fileparts(spmc(sc).filename), '\', '/');
                spmodir = spmc(sc).swd;

                % re-location not possible
                if strcmp(spmodir, spmndir)
                    continue;
                end
                ncomp = min(numel(spmndir), numel(spmodir));

                % to re-locate with respect to swd in SPM, find consistency
                lmatch = findfirst( ...
                    spmodir(end+1-ncomp:end) ~= spmndir(end+1-ncomp:end), -1);

                % no inconsistency or no consitency in same length part?
                if isempty(lmatch) || ...
                    lmatch == ncomp

                    % can't do it...
                    continue;
                end

                % parts that match
                nsame = ncomp - lmatch;

                % try to replace this part
                inames = strrep(inames, ...
                    spmodir(1:end-nsame), spmndir(1:end-nsame));
                tnames = unique(regexprep(inames, '\,\d+$', ''));

                % check each file!
                allfexist = true;
                for pcf = 1:numel(tnames)
                    if exist(tnames{pcf}, 'file') ~= 2
                        allfexist = false;
                        break;
                    end
                end

                % if some are missing
                if ~allfexist
                    continue;
                end
            end

            % with error handling
            try

                % attempt import
                vtcs{pc} = importvtcfromanalyze( ...
                    inames, opts.bbox, opts.res, opts.imeth);

                % make settings
                if ~isempty(prts)
                    vtcs{pc}.NrOfLinkedPRTs = 1;
                    vtcs{pc}.NameOfLinkedPRT = prts{pc};
                    vtcs{pc}.NrOfCurrentPRT = 1;
                end
                vtcs{pc}.TR = spmc(pc).TR;

                % save as...
                tfile = sprintf([impvtc '.vtc'], pc);
                if exist(fileparts(tfile), 'dir') ~= 7
                    mkadir(fileparts(tfile));
                    if exist(fileparts(tfile), 'dir') ~= 7
                        error('DIR_NOT_CREATED');
                    end
                end
                vtcs{pc}.SaveAs(tfile);
                if sdmn(pc) == 0
                    sdmn(pc) = vtcs{pc}.NrOfVolumes;
                end

                % then reload with transio (for filename)
                vtcs{pc}.ClearObject;
                vtcs{pc} = sprintf([impvtc '.vtc'], pc);

            % handle errors
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                clearxffobjects(vtcs(pc));
                vtcs{pc} = [];
                continue;
            end

            % keep files flag
            kprt(pc) = true;
        end

        % keep values
        vtcs = vtcs(kprt);
        if ~isempty(prts)
            prts = prts(kprt);
        else
            prts = repmat({''}, size(vtcs));
        end
        if ~isempty(sdms)
            sdms = sdms(kprt);
            sdmc = sdmc(kprt);
        else
            sdms = repmat({'interactive'}, size(vtcs));
            sdmc = ones(size(vtcs));
        end
        sdmn = sdmn(kprt);
        studyxtcs = [studyxtcs(:); vtcs(:)];
        studyprts = [studyprts(:); prts(:)];
        studysdms = [studysdms(:); sdms(:)];
        studysdmc = [studysdmc(:); sdmc(:)];
        studysdmn = [studysdmn(:); sdmn(:)];
    end

    % temporarily import beta maps into VMP
    bvmp = importvmpfromspms(spmc(sc).Vbeta, 'b', ...
        opts.bbox, opts.res, opts.imeth);

    % GLM settings
    if sc == 1
        glm.XStart = bvmp.XStart;
        glm.XEnd = bvmp.XEnd;
        glm.YStart = bvmp.YStart;
        glm.YEnd = bvmp.YEnd;
        glm.ZStart = bvmp.ZStart;
        glm.ZEnd = bvmp.ZEnd;
        % account for VMP with 1mm resolution (STILL NEEDS BVQX CHECKING!!)
        if glm.Resolution == 1
            glm.XEnd = glm.XEnd + 1;
            glm.YEnd = glm.YEnd + 1;
            glm.ZEnd = glm.ZEnd + 1;
        end
        glm.GLMData.RFXGlobalMap = single(zeros(size(bvmp.Map(1).VMPData)));
        glm.GLMData.Subject(1).BetaMaps = ...
            single(zeros([size(bvmp.Map(1).VMPData), glm.NrOfSubjectPredictors]));
        glm.GLMData.Subject(2:numel(spmc)) = glm.GLMData.Subject(1);
        omask = zeros(size(bvmp.Map(1).VMPData));
        odsz = size(omask);
    end

    % update progress bar for import
    if ~isempty(pbar)
        pbar.Progress((sc - 0.4) / numel(spmc), sprintf( ...
            '%sCombining runs and transforming data (%d/%d)...', ...
            pbarn, sc, numel(spmc)));
    end

    % set study & predictors
    glm.Study(sc).NrOfTimePoints = timepoints(sc);
    glm.Study(sc).NameOfAnalyzedFile = [spmc(sc).filepath, ...
        mstrrep(opts.subjids{sc}, {'/','_'}, {'', ''}), '_SPM.vtc'];
    glm.Study(sc).NameOfSDMFile = '<SPM.xX.X>';
    for pc = 1:(numel(opts.cond) - 1)
        glm.Predictor(prc).Name1 = sprintf('Predictor: %d', prc);
        glm.Predictor(prc).Name2 = ...
            sprintf('Subject %s: %s', opts.subjids{sc}, opts.cond(pc).bvname);
        glm.Predictor(prc).RGB = [opts.cond(pc).color; zeros(3, 3)];
        prc = prc + 1;
    end

    % put VMPData together
    vmpd = bvmp.Map;
    vmpd = {vmpd(:).VMPData};
    vmpd = cat(4, vmpd{:});
    bvmp.ClearObject;

    % update mask Inf/NaN
    mask = (~any(isinf(vmpd), 4) & ~any(isnan(vmpd), 4));
    vmpd(isinf(vmpd) | isnan(vmpd)) = 0;

    % update mask with mean-masking factor
    mnm = (1 ./ spmc(sc).nSess) .* sum(vmpd(:, :, :, end+1-spmc(sc).nSess:end), 4);
    mmm = minmaxmean(mnm);
    mnm = mnm < (opts.mmaskfac * mmm(3));
    vmpd(mnm(:, :, :, ones(1, size(vmpd, 4)))) = 0;

    % update mask general and across subjects
    mask = mask & (~all(vmpd == 0, 4));
    omask = omask + double(mask);

    % PSC-transform ?
    if strcmpi(opts.trans, 'psc')
        for pc = 1:size(vmpd, 4)
            vmpd(:, :, :, pc) = 100 .* (vmpd(:, :, :, pc) ./ ...
                vmpd(:, :, :, end - spmc(sc).nSess + str2double(spmc(sc).xXsess{pc})));
        end
        vmpd(isinf(vmpd) | isnan(vmpd)) = 0;
    end

    % combine beta maps (plain average for now)
    bmaps = glm.GLMData.Subject(sc).BetaMaps;
    for pc = 1:numel(opts.cond)

        % find beta maps of condition
        if opts.keepbf
            uidx = find(strcmpi(spmc(sc).xXuname(:), opts.cond(pc).spmname(2:end-1)));
        else
            uidx = find(~cellfun('isempty', regexpi(spmc(sc).xXuname, opts.cond(pc).spmname)));
        end
        bmci = [];
        for uic = 1:numel(uidx)
            bmca = find(spmc(sc).xXuidx == uidx(uic));
            bmci = [bmci, bmca(:)'];
        end

        % combine maps as average
        if opts.vweight
            ixxdiag = repmat(shiftdim(lsqueeze(spmc(sc).iXXdiag(bmci)), -3), odsz);
            ixxdiag(vmpd(:, :, :, bmci) == 0) = 0;
            bmaps(:, :, :, pc) = ...
                sum(vmpd(:, :, :, bmci), 4) ./ sum(ixxdiag, 4);
        else
            bmaps(:, :, :, pc) = ...
                (1 / numel(bmci)) .* sum(vmpd(:, :, :, bmci), 4);
        end
    end

    % eliminate invalid content
    bmaps(isinf(bmaps) | isnan(bmaps)) = 0;
    glm.GLMData.Subject(sc).BetaMaps = bmaps;
end

% apply final masking
glm.GLMData.RFXGlobalMap = single(omask >= (0.75 * numel(spmc)));
glm.NrOfVoxelsForBonfCorrection = ...
    double(sum(glm.GLMData.RFXGlobalMap(:)));
omask = repmat(omask < (0.75 * numel(spmc)), [1, 1, 1, numel(opts.cond)]);
for sc = 1:numel(spmc)
    glm.GLMData.Subject(sc).BetaMaps(omask) = 0;
end

% add last predictors
for sc = 1:numel(spmc)
    glm.Predictor(prc).Name1 = sprintf('Predictor: %d', prc);
    glm.Predictor(prc).Name2 = ...
        sprintf('Subject %s: constant', opts.subjids{sc});
    glm.Predictor(prc).RGB = [opts.cond(end).color; zeros(3, 3)];
    prc = prc + 1;
end

% replace study information
if ~isempty(studyxtcs)
    glm.NrOfStudies = numel(studyxtcs);
    glm.NrOfConfounds = sum(studysdmc);
    glm.NrOfStudiesWithConfounds = numel(studyxtcs);
    glm.NrOfConfoundsPerStudy = studysdmc(:)';
    glm.Study = glm.Study(1, ones(1, numel(studyxtcs)));
    for sc = 1:numel(glm.Study)
        glm.Study(sc).NrOfTimePoints = studysdmn(sc);
        glm.Study(sc).NameOfAnalyzedFile = studyxtcs{sc};
        glm.Study(sc).NameOfSDMFile = studysdms{sc};
    end
end

% close progress bar
if ~isempty(pbar) && ...
    isempty(opts.pbar)
    closebar(pbar);
end

% save ?
if ~isempty(opts.filename)
    try
        glm.SaveAs(opts.filename(:)');
    catch ne_eo;
        warning( ...
            'neuroelf:SaveError', ...
            'Error saving GLM file (%s). Please do so manually.', ...
            ne_eo.message ...
        );
    end
end



% internal functions
function clist = ig_getclist(spmfile, keepbf)
% try loading
try
    tSPM = load(spmfile);
    tSPM = tSPM.SPM;
    if ~isfield(tSPM, 'Sess') || ...
       ~isfield(tSPM, 'xX')
        error('BAD_SPM_FILE');
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    uiwait(warndlg('Invalid SPM.mat file selected.', ...
        'Import RFX-GLM from SPM: Warning', 'modal'));
    return;
end

% look for conditions?
if keepbf
    clist = uunion(regexprep(tSPM.xX.name, '^.*Sn\(\d+\)\s*', ''), {});
else
    clist = uunion(regexprep(regexprep(tSPM.xX.name, ...
        '^.*Sn\(\d+\)\s*', ''), '\s*\*bf\(\d+\).*$', ''), {});
end



% UI functions
function ig_closeui(varargin)
global ne_ui;
hFig = ne_ui.importrfxglmfromspms.hFig;
hFig.CloseRequestFcn = '';
hFig.Delete;

function ig_btbrowse(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;

% disable import button
hTag.BT_importrfxglm_import.Enable = 'off';

% request one SPM.mat file
[spmfile, spmpath] = uigetfile( ...
    {'SPM.mat', 'SPM design matrices (SPM.mat)'}, ...
    'Please select one of the SPM.mat files to import...', ...
    'MultiSelect', 'off');
if isequal(spmfile, 0) || ...
    isequal(spmpath, 0)
    return;
end
spmfile = strrep([strrep(spmpath, '\', '/') '/' spmfile], '//', '/');

% get condition list
hTag.CB_importrfxglm_keepbf.Enable = 'off';
if hTag.CB_importrfxglm_autcond.Value > 0
    try
        clist = ig_getclist(spmfile, hTag.CB_importrfxglm_keepbf.Value > 0);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        hTag.LB_importrfxglm_spmcond.Value = [];
        hTag.LB_importrfxglm_spmcond.String = {};
        hTag.LB_importrfxglm_glmcond.Value = [];
        hTag.LB_importrfxglm_glmcond.String = {};
        ne_ui.importrfxglmfromspms.ccol = zeros(0, 3);
        uiwait(warndlg('Invalid SPM.mat file selected.', ...
            'Import RFX-GLM from SPM: Warning', 'modal'));
        return;
    end
    hTag.LB_importrfxglm_spmcond.Value = [];
    hTag.LB_importrfxglm_spmcond.String = clist;
    hTag.LB_importrfxglm_glmcond.Value = [];
    hTag.LB_importrfxglm_glmcond.String = clist;
    ne_ui.importrfxglmfromspms.ccol = ...
        [floor(255.999 .* rand(numel(clist) - 1, 3)); 255 .* ones(1, 3)];
end

% put the filename into the edit field
hTag.ED_importrfxglm_folder.String = spmfile;
hTag.ED_importrfxglm_folder.Enable = 'on';
hTag.BT_importrfxglm_search.Enable = 'on';

% and also set filename in list (for access of condition list)
[spmsuppath, spmpath] = fileparts(fileparts(spmfile));
hTag.LB_importrfxglm_subids.Value = [];
hTag.LB_importrfxglm_subids.String = {spmpath};
hTag.LB_importrfxglm_files.Value = [];
hTag.LB_importrfxglm_files.String = {spmfile};
hTag.BT_importrfxglm_import.Enable = 'off';


function ig_btsearch(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;
keepbf = (hTag.CB_importrfxglm_keepbf.Value > 0);

% disable import button
hTag.BT_importrfxglm_import.Enable = 'off';

% perform findfiles
pattf = fileparts(hTag.ED_importrfxglm_folder.String);
try
    spms = findfiles(pattf, 'SPM.mat', 'depth=1');
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    uiwait(warndlg('Invalid search pattern given.', ...
        'Import RFX-GLM from SPM: Warning', 'modal'));
    return;
end

% get condition list
if hTag.CB_importrfxglm_autcond.Value > 0
    try
        clist = ig_getclist(spms{1}, keepbf);
    catch ne_eo;
        hTag.BT_importrfxglm_import.Enable = 'off';
        hTag.LB_importrfxglm_spmcond.Value = [];
        hTag.LB_importrfxglm_spmcond.String = {};
        hTag.LB_importrfxglm_glmcond.Value = [];
        hTag.LB_importrfxglm_glmcond.String = {};
        ne_ui.importrfxglmfromspms.ccol = zeros(0, 3);
        neuroelf_lasterr(ne_eo);
        uiwait(warndlg('Invalid SPM.mat file selected.', ...
            'Import RFX-GLM from SPM: Warning', 'modal'));
        return;
    end
    hTag.LB_importrfxglm_spmcond.Value = [];
    hTag.LB_importrfxglm_spmcond.String = clist;
    hTag.LB_importrfxglm_glmcond.Value = [];
    hTag.LB_importrfxglm_glmcond.String = clist;
    ne_ui.importrfxglmfromspms.ccol = ...
        [floor(255.999 .* rand(numel(clist) - 1, 3)); 255 .* ones(1, 3)];
end

% detect subject IDs
spms = strrep(spms, '\', '/');
spmp = spms;
spmpn = zeros(numel(spmp), 1);
for sc = 1:numel(spmp)
    spmp{sc} = splittocellc(fileparts(spmp{sc}), '/');
    spmpn(sc) = numel(spmp{sc});
end
spmpp = cell(numel(spmp), max(spmpn));
for sc = 1:numel(spmp)
    spmpp(sc, 1:numel(spmp{sc})) = spmp{sc};
end
for pc = 1:min(spmpn)
    if all(strcmp(spmpp{1, 1}, spmpp(:, 1)))
        spmpp(:, 1) = [];
        spmpn = spmpn - 1;
    else
        break;
    end
end
if ~all(spmpn == spmpn(1))
    maxspmpn = max(spmpn);
    for sc = 1:numel(spms)
        if spmpn(sc) < maxspmpn
            spmpp(sc, :) = [repmat({''}, 1, maxspmpn - spmpn(sc)), ...
                spmpp(sc, 1:spmpn(sc))];
        end
    end
    for pc = 1:min(spmpn)
        if all(strcmp(spmpp{1, end}, spmpp(:, end)))
            spmpp(:, end) = [];
            spmpn = spmpn - 1;
        else
            break;
        end
    end
    maxspmpn = max(spmpn);
    for sc = 1:numel(spms)
        if spmpn(sc) < maxspmpn
            spmpp(sc, :) = [spmpp(sc, (maxspmpn+1-spmpn(sc)):end), ...
                repmat({''}, 1, maxspmpn - spmpn(sc))];
        end
    end
end
for sc = 1:numel(spms)
    spmp{sc} = gluetostringc(spmpp(sc, 1:spmpn(sc)), '_');
end

% put into list of SPMs
hTag.LB_importrfxglm_subids.Value = [];
hTag.LB_importrfxglm_subids.String = spmp;
hTag.LB_importrfxglm_files.Value = [];
hTag.LB_importrfxglm_files.String = spms;
if numel(spms) > 2
    hTag.BT_importrfxglm_import.Enable = 'on';
end

function ig_subidselect(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;
try
    hTag.LB_importrfxglm_files.Value = hTag.LB_importrfxglm_subids.Value;
    hTag.LB_importrfxglm_files.ListBoxTop = hTag.LB_importrfxglm_subids.ListBoxTop;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

function ig_spmmatselect(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;
try
    hTag.LB_importrfxglm_subids.Value = hTag.LB_importrfxglm_files.Value;
    hTag.LB_importrfxglm_subids.ListBoxTop = hTag.LB_importrfxglm_files.ListBoxTop;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

function ig_btfileup(varargin)
global ne_ui;
ne_ui.importrfxglmfromspms.hTag.LB_importrfxglm_subids.MMove(-1);
ne_ui.importrfxglmfromspms.hTag.LB_importrfxglm_files.MMove(-1);

function ig_btfiledown(varargin)
global ne_ui;
ne_ui.importrfxglmfromspms.hTag.LB_importrfxglm_subids.MMove(1);
ne_ui.importrfxglmfromspms.hTag.LB_importrfxglm_files.MMove(1);

function ig_btsubjprop(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;

% none selected
idsel = hTag.LB_importrfxglm_subids.Value;
if isempty(idsel)
    return;
end

% get string
subjids = hTag.LB_importrfxglm_subids.String;
if ~iscell(subjids)
    subjids = cellstr(subjids);
end

% single ID
if numel(idsel) == 1

    % ask for replacement ID
    newid = inputdlg({'Subject ID for selected subject:'}, ...
        'Import GLM from SPM.mat files - Request', 1, ...
        subjids(idsel));
    if ~iscell(newid) || ...
        numel(newid) ~= 1 || ...
        isempty(newid{1}) || ...
        strcmpi(newid{1}, subjids{idsel})
        return;
    end
    if any(strcmpi(newid{1}, subjids))
        uiwait(warndlg('Subject ID already in the list!', ...
            'Import RFX-GLM from SPM: Warning', 'modal'));
        return;
    end
    subjids(idsel) = newid;
    hTag.LB_importrfxglm_subids.String = subjids;

% multiple IDs
else

    % ask for replacement strategy
    newid = inputdlg({'Subject ID replace-from pattern:', ...
        'Subject ID replace-to pattern'}, ...
        'Import GLM from SPM.mat files - Request', 1, ...
        {'', ''});
    if ~iscell(newid) || ...
        numel(newid) ~= 2 || ...
        isempty(newid{1}) || ...
        strcmpi(newid{1}, newid{2})
        return;
    end

    % perform replacement
    newsubjids = subjids;
    newsubjids(idsel) = regexprep(subjids(idsel), newid{:}, 'preservecase');
    if ~all(strcmp(subjids(idsel), newsubjids(idsel))) && ...
        numel(newsubjids) == numel(unique(newsubjids))
        hTag.LB_importrfxglm_subids.String = newsubjids;
    end
end

function ig_btfileadd(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;
keepbf = (hTag.CB_importrfxglm_keepbf.Value > 0);

% request one SPM.mat file
[spmfile, spmpath] = uigetfile( ...
    {'SPM.mat', 'SPM design matrices (SPM.mat)'}, ...
    'Please select an SPM.mat file to import...', ...
    'MultiSelect', 'off');
if isequal(spmfile, 0) || ...
    isequal(spmpath, 0)
    return;
end

% try loading
spmfile = [spmpath '/' spmfile];
try
    tSPM = load(spmfile);
    tSPM = tSPM.SPM;
    if ~isfield(tSPM, 'Sess') || ...
       ~isfield(tSPM, 'xX')
        error('BAD_SPM_MAT');
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    uiwait(warndlg('Invalid SPM.mat file selected.', ...
        'Import RFX-GLM from SPM: Warning', 'modal'));
    return;
end

% get subject ID
subjids = hTag.LB_importrfxglm_subids.String;
if ~iscell(subjids)
    subjids = cellstr(subjids);
end
subjid = inputdlg({'Please enter a subject ID:'}, ...
    'Import GLM from SPM.mat files - Request', 1, ...
    {sprintf('subj%02d', numel(subjids) + 1)});
if ~iscell(subjid) || ...
    numel(subjid) ~= 1 || ...
    isempty(subjid{1}) || ...
    any(strcmpi(subjid{1}, subjids))
    return;
end

% put the filename into the edit field
hTag.LB_importrfxglm_files.AddString(subjid{1});
hTag.LB_importrfxglm_files.AddString(spmfile);
if hTag.LB_importrfxglm_files.MSize > 2
    hTag.BT_importrfxglm_import.Enable = 'on';
elseif hTag.CB_importrfxglm_autcond.Value > 0
    try
        clist = ig_getclist(spmfile, keepbf);
    catch ne_eo;
        hTag.LB_importrfxglm_spmcond.Value = [];
        hTag.LB_importrfxglm_spmcond.String = {};
        hTag.LB_importrfxglm_glmcond.Value = [];
        hTag.LB_importrfxglm_glmcond.String = {};
        ne_ui.importrfxglmfromspms.ccol = zeros(0, 3);
        neuroelf_lasterr(ne_eo);
        uiwait(warndlg('Invalid SPM.mat file selected.', ...
            'Import RFX-GLM from SPM: Warning', 'modal'));
        return;
    end
    hTag.LB_importrfxglm_spmcond.Value = [];
    hTag.LB_importrfxglm_spmcond.String = clist;
    hTag.LB_importrfxglm_glmcond.Value = [];
    hTag.LB_importrfxglm_glmcond.String = clist;
    ne_ui.importrfxglmfromspms.ccol = ...
        [floor(255.999 .* rand(numel(clist) - 1, 3)); 255 .* ones(1, 3)];
end

function ig_btfileremove(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;

% remove any selected values
sidx = hTag.LB_importrfxglm_files.Value;
if ~isempty(sidx)
    subjids = hTag.LB_importrfxglm_subids.String;
    if ~iscell(subjids)
        subjids = cellstr(subjids);
    end
    files = hTag.LB_importrfxglm_files.String;
    if ~iscell(files)
        files = cellstr(files);
    end
    subjids(sidx) = [];
    files(sidx) = [];
    hTag.LB_importrfxglm_subids.String = subjids;
    hTag.LB_importrfxglm_subids.Value = [];
    hTag.LB_importrfxglm_files.String = files;
    hTag.LB_importrfxglm_files.Value = [];
    if numel(files) < 3
        hTag.BT_importrfxglm_import.Enable = 'off';
    end
end

function ig_lbspmcselect(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;
try
    hTag.LB_importrfxglm_glmcond.Value = hTag.LB_importrfxglm_spmcond.Value;
    hTag.LB_importrfxglm_glmcond.ListBoxTop = hTag.LB_importrfxglm_spmcond.ListBoxTop;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

function ig_lbglmcselect(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;
try
    hTag.LB_importrfxglm_spmcond.Value = hTag.LB_importrfxglm_glmcond.Value;
    hTag.LB_importrfxglm_spmcond.ListBoxTop = hTag.LB_importrfxglm_glmcond.ListBoxTop;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

function ig_btcondup(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;
mvidx = hTag.LB_importrfxglm_spmcond.Value;
mvstr = hTag.LB_importrfxglm_spmcond.String;
glmstr = hTag.LB_importrfxglm_glmcond.String;
if ~iscell(mvstr)
    mvstr = cellstr(mvstr);
end
if ~iscell(glmstr)
    glmstr = cellstr(glmstr);
end
if isempty(mvidx) || ...
    isequal(mvidx, 1)
    return;
end
[reidx, newidx] = moveinlist(numel(mvstr), mvidx, -1, true);
hTag.LB_importrfxglm_spmcond.String = mvstr(reidx);
hTag.LB_importrfxglm_spmcond.Value = newidx;
hTag.LB_importrfxglm_glmcond.String = glmstr(reidx);
hTag.LB_importrfxglm_glmcond.Value = newidx;
ne_ui.importrfxglmfromspms.ccol = ne_ui.importrfxglmfromspms.ccol(reidx, :);

function ig_btconddown(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;
mvidx = hTag.LB_importrfxglm_spmcond.Value;
mvstr = hTag.LB_importrfxglm_spmcond.String;
glmstr = hTag.LB_importrfxglm_glmcond.String;
if ~iscell(mvstr)
    mvstr = cellstr(mvstr);
end
if ~iscell(glmstr)
    glmstr = cellstr(glmstr);
end
if isempty(mvidx) || ...
    isequal(mvidx, numel(mvstr))
    return;
end
[reidx, newidx] = moveinlist(numel(mvstr), mvidx, 1, true);
hTag.LB_importrfxglm_spmcond.String = mvstr(reidx);
hTag.LB_importrfxglm_spmcond.Value = newidx;
hTag.LB_importrfxglm_glmcond.String = glmstr(reidx);
hTag.LB_importrfxglm_glmcond.Value = newidx;
ne_ui.importrfxglmfromspms.ccol = ne_ui.importrfxglmfromspms.ccol(reidx, :);

function ig_btcondprop(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;

% none selected
idsel = hTag.LB_importrfxglm_glmcond.Value;
if isempty(idsel)
    return;
end

% get strings
spmconds = hTag.LB_importrfxglm_spmcond.String;
glmconds = hTag.LB_importrfxglm_glmcond.String;
if ~iscell(spmconds)
    spmconds = cellstr(spmconds);
end
if ~iscell(glmconds)
    glmconds = cellstr(glmconds);
end

% single selection
if numel(idsel) == 1

    % ask for replacement ID
    newname = inputdlg({ ...
        ['BV GLM Condition name for SPM condition ' spmconds{idsel} ':']}, ...
        'Import GLM from SPM.mat files - Request', 1, ...
        glmconds(idsel));
    if ~iscell(newname) || ...
        numel(newname) ~= 1 || ...
        isempty(newname{1}) || ...
        any(strcmpi(newname{1}, {glmconds{idsel}, 'constant'}))
        return;
    end
    if any(strcmpi(newname{1}, glmconds))
        uiwait(warndlg('Condition name already in the list!', ...
            'Import RFX-GLM from SPM: Warning', 'modal'));
        return;
    end
    glmconds(idsel) = newname;
    hTag.LB_importrfxglm_glmcond.String = glmconds;

% multiple IDs
else

    % ask for replacement strategy
    newname = inputdlg({'BV GLM condition name replace-from pattern:', ...
        'BV GLM condition name replace-to pattern'}, ...
        'Import GLM from SPM.mat files - Request', 1, ...
        {'', ''});
    if ~iscell(newname) || ...
        numel(newname) ~= 2 || ...
        isempty(newname{1}) || ...
        strcmpi(newname{1}, newname{2})
        return;
    end

    % perform replacement
    newglmconds = glmconds;
    newglmconds(idsel) = regexprep(glmconds(idsel), newname{:}, 'preservecase');
    if ~all(strcmp(glmconds(idsel), newglmconds(idsel))) && ...
        numel(newglmconds) == numel(unique(newglmconds))
        hTag.LB_importrfxglm_glmcond.String = newglmconds;
    end
end

function ig_btcondcol(varargin)
global ne_ui;
condstr = ne_ui.importrfxglmfromspms.hTag.LB_importrfxglm_spmcond.String;
if ~iscell(condstr)
    condstr = cellstr(condstr);
end
ne_ui.importrfxglmfromspms.ccol = ...
    colorpicker(ne_ui.importrfxglmfromspms.ccol, condstr);

function ig_btcondadd(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;

% request two names
cnames = inputdlg({'Please enter the name as it occurs in the SPM.mat files:', ...
    'Please enter the name and color as it is to appear in the GLM file:'}, ...
    'Import GLM from SPM.mat files - Request', 1, ...
    {'spmcond', 'glmcond', sprintf('  %d', floor(255.999 .* rand(1, 3)))}, ...
    struct('Resize', 'off', 'WindowStyle', 'modal'));
if ~iscell(cnames) || ...
    numel(cnames) ~= 3 || ...
    isempty(regexpi(cnames{3}, '^\s*\d+\s+\d+\s+\d+\s*$'))
    return;
end

% test that the names are unique
spmc = hTag.LB_importrfxglm_spmcond.String;
if ~iscell(spmc)
    spmc = cellstr(spmc);
end
glmc = hTag.LB_importrfxglm_glmcond.String;
if ~iscell(glmc)
    glmc = cellstr(glmc);
end
if any(strcmpi(cnames{1}, spmc)) || ...
    any(strcmpi(cnames{2}, glmc))
    uiwait(warndlg('Condition names must be unique.', ...
        'Import RFX-GLM from SPM: Warning', 'modal'));
    return;
end
hTag.LB_importrfxglm_spmcond.AddString(cnames{1});
hTag.LB_importrfxglm_glmcond.AddString(cnames{2});
ne_ui.importrfxglmfromspms.ccol(end+1, :) = ...
    min(255, max(0, eval(['[' cnames{3} ']'])));

function ig_btcondremove(varargin)
global ne_ui;
hTag = ne_ui.importrfxglmfromspms.hTag;

% remove any selected values
sidx = hTag.LB_importrfxglm_spmcond.Value;
sid2 = hTag.LB_importrfxglm_glmcond.Value;
if ~isempty(sidx) && ...
    isequal(sidx, sid2)
    spmc = hTag.LB_importrfxglm_spmcond.String;
    if ~iscell(spmc)
        spmc = cellstr(spmc);
    end
    glmc = hTag.LB_importrfxglm_glmcond.String;
    if ~iscell(glmc)
        glmc = cellstr(glmc);
    end
    spmc(sidx) = [];
    glmc(sidx) = [];
    ne_ui.importrfxglmfromspms.ccol(sidx, :) = [];
    hTag.LB_importrfxglm_spmcond.String = spmc;
    hTag.LB_importrfxglm_glmcond.String = glmc;
    hTag.LB_importrfxglm_spmcond.Value = [];
    hTag.LB_importrfxglm_glmcond.Value = [];
end

function ig_toggleautocond(varargin)
global ne_ui;
st = ne_ui.importrfxglmfromspms;
if st.hTag.CB_importrfxglm_autcond.Value > 0
    st.hFig.SetGroupEnabled('MCond', 'off');
else
    st.hFig.SetGroupEnabled('MCond', 'on');
end

function ig_toggleautores(varargin)
global ne_ui;
st = ne_ui.importrfxglmfromspms;
if st.hTag.CB_importrfxglm_autores.Value > 0
    st.hFig.SetGroupEnabled('MResBox', 'off');
else
    st.hFig.SetGroupEnabled('MResBox', 'on');
end

function ig_toggleimpvtc(varargin)
global ne_ui;
st = ne_ui.importrfxglmfromspms;
if st.hTag.CB_importrfxglm_impvtc.Value > 0
    st.hFig.SetGroupEnabled('ImpVTC', 'on');
else
    st.hFig.SetGroupEnabled('ImpVTC', 'off');
end

function ig_btsaveas(varargin)
global ne_ui;
[glmfile, glmpath] = uiputfile( ...
    {'*.glm', 'BrainVoyager QX RFX-GLM file'}, 'Save imported GLM as');
if isequal(glmfile, 0) || ...
    isequal(glmpath, 0)
    return;
end
ne_ui.importrfxglmfromspms.hTag.ED_importrfxglm_target.String = ...
    strrep([strrep(glmpath, '\', '/') '/' glmfile], '//', '/');

function ig_btimport(varargin)
global ne_ui;
st = ne_ui.importrfxglmfromspms;
ccol = st.ccol;
hTag = st.hTag;

% get settings
spms = hTag.LB_importrfxglm_files.String;
if ~iscell(spms)
    spms = cellstr(spms);
end
if numel(spms) < 3
    uiwait(warndlg('Too few SPM.mat files to import.', ...
        'Import RFX-GLM from SPM: Warning', 'modal'));
    return;
end
opts = struct;
imeths = {'linear', 'cubic', 'lanczos3'};
opts.imeth = imeths{hTag.DD_importrfxglm_imeth.Value};
opts.pbar = st.pbar;
subjids = hTag.LB_importrfxglm_subids.String;
if ~iscell(subjids)
    subjids = cellstr(subjids);
end
opts.subjids = subjids;
spmc = hTag.LB_importrfxglm_spmcond.String;
glmc = hTag.LB_importrfxglm_glmcond.String;
if ~iscell(spmc)
    spmc = cellstr(spmc);
end
if ~iscell(glmc)
    glmc = cellstr(glmc);
end
if numel(spmc) ~= numel(glmc) || ...
    numel(unique(spmc)) ~= numel(spmc) || ...
    numel(unique(spmc)) ~= numel(spmc) || ...
    isempty(spmc)
    uiwait(warndlg('Invalid condition name configuration.', ...
        'Import RFX-GLM from SPM: Warning', 'modal'));
    return;
end
opts.vweight = (hTag.CB_importrfxglm_vweight.Value > 0);
if hTag.CB_importrfxglm_autores.Value > 0
    opts.bbox = [];
    opts.res = [];
else
    try
        opts.bbox = ...
            reshape(eval(['[' hTag.ED_importrfxglm_bbox.String ']']), ...
            2, 3);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        uiwait(warndlg('Invalid bounding box.', ...
            'Import RFX-GLM from SPM: Warning', 'modal'));
        return;
    end
    try
        opts.res = round(str2double(hTag.ED_importrfxglm_res.String));
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        uiwait(warndlg('Invalid resolution value.', ...
            'Import RFX-GLM from SPM: Warning', 'modal'));
        return;
    end
end
if hTag.CB_importrfxglm_impvtc.Value > 0
    impvtc = ddeblank(hTag.ED_importrfxglm_impvtc.String);
    if ~isempty(impvtc) && ...
        any([1, 2] == sum(impvtc == '#')) && ...
        sum(impvtc == '%') == 1
        opts.impvtc = impvtc;
    else
        uiwait(warndlg('Invalid VTC import filename pattern.', ...
            'Import RFX-GLM from SPM: Warning', 'modal'));
        return;
    end
end
opts.cond = cell2struct([glmc(:), spmc(:)], {'bvname', 'spmname'}, 2);
for cc = 1:numel(opts.cond)
    opts.cond(cc).color = ccol(cc, :);
end
fname = hTag.ED_importrfxglm_target.String;
fname = fname(:)';
if ~any(fname == '<' & fname == '>') && ...
    numel(fname) > 4 && ...
    strcmpi(fname(end-3:end), '.glm')
    opts.filename = fname;
end
opts.keepbf = (hTag.CB_importrfxglm_keepbf.Value > 0);
opts.pbar = st.pbar;
if hTag.CB_importrfxglm_psctr.Value > 0
    opts.trans = 'psc';
else
    opts.trans = 'none';
end
try
    st.hFig.Visible = 'off';
    drawnow;
    ne_ui.importrfxglmfromspms.glm = importrfxglmfromspms(spms, opts);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    st.hFig.Visible = 'on';
    uiwait(warndlg(sprintf('An error occurred during the import:\n\n%s', ...
        ne_eo.message), 'Import RFX-GLM from SPM: Warning', 'modal'));
    return;
end
st.hFig.CloseRequestFcn = '';
st.hFig.Delete;
