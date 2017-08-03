function glm = importrfxglmfromafni(heads, opts)
% importrfxglmfromafni  - import AFNI beta maps to a BrainVoyager GLM file
%
% FORMAT:       glm = importrfxglmfromafni(heads, opts)
%
% Input fields:
%
%       heads       list of AFNI GLM *.HEAD filenames to use for import
%       opts        options for the import
%        .bbox      bounding box to use (in BVS notation)
%        .cond      1xC structure with fields (default: auto detect)
%         .afniname AFNI based name (or pattern)
%         .bvname   condition name for the BV file
%         .color    1x3 RGB color for condition
%        .filename  output GLM filename (otherwise unsaved)
%        .imeth     interpolation 'cubic', 'lanczos3', {'linear'}, 'nearest'
%        .pbar      either xprogress or xfigure:XProgress object
%        .res       resolution (default: floor of AFNI resolution)
%        .subjids   subject identifiers (default: auto detect)
%
% Output fields:
%
%       glm         GLM object (saved if .filename is given)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:21 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% argument check
if nargin < 1 || ...
   ~iscell(heads) || ...
    numel(heads) < 3
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
   ~isfield(opts.cond, 'afniname') || ...
   ~isfield(opts.cond, 'bvname') || ...
   ~isfield(opts.cond, 'color')
    opts.cond = emptystruct({'afniname', 'bvname', 'color'});
else
    afninames = {opts.cond.afniname};
    glmnames = {opts.cond.bvname};
    if any(cellfun('isempty', afninames)) || ...
        numel(unique(lower(afninames))) ~= numel(afninames) || ...
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
    numel(opts.subjids) ~= numel(heads)
    opts.subjids = {};
end
if numel(heads) == 1 && ...
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
    if numel(opts.subjids) ~= numel(heads) || ...
        numel(unique(lower(opts.subjids))) ~= numel(heads)
        opts.subjids = {};
    end
end

% try load HEADs and store essential information
hx = cell(numel(heads), 1);
hxm = hx;
hxn = zeros(numel(heads), 2);
for sc = 1:numel(heads)
    try
        if ~ischar(heads{sc}) || ...
            isempty(heads{sc})
            heads{sc} = sprintf('No/empty string at position %d.', sc);
            error('BADSTRING');
        end
        heads{sc} = heads{sc}(:)';
        hx{sc} = xff(heads{sc});
        s = hx{sc};
        if ~isxff(s, 'head') || ...
            all(cellfun('isempty', regexp(s.MapNames, '\#0_Coef$')))
            error('BADHEADXFF');
        end
        hxm{sc} = s.MapNames;
        m = hxm{sc};
        ffs = findfirst(strcmpi(m, 'full_fstat'));
        if isempty(ffs)
            error('NOFULLFSTAT');
        end
        hxn(sc, :) = s.Brick(ffs).FuncParams(1:2);
    catch ne_eo;
        error( ...
            'neuroelf:BadAFNIFile', ...
            'Invalid AFNI GLM file: %s (%s).', ...
            heads{sc}, ne_eo.message ...
        );
    end
end

% subject ID's ?
if isempty(opts.subjids)
    subjids = char(heads);
    subjids(subjids(:) == '.' | subjids(:) == ' ') = '_';
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
    opts.subjids = strrep(subjids, filesep, '_');
end

% time points
timepoints = sum(hxn, 2);
nroftotaltimepoints = sum(timepoints);

% check spatial layout/dimension
VResDims = zeros(numel(hx), 19);
for sc = 1:numel(hx)
    VResDims(sc, :) = [lsqueeze(hx{sc}.CoordinateFrame.Trf(:))', ...
        hx{sc}.CoordinateFrame.Dimensions(1:3)];
end
if any(any(diff(VResDims, 1, 1)))
    warning( ...
        'neuroelf:warning', ...
        'Incompatible map sizes. Output must be taken carefully!!' ...
    );
end

% create list of regressors from first subject (in order of appearance)
uname = hxm{1};
sname = uname(~cellfun('isempty', regexp(uname, '\#0_Coef$')));
if numel(sname) < hxn(1)
    clearxffobjects(hx);
    error( ...
        'neuroelf:BadGLMObject', ...
        'Too few NAME#0_Coef maps in GLMs.' ...
    );
end
sname = sname(1:hxn(1));

% and ensure regressors are available in all other subjects!
for sc = 2:numel(hxm)
    sname(multimatch(sname, hxm{sc}) < 1) = [];
    if isempty(sname)
        break;
    end
end
sname = strrep(sname, '#0_Coef', '');

% no common regressors?
if isempty(sname)
    clearxffobjects(hx);
    error( ...
        'neuroelf:BadArgument', ...
        'Condition names don''t match between subject.' ...
    );
end

% no conditions named
if isempty(opts.cond)

    % put into cond list
    opts.cond(numel(sname)).afniname = '';
    for cc = 1:numel(sname)
        opts.cond(cc).afniname = sname{cc};
        opts.cond(cc).bvname = sname{cc};
        opts.cond(cc).color = floor(255.999 * rand(1, 3));
    end
    opts.cond(end).color = [255, 255, 255];

% for named conditions
else

    % check they exist in all subjects
    for cc = 1:numel(opts.cond)
        if ~any(strcmpi(opts.cond(cc).afniname, sname))
            clearxffobjects(hx);
            error( ...
                'neuroelf:BadArgument', ...
                'Requested condition (%s) not found in HEAD files.', ...
                opts.cond(cc).afniname ...
            );
        end
    end
end

% get HEAD condition names
afnicn = {opts.cond(:).afniname};
consti = find(strcmpi(afnicn, 'constant') | strcmpi(afnicn, 'rest'));
if ~isempty(consti)
    opts.cond(consti) = [];
end
afnicn = {opts.cond(:).afniname};
glmcn = {opts.cond(:).bvname};
glmcc = cat(1, opts.cond(:).color);
opts.cond(end + 1).afniname = 'Rest';
opts.cond(end).bvname = 'Constant';
opts.cond(end).color = [255, 255, 255];

% check bounding box and resolution
if isempty(opts.bbox)
    hxtrf = hx{1}.CoordinateFrame.Trf;
    pos1 = hxtrf * ones(4, 1);
    posn = hxtrf * [(1 + size(hx{1}.Brick(1).Data)'); 1];
    bbox = 128 - [max([pos1, posn], [], 2), min([pos1, posn], [], 2)]';
    opts.bbox = round([bbox(1, [2, 3, 1]); bbox(2, [2, 3, 1]) - 1]);
end
if isempty(opts.res)
    opts.res = floor(0.05 + mean(sqrt(sum( ...
        hx{1}.CoordinateFrame.Trf(1:3, 1:3) .^ 2))));
end

% create GLM structure
glm = xff('new:glm');

% make some initial settings
glm.ProjectType = 1;
glm.ProjectTypeRFX = 1;
glm.NrOfSubjects = numel(hx);
glm.NrOfSubjectPredictors = numel(opts.cond);
glm.NrOfTimePoints = nroftotaltimepoints;
glm.NrOfPredictors = numel(hx) * glm.NrOfSubjectPredictors;
glm.NrOfConfounds = numel(hx);
glm.NrOfStudies = numel(hx);
glm.NrOfStudiesWithConfounds = numel(hx);
glm.NrOfConfoundsPerStudy = ones(1, numel(hx));
glm.SeparatePredictors = 2;
glm.Resolution = opts.res;
glm.SerialCorrelation = 0;
glm.Study(numel(hx)).NrOfTimePoints = 0;
glm.Predictor(glm.NrOfPredictors).Name1 = '';
prc = 1;
studyxtcs = {};
studysdms = {};
studysdmc = [];
studysdmn = [];

% initialize progress bar
if isempty(opts.pbar)
    try
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 264, 640, 36]);
        xprogress(pbar, 'settitle', ...
            sprintf('Importing %d subjects'' HEAD files to RFX-GLM...', numel(hx)));
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
        sprintf('importrfxglmfromafni: Importing subject 1/%d...', numel(hx)));
    pbarn = 'importrfxglmfromafni: ';
end

% iterate over subjects
for sc = 1:numel(hx)

    % update progress bar
    if ~isempty(pbar)
        pbar.Progress((sc - 1) / numel(hx), ...
            sprintf('%sImporting subject ''%s'' (%d/%d)...', pbarn, ...
            opts.subjids{sc}, sc, numel(hx)));
    end

    % temporarily import beta maps into VMP
    bvmp = importvmpfromspms(hx(sc), 'b', ...
        opts.bbox, opts.res, opts.imeth);

    % GLM settings
    if sc == 1
        glm.XStart = bvmp.XStart;
        glm.XEnd = bvmp.XEnd;
        glm.YStart = bvmp.YStart;
        glm.YEnd = bvmp.YEnd;
        glm.ZStart = bvmp.ZStart;
        glm.ZEnd = bvmp.ZEnd;
        glm.GLMData.RFXGlobalMap = single(ones(size(bvmp.Map(1).VMPData)));
        glm.NrOfVoxelsForBonfCorrection = numel(glm.GLMData.RFXGlobalMap);
        glm.GLMData.Subject(1).BetaMaps = ...
            single(zeros([size(bvmp.Map(1).VMPData), glm.NrOfSubjectPredictors]));
        glm.GLMData.Subject(2:numel(hx)) = glm.GLMData.Subject(1);
        omask = zeros(size(bvmp.Map(1).VMPData));
        odsz = size(omask);
    end

    % set study & predictors
    glm.Study(sc).NrOfTimePoints = timepoints(sc);
    glm.Study(sc).NameOfAnalyzedFile = [opts.subjids{sc} '_AFNI.vtc'];
    glm.Study(sc).NameOfSDMFile = '<AFNI-GLM>';
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

    % combine beta maps (plain average for now)
    bmaps = glm.GLMData.Subject(sc).BetaMaps;
    for pc = 1:numel(opts.cond)

        % find beta maps of condition
        uidx = findfirst(strcmpi(hxm{sc}, [opts.cond(pc).afniname '#0_coef']));
        if ~isempty(uidx)
            bmaps(:, :, :, pc) = vmpd(:, :, :, uidx);
        end
    end

    % eliminate invalid content
    bmaps(isinf(bmaps) | isnan(bmaps)) = 0;
    glm.GLMData.Subject(sc).BetaMaps = bmaps;
end

% add last predictors
for sc = 1:numel(hx)
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
