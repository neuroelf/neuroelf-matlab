function glm = importrfxglmfromfsl(opts)
% importrfxglmfromfsl  - import FSL stats maps to a BrainVoyager GLM file
%
% FORMAT:       glm = importrfxglmfromfsl(opts)
%
% Input fields:
%
%       opts        options for the import
%        .bbox      bounding box to use (in BVS notation)
%        .cond      1xC structure with fields
%         .bvname   condition name for the BV file
%         .color    1x3 RGB color for condition (if unset, randomly chosen)
%         .pathpat  pattern for folder to search for conditions
%        .filename  output GLM filename (otherwise unsaved)
%        .imeth     interpolation 'cubic', 'lanczos3', {'linear'}, 'nearest'
%        .pbar      either xprogress or xfigure:XProgress object
%        .statspath FSL filename to be imported (per folder, 'pe1.nii.gz')
%        .trspersub TRs (volumes) per subject (if not given, set to 500)
%
% Output fields:
%
%       glm         GLM object (saved if .filename is given)
%

% Version:  v1.1
% Build:    21092811
% Date:     Sep-28 2021, 11:51 AM EST
% Author:   Jochen Weber, Memorial Sloan Kettering Cancer Center, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2021, Jochen Weber
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
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    error( ...
        'neuroelf:BadArgument', ...
        'Argument opts must be set and 1x1 struct for this function.' ...
    );
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
   ~isfield(opts.cond, 'pathpat')
    error( ...
        'neuroelf:BadArgument', ...
        'Argument opts must have a valid cond and pathpat fields.' ...
    );
else
    glmnames = {opts.cond.bvname};
    pathpats = {opts.cond.pathpat};
    if any(cellfun('isempty', glmnames)) || ...
        any(cellfun('isempty', pathpats)) || ...
        numel(unique(lower(glmnames))) ~= numel(glmnames)
        error( ...
            'neuroelf:BadArgument', ...
            'Bad condition names or path patterns requested.' ...
        );
    end
end
if ~isfield(opts.cond, 'color')
    opts.cond(1).color = [];
end
numconds = numel(opts.cond);
condfolds = pathpats;
condnums = zeros(numconds, 1);
for cc = 1:numconds
    [stf, stp, ste] = fileparts(pathpats{cc});
    condfolds{cc} = findfiles(stf, [stp, ste], '-d1D');
    condnums(cc) = numel(condfolds{cc});
    if ~isa(opts.cond(cc).color, 'double') || ...
        numel(opts.cond(cc).color) ~= 3 || ...
        any(opts.cond(cc).color < 0 | opts.cond(cc).color > 255)
        opts.cond(cc).color = floor(255.999 .* rand(1, 3));
    else
        opts.cond(cc).color = round(opts.cond(cc).color(:)');
    end
end
[numsubs, maxp] = max(condnums);
if min(condnums) < 3
    error( ...
        'neuroelf:importrfxglmfromfsl:tooFewSubjects', ...
        'Less than 3 subjects detected across conditions.' ...
    );
end
if ~all(condnums == numsubs)
    warning( ...
        'neuroelf:importrfxglmfromfsl:subjectsMismatch', ...
        'Number of subject folders found mismatch across conditions.' ...
    );
end
pathparts = condfolds{maxp};
for sc = 1:numsubs
    pathparts{sc} = splittocellc(pathparts{sc}, filesep);
end
maxparts = min(cellfun('prodofsize', pathparts));
partidx = 0;
for pc = 1:maxparts
    subparts = cell(numsubs, 1);
    for sc = 1:numsubs
        subparts{sc} = pathparts{sc}{pc};
    end
    if numel(unique(subparts)) == numsubs
        partidx = pc;
        break;
    end
end
if partidx < 1
    error( ...
        'neuroelf:importrfxglmfromfsl:detectionError', ...
        'Could not detect subject IDs.' ...
    );
end
subsparts = subparts;
for sc = 1:numsubs
    subsparts{sc} = splittocellc(subparts{sc}, ' -_', true, true);
end
maxsparts = min(cellfun('prodofsize', subsparts));
spartidx = 0;
for pc = 1:maxsparts
    subssparts = cell(numsubs, 1);
    for sc = 1:numsubs
        subssparts{sc} = subsparts{sc}{pc};
    end
    if numel(unique(subssparts)) == numsubs
        spartidx = pc;
        break;
    end
end
if spartidx == 0
    opts.subjids = subparts;
else
    opts.subjids = subssparts;
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
if ~isfield(opts, 'statspath') || ...
   ~ischar(opts.statspath) || ...
    isempty(opts.statspath)
    opts.statspath = 'pe1.nii.gz';
end
if ~isfield(opts, 'trspersub') || ...
   ~isa(opts.trspersub, 'double') || ...
    numel(opts.trspersub) ~= 1 || ...
    isinf(opts.trspersub) || ...
    isnan(opts.trspersub) || ...
    opts.trspersub < 80
    opts.trspersub = 500;
else
    opts.trspersub = round(opts.trspersub);
end

% try to locate FSL stats files and store essential information
condfiles = cell(numsubs, numconds);
for cc = 1:numconds
    for sc = 1:numsubs
        sci = find(multimatch(condfolds{cc}, opts.subjids(sc), true));
        if numel(sci) > 1
            error( ...
                'neuroelf:importrfxglmfromfsl:invalidSubjectId', ...
                'Invalid subject ID %s for condition %d (%s).', ...
                opts.subjids{sc}, cc, glmnames{cc} ...
            );
        end
        if isempty(sci)
            continue;
        end
        scondfile = findfiles(condfolds{cc}{sci}, opts.statspath, '-d1');
        if numel(scondfile) ~= 1
            error( ...
                'neuroelf:importrfxglmfromfsl:missingStatsFile', ...
                'Missing stats file for subject ID %s, condition %d (%s).', ...
                opts.subjids{sc}, cc, glmnames{cc} ...
            );
        end
        condfiles(sc, cc) = scondfile;
    end
end
ucondfiles = (cellfun('prodofsize', condfiles) > 0);
fcondfile = findfirst(ucondfiles(:));

% time points
nroftotaltimepoints = numsubs * opts.trspersub;

% check spatial layout/dimension
Vbetaobj = {[]};
try
    Vbetaobj{1} = xff(condfiles{fcondfile});
catch ne_eo;
    clearxffobjects(Vbetaobj);
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid Analyze/NIftI file: %s', ...
        condfiles{fcondfile}, ne_eo.message ...
    );
end

% add constant condition
opts.cond(end + 1).bvname = 'Constant';
opts.cond(end).color = [255, 255, 255];

% get resolution
opts.res = floor(0.05 + mean(sqrt(sum( ...
    Vbetaobj{1}.CoordinateFrame.Trf(1:3, 1:3) .^ 2))));

% temporary clean-up
clearxffobjects(Vbetaobj);

% create GLM structure
glm = xff('new:glm');

% make some initial settings
glm.ProjectType = 1;
glm.ProjectTypeRFX = 1;
glm.NrOfSubjects = numsubs;
glm.NrOfSubjectPredictors = numel(opts.cond);
glm.NrOfTimePoints = nroftotaltimepoints;
glm.NrOfPredictors = numsubs * glm.NrOfSubjectPredictors;
glm.NrOfConfounds = numsubs;
glm.NrOfStudies = numsubs;
glm.NrOfStudiesWithConfounds = numsubs;
glm.NrOfConfoundsPerStudy = numsubs;
glm.SeparatePredictors = 2;
glm.Resolution = opts.res;
glm.SerialCorrelation = 0;
glm.Study(numsubs).NrOfTimePoints = 0;
glm.Predictor(glm.NrOfPredictors).Name1 = '';
prc = 1;

% initialize progress bar
if isempty(opts.pbar)
    try
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 264, 640, 36]);
        xprogress(pbar, 'settitle', ...
            sprintf('Importing %d subjects'' FSL stats to RFX-GLM...', numsubs));
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
        sprintf('importrfxglmfromfsl: Importing subject 1/%d...', numsubs));
    pbarn = 'importrfxglmfromfsl: ';
end

% iterate over subjects
for sc = 1:numsubs

    % update progress bar
    if ~isempty(pbar)
        pbar.Progress((sc - 1) / numsubs, ...
            sprintf('%sImporting subject ''%s'' (%d/%d)...', pbarn, ...
            opts.subjids{sc}, sc, numsubs));
    end

    % temporarily import beta maps into VMP
    umaps = find(ucondfiles(sc, :));
    bvmp = importvmpfromspms(condfiles(sc, umaps)', 'b', ...
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
        glm.GLMData.Subject(2:numsubs) = glm.GLMData.Subject(1);
        omask = zeros(size(bvmp.Map(1).VMPData));
    end

    % set study & predictors
    glm.Study(sc).NrOfTimePoints = opts.trspersub;
    glm.Study(sc).NameOfAnalyzedFile = sprintf('<subject %s>', opts.subjids{sc});
    glm.Study(sc).NameOfSDMFile = '<model.txt>';
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
    mask = mask & (~all(vmpd == 0, 4));
    omask = omask + double(mask);

    % combine beta maps (plain average for now)
    bmaps = glm.GLMData.Subject(sc).BetaMaps;
    for pc = 1:numel(umaps)
        bmaps(:, :, :, umaps(pc)) = vmpd(:, :, :, pc);
    end
    glm.GLMData.Subject(sc).BetaMaps = bmaps;
end

% apply final masking
glm.GLMData.RFXGlobalMap = single(omask >= (0.75 * numsubs));
glm.NrOfVoxelsForBonfCorrection = ...
    double(sum(glm.GLMData.RFXGlobalMap(:)));

% add last predictors
for sc = 1:numsubs
    glm.Predictor(prc).Name1 = sprintf('Predictor: %d', prc);
    glm.Predictor(prc).Name2 = ...
        sprintf('Subject %s: constant', opts.subjids{sc});
    glm.Predictor(prc).RGB = [opts.cond(end).color; zeros(3, 3)];
    prc = prc + 1;
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
        glm.SaveRunTimeVars;
    catch ne_eo;
        warning( ...
            'neuroelf:SaveError', ...
            'Error saving GLM file (%s). Please do so manually.', ...
            ne_eo.message ...
        );
    end
end
