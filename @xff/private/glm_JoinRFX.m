function xo3 = glm_JoinRFX(xo, xo2, prds, opts)
% GLM::JoinRFX  - join two random effects GLMs
%
% FORMAT:       combined = glm1.JoinRFX(glm2 [, prds [, opts]]);
%
% Input fields:
%
%       glm2        second RFX GLM to be added to glm1
%       prds        list of predictors (if empty, use first GLM's)
%       opts        additional options
%        .imissing  ignore missing conditions (default: false)
%
% Joins the Random Effects GLM results of glm1 and glm2 to one
% single structure.
%
% Using: catstruct.

% Version:  v1.1
% Build:    20030419
% Date:     Mar-04 2020, 7:50 PM EST
% Author:   Jochen Weber, NeuroElf.net, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010 - 2020, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || numel(xo2) ~= 1 || ...
   ~xffisobject(xo, true, 'glm') || ~xffisobject(xo2, true, 'glm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc1 = xo.C;
bc2 = xo2.C;
rtv1 = bc1.RunTimeVars;
rtv2 = bc2.RunTimeVars;
if bc1.ProjectTypeRFX ~= 1 || bc2.ProjectTypeRFX ~= 1 || bc1.SeparatePredictors < 1
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
prds1 = glm_SubjectPredictors(xo);
prds2 = glm_SubjectPredictors(xo2);
if nargin < 3 || ~iscell(prds) || isempty(prds) || ~all(cellfun(@ischar, prds(:)))
    prds = prds1;
else
    prds = prds(:);
    if any(strcmpi(prds, 'constant'))
        prds(strcmpi(prds, 'constant')) = [];
    end
    prds = [prds; {'Constant'}];
end
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'imissing') || ~islogical(opts.imissing) || numel(opts.imissing) ~= 1
    opts.imissing = false;
end
if ~opts.imissing && ...
   (bc1.NrOfSubjectPredictors ~= bc2.NrOfSubjectPredictors || ...
    numel(prds1) ~= numel(prds2) || ...
    ~all(strcmpi(sort(prds1), sort(prds2))))
    error('neuroelf:xff:badObject', ...
        'Invalid object(s) given in call. Condition names don''t match.');
end

% check type and sizes of GLMs
if ~isfield(bc1.GLMData, 'RFXGlobalMap') || ~isfield(bc2.GLMData, 'RFXGlobalMap') || ...
    bc1.ProjectType ~= bc2.ProjectType || bc1.FileVersion ~= bc2.FileVersion || ...
    bc1.SeparatePredictors ~= bc2.SeparatePredictors || ...
    bc1.TransformationType ~= bc2.TransformationType || ...
    bc1.Resolution ~= bc2.Resolution || bc1.SerialCorrelation ~= bc2.SerialCorrelation || ...
    bc1.XStart ~= bc2.XStart || bc1.XEnd ~= bc2.XEnd || ...
    bc1.YStart ~= bc2.YStart || bc1.YEnd ~= bc2.YEnd || ...
    bc1.ZStart ~= bc2.ZStart || bc1.ZEnd ~= bc2.ZEnd || ...
    bc1.CortexBasedStatistics ~= bc2.CortexBasedStatistics || ...
   ~isequal(bc1.RunTimeVars.TrfPlus, bc2.RunTimeVars.TrfPlus)
    error('neuroelf:xff:badObject', ...
        'Invalid object(s) given in call. Crucial fields mismatch.');
end

% get two subject ID lists and don't allow any overlap
id1 = glm_Subjects(xo);
id2 = glm_Subjects(xo2);
idc = intersect(id1, id2);
if ~isempty(idc)
    idc = sprintf('%s, ', idc{:});
    error('neuroelf:xff:badArgument', ...
        'Subjects %s are in both GLMs; invalid input.', idc(1:end-2));
end

% build new file
NS1 = bc1.NrOfStudies;
NS2 = bc2.NrOfStudies;
NSS = NS1 + NS2;
NP1 = bc1.NrOfPredictors;
NP2 = bc2.NrOfPredictors;
NPS = NP1 + NP2;
NSB = bc1.NrOfSubjects + bc2.NrOfSubjects;
NTP = bc1.NrOfTimePoints + bc2.NrOfTimePoints;
NC1 = bc1.NrOfConfoundsPerStudy;
NC2 = bc2.NrOfConfoundsPerStudy;
SC1 = numel(bc1.GLMData.Subject);
SC2 = numel(bc2.GLMData.Subject);
xo3 = aft_CopyObject(xo);
xo3.F = '';
bc3 = xo3.C;
bc3.NrOfSubjects = NSB;
bc3.NrOfTimePoints = NTP;
bc3.NrOfConfounds = bc3.NrOfConfounds + bc2.NrOfConfounds;
bc3.NrOfStudies = NSS;
bc3.NrOfStudiesWithConfounds = bc3.NrOfStudiesWithConfounds + bc2.NrOfStudiesWithConfounds;
bc3.NrOfConfoundsPerStudy = [NC1(:); NC2(:)]';
bc3.NrOfVoxelsForBonfCorrection = ...
    max(bc3.NrOfVoxelsForBonfCorrection, bc2.NrOfVoxelsForBonfCorrection);
bc3.Study = ne_methods.catstruct(bc3.Study(:), bc2.Study(:))';

% global map
gm1 = bc1.GLMData.RFXGlobalMap;
gm1 = reshape(gm1(:), size(gm1));
gm2 = bc2.GLMData.RFXGlobalMap;
gm2 = reshape(gm2(:), size(gm2));
bc3.GLMData.RFXGlobalMap = max(gm1, gm2);

% short route if all conditions are the same
if isequal(prds, prds1) && isequal(prds, prds2)
    
    % set fields that change
    bc3.NrOfPredictors = NPS;
    bc3.Predictor = ne_methods.catstruct(bc3.Predictor(:), bc2.Predictor(:))';

    % reorder predictors
    bc3.Predictor = bc3.Predictor([ ...
                            1:(NP1 - SC1), ...
                    (NP1 + 1):(NP1 + NP2 - SC2), ...
              (1 + NP1 - SC1):NP1, ...
        (1 + NP1 + NP2 - SC2):(NP1 + NP2)]);

    % rename predictors
    for pc = 1:(NP1 + NP2)
        bc3.Predictor(pc).Name1 = sprintf('Predictor: %d', pc);
    end

    % combine GLMData Subjects
    bc3.GLMData.Subject = ne_methods.catstruct(bc3.GLMData.Subject(:), bc2.GLMData.Subject(:))';

% more work if conditions don't match
else
    
    % set fields that change
    bc3.NrOfSubjectPredictors = numel(prds);
    bc3.NrOfPredictors = numel(prds) * bc3.NrOfSubjects;
    pnames1 = {bc1.Predictor.Name2};
    pnames2 = {bc2.Predictor.Name2};
    bc3.Predictor = repmat(bc3.Predictor(1), 1, bc3.NrOfPredictors);
    tp = 1;
    for sc = 1:numel(id1)
        for pc = 1:(numel(prds)-1)
            pstring = sprintf('Subject %s: %s', id1{sc}, prds{pc});
            pmatch = find(strcmpi(pstring, pnames1));
            if ~isempty(pmatch)
                bc3.Predictor(tp) = bc1.Predictor(pmatch(1));
            else
                bc3.Predictor(tp).Name2 = pstring;
                bc3.Predictor(tp).RGB = [floor(255.999 * rand(1, 3)); zeros(3, 3)];
            end
            bc3.Predictor(tp).Name1 = sprintf('Predictor: %d', tp);
            tp = tp + 1;
        end
    end
    for sc = 1:numel(id2)
        for pc = 1:(numel(prds)-1)
            pstring = sprintf('Subject %s: %s', id2{sc}, prds{pc});
            pmatch = find(strcmpi(pstring, pnames2));
            if ~isempty(pmatch)
                bc3.Predictor(tp) = bc2.Predictor(pmatch(1));
            else
                bc3.Predictor(tp).Name2 = pstring;
                bc3.Predictor(tp).RGB = [floor(255.999 * rand(1, 3)); zeros(3, 3)];
            end
            bc3.Predictor(tp).Name1 = sprintf('Predictor: %d', tp);
            tp = tp + 1;
        end
    end
    for sc = 1:numel(id1)
        bc3.Predictor(tp).Name1 = sprintf('Predictor: %d', tp);
        bc3.Predictor(tp).Name2 = sprintf('Subject %s: Constant', id1{sc});
        bc3.Predictor(tp).RGB = [255, 255, 255; zeros(3, 3)];
        tp = tp + 1;
    end
    for sc = 1:numel(id2)
        bc3.Predictor(tp).Name1 = sprintf('Predictor: %d', tp);
        bc3.Predictor(tp).Name2 = sprintf('Subject %s: Constant', id2{sc});
        bc3.Predictor(tp).RGB = [255, 255, 255; zeros(3, 3)];
        tp = tp + 1;
    end
    
    % create new subjects maps data
    if ndims(gm1) < 3
        m2d = true;
        bc3.GLMData.Subject = repmat(struct('BetaMaps',  ...
            single(NaN * zeros([numel(gm1), numel(prds)]))), 1, bc3.NrOfSubjects);
    else
        m2d = false;
        gs1 = size(gm1);
        gs1(end+1) = numel(prds);
        bc3.GLMData.Subject = repmat(struct('BetaMaps',  ...
            single(NaN * zeros(gs1))), 1, bc3.NrOfSubjects);
    end
    
    % assign beta content
    pmatch = ne_methods.multimatch(prds, prds1);
    pmatchm = find(pmatch > 0);
    for sc1 = 1:numel(id1)
        if m2d
            bc3.GLMData.Subject(sc1).BetaMaps(:, pmatchm) = ...
                bc1.GLMData.Subject(sc1).BetaMaps(:, pmatch(pmatchm));
        else
            bc3.GLMData.Subject(sc1).BetaMaps(:, :, :, pmatchm) = ...
                bc1.GLMData.Subject(sc1).BetaMaps(:, :, :, pmatch(pmatchm));
        end
    end
    pmatch = ne_methods.multimatch(prds, prds2);
    pmatchm = find(pmatch > 0);
    for sc2 = 1:numel(id2)
        if m2d
            bc3.GLMData.Subject(sc1+sc2).BetaMaps(:, pmatchm) = ...
                bc2.GLMData.Subject(sc2).BetaMaps(:, pmatch(pmatchm));
        else
            bc3.GLMData.Subject(sc1+sc2).BetaMaps(:, :, :, pmatchm) = ...
                bc2.GLMData.Subject(sc2).BetaMaps(:, :, :, pmatch(pmatchm));
        end
    end
end
if isfield(bc1.GLMData, 'RunTimeVars') && isfield(bc2.GLMData, 'RunTimeVars')
    if isfield(bc1.GLMData.RunTimeVars, 'Subject') && isfield(bc2.GLMData.RunTimeVars, 'Subject')
        bc3.GLMData.RunTimeVars.Subject = ne_methods.catstruct( ...
            bc1.GLMData.RunTimeVars.Subject(:), ...
            bc2.GLMData.RunTimeVars.Subject(:));
    else
        if isfield(bc3.GLMData.RunTimeVars, 'Subject')
            bc3.GLMData.RunTimeVars = rmfield(bc3.GLMData.RunTimeVars, 'Subject');
        end
    end
else
    if isfield(bc3.GLMData, 'RunTimeVars')
        bc3.GLMData = rmfield(bc3.GLMData, 'RunTimeVars');
    end
end

% clear RunTimeVars variables that are unsafe
if opts.imissing
    bc3.RunTimeVars.Contrasts = cell(0, 2);
    bc3.RunTimeVars.ContrastColors = zeros(0, 3);
end
bc3.RunTimeVars.CovariatesColors = cell(0, 2);
bc3.RunTimeVars.CovariatesData = [];
bc3.RunTimeVars.CovariatesNames = {};
bc3.RunTimeVars.Groups = cell(0, 2);
bc3.RunTimeVars.Map(:) = [];
bc3.RunTimeVars.MapSelection = {{}, []};
if ~isempty(rtv1.MotionParameters) && ~isempty(rtv2.MotionParameters)
    bc3.RunTimeVars.MotionParameters = [rtv1.MotionParameters(:); rtv2.MotionParameters(:)];
else
    bc3.RunTimeVars.MotionParameters = {};
end
if isfield(rtv1, 'PerSubjectGLMs') && ~isempty(rtv1.PerSubjectGLMs) && ...
    isfield(rtv2, 'PerSubjectGLMs') && ~isempty(rtv2.PerSubjectGLMs)
    bc3.RunTimeVars.PerSubjectGLMs = [rtv1.PerSubjectGLMs(:)', rtv2.PerSubjectGLMs(:)'];
else
    bc3.RunTimeVars.PerSubjectGLMs = ne_methods.emptystruct({ ...
        'SourceFile', 'SubjectID', 'Predictors', 'iXX', 'DF1', 'SEMap', 'ARLag'});
end
bc3.RunTimeVars.PredictorColors = [];
try
    snfc = fieldnames(rtv2.SubjectSPMsn);
catch
    snfc = {};
end
if ~isstruct(rtv1.SubjectSPMsn)
    bc3.RunTimeVars.SubjectSPMsn = struct;
end
for fc = snfc(:)'
    bc3.RunTimeVars.SubjectSPMsn(1).(fc{1}) = rtv2.SubjectSPMsn.(fc{1});
end
try
    snfc = fieldnames(rtv2.SubjectTrfPlus);
catch
    snfc = {};
end
if ~isstruct(rtv1.SubjectTrfPlus)
    bc3.RunTimeVars.SubjectTrfPlus = struct;
end
for fc = snfc(:)'
    bc3.RunTimeVars.SubjectTrfPlus(1).(fc{1}) = rtv2.SubjectTrfPlus.(fc{1});
end
bc3.RunTimeVars.SubSels = {'default', {}};
bc3.RunTimeVars.SubSelsSel = 'default';
bc3.RunTimeVars.SubSel = {};
bc3.RunTimeVars.UseGroups = false;

% put back
xo3.C = bc3;
