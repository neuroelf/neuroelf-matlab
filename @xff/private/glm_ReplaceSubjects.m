function xo = glm_ReplaceSubjects(xo, xon, sid)
% GLM::ReplaceSubject  - replace subject(s) in one GLM by data from another
%
% FORMAT:       [glm = ] glm.ReplaceSubject(source, sid);
%
% Input fields:
%
%       source      second (source) GLM from which to take data to replace
%       sid         either single subject (char) or list of subjects (cell)
%
% Output fields:
%
%       glm         altered GLM object
%
% Using: multimatch.

% Version:  v1.1
% Build:    16032621
% Date:     Mar-26 2016, 9:36 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% check arguments
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm') || ...
    numel(xon) ~= 1 || ~xffisobject(xon, true, 'glm') || xo.C.ProjectTypeRFX ~= 1 || ...
    xo.C.ProjectType ~= xon.C.ProjectType || xon.C.ProjectTypeRFX ~= 1 || ...
   ~isequal(aft_Layout(xo), aft_Layout(xon)) || ...
   ((~ischar(sid) || isempty(sid) || numel(sid) ~= size(sid, 2)) && ...
    (~iscell(sid) || isempty(sid) || ~ischar(sid{1}) || isempty(sid{1})))
    error('neuroelf:xff:badArgument', 'Bad object or argument in call.');
end

% prepare array for subjects to keep
sl = glm_Subjects(xo);
ssl = glm_Subjects(xon);
spred = glm_SubjectPredictors(xo);
sspred = glm_SubjectPredictors(xon);
if ~isequal(sort(spred), sort(sspred))
    error('neuroelf:xff:badArgument', 'GLMs must match in predictors.');
end
predm = ne_methods.multimatch(spred, sspred);
if any(predm < 1 ) || numel(unique(predm)) ~= numel(predm)
    error('neuroelf:xff:badArgument', 'GLMs must have unique predictors.');
end
samepred = isequal(predm, (1:numel(predm))');

% ensure sid is in a cell
if ~iscell(sid)
    sid = {sid};
end
sid = sid(:);
if any(cellfun('isempty', sid) | ~cellfun(@ischar, sid))
    error('neuroelf:xff:badArgument', 'Bad object or argument in call.');
end
sid = unique(sid);

% make sure sid exists in both
msource = ne_methods.multimatch(sid, ssl);
mtarget = ne_methods.multimatch(sid, sl);
if any(msource == 0 | mtarget == 0)
    error('neuroelf:xff:badArgument', 'Subjects must be in both GLMs.');
end

% get data for temporary, non-overwriting access
bc = xo.C;
rtv = bc.RunTimeVars;
sbc = xon.C;
srtv = sbc.RunTimeVars;
fsl = glm_Subjects(xo, true);
fssl = glm_Subjects(xon, true);

% replace covariates?
repc = false;
if isfield(rtv, 'CovariatesData') && ~isempty(rtv.CovariatesData) && ...
    isfield(rtv, 'CovariatesNames') && ~isempty(rtv.CovariatesNames) && ...
    isfield(srtv, 'CovariatesData') && ~isempty(srtv.CovariatesData) && ...
    isfield(srtv, 'CovariatesNames') && ~isempty(srtv.CovariatesNames) && ...
    numel(rtv.CovariatesNames) == numel(srtv.CovariatesNames)
    covm = ne_methods.multimatch(rtv.CovariatesNames(:), srtv.CovariatesNames(:));
    if ~any(covm < 1) && numel(unique(covm)) == numel(covm)
        repc = true;
    end
end

% iterate over subjects
for sc = 1:numel(sid)

    % find studies
    sst = find(strcmp(fssl, sid{sc}));
    st = find(strcmp(fsl, sid{sc}));
    kst = 1:numel(bc.Study);
    kst(st) = [];

    % replace data
    bc.NrOfTimePoints = bc.NrOfTimePoints - sum(cat(1, bc.Study(st).NrOfTimePoints)) + ...
        sum(cat(1, sbc.Study(sst).NrOfTimePoints));
    bc.NrOfStudies = bc.NrOfStudies - numel(st) + numel(sst);
    bc.NrOfStudiesWithConfounds = bc.NrOfStudiesWithConfounds - numel(st) + numel(sst);
    bc.NrOfConfoundsPerStudy = [bc.NrOfConfoundsPerStudy(1, kst(1:st(1)-1)), ...
        sbc.NrOfConfoundsPerStudy(1, sst), bc.NrOfConfoundsPerStudy(1, kst(st(1):end))];
    bc.Study = [bc.Study(1, kst(1:st(1)-1)), sbc.Study(1, sst), bc.Study(1, kst(st(1):end))];
    if samepred
        bc.GLMData.Subject(mtarget(sc)).BetaMaps = ...
            sbc.GLMData.Subject(msource(sc)).BetaMaps;
    elseif bc.ProjectType < 2
        bc.GLMData.Subject(mtarget(sc)).BetaMaps = ...
            sbc.GLMData.Subject(msource(sc)).BetaMaps(:, :, :, predm);
    else
        bc.GLMData.Subject(mtarget(sc)).BetaMaps = ...
            sbc.GLMData.Subject(msource(sc)).BetaMaps(:, predm);
    end

    % covariates in new GLM
    if repc
        bc.RunTimeVars.CovariatesData(mtarget(sc), :) = ...
            srtv.CovariatesData(msource(sc), covm);
    end

    % motion parameters
    if isfield(rtv, 'MotionParameters') && iscell(rtv.MotionParameters) && ...
        numel(rtv.MotionParameters) == numel(bc.Study) && ...
        isfield(srtv, 'MotionParameters') && iscell(srtv.MotionParameters) && ...
        numel(srtv.MotionParameters) == numel(sbc.Study)
        bc.RunTimeVars.MotionParameters = ...
            [bc.RunTimeVars.MotionParameters(kst(1:st(1)-1)); ...
             sbc.RunTimeVars.MotionParameters(1, sst), ...
             bc.RunTimeVars.MotionParameters(1, kst(st(1):end))];
    end

    % per-subjects GLM
    if isfield(rtv, 'PerSubjectGLMs') && isstruct(rtv.PerSubjectGLMs) && ...
        numel(rtv.PerSubjectGLMs) == 1 && isfield(rtv.PerSubjectGLMs, sid{sc}) && ...
        isfield(srtv, 'PerSubjectGLMs') && isstruct(srtv.PerSubjectGLMs) && ...
        numel(srtv.PerSubjectGLMs) == 1 && isfield(srtv.PerSubjectGLMs, sid{sc})
        bc.RunTimeVars.PerSubjectGLMs.(sid{sc}) = srtv.PerSubjectGLMs.(sid{sc});
    end

    % make sure list is refreshed (in target!)
    fsl = [fsl(kst(1:st(1)-1)); fssl(sst); fsl(kst(st(1):end))];
end

% set content
xo.C = bc;
