function xo3 = glm_JoinRFX(xo, xo2)
% GLM::JoinRFX  - join two random effects GLMs
%
% FORMAT:       combined = glm1.JoinRFX(glm2);
%
% Input fields:
%
%       glm2        second RFX GLM to be added to glm1
%
% Joins the Random Effects GLM results of glm1 and glm2 to one
% single structure.
%
% Using: catstruct.

% Version:  v1.1
% Build:    16020312
% Date:     Feb-03 2016, 12:01 PM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || numel(xo2) ~= 1 || ...
   ~xffisobject(xo, true, 'glm') || ~xffisobject(xo2, true, 'glm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc1 = xo.C;
bc2 = xo2.C;
if bc1.ProjectTypeRFX ~= 1 || bc2.ProjectTypeRFX ~= 1 || bc1.SeparatePredictors < 1
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
% check type and sizes of GLMs
if ~isfield(bc1.GLMData, 'RFXGlobalMap') || ~isfield(bc2.GLMData, 'RFXGlobalMap') || ...
    bc1.ProjectType ~= bc2.ProjectType || bc1.FileVersion ~= bc2.FileVersion || ...
    bc1.NrOfSubjectPredictors ~= bc2.NrOfSubjectPredictors || ...
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

% set fields that change
bc3.NrOfSubjects = NSB;
bc3.NrOfTimePoints = NTP;
bc3.NrOfPredictors = NPS;
bc3.NrOfConfounds = bc3.NrOfConfounds + bc2.NrOfConfounds;
bc3.NrOfStudies = NSS;
bc3.NrOfStudiesWithConfounds = bc3.NrOfStudiesWithConfounds + bc2.NrOfStudiesWithConfounds;
bc3.NrOfConfoundsPerStudy = [NC1(:); NC2(:)]';
bc3.NrOfVoxelsForBonfCorrection = ...
    max(bc3.NrOfVoxelsForBonfCorrection, bc2.NrOfVoxelsForBonfCorrection);
bc3.Study = ne_methods.catstruct(bc3.Study(:), bc2.Study(:))';
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

% put back
xo3.C = bc3;
