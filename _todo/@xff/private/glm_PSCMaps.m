function xo2 = glm_PSCMaps(xo)
% GLM::PSCMaps  - calculate PSC maps for single study GLM
%
% FORMAT:       pscmap = glm.PSCMaps;
%
% Output fields:
%
%       pscmap      VMP object
%
% Using: newnatresvmp.

% Version:  v1.1
% Build:    16020312
% Date:     Feb-03 2016, 12:00 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if bc.NrOfStudies ~= 1 || bc.ProjectType ~= 1 || bc.TransformationType == 1
    error('neuroelf:xff:badObject', 'Invalid call to %s.', mfilename);
end
bbox = aft_BoundingBox(xo);

% create VMP
xo2 = ne_methods.newnatresvmp(bbox.BBox, bc.Resolution, 11 * ones(1, bc.NrOfPredictors - 1));
bc2 = xo2.C;

% update TrfPlus
bc2.RunTimeVars.TrfPlus = bc.RunTimeVars.TrfPlus;

% calculate maps
for pc = 1:(bc.NrOfPredictors - 1)
    bc2.Map(pc).LowerThreshold = 0.5;
    bc2.Map(pc).UpperThreshold = 5;
    bc2.Map(pc).DF1 = bc.NrOfTimePoints - bc.NrOfPredictors;
    bc2.Map(pc).DF2 = 0;
    bc2.Map(pc).BonferroniValue = bc.NrOfVoxelsForBonfCorrection;
    bc2.Map(pc).Name = sprintf('PSC map for %s', bc.Predictor(pc).Name2);
    mapv = 100 * bc.GLMData.BetaMaps(:, :, :, pc) ./ bc.GLMData.BetaMaps(:, :, :, end);
    mapv(isinf(mapv) | isnan(mapv)) = 0;
    bc2.Map(pc).VMPData = mapv;
end

% set back
xo2.C = bc2;
