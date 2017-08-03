function xo = vtc_JoinAvgTCs(xo, xo2)
% VTC::JoinAvgTCs  - join two AvgVTC objects
%
% FORMAT:       [vtc = ] vtc.JoinAvgTCs(vtc2)
%
% Input fields:
%
%       vtc2        second AvgVTC object (must match in sampling parameters)
%
% Output fields:
%
%       vtc         (altered) VTC with joined (combined) time course data

% Version:  v1.1
% Build:    16021320
% Date:     Feb-13 2016, 8:56 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2013, 2014, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc') || ...
    numel(xo2) ~= 1 || ~xffisobject(xo2, true)
    error('neuroelf:xff:badArgument', 'Invalid object in call.');
end
bc = xo.C;
bc2 = xo2.C;
if ~isfield(bc.RunTimeVars, 'AvgVTC') || ~islogical(bc.RunTimeVars.AvgVTC) || ...
    numel(bc.RunTimeVars.AvgVTC) ~= 1 || ~bc.RunTimeVars.AvgVTC || ...
   ~isfield(bc2.RunTimeVars, 'AvgVTC') || ~islogical(bc2.RunTimeVars.AvgVTC) || ...
    numel(bc2.RunTimeVars.AvgVTC) ~= 1 || ~bc2.RunTimeVars.AvgVTC
    error('neuroelf:xff:badArgument', 'Only valid for Average-VTCs.');
end
if bc.Resolution ~= bc2.Resolution || bc.XStart ~= bc2.XStart || bc.XEnd ~= bc2.XEnd || ...
    bc.YStart ~= bc2.YStart || bc.YEnd ~= bc2.YEnd || bc.ZStart ~= bc2.ZStart || ...
    bc.ZEnd ~= bc2.ZEnd || bc.TR ~= bc2.TR || ...
    bc.RunTimeVars.AvgTransformationType ~= bc2.RunTimeVars.AvgTransformationType || ...
    bc.RunTimeVars.AvgWindowFrom ~= bc2.RunTimeVars.AvgWindowFrom || ...
    bc.RunTimeVars.AvgWindowStep ~= bc2.RunTimeVars.AvgWindowStep || ...
    bc.RunTimeVars.AvgWindowTo ~= bc2.RunTimeVars.AvgWindowTo || ...
    bc.RunTimeVars.BaseWindowFrom ~= bc2.RunTimeVars.BaseWindowFrom || ...
    bc.RunTimeVars.BaseWindowStep ~= bc2.RunTimeVars.BaseWindowStep || ...
    bc.RunTimeVars.BaseWindowTo ~= bc2.RunTimeVars.BaseWindowTo || ...
    bc.RunTimeVars.NrOfTCsPerCondition ~= bc2.RunTimeVars.NrOfTCsPerCondition || ...
    bc.RunTimeVars.NrOfVolumesPerTC ~= bc2.RunTimeVars.NrOfVolumesPerTC || ...
   ~all(strcmpi(bc.RunTimeVars.TCNames(:), bc2.RunTimeVars.TCNames(:)))
    error('neuroelf:xff:badArgument', 'VTCs must match in settings and spatial layout.');
end

% set and extend data
if ~isempty(bc.NameOfSourceFMR) && ~isempty(bc2.NameOfSourceFMR)
    bc.NameOfSourceFMR = [bc.NameOfSourceFMR '+' bc2.NameOfSourceFMR];
end
bc.NrOfVolumes = bc.NrOfVolumes + bc2.NrOfVolumes;
vd = single(0);
vd(bc.NrOfVolumes, size(bc.VTCData, 2), size(bc.VTCData, 3), size(bc.VTCData, 4)) = 0;
for sc = 1:size(bc.VTCData, 4)
    vd(:, :, :, sc) = cat(1, bc.VTCData(:, :, :, sc), bc2.VTCData(:, :, :, sc));
end
bc.VTCData = vd;
bc.RunTimeVars.Map = [bc.RunTimeVars.Map(:)', bc2.RunTimeVars.Map(:)'];
bc.RunTimeVars.NrOfConditions = bc.RunTimeVars.NrOfConditions + bc2.RunTimeVars.NrOfConditions;
bc.RunTimeVars.NrOfConditionOnsets = [bc.RunTimeVars.NrOfConditionOnsets, bc2.RunTimeVars.NrOfConditionOnsets];
bc.RunTimeVars.NrOfSourceVTCs = bc.RunTimeVars.NrOfSourceVTCs + bc2.RunTimeVars.NrOfSourceVTCs;
bc.RunTimeVars.NrOfSubjects = bc.RunTimeVars.NrOfSubjects + bc2.RunTimeVars.NrOfSubjects;
bc.RunTimeVars.ConditionColors = [bc.RunTimeVars.ConditionColors; bc2.RunTimeVars.ConditionColors];
bc.RunTimeVars.ConditionNames = [bc.RunTimeVars.ConditionNames(:); bc2.RunTimeVars.ConditionNames(:)];
co1 = bc.RunTimeVars.ConditionOnsets;
co2 = bc2.RunTimeVars.ConditionOnsets;
co = repmat({zeros(0, 2)}, size(co1) + size(co2));
co(1:size(co1, 1), 1:size(co1, 2)) = co1;
co(end+1-size(co2, 1):end, end+1-size(co2, 2):end) = co2;
bc.RunTimeVars.ConditionOnsets = co;
bc.RunTimeVars.ConditionThresholds = cat(1, ...
    bc.RunTimeVars.ConditionThresholds, bc2.RunTimeVars.ConditionThresholds);
bc.RunTimeVars.SubjectNames = [bc.RunTimeVars.SubjectNames(:); bc2.RunTimeVars.SubjectNames(:)];
bc.RunTimeVars.SourcePRTs = [bc.RunTimeVars.SourcePRTs(:); bc2.RunTimeVars.SourcePRTs(:)];
bc.RunTimeVars.SourceVTCs = [bc.RunTimeVars.SourceVTCs(:); bc2.RunTimeVars.SourceVTCs(:)];
bc.RunTimeVars.TCOnsetWeights = [bc.RunTimeVars.TCOnsetWeights(:); bc2.RunTimeVars.TCOnsetWeights(:)];

% set data
xo.C = bc;
