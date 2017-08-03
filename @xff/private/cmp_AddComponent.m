function xo = cmp_AddComponent(xo, cmap)
% CMP::AddComponent  - add a component map to a CMP/ICA file
%
% FORMAT:       cmp.AddComponent([cmap]);
%
% Input fields:
%
%       cmap        1x1 struct with optional fields (as in Map)
%        .Type               component type (12, ...)
%        .LowerThreshhold    lower threshold (default: 0)
%        .UpperThreshhold    upper threshold (default: 1)
%        .Name               component name (default: 'component %d')
%        .RGBLowerThreshPos  RGB color code for positive lower thresh
%        .RGBUpperThreshPos  RGB color code for positive upper thresh
%        .RGBLowerThreshNeg  RGB color code for negative lower thresh
%        .RGBUpperThreshNeg  RGB color code for negative upper thresh
%        .UseRGBColor        whether or not use CMP color (default: 1)
%        .TransColorFactor   transparent coloring factor (default: 1)
%        .TimePointData      timecourse (default: zeros(NrOfTimePoints,1))
%        .CMPData            XxYxZ component map (dims must match)
%                   additionally the parameter values can be given in
%        .ParamValues        values will go into .MapParameters
%
% No output fields.
%
% Note: if no input argument is given, the dimensions will be taken
%       from the main structure

% Version:  v1.1
% Build:    16021413
% Date:     Feb-14 2016, 1:22 PM EST
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

% only valid for single file
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'cmp')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get current X, Y, Z dims
bc = xo.C;
pres = bc.Resolution;
dimx = (bc.XEnd - bc.XStart) / pres;
dimy = (bc.YEnd - bc.YStart) / pres;
dimz = (bc.ZEnd - bc.ZStart) / pres;

% no component struct
if nargin < 2 || ~isstruct(cmap) || isempty(cmap)
    cmap = struct('CMPData', zeros([dimx, dimy, dimz]));

% no or non numeric CMPData
elseif ~isfield(cmap, 'CMPData') || ~isnumeric(cmap.CMPData)
    cmap.CMPData = zeros([dimx, dimy, dimz]);

% check CMPData size
elseif numel(size(cmap.CMPData)) ~= 3 || any(size(cmap.CMPData) ~= [dimx, dimy, dimz])
    error('neuroelf:xff:badArgument', 'Bad dimension of component Map.');
end

% is this the first component?
nc = numel(bc.Map);
if nc < 1
    defvals = struct('Type', 1, 'LowerThreshold', 2.5, 'UpperThreshold', 8, ...
        'Name', 'component 1', 'RGBLowerThreshPos', [255, 0, 0], 'RGBUpperThreshPos', [255, 192, 0], ...
        'RGBLowerThreshNeg', [0, 0, 128], 'RGBUpperThreshNeg', [0, 128, 255], ...
        'UseRGBColor', 1, 'TransColorFactor', 1, 'TimePointData', zeros(bc.NrOfTimePoints, 1), ...
        'ParamValues', zeros(1, bc.NrOfMapParameters));

% adding another component
else
    defvals = bc.Map(1);
    defvals.Name = 'component %d';
    defvals.TimePointData = zeros(bc.NrOfTimePoints, 1);
    defvals.ParamValues = zeros(1, bc.NrOfMapParameters);
end

% add missing fields to struct
mfields = fieldnames(defvals);
for fc = 1:numel(mfields)
    if ~isfield(cmap, mfields{fc})
        for cc = 1:numel(cmap)
            cmap(cc).(mfields{fc}) = defvals.(mfields{fc});
        end
    end
end

% add component(s)
for cc = 1:numel(cmap)

    % increase counter
    nc = nc + 1;

    % make settings
    bc.Map(nc).Type = cmap(cc).Type;
    bc.Map(nc).LowerThreshold = cmap(cc).LowerThreshold;
    bc.Map(nc).UpperThreshold = cmap(cc).UpperThreshold;
    if sum(cmap(cc).Name == '%') ~= 1
        bc.Map(nc).Name = cmap(cc).Name;
    else
        bc.Map(nc).Name = sprintf(cmap(cc).Name, nc);
    end
    bc.Map(nc).RGBLowerThreshPos = cmap(cc).RGBLowerThreshPos;
    bc.Map(nc).RGBUpperThreshPos = cmap(cc).RGBUpperThreshPos;
    bc.Map(nc).RGBLowerThreshNeg = cmap(cc).RGBLowerThreshNeg;
    bc.Map(nc).RGBUpperThreshNeg = cmap(cc).RGBUpperThreshNeg;
    bc.Map(nc).UseRGBColor = cmap(cc).UseRGBColor;
    bc.Map(nc).TransColorFactor = cmap(cc).TransColorFactor;
    bc.Map(nc).TimePointData = cmap(cc).TimePointData;

    % set CMPData
    bc.Map(nc).CMPData = cmap.CMPData;

    % set ParameterValue
    for pc = 1:numel(cmap(cc).ParamValues)
        bc.MapParameter(pc).Values(nc) = cmap.ParamValues(pc);
    end
end

% set correct number
bc.NrOfMaps = nc;
bc.NrOfMapParameters = numel(bc.MapParameter);

% set into storage
xo.C = bc;
