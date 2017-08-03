% FUNCTION ne_maskstatswithcls: mask the current statsvar with selected VOIs
function ne_maskstatswithcls(varargin)

% Version:  v0.9b
% Build:    11050921
% Date:     Apr-09 2011, 11:48 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg;
cv = ne_gcfg.h.Clusters.Value;
mv = ne_gcfg.h.StatsVarMaps.Value;

% only valid for certain combinations
if numel(cc.StatsVar) ~= 1 || ...
   ~isxff(cc.StatsVar, {'glm', 'vmp'}) || ...
   (isxff(cc.StatsVar, 'glm') && ...
    cc.StatsVar.ProjectType ~= 1) || ...
    numel(ne_gcfg.voi) ~= 1 || ...
   ~isxff(ne_gcfg.voi, 'voi') || ...
    isempty(ne_gcfg.voi.VOI) || ...
    isempty(cv) || ...
    isempty(mv)
    return;
end

% block further calls
ne_gcfg.c.incb = true;
ne_gcfg.h.MainFig.Pointer = 'watch';

% get VOI names and number of maps
if numel(cv) == 1
    vn = ne_gcfg.voi.VOINames;
    vn = deblank(strrep(vn{cv}, '_', ' '));
else
    vn = sprintf('%d clusters', numel(cv));
end
nummaps = numel(cc.StatsVar.Map);

% try the rest
obj = {[], []};
try

    % create MSK object (of selected clusters)
    obj{1} = ne_gcfg.voi.CreateMSK(cc.StatsVar, cv);

    % dilate mask by one voxel
    obj{1}.Mask = uint8(dilate3d(obj{1}.Mask > 0));

    % create mask version (as a copy)
    obj{2} = obj{1}.ApplyTo(cc.StatsVar, true);

    % copy selected maps at end of original
    cc.StatsVar.Map = [cc.StatsVar.Map(:); lsqueeze(obj{2}.Map(mv))];
    cc.StatsVar.NrOfMaps = numel(cc.StatsVar.Map);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% clear objects
clearxffobjects(obj);

% update map names
for mc = 1:numel(mv)
    cc.StatsVar.Map(mc + nummaps).Name = ...
        sprintf('%s (masked with %s)', cc.StatsVar.Map(mv(mc)).Name, vn);
end

% set new selection
m = cc.StatsVar.Map;
cc.StatsVar.RunTimeVars.MapSelection = ...
    {{m((1:numel(mv)) + nummaps).Name}, (1:numel(mv)) + nummaps};

% re-allow calls
ne_gcfg.c.incb = false;
ne_gcfg.h.MainFig.Pointer = 'arrow';

% update display
ne_openfile(0, 0, cc.StatsVar);
