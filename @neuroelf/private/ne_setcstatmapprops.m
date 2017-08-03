% FUNCTION ne_setcstatmapprops: set map properties of current maps
function varargout = ne_setcstatmapprops(varargin)

% Version:  v1.1
% Build:    16031523
% Date:     Mar-15 2016, 11:53 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2016, Jochen Weber
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
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% only valid if a single VMP is selected
if nargin > 2 && ...
    ischar(varargin{3}) && ...
   ~isempty(varargin{3}) && ...
    isfield(cc, varargin{3}(:)')
    stvar = cc.(varargin{3}(:)');
else
    stvar = cc.StatsVar;
end
if ~isxff(stvar, {'cmp', 'fsmf', 'smp', 'vmp'})
    return;
end

% index given
if nargin > 2 && ...
    isa(varargin{3}, 'double') && ...
    numel(varargin{3}) == 1 && ...
    any(1:numel(stvar.Map) == varargin{3})

    % take index
    stvix = varargin{3};

% otherwise
else

    % get from control
    if isxff(stvar, {'fsmf', 'smp'})
        surfs = true;
        stvix = ch.SurfStatsVarMaps.Value;
    else
        surfs = false;
        stvix = ch.StatsVarMaps.Value;
    end
end

% only allow singular index
if numel(stvix) ~= 1
    return;
end

% compile settings
map = stvar.Map(stvix);
maptypes = {'t', 'r', 'CC', 'F', '', '', '', '', 'm', '', ...
    '%', 'z_ica', 'TH', '', 'beta', 'prob', '', '', '', ...
    'MD', 'FA'};
try
    maptype = maptypes{map.Type};
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    maptype = '';
end
if isempty(maptype)
    maptype = sprintf('<unknown %d>', map.Type);
end
mapset = {map.Name, ['  ' maptype], sprintf('  %d', map.DF1), ...
    sprintf('  %d', map.DF2), sprintf('  %d  %d  %d', map.NrOfLags, ...
    map.MinLag, map.MaxLag), sprintf('  %d', map.BonferroniValue)};

% request updated settings
mapset = inputdlg({'Map name:', 'Map type (either of t, r, CC, or F):', ...
    'DF1:', 'DF2:', 'CC lag information:', ...
    'Value (number of independent tests) for Bonferroni correction:'}, ...
    'NeuroElf GUI - Map settings', 1, mapset);
if ~iscell(mapset) || ...
    numel(mapset) ~= 6
    return;
end

% check and update return values (separately)
if ~isempty(mapset{1}) && ...
    ischar(mapset{1})
    map.Name = mapset{1};

    % replace in listbox
    if surfs
        mapnames = ch.SurfStatsVarMaps.String;
    else
        mapnames = ch.StatsVarMaps.String;
    end
    if ~iscell(mapnames)
        mapnames = cellstr(mapnames);
    end
    mapnames{stvix} = mapset{1};
    if surfs
        ch.SurfStatsVarMaps.String = mapnames;
    else
        ch.StatsVarMaps.String = mapnames;
    end
end
if ~isempty(mapset{2}) && ...
    any(strcmpi(maptypes, mapset{2}))
    map.Type = findfirst(strcmpi(maptypes, mapset{2}));
end
if ~isempty(mapset{3}) && ...
   ~isempty(regexpi(mapset{3}, '^\d+$')) && ...
    mapset{3}(1) > '0'
    map.DF1 = str2double(mapset{3});
end
if ~isempty(mapset{4}) && ...
   ~isempty(regexpi(mapset{4}, '^\d+$')) && ...
    mapset{4}(1) > '0'
    map.DF2 = str2double(mapset{4});
end
if surfs
    ne_gcfg.fcfg.SurfStatsVarPar = ...
        {maptypes{map.Type}, map.DF1, map.DF2};
else
    ne_gcfg.fcfg.StatsVarPar = ...
        {maptypes{map.Type}, map.DF1, map.DF2};
end

% put back into map object
stvar.Map(stvix) = map;
if ne_gcfg.c.extmapnames
    if surfs
        ch.SurfStatsVarMaps.String = stvar.MapNames(true);
    else
        ch.StatsVarMaps.String = stvar.MapNames(true);
    end
end

% update viewer
if ~surfs
    ne_setcstatmap;
end
