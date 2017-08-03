% FUNCTION ne_statsmenu: populate Statistics menu and handle events
function ne_statsmenu(varargin)

% Version:  v1.1
% Build:    16031522
% Date:     Mar-15 2016, 10:26 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, 2016, Jochen Weber
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

% get handles for menu, voi, and voi index object
ch = ne_gcfg.h.Menu.Stats;
mh = ch.MLHandle;

% depending on page
if ne_gcfg.fcfg.page == 3
    stvar = ne_gcfg.fcfg.SurfStatsVar;
    stvarh = ne_gcfg.h.SurfStatsVar;
    stmaph = ne_gcfg.h.SurfStatsVarMaps;
    closearg = 'SurfStatsVar';
    setstatsfun = @ne_setcsrfstats;
    setmapfun = @ne_setcsrfstatmap;
else
    stvar = ne_gcfg.fcfg.StatsVar;
    stvarh = ne_gcfg.h.StatsVar;
    stmaph = ne_gcfg.h.StatsVarMaps;
    closearg = 'StatsVar';
    setstatsfun = @ne_setcstats;
    setmapfun = @ne_setcstatmap;
end
vidx = stvarh.Value;
midx = stmaph.Value;
midxn = numel(midx);

% for no further arguments
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3})

    % remove all children first
    delete(ch.Children);

    % only once at a time
    if ne_gcfg.c.incb
        return;
    end
    ne_gcfg.c.incb = true;

    % no or empty stvar
    if numel(stvar) ~= 1 || ...
       ~isxff(stvar, {'cmp', 'fsmf', 'glm', 'hdr', 'head', 'smp', 'vmp'})

        % simply inform user...
        uimenu(mh, 'Label', 'No valid statistics object loaded', 'Enable', 'off');
        ne_gcfg.c.incb = false;
        return;
    end
    stvartype = lower(stvar.Filetype);
    stvarmaps = stvar.MapNames(ne_gcfg.c.extmapnames);
    stvarfn = stvar.FilenameOnDisk(true);
    if ~isempty(stvarfn)
        [stvarfp, stvarfn, stvarfe] = fileparts(stvarfn);
        stvarfn = [stvarfn, stvarfe];
        if numel(stvarfn) > 34
            stvarfn = [stvarfn(1:18) '...' stvarfn(end-10:end)];
        end
    else
        stvarfn = sprintf('%s (xff #%d)', upper(stvartype), stvar.Filenumber);
    end

    % first add tools entries
    sep = 'on';
    uimenu(mh, 'Label', ['Close ' stvarfn], 'Callback', {@ne_closefile, closearg});
    if midxn == 1 && ...
        strcmp(stvartype, 'vmp')
        uimenu(mh, 'Label', 'Map properties', 'Callback', @ne_setcstatmapprops);
    end
    if ~isempty(stvarmaps) && ...
        strcmp(stvartype, 'vmp')
        uimenu(mh, 'Label', 'Compute formula', 'Callback', @ne_setcstatmapformula);
    end
    if midxn > 0 && ...
       ~strcmp(stvartype, 'smp') && ...
       ~isempty(ne_gcfg.h.Scenery.Value)
        uimenu(mh, 'Label', 'Create SMP from selected maps', ...
            'Callback', {@ne_vmp_createsmp, 'StatsVar', 'cursel'});
    end
    if midxn == 1 && ...
        any(strcmp(stvartype, {'fsmf', 'smp', 'vmp'}))
        uimenu(mh, 'Label', 'Clustertable', 'Separator', 'on', ...
            'Callback', {@ne_clustertable, stvartype});
    end

    % select statistics object menu
    stvars = stvarh.String;
    if ~iscell(stvars)
        stvars = cellstr(stvars);
    end
    if numel(stvars) > 1
        sch = uimenu(mh, 'Label', 'Switch to stats object', 'Separator', sep);
        for vc = 1:numel(stvars)
            mmh = uimenu(sch, 'Label', stvars{vc}, ...
                'Callback', {setstatsfun, vc});
            if vc == vidx
                set(mmh, 'Checked', 'on');
            end
            if mod(vc, 30) == 0 && ...
                numel(stvars) > vc
                sch = uimenu(mh, 'Label', 'Switch to stats (cont''d)');
            end
        end
        sep = 'on';
    end
    if numel(stvarmaps) > 0
        sch = uimenu(mh, 'Label', 'Select map', 'Separator', sep);
        for vc = 1:numel(stvarmaps)
            mmh = uimenu(sch, 'Label', stvarmaps{vc}, ...
                'Callback', {setmapfun, vc});
            if any(midx == vc)
                set(mmh, 'Checked', 'on')
            end
            if mod(vc, 30) == 0 && ...
                numel(stvarmaps) > vc
                sch = uimenu(mh, 'Label', 'Select map (cont''d');
            end
        end
        if numel(stvarmaps) > 1
            uimenu(mh, 'Label', 'Multi-select maps', ...
                'Callback', {@ne_statsmenu, 'multiselect'});
        end
        if strcmp(stvartype, 'glm') && ...
            stvar.ProjectTypeRFX > 0 && ...
            midxn == 1
            mapname = regexprep(stvarmaps{midx}, '^.*\:\s+', '');
            uimenu(mh, 'Label', ['Select ' mapname ' for all subjects'], ...
                'Callback', @ne_setcstatmapproj);
        end
        if strcmp(stvartype, 'vmp')
            uimenu(mh, 'Label', 'Delete selected maps', 'Separator' , 'on', ...
                'Callback', {@ne_movestatmap, 0});
        end
    end

    % remove only-once flag
    ne_gcfg.c.incb = false;
    return;
end

% no or empty stvar
if numel(stvar) ~= 1 || ...
   ~isxff(stvar, {'cmp', 'fsmf', 'glm', 'hdr', 'head', 'smp', 'vmp'})
    return;
end
stvarmaps = stvar.MapNames;

% only once at a time
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

% action
action = lower(varargin{3}(:)');
switch (action)

    % select a different cluster
    case {'multiselect'}

        % request new selection
        [newsel, selok] = listdlg( ...
            'PromptString',  'Please select the maps to visualize...', ...
            'SelectionMode', 'multiple', ...
            'ListString',    stvarmaps(:), ...
            'ListSize',      [520, 320], ...
            'InitialValue',  midx, ...
            'Name',          'NeuroElf - user input');
        ne_gcfg.c.incb = false;
        if ~isequal(selok, 1)
            return;
        end

        % set maps
        feval(setmapfun, 0, 0, newsel(:));
        drawnow;
end

% remove only-once flag
ne_gcfg.c.incb = false;
