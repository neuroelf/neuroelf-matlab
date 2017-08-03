% FUNCTION ne_voimenu: populate VOI menu and handle events
function ne_voimenu(varargin)

% Version:  v0.9c
% Build:    11050917
% Date:     Apr-29 2011, 8:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
ch = ne_gcfg.h.Menu.VOI;
mh = ch.MLHandle;
voi = ne_gcfg.voi;
vidxh = ne_gcfg.h.Clusters;
vidx = vidxh.Value;
vidxn = numel(vidx);

% for no further arguments
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3})

    % remove all children first
    delete(ch.Children);

    % currently only for VOI (not POI)
    if ne_gcfg.fcfg.page == 3
        return;
    end

    % only once at a time
    if ne_gcfg.c.incb
        return;
    end
    ne_gcfg.c.incb = true;

    % no or empty VOI
    if numel(voi) ~= 1 || ...
       ~isxff(voi, 'voi') || ...
        isempty(voi.VOI)

        % simply inform user...
        uimenu(mh, 'Label', 'No clusters defined', 'Enable', 'off');
        ne_gcfg.c.incb = false;
        return;
    end
    numvois = numel(voi.VOI);

    % first add tools entries
    sep = 'off';
    if vidxn > 1
        uimenu(mh, 'Label', 'Combine selected clusters (OR)', ...
            'Callback', {@ne_voimenu, 'combine', 'or'});
    end
    if numvois > 0
        uimenu(mh, 'Label', 'Go to nearest cluster (peak)', ...
            'Callback', {@ne_voimenu, 'select', 'nearest'});
        sep = 'on';
    end

    % select cluster sub-menu
    sch = uimenu(mh, 'Label', 'Select cluster...', 'Separator', sep);
    vnames = voi.VOINames;
    for vc = 1:numel(vnames)
        uimenu(sch, 'Label', vnames{vc}, ...
            'Callback', {@ne_voimenu, 'select', 'set', vc});
        if mod(vc, 20) == 0 && ...
            numel(vnames) > vc
            sch = uimenu(mh, 'Label', 'select cluster (cont''d)');
        end
    end

    % remove only-once flag
    ne_gcfg.c.incb = false;
    return;
end

% only once at a time
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

% action
action = lower(varargin{3}(:)');
switch (action)

    % select a different cluster
    case {'select'}

        % pass on
        try
            ne_setcluster(0, 0, varargin{4:end});
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

    % combine clusters
    case {'combine'}

        % sub-action
        if nargin < 4 || ...
           ~ischar(varargin{4}) || ...
           ~any(strcmpi(varargin{4}(:)', {'and', 'avg', 'or', 'xor'}))
            ne_gcfg.c.incb = false;
            return;
        end
end

% remove only-once flag
ne_gcfg.c.incb = false;
