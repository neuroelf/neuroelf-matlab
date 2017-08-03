% FUNCTION ne_glmplotbetaskyp: handle keyboard (if no control has focus)
function ne_glmplotbetaskyp(src, ke, varargin)

% Version:  v0.9d
% Build:    14062115
% Date:     Jun-21 2014, 3:57 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

% only valid for valid Tag string
if nargin < 3 || ...
   ~ischar(varargin{1}) || ...
    numel(varargin{1}) ~= 8 || ...
   ~isfield(ne_gcfg.cc, varargin{1}(:)')
    return;
end
tsat = varargin{1}(:)';

% get Key and Modifier from keyboard event (see Matlab docu!)
kk = ke.Key;
mn = ke.Modifier;
ne_gcfg.cc.(tsat).Config.mods = mn;

% get configuration and handles
ch = ne_gcfg.cc.(tsat);
hFig = ch.Satellite;
ch = ch.Tags;

% determine which modifiers are pressed
km = false(1, 4);
if ~isempty(mn)
    try
        km = [ ...
            any(strcmpi('alt', mn)), ...
            any(strcmpi('control', mn)), ...
            any(strcmpi('shift', mn)), ...
            any(strcmpi('command', mn))];
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end

% handle events without modifiers
if ~any(km)
    switch (lower(kk))

        % save axes screenshot
        case {'a'}
            ne_screenshot(0, 0, ch.Axes.MLHandle);

        % background color
        case {'b'}
            ne_glmplotbetasgui(0, 0, tsat, 'BackColor');

        % set condition/contrast colors
        case {'c'}
            if isempty(ch.Conditions.Value)
                ne_glmplotbetasgui(0, 0, tsat, 'ContColors');
            else
                ne_glmplotbetasgui(0, 0, tsat, 'CondColors');
            end

        % copy data
        case {'d'}
            ne_glmplotbetasgui(0, 0, 'CopyData');

        % dis/enable error bars
        case {'e'}
            ch.EBars.Value = 1 - double(ch.EBars.Value > 0);
            ne_glmplotbetasgui(0, 0, tsat, 'EBars');

        % save figure screenshot
        case {'f'}
            ne_screenshot(0, 0, hFig.MLHandle, '', 'high-q');

        % dis/enable grouping
        case {'g'}
            if strcmpi(ch.DoGroups.Enable, 'on')
                ch.DoGroups.Value = 1 - double(ch.DoGroups.Value > 0);
                ne_glmplotbetasgui(0, 0, tsat, 'DoGroups');
            end

        % individual subject indicators
        case {'i'}
            ch.PPts.Value = 1 - double(ch.PPts.Value > 0);
            ne_glmplotbetasgui(0, 0, tsat, 'PPts');

        % subject profile lines
        case {'l'}
            ne_glmplotbetasgui(0, 0, tsat, 'SubLines');

        % toggle options visible
        case {'o'}
            ch.OptVis.Value = 1 - double(ch.OptVis.Value > 0);
            ne_glmplotbetasrsz(0, 0, tsat, 'OptVis');

        % copy plot
        case {'p'}
            ne_glmplotbetasgui(0, 0, tsat, 'CopyPlot');

        % toggle robust stats
        case {'r'}
            ne_glmplotbetasgui(0, 0, tsat, 'Robust');

        % subject labels
        case {'s'}
            ne_glmplotbetasgui(0, 0, tsat, 'SubLabels');

        % plot type
        case {'t'}
            if strcmpi(ch.Covariates.Enable, 'on')
                ch.Type.Value = 3 - ch.Type.Value;
                ne_glmplotbetasgui(0, 0, tsat, 'Type');
            end

        % toggle update flag
        case {'u'}
            ch.DoUpdate.Value = 1 - double(ch.DoUpdate.Value > 0);
            ne_glmplotbetasgui(0, 0, tsat, 'DoUpdate');

        % remove 0 values
        case {'0'}
            ne_glmplotbetasgui(0, 0, tsat, 'Remove0s');
    end

% handle events with ALT
elseif km(1)

% handle events with CTRL
elseif km(2)
    switch (lower(kk))

        % marker type
        case {'a', 'c', 'd', 'o', 'p', 's', 'x'}
            if strcmpi(get(ch.ScMarkAsterisk.Parent, 'Enable'), 'on')
                kk = lower(kk);
                if kk == 'a'
                    kk = '*';
                elseif kk == 'c'
                    kk = '.';
                elseif kk == 'p'
                    kk = '+';
                end
                ne_glmplotbetasgui(0, 0, tsat, 'ScatterMarker', kk);
            end

        % confidence ellipse
        case {'e'}
            if strcmpi(ch.ScConfEllipse.Enable, 'on')
                ne_glmplotbetasgui(0, 0, tsat, 'ScatterConfEllipse');
            end

        % filled markers
        case {'f'}
            if strcmpi(ch.ScFilledMarkers.Enable, 'on');
                ne_glmplotbetasgui(0, 0, tsat, 'ScatterFilledMarkers');
            end

        % scatter groups separately
        case {'g'}
            ne_glmplotbetasgui(0, 0, tsat, 'ScatterGroups');

        % best linear fit
        case {'l'}
            if strcmpi(ch.ScLine.Enable, 'on')
                ne_glmplotbetasgui(0, 0, tsat, 'ScatterLine');
            end

        % best quadratic fit
        case {'q'}
            if strcmpi(ch.ScQuadratic.Enable, 'on')
                ne_glmplotbetasgui(0, 0, tsat, 'ScatterQuadratic');
            end

        % rank transform
        case {'r'}
            if strcmpi(ch.RankTrans.Enable, 'on')
                ne_glmplotbetasgui(0, 0, tsat, 'RankTrans');
            end

        % show stats output
        case {'t'}
            if strcmpi(ch.StatsOut.Enable, 'on')
                ne_glmplotbetasgui(0, 0, tsat, 'StatsOut');
            end

        % mark one group (up to 8) or all groups (0)
        case {'0', '1', '2', '3', '4', '5', '6', '7', '8'}
            if ch.DoGroups.Value > 0 && ...
                strcmpi(ch.DoGroups.Enable, 'on')
                if kk(1) == '0'
                    ag = ch.Groups.String;
                    if ~iscell(ag)
                        ag = cellstr(ag);
                    end
                    ch.Groups.Value = (1:numel(ag))';
                else
                    ch.Groups.Value = str2double(kk);
                end
                ne_glmplotbetasgui(0, 0, tsat, 'Groups');
            end
    end

% handle events with SHIFT
elseif km(3)
    switch (lower(kk))

        % full axes screen shot
        case {'a', 'f'}
            ne_glmplotbetasgui(0, 0, tsat, 'SaveAxes');

        % error bar colors
        case {'e'}
            ne_glmplotbetasgui(0, 0, tsat, 'EBarColors');

        % legend
        case {'l'}
            if ch.Type.Value == 1
                ne_glmplotbetasgui(0, 0, tsat, 'LegendBars');
            else
                ne_glmplotbetasgui(0, 0, tsat, 'LegendScat');
            end

        % subject selection
        case {'s'}
            ne_glmplotbetasgui(0, 0, tsat, 'SubSel');

        % close window
        case {'q'}
            ne_closesatwindow(0, 0, tsat);

        % legend position
        case {'1'}
            ne_glmplotbetasgui(0, 0, tsat, 'LegendPos', 'NorthWest');
        case {'2'}
            ne_glmplotbetasgui(0, 0, tsat, 'LegendPos', 'NorthEast');
        case {'3'}
            ne_glmplotbetasgui(0, 0, tsat, 'LegendPos', 'SouthEast');
        case {'4'}
            ne_glmplotbetasgui(0, 0, tsat, 'LegendPos', 'SouthWest');
    end
end
