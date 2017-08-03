% FUNCTION ne_btmove: update window with pressed button
function ne_btmove(varargin)

% Version:  v1.1
% Build:    16050816
% Date:     May-08 2016, 4:59 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% only once
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

% hide ToolTip (menu)
ne_gcfg.h.ToolTip.Visible = 'off';

% if the mouse button isn't currently pressed (while in this figure)
if ~any(ne_gcfg.c.btdown == gcbf)

    % width enough for tooltip menu
    wpos = get(ne_gcfg.h.MainFigMLH, 'Position');
    if wpos(3) > 920

        % get relative mouse cursor position and children positions
        cpos = get(ne_gcfg.h.MainFigMLH, 'CurrentPoint');
        chpos = ne_gcfg.fcfg.chpos{3};

        % control(s) we're hovering over
        hpos = find(chpos(:, 1) <= cpos(1) & chpos(:, 2) <= cpos(2) & ...
            chpos(:, 3) >= cpos(1) & chpos(:, 4) >= cpos(2));

        % remove invisible ones
        if ~isempty(hpos)
            if isa(ne_gcfg.fcfg.chpos{2}, 'double')
                hpos(~ishandle(ne_gcfg.fcfg.chpos{2}(hpos))) = [];
            else
                hpos(~isvalid(ne_gcfg.fcfg.chpos{2}(hpos))) = [];
            end
            if ~isempty(hpos)
                hpos(~strcmpi(get(ne_gcfg.fcfg.chpos{2}(hpos), 'Visible'), 'on')) = [];
            end
        end

        % still a control we're hovering over?
        if ~isempty(hpos)

            % get component xfigure's ToolTip
            try
                xh = xfigure(ne_gcfg.fcfg.chpos{1}(hpos(1)));
                tt = xh.LoadProps.ToolTip;
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
                tt = '';
            end

            % set
            if ~isempty(tt)
                ne_gcfg.h.ToolTip.Label = tt;
                ne_gcfg.h.ToolTip.Visible = 'on';
            end
        end
    end

    % then don't do anything (for now)
    ne_gcfg.c.incb = false;
    return;
end

% main figure
if ne_gcfg.c.btdown == ne_gcfg.h.MainFigMLH

    % do nothing if same position
    last = get(ne_gcfg.h.MainFigMLH, 'CurrentPoint');
    if isequal(last, ne_gcfg.fcfg.mpos.last)
        ne_gcfg.c.incb = false;
        return;
    end
    ne_gcfg.fcfg.mpos.last = last;

    % figure configuration
    cc = ne_gcfg.fcfg;

    % slicing pages
    if cc.page < 3

        % no modifiers
        if isempty(cc.mpos.mods) || ...
            any(strcmpi(cc.mpos.mods, 'command'))

            % update the window regularly
            fPos = [cc.slicepos; cc.zslicepos(1, :); cc.tcpos; cc.histpos];
            oPos = cc.mpos.down;
            if cc.page == 1 && ...
                any(all(oPos([1, 1, 1], :) >= (fPos(1:3, 1:2) - 1), 2) & all(oPos([1, 1, 1], :) <= (fPos(1:3, 3:4) + 1), 2))
                ne_setslicepos(0, 0, [], 'OnMouse');
            elseif cc.page == 2 && ...
                all(oPos >= (fPos(4, 1:2) - 1)) && ...
                all(oPos <= (fPos(4, 3:4) + 1))
                ne_setslicepos(0, 0, [], 'OnMouse');
            elseif strcmpi(get(ne_gcfg.h.TCPlot.MLHandle, 'Visible'), 'on') && ...
                all(oPos >= (fPos(5, 1:2) - 1)) && ...
                all(oPos <= (fPos(5, 3:4) + 1))
                ne_setslicepos;
            elseif all(oPos >= (fPos(6, 1:2) - 1)) && ...
                all(oPos <= (fPos(6, 3:4) + 1))
                ne_setslicepos;
            end

        % otherwise
        else

            % update global rotation, etc.
            ne_setglobaltrf;

            % update with current coordinate
            ne_setslicepos(0, 0, cc.cpos);
        end

    % scenery view
    elseif ne_gcfg.fcfg.page == 3

        % update surface view
        ne_setsurfpos;

    % render view
    elseif ne_gcfg.fcfg.page == 4

        % update render view
        ne_render_setview(0, 0, [], 'preview');
    end

% satellite figure
else

    % test available figures
    sats = fieldnames(ne_gcfg.cc);
    try
        last = get(ne_gcfg.c.btdown, 'CurrentPoint');
        for sc = 1:numel(sats)
            if ne_gcfg.c.btdown == ne_gcfg.cc.(sats{sc}).SatelliteMLH
                if isequal(last, ne_gcfg.cc.(sats{sc}).Config.mpos.last)
                    ne_gcfg.c.incb = false;
                    return;
                end
                ne_gcfg.cc.(sats{sc}).Config.mpos.last = last;

                switch (ne_gcfg.cc.(sats{sc}).Config.sattype)
                    case {'slice'}
                        ne_setsatslicepos(0, 0, ne_gcfg.cc.(sats{sc}).Config.sattag);
                        break;
                    case {'surf'}
                        ne_setsurfpos(0, 0, ne_gcfg.cc.(sats{sc}).Config.sattag);
                        break;
                    case {'render'}
                        ne_render_setview(0, 0, ne_gcfg.cc.(sats{sc}).Config.sattag, 'preview');
                end
            end
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end

% lift restriction
ne_gcfg.c.incb = false;

% check if it was up inbetween
if ne_gcfg.c.btdoup
    ne_gcfg.c.btdown = [];
    ne_gcfg.c.btdoup = false;
    ne_gcfg.fcfg.mpos.ddat = {};
end
