% FUNCTION ne_mdm_acond: VTCCondAverage sub-interface
function ne_mdm_acond(varargin)

% Version:  v1.1
% Build:    16052411
% Date:     May-24 2016, 11:01 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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

% check preliminaries
if nargin < 3 || ~ischar(varargin{3}) || isempty(varargin{3}) || ...
   ~isfield(ne_gcfg.h, 'MDM') || ~isstruct(ne_gcfg.h.MDM) || ~isfield(ne_gcfg.h.MDM, 'MDMFig') || ...
   ~isxfigure(ne_gcfg.h.MDM.MDMFig, true) || ~isfield(ne_gcfg.fcfg, 'MDM') || ...
   ~isstruct(ne_gcfg.fcfg.MDM) || ~isfield(ne_gcfg.fcfg.MDM, 'mdm') || ~isxff(ne_gcfg.fcfg.MDM.mdm, 'mdm')
    return;
end
hFig = ne_gcfg.h.MDM.MDMFig;
cc = ne_gcfg.fcfg.MDM;
ch = ne_gcfg.h.MDM.h;
mdm = cc.mdm;
action = lower(varargin{3}(:)');

% what to do
switch (action)

    % add conditions to averaging
    case 'addcoll'

        % get list to add
        addlist = ch.PRTConds.String(ch.PRTConds.Value);
        if isempty(addlist)
            return;
        end
        if ~iscell(addlist)
            addlist = cellstr(addlist);
        end
        addlist = addlist(:);
        if ~isempty(ch.VTCAvgConds.UserData) && ...
            any(multimatch(addlist, ch.VTCAvgConds.UserData(:, 1)) > 0)
            uiwait(warndlg('Condition(s) already selected.', 'NeuroElf - info', 'modal'));
            return;
        end

        % multiple conditions
        if numel(addlist) == 1
            ne_mdm_acond(0, 0, 'addcond');
            return;
        end

        % combine
        cans = inputdlg({'Please enter the collapsed condition name:'}, ...
            'NeuroElf - user input', 1, {'  '});
        if ~iscell(cans) || numel(cans) ~= 1 || isempty(ddeblank(cans{1}))
            return;
        end
        cans = ddeblank(cans{1});
        colstring = sprintf('%s|', addlist{:});
        colstring = sprintf('^(%s)$', colstring(1:end-1));
        addlist = [addlist(:), repmat({colstring, cans}, numel(addlist), 1)];
        addstring = {cans};

        % add to list
        if isempty(ch.VTCAvgConds.UserData)
            ch.VTCAvgConds.String = addstring;
            ch.VTCAvgConds.UserData = addlist;
        else
            ch.VTCAvgConds.String = [ch.VTCAvgConds.String; addstring];
            ch.VTCAvgConds.UserData = [ch.VTCAvgConds.UserData; addlist];
        end
        ch.VTCAvgRun.Enable = 'on';

    % add conditions to averaging
    case 'addconds'

        % get list to add
        addlist = ch.PRTConds.String(ch.PRTConds.Value);
        if isempty(addlist)
            return;
        end
        if ~iscell(addlist)
            addlist = cellstr(addlist);
        end
        addlist = addlist(:);
        if ~isempty(ch.VTCAvgConds.UserData) && ...
            any(multimatch(addlist, ch.VTCAvgConds.UserData(:, 1)) > 0)
            uiwait(warndlg('Condition(s) already selected.', 'NeuroElf - info', 'modal'));
            return;
        end
        addstring = addlist;
        addlist = [addlist, addlist, addlist];

        % add to list
        if isempty(ch.VTCAvgConds.UserData)
            ch.VTCAvgConds.String = addstring;
            ch.VTCAvgConds.UserData = addlist;
        else
            ch.VTCAvgConds.String = [ch.VTCAvgConds.String; addstring];
            ch.VTCAvgConds.UserData = [ch.VTCAvgConds.UserData; addlist];
        end
        ch.VTCAvgRun.Enable = 'on';

    % remove conditions
    case 'delconds'

        % get indices to remove
        rmidx = ch.VTCAvgConds.Value;
        if isempty(rmidx)
            return;
        end

        % something went wrong earlier
        if isempty(ch.VTCAvgConds.UserData)
            ch.VTCAvgConds.UserData = cell(0, 3);
            ch.VTCAvgConds.String = cell(0, 1);
            ch.VTCAvgConds.Value = [];
            ch.VTCAvgConds.ListboxTop = 1;
            return;
        end

        % select which rows from UserData to remove
        fulllist = ch.VTCAvgConds.String;
        if ~iscell(fulllist)
            fulllist = cellstr(fulllist);
        end
        fulllist = fulllist(:);
        remlist = ch.VTCAvgConds.String(rmidx);
        if ~iscell(remlist)
            remlist = cellstr(remlist);
        end
        remlist = remlist(:);
        rmuidx = (multimatch(ch.VTCAvgConds.UserData(:, 3), remlist) > 0);
        ch.VTCAvgConds.UserData(rmuidx, :) = [];
        fulllist(rmidx) = [];
        if isempty(fulllist)
            ch.VTCAvgRun.Enable = 'off';
            fulllist = {'<no conditions selected>'};
        end
        ch.VTCAvgConds.String = fulllist;
        ch.VTCAvgConds.Value = [];
        ch.VTCAvgConds.ListboxTop = 1;

    % move down
    case 'movedown'

    % move up
    case 'moveup'

    % naive mode
    case 'naive'
        if ch.VTCAvgNaive.Value > 0
            gonoff = 'off';
        else
            gonoff = 'on';
        end
        hFig.SetGroupEnabled('AVNaive', gonoff);

    % average
    case 'average'

        % collect options
        conds = ch.VTCAvgConds.String;
        if ~iscell(conds)
            conds = cellstr(conds);
        end
        conds = conds(:);
        avgwin = round(1000 * str2double(ddeblank(ch.VTCAvgWin.String)));
        if isinf(avgwin) || isnan(avgwin) || avgwin < 5000 || avgwin > 180000
            uiwait(warndlg('Invalid averaging window.', 'NeuroElf - error', 'modal'));
            return;
        end
        basewin = ch.VTCAvgBaseWin.String;
        if sum(basewin == ':') ~= 2
            uiwait(warndlg('Invalid baseline window.', 'NeuroElf - error', 'modal'));
            return;
        end
        basewin = splittocell(ddeblank(basewin), ':');
        basewin = str2double(basewin(:));
        if any(isinf(basewin) | isnan(basewin)) || basewin(2) <= 0
            uiwait(warndlg('Invalid baseline window.', 'NeuroElf - error', 'modal'));
            return;
        end
        basewin = round(1000 .* (basewin(1):basewin(2):basewin(3)));
        if numel(basewin) > 10
            uiwait(warndlg('Invalid baseline window.', 'NeuroElf - error', 'modal'));
            return;
        end
        collapsed = ch.VTCAvgConds.UserData(:, [2, 3]);
        [ucl, uci] = unique(collapsed(:, 2));
        collapsed = collapsed(uci, :);
        doffx = (ch.VTCAvgFFX.Value > 0);
        dorfx = (ch.VTCAvgRFX.Value > 0);
        dowrfx = (ch.VTCAvgWRFX.Value > 0);
        samptr = round(1000 * str2double(ddeblank(ch.VTCAvgSampleTR.String)));
        if isinf(samptr) || isnan(samptr) || samptr < 100 || samptr > 2500
            uiwait(warndlg('Invalid sampling rate.', 'NeuroElf - error', 'modal'));
            return;
        end
        if ch.TransPSC.Value > 0
            tctrans = 'psc';
        elseif ch.Transz.Value > 0
            tctrans = 'z';
        else
            tctrans = 'none';
        end

        % compile to struct
        opts = struct( ...
            'avgwin',   avgwin, ...
            'basewin',  basewin, ...
            'collapse', {collapsed}, ...
            'ffx',      doffx, ...
            'naive',    (ch.VTCAvgNaive.Value > 0), ...
            'remgsig',  (ch.VTCAvgRemGS.Value > 0), ...
            'rfx',      dorfx, ...
            'robust',   (ch.VTCAvgRobust.Value > 0), ...
            'rsngtrial', (ch.VTCAvgSngTrl.Value > 0), ...
            'samptr',   samptr, ...
            'trans',    tctrans, ...
            'wrfx',     dowrfx);

        % compute
        vtcs = cell(1, sum(double([doffx, dorfx, dowrfx])));
        try
            hFig.Visible = 'off';
            drawnow;
            [vtcs{:}] = mdm.VTCCondAverage(conds, opts);
        catch ne_eo;
            uiwait(warndlg(sprintf('Error averaging VTCs: %s.', ne_eo.message), ...
                'NeuroElf - error', 'modal'));
            if isxfigure(hFig, true)
                hFig.Visible = 'on';
                drawnow;
            end
            return;
        end

        % add to interface
        for vc = 1:numel(vtcs)
            if isxff(vtcs{vc}, 'vtc')
                ne_openfile(0, 0, vtcs{vc}, true);
            end
        end
        if isxfigure(hFig, true)
            hFig.Visible = 'on';
            drawnow;
        end

    otherwise
        uiwait(warndlg(sprintf('Unknown action: %s.', action), 'NeuroElf - error', 'modal'));
end
