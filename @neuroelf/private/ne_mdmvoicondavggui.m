% FUNCTION ne_mdmvoicondavggui: GUI callback (updating the config)
function ne_mdmvoicondavggui(varargin)

% Version:  v1.1
% Build:    16021714
% Date:     Feb-17 2016, 2:18 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2015, 2016, Jochen Weber
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

% check argument
if nargin < 4 || ...
   ~ischar(varargin{3}) || ...
   ~isfield(ne_gcfg.cc, varargin{3}(:)') || ...
   ~ischar(varargin{4})
    return;
end
tstr = varargin{3}(:)';
cc = ne_gcfg.cc.(tstr);
hFig = cc.Satellite;
tags = cc.Tags;
cc = cc.Config;
ax = cc.ax;
bands = cc.bands;

% depending on where the call comes from
action = lower(varargin{4}(:)');
switch (action)

    % make sure AvgTR is valid
    case {'avgtr'}
        try
            avgtr = min(3000, max(10, fix(1000 .* str2double(tags.AvgTR.String))));
            if isinf(avgtr) || ...
                isnan(avgtr)
                error('BAD_AVGTR');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            avgtr = cc.avgopt.samptr;
        end
        cc.avgopt.samptr = avgtr;
        tags.AvgTR.String = sprintf('%.2f', 0.001 * avgtr);

        % update timcourses
        ne_gcfg.cc.(tstr).Config = ne_mdmvoicondavgrex(cc);
        ne_mdmvoicondavgup(0, 0, tstr);

    % make sure AvgWin is valid
    case {'avgwin'}
        try
            avgbwin = splittocellc(tags.AvgWin.String, ':');
            if numel(avgbwin) ~= 2
                error('BAD_AVGWIN');
            end
            avgbegin = min(30000, max(-30000, fix(1000 .* str2double(avgbwin{1}))));
            avgwin = min(120000, max(5000, fix(1000 .* str2double(avgbwin{2})))) - avgbegin;
            if isinf(avgbegin) || ...
                isnan(avgbegin) || ...
                isinf(avgwin) || ...
                isnan(avgwin)
                error('BAD_AVGWIN');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            avgbegin = cc.avgopt.avgbegin;
            avgwin = cc.avgopt.avgwin;
        end
        cc.avgopt.avgbegin = avgbegin;
        cc.avgopt.avgwin = avgwin;
        tags.AvgWin.String = sprintf('%.1f:%.1f', 0.001 * [avgbegin, avgwin + avgbegin]);

        % update timcourses
        ne_gcfg.cc.(tstr).Config = ne_mdmvoicondavgrex(cc);
        ne_mdmvoicondavgup(0, 0, tstr);

    % set new axes background color
    case {'backcolor'}

        % pick new color
        try
            nc = colorpicker(255 .* get(tags.Axes.MLHandle, 'Color'));
            set(ax, 'Color', (1 / 255) .* nc);
            ne_gcfg.cc.(tstr).Config.axcol = nc;
            tags.BackColor.CData = repmat(uint8(reshape(nc, [1, 1, 3])), [12, 18]);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

        % update to correctly alter legend
        ne_mdmvoicondavgup(0, 0, tstr);

    % make sure BaseWin is valid
    case {'basewin'}
        try
            basewin = tags.BaseWin.String;
            basewin = basewin(:)';
            if sum(basewin == ':') ~= 2
                error('BAD_BASEWIN');
            end
            basewin = splittocellc(basewin, ':');
            if numel(basewin) ~= 3
                error('BAD_BASEWIN');
            end
            basewin{1} = fix(1000 .* str2double(basewin{1}));
            basewin{2} = fix(1000 .* str2double(basewin{2}));
            basewin{3} = fix(1000 .* str2double(basewin{3}));
            if numel(basewin{1}) ~= 1 || ...
                isinf(basewin{1}) || ...
                isnan(basewin{1}) || ...
                numel(basewin{2}) ~= 1 || ...
                isinf(basewin{2}) || ...
                isnan(basewin{2}) || ...
                basewin{2} <= 0 || ...
                numel(basewin{3}) ~= 1 || ...
                isinf(basewin{3}) || ...
                isnan(basewin{3})
                error('BAD_BASEWIN');
            end
            while numel(basewin{1}:basewin{2}:basewin{3}) > 20
                basewin{2} = 2 * basewin{2};
            end
            basewin = basewin{1}:basewin{2}:basewin{3};
            if isempty(basewin)
                error('BAD_BASEWIN');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            basewin = cc.avgopt.basewin;
        end
        cc.avgopt.basewin = basewin;
        basewin = 0.001 .* basewin;
        if numel(basewin) > 1
            tags.BaseWin.String = sprintf('%g:%g:%g', basewin(1),  ...
                basewin(2) - basewin(1), basewin(end));
        else
            tags.BaseWin.String = sprintf('%g:1:%g', basewin(1), basewin(1));
        end

        % update timcourses
        ne_gcfg.cc.(tstr).Config = ne_mdmvoicondavgrex(cc);
        ne_mdmvoicondavgup(0, 0, tstr);

    % add collapsing
    case {'colladd'}

        % get indices of plot curves
        cval = tags.Conditions.Value;

        % if empty, return
        if isempty(cval)
            return;
        end

        % get name for collapsing
        cvalname = sprintf('%s, ', cc.conds{cval});
        collname = inputdlg({'Please enter a name for collapsing:'}, ...
            'NeuroElf - input', 1, {cvalname(1:end-2)});
        if isempty(collname)
            return;
        end
        collname = collname{1};

        % add to collapsed
        ne_gcfg.cc.(tstr).Config.colls(end+1, :) = ...
            {collname, cval(:)', max(0, min(255, round(mean(cc.condcol(cval, :), 1))))};

        % and to control
        collnames = tags.Collapsed.String;
        if ~iscell(collnames)
            collnames = cellstr(collnames);
        end
        if numel(collnames) == 1 && ...
            strcmpi(collnames, 'none defined')
            collnames = cell(0, 1);
        end
        collnames{end+1} = collname;
        tags.Collapsed.String = collnames;
        tags.Collapsed.Value = numel(collnames);

        % unset conditions
        tags.Conditions.Value = [];

        % enable group
        hFig.SetGroupEnabled('CGrp', 'on');

        % update display and return
        ne_mdmvoicondavgup(0, 0, tstr);

    % select collapsed colors
    case {'collcolors'}

        % update colors and selection (if conditions selected)
        if ~isempty(cc.colls)
            collcols = colorpicker(cat(1, cc.colls{:, 3}), cc.colls(:, 1));
            for c = 1:size(collcols, 1)
                ne_gcfg.cc.(tstr).Config.colls{c, 3} = collcols(c, :);
            end
        end
        ne_mdmvoicondavgup(0, 0, tstr);

    % delete collapsing
    case {'colldelete'}

        % get value to delete
        clval = tags.Collapsed.Value;

        % empty, return
        if isempty(clval)
            return;
        end

        % delete from lists
        ne_gcfg.cc.(tstr).Config.colls(clval, :) = [];
        collnames = tags.Collapsed.String;
        if ~iscell(collnames)
            collnames = cellstr(collnames);
        end
        collnames(clval) = [];
        if isempty(collnames)
            collnames = {'none defined'};
            hFig.SetGroupEnabled('CGrp', 'off');
        end
        tags.Collapsed.String = collnames;
        tags.Collapsed.Value = [];

        % update display and return
        ne_mdmvoicondavgup(0, 0, tstr);

    % select condition colors
    case {'condcolors'}

        % update colors and selection (if conditions selected)
        ne_gcfg.cc.(tstr).Config.condcol = colorpicker(cc.condcol, cc.conds);
        ne_mdmvoicondavgup(0, 0, tstr);

    % set data in base workspace
    case {'copydata'}

        % get and copy data
        if ~isempty(cc.lines)
            cdata = get(cc.lines, 'YData');
            if ~iscell(cdata)
                cdata = {cdata};
            end
            for c = 1:numel(cdata)
                cdata{c} = lsqueeze(cdata{c});
            end
            assignin('base', 'lines', cdata);
        end
        if ~isempty(cc.bands)
            bdata = get(cc.bands, 'YData');
            if ~iscell(bdata)
                bdata = {bdata};
            end
            for c = 1:numel(bdata)
                bdata{c} = lsqueeze(bdata{c}(1:floor(0.5 * numel(bdata{c})))) - cdata{c};
            end
            assignin('base', 'bands', bdata);
        end
        if ~isempty(cc.curves)
            cdata = get(cc.curves, 'YData');
            if ~iscell(cdata)
                cdata = {cdata};
            end
            for c = 1:numel(cdata)
                cdata{c} = lsqueeze(cdata{c});
            end
            assignin('base', 'curves', cdata);
        end

    % copy plot to new figure
    case {'copyplot'}

        % first use saveaxes feature
        tfigname = [tempname '.fig'];
        ne_mdmvoicondavggui(0, 0, tstr, 'saveaxes', tfigname);

        % then open again
        try
            tfig = openfig(tfigname, 'new', 'invisible');
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end

        % and delete file
        try
            delete(tfigname);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

        % then set the figure position and show
        try
            set(tfig, 'Position', get(tfig, 'Position') + [32, -32, 0, 0]);
            set(tfig, 'Units', 'normalized');
            set(get(tfig, 'Children'), 'Units', 'normalized');
            set(tfig, 'Resize', 'on');
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        try
            set(tfig, 'Visible', 'on');
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

    % set raw data in base workspace
    case {'copyrawdata'}

        % copy important raw data elements
        assignin('base', 'mvc_avgopt', cc.avgopt);
        assignin('base', 'mvc_baseline', cc.btc);
        assignin('base', 'mvc_extracts', cc.mtc);
        assignin('base', 'mvc_mdm', cc.mdm);
        assignin('base', 'mvc_prtconds', cc.conds);
        assignin('base', 'mvc_prtonsets', cc.onsets);
        assignin('base', 'mvc_tcdata', cc.rtc(:));
        assignin('base', 'mvc_tcfiles', cc.tcf);
        assignin('base', 'mvc_tctr', cc.tr);
        assignin('base', 'mvc_voi', cc.voi);

    % set diff mode
    case {'diff1', 'diff2'}
        diffmode = double(action(end)) - 48;
        if cc.diffmode == diffmode
            diffmode = 0;
        end
        ne_gcfg.cc.(tstr).Config.diffmode = diffmode;
        if diffmode > 0 && ...
            tags.Type.Value == 1
            tags.Type.Value = 2;
            drawnow;
        end
        ne_mdmvoicondavgup(0, 0, tstr);

    % set error bar color(s)
    case {'ebandalpha'}

        % get value
        try
            sealpha = min(1, max(0, str2double(tags.EBandAlpha.String)));
            if numel(sealpha) ~= 1 || ...
                isinf(sealpha) || ...
                isnan(sealpha)
                error('BAD_SEALPHA');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            sealpha = cc.sealpha;
        end
        ne_gcfg.cc.(tstr).Config.sealpha = sealpha;
        tags.EBandAlpha.String = sprintf('%.2f', sealpha);
        ne_mdmvoicondavgup(0, 0, tstr);

    % en/disable error bars
    case {'ebands'}

        % set visibility
        if tags.EBands.Value > 0
            visflag = 'on';
        else
            visflag = 'off';
        end
        set(bands, 'Visible', visflag);
        drawnow;

    % set groups
    case {'groups'}

        % set groups
        ne_gcfg.cc.(tstr).Config.groups = ...
            tags.Groups.Value;

        % and update
        ne_mdmvoicondavgup(0, 0, tstr);

    % set legend location
    case {'legloc'}

        % check input
        if nargin < 5 || ...
           ~ischar(varargin{5}) || ...
           ~any(strcmp(varargin{5}(:)', ...
                {'NorthEast', 'NorthWest', 'SouthEast', 'SouthWest'}))
            return;
        end

        % update setting
        legloc = varargin{5}(:)';
        ne_gcfg.cc.(tstr).Config.legloc = legloc;

        % update menu string
        tags.LLNE.Checked = 'off';
        tags.LLNW.Checked = 'off';
        tags.LLSE.Checked = 'off';
        tags.LLSW.Checked = 'off';
        switch (legloc)
            case {'NorthEast'}
                tags.LLNE.Checked = 'on';
            case {'NorthWest'}
                tags.LLNW.Checked = 'on';
            case {'SouthEast'}
                tags.LLSE.Checked = 'on';
            case {'SouthWest'}
                tags.LLSW.Checked = 'on';
        end

        % update graphics
        ne_mdmvoicondavgup(0, 0, tstr);

    % toggle p<0.05 flag
    case {'p05'}
        ne_gcfg.cc.(tstr).Config.p05 = ~cc.p05;
        if ne_gcfg.cc.(tstr).Config.p05
            tags.P05.Checked = 'on';
            if tags.SDSE.Value == 0
                tags.SDSE.Value = 1;
                drawnow;
            end
        else
            tags.P05.Checked = 'off';
        end

        % update
        ne_mdmvoicondavgup(0, 0, tstr);

    % robust computation
    case {'robust'}

        % update flag
        cc.avgopt.robust = (tags.Robust.Value > 0);
        ne_gcfg.c.ini.Statistics.RobustVOICondAvg = cc.avgopt.robust;

        % update timcourses
        ne_gcfg.cc.(tstr).Config = ne_mdmvoicondavgrex(cc);

        % for now, just redo graphs
        ne_mdmvoicondavgup(0, 0, tstr);

    % save screenshot
    case {'saveaxes'}

        % get optvisible setting
        ovis = tags.OptVis.Value;

        % make setting
        if ovis > 0
            tags.OptVis.Value = 0;
            ne_mdmvoicondavgrsz(0, 0, tstr, 'OptVis');
        end
        tags.OptVis.Visible = 'off';
        set(tags.OptVis.Parent, 'Resize', 'off');
        drawnow;

        % try saving
        try
            if nargin > 4
                ne_screenshot(0, 0, tags.OptVis.Parent, varargin{5}, 'high-q');
            else
                ne_screenshot(0, 0, tags.OptVis.Parent, '', 'high-q');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

        % restore settings
        set(tags.OptVis.Parent, 'Resize', 'on');
        tags.OptVis.Visible = 'on';
        if ovis
            tags.OptVis.Value = 1;
            ne_mdmvoicondavgrsz(0, 0, tstr, 'OptVis');
        end
        drawnow;

    % save data file
    case {'savedata'}

        % get filename
        [savefile, savepath] = uiputfile( ...
            {'*.mat', 'MDM voi-extracted MAT-files (*.mat)'}, ...
            'Save extracted data as...');
        if isequal(savefile, 0) || ...
            isequal(savepath, 0) || ...
            isempty(savefile)
            return;
        end
        if isempty(savepath)
            savepath = pwd;
        end

        % compile data
        mvc_data = struct;
        mvc_data.btc = cc.btc;
        mvc_data.mtc = cc.mtc;
        mvc_data.rtc = cc.rtc;
        mvc_data.avgopt = rmfield(cc.avgopt, 'pbar');
        mvc_data.colls = cc.colls;
        mvc_data.condcol = cc.condcol;
        mvc_data.conds = cc.conds;
        mvc_data.mdm = getcont(cc.mdm);
        mvc_data.onsets = cc.onsets;
        mvc_data.sdms = cc.sdms;
        mvc_data.subsel = cc.subsel;
        mvc_data.tcf = cc.tcf;
        mvc_data.tr = cc.tr;
        mvc_data.voi = getcont(cc.voi);

        % try saving
        try
            save([savepath '/' savefile], 'mvc_data');
        catch ne_eo;
            uiwait(warndlg(sprintf('Error saving file: %s', ne_eo.message), ...
                'NeuroElf - warning', 'modal'));
        end

    % subject selection
    case {'subsel'}

        % get new subject selection
        [subsel, ssok] = listdlg( ...
            'PromptString', 'Please select subjects to use...', ...
            'ListString',   cc.mdm.Subjects, ...
            'InitialValue', cc.subsel, ...
            'Name',         'NeuroElf - MDM VOI cond average request', ...
            'ListSize',     [min(200, max(640, 10 .* size(char(cc.mdm.Subjects), 2))), 300]);
        if ~isequal(ssok, 1)
            return;
        end

        % don't allow empty selection
        if isempty(subsel)
            uiwait(warndlg('The subject selection must not be empty!', ...
                'NeuroElf - warning', 'modal'));
            return;
        end

        % set new subject selection and re-run
        ne_gcfg.cc.(tstr).Config.subsel = subsel(:);
        if numel(subsel) < 2
            cc.typeuic.Value = 1;
            cc.typeuic.Enable = 'off';
        else
            cc.typeuic.Enable = 'on';
        end
        ne_mdmvoicondavgup(0, 0, tstr);

    % temporal smoothing
    case {'tsmooth'}

        % request new value
        tsmval = inputdlg({'Temporal smoothing kernel:'}, ...
            'NeuroElf - user input', 1, {sprintf(' %g', cc.tsmooth)});
        if ~iscell(tsmval) || ...
            numel(tsmval) ~= 1 || ...
           ~ischar(tsmval{1}) || ...
            isempty(ddeblank(tsmval{1}))
            return;
        end
        try
            tsmval = str2double(ddeblank(tsmval{1}));
            if numel(tsmval) ~= 1 || ...
               ~isa(tsmval, 'double') || ...
                isinf(tsmval) || ...
                isnan(tsmval) || ...
                tsmval < 0
                return;
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end
        tsmval = 0.01 * round(100 * min(24, tsmval));
        ne_gcfg.cc.(tstr).Config.tsmooth = tsmval;
        tags.TSmooth.Label = sprintf('Temporal smoothing (%gs)', tsmval);
        ne_mdmvoicondavgup(0, 0, tstr);

    % switch VOI
    case {'voi'}

        % get index
        if nargin < 5 || ...
           ~isa(varargin{5}, 'double') || ...
            numel(varargin{5}) ~= 1 || ...
            isinf(varargin{5}) || ...
            isnan(varargin{5}) || ...
            varargin{5} < 1 || ...
            varargin{5} ~= fix(varargin{5})
            voiidx = tags.VOI.Value;
        else
            voiidx = min(varargin{5}, numel(cc.voi.VOI));
            tags.VOI.Value = voiidx;
        end

        % update UI
        ne_mdmvoicondavgup(0, 0, tstr);

        % still the same object
        if tags.VOIUpdate.Value > 0 && isxff(ne_gcfg.voi, 'voi') && ...
            ne_gcfg.voi == cc.voi.Handles.SourceObject

            % set cluster
            ne_gcfg.h.Clusters.Value = voiidx;
            ne_setcluster(0, 0, 'set');

        % still update the position (at least)
        else

            % set to first coordinate of VOI
            if ~isempty(cc.voi.VOI(voiidx).Voxels)
                ne_setslicepos(0, 0, cc.voi.VOI(voiidx).Voxels(1, :));
            end
        end

    % update range
    case {'xrange', 'yrange'}

        % get current and new range values
        cxlim = get(ax, 'XLim');
        cylim = get(ax, 'YLim');

        % update from UI
        try
            cxlim(1) = str2double(tags.XFrom.String);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        try
            cxlim(2) = str2double(tags.XTo.String);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        try
            cylim(1) = str2double(tags.YFrom.String);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        try
            cylim(2) = str2double(tags.YTo.String);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

        % checks
        if any(isinf(cxlim) | isnan(cxlim) | isinf(cylim) | isnan(cylim)) || ...
            cxlim(2) <= cxlim(1) || ...
            cylim(2) <= cylim(1)
            cxlim = get(ax, 'XLim');
            cylim = get(ax, 'YLim');
            if strcmpi(varargin{4}(:)', 'xrange')
                set(ax, 'XLimMode', 'auto');
                cxlim = get(ax, 'XLim');
            else
                set(ax, 'YLimMode', 'auto');
                cylim = get(ax, 'YLim');
            end
        else
            if strcmpi(varargin{4}(:)', 'xrange')
                set(ax, 'XLim', cxlim);
            else
                set(ax, 'YLim', cylim);
            end
        end

        % update text boxes
        tags.XFrom.String = sprintf('%.2f', cxlim(1));
        tags.XTo.String = sprintf('%.2f', cxlim(2));
        tags.YFrom.String = sprintf('%.2f', cylim(1));
        tags.YTo.String = sprintf('%.2f', cylim(2));
end



% % % SUB-FUNCTIONS % % %

function cc = ne_mdmvoicondavgrex(cc)
cc.fig.SetGroupEnabled('All', 'off');
cc.fig.Pointer = 'watch';
drawnow;
try
    [cc.mtc, mtcse, cc.btc] = cc.mdm.VOICondAverage( ...
        cc.voi, cc.conds, cc.avgopt, cc.rtc, cc.tcf, cc.tr, cc.onsets, cc.sdms);
    cc.btc = meannoinfnan(cc.btc, 1);
    cc.btc(cc.btc == 0) = NaN;
catch ne_eo;
    uiwait(warndlg(['Error updating timecourses: ' ne_eo.message], ...
        'NeuroElf - warning', 'modal'));
end
cc.fig.SetGroupEnabled('All', 'on');
if isempty(cc.colls)
    cc.fig.SetGroupEnabled('CGrp', 'off');
end
if isempty(cc.mdm.RunTimeVars.Groups)
    cc.fig.SetGroupEnabled('UGrp', 'off');
end
if numel(cc.subsel) < 2
    cc.typeuic.Value = 1;
    cc.typeuic.Enable = 'off';
else
    cc.typeuic.Enable = 'on';
end
cc.fig.Pointer = 'arrow';
drawnow;
