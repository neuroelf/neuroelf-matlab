% FUNCTION ne_glmplotbetasgui: GUI callback (updating the config)
function ne_glmplotbetasgui(varargin)

% Version:  v1.1
% Build:    16041810
% Date:     Apr-18 2016, 10:35 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2015, 2016, Jochen Weber
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
cini = ne_gcfg.c.ini;

% check argument
if nargin < 4 || ...
   ~ischar(varargin{3}) || ...
   ~isfield(ne_gcfg.cc, varargin{3}(:)') || ...
   ~ischar(varargin{4})
    return;
end
tstr = varargin{3}(:)';
cc = ne_gcfg.cc.(tstr);
if ~isfield(cc, 'Config') || ...
   ~isstruct(cc.Config) || ...
   ~isfield(cc.Config, 'glm') || ...
   ~isxff(cc.Config.glm, 'glm')
    return;
end
glm = cc.Config.glm;
grp = cc.Config.groups;
glmh = handles(glm);
if ~isfield(glmh, 'PlotFig') || ...
   ~iscell(glmh.PlotFig) || ...
    isempty(glmh.PlotFig)
    return;
end
pfi = [];
for fc = 1:numel(glmh.PlotFig)
    if strcmpi(tstr, glmh.PlotFig{fc}.Tag(1:8))
        pfi = fc;
        break;
    end
end
if isempty(pfi)
    return;
end
ch = glmh.PlotHnd{pfi};
hFig = cc.Satellite;
rtv = glm.RunTimeVars;
tags = cc.Tags;
ax = tags.Axes.MLHandle;

% depending on where the call comes from
switch(lower(varargin{4}(:)'))

    % set new axes background color
    case {'backcolor'}

        % pick new color
        try
            nc = colorpicker(255 .* get(ax, 'Color'));
            set(ax, 'Color', (1 / 255) .* nc);
            ne_gcfg.cc.(tstr).Config.axcol = nc;
            cini.BetaPlot.BackgroundColor = nc;
            tags.BackColor.CData = repmat(uint8(reshape(nc, [1, 1, 3])), [12, 18]);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        return;

    % clone plot window
    case {'cloneplot'}

        % open another beta plotter on same GLM
        miflag = cini.BetaPlot.MultiInstance;
        if ~miflag
            cini.BetaPlot.MultiInstance = true;
        end
        [newaxes, newtstr] = ne_glmplotbetas(0, 0, glm);
        if ~miflag
            cini.BetaPlot.MultiInstance = false;
        end

        % copy configuration
        newhFig = ne_gcfg.cc.(newtstr).Satellite;
        newtags = ne_gcfg.cc.(newtstr).Tags;
        for cpfields = {'axcol', 'barcolors', 'condcol', 'contcol', 'contrasts', ...
            'covariate', 'covcol', 'covsize', 'drawse', 'groups', 'legbars', ...
            'legpos', 'legscat', 'linecolor', 'optvis', 'plotlines', 'plotpts', ...
            'quadratic', 'radius', 'radvox', 'ranktrans', 'regress', 'remove0s', ...
            'robust', 'scellipse', 'scfmarker', 'scgroups', 'scline', 'scmarker', ...
            'scrange', 'scrangea', 'scstats', 'sublabels', 'subsel', 'title', ...
            'type', 'upplot', 'xrange'}
            ne_gcfg.cc.(newtstr).Config.(cpfields{1}) = cc.Config.(cpfields{1});
        end
        set(newaxes, 'Color', (1 / 255) .* ne_gcfg.cc.(newtstr).Config.axcol);
        newtags.BackColor.CData = ...
            repmat(uint8(reshape(cc.Config.axcol, [1, 1, 3])), [12, 18]);
        for cpfields = {'Conditions', 'Contrasts', 'Covariates', 'DoGroups', ...
            'DoUpdate', 'EBars', 'Groups', 'OptVis', 'PPts', 'Type'}
            newtags.(cpfields{1}).Value = tags.(cpfields{1}).Value;
        end
        for cpfields = {'Conditions', 'Contrasts', 'Covariates', 'Groups'}
            newtags.(cpfields{1}).ListboxTop = tags.(cpfields{1}).ListboxTop;
        end
        for cpfields = {'LegendBars', 'LegendScat', 'LegPosNE', 'LegPosNW', ...
            'LegPosSE', 'LegPosSW', 'RankTrans', 'Remove0s', 'Robust', ...
            'ScConfEllipse', 'ScFilledMarkers', 'ScGroups', 'ScLine', ...
            'ScMarkAsterisk', 'ScMarkCircle', 'ScMarkDiamond', ...
            'ScMarkPeriod', 'ScMarkPlus', 'ScMarkSquare', 'ScMarkX', ...
            'ScQuadratic', 'ShowGenInfo', 'StatsOut', 'SubLabels', ...
            'SubLines'}
            newtags.(cpfields{1}).Checked = tags.(cpfields{1}).Checked;
            newtags.(cpfields{1}).Enable = tags.(cpfields{1}).Enable;
        end
        newtags.EBarColors.CData = repmat(uint8(reshape( ...
            ne_gcfg.cc.(newtstr).Config.linecolor, [1, 1, 3])), [12, 18]);
        if numel(newtags.Groups.Value) < 2 || ...
           ~ne_gcfg.cc.(newtstr).Config.scgroups
            newhFig.SetGroupEnabled('SGrps', 'on');
        else
            newhFig.SetGroupEnabled('SGrps', 'off');
        end
        
        % update
        ne_glmplotbetasgui(0, 0, newtstr, 'type');
        ne_glmplotbetasup(0, 0, glm, 'fromdata', newtstr);
        ne_glmplotbetasrsz(0, 0, newtstr);

        % final touches
        for cpfields = {'Radius', 'XFrom', 'XTo', 'YFrom', 'YTo'}
            newtags.(cpfields{1}).String = tags.(cpfields{1}).String;
        end
        for cpfields = {'scrange', 'scrangea', 'xrange'}
            ne_gcfg.cc.(newtstr).Config.(cpfields{1}) = cc.Config.(cpfields{1});
        end

        % return early
        return;

    % close window
    case {'close'}
        glmh.PlotFig(pfi) = [];
        glm.SetHandle('PlotFig', glmh.PlotFig);
        glmh.PlotHnd(pfi) = [];
        glm.SetHandle('PlotHnd', glmh.PlotHnd);
        figpos = hFig.Position;
        try
            hFig.Delete;
            if isfield(ne_gcfg.cc, tstr)
                ne_gcfg.cc = rmfield(ne_gcfg.cc, tstr);
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        cini.BetaPlot.WindowSize = figpos(3:4);
        cini.Children.BetaPlotPosition = figpos(1:2);
        return;

    % color by covariate
    case {'colorbycov'}

        % set in config
        colval = tags.ColorByCov.Value;
        if colval > 1
            colrange = rtv.CovariatesColors{colval - 1, 1};
            colrange = flexinterpn(colrange, [Inf, Inf; 1, 1; 0.01, 1; size(colrange)]);
            ne_gcfg.cc.(tstr).Config.covcol = ...
                {colval - 1, (1 / 255) .* double(colrange), rtv.CovariatesColors{colval - 1, 2}};
        else
            ne_gcfg.cc.(tstr).Config.covcol = {0, zeros(0, 3), [NaN, NaN]};
        end

    % select condition colors
    case {'condcolors'}

        % update colors and selection (if conditions selected)
        pnames = glm.SubjectPredictors;
        if numel(pnames) > size(ne_gcfg.cc.(tstr).Config.condcol, 1)
            pnames(end) = [];
        end
        newcols = colorpicker(cc.Config.condcol, pnames);
        ne_gcfg.cc.(tstr).Config.condcol = newcols;
        glm.RunTimeVars.PredictorColors = newcols;
        if ~isempty(tags.Conditions.Value)
            ne_gcfg.cc.(tstr).Config.barcolors = newcols(tags.Conditions.Value, :);
        else
            return;
        end

    % select conditions
    case {'conditions'}

        % unselect contrasts
        tags.Contrasts.Value = [];

        % none -> all
        subpreds = glm.SubjectPredictors;
        if ~isempty(subpreds) && ...
            strcmpi(subpreds{end}, 'constant')
            subpreds(end) = [];
        end
        nsubpreds = numel(subpreds);
        if isempty(tags.Conditions.Value)
            if nsubpreds <= 18
                tags.Conditions.Value = (1:nsubpreds)';
            else
                tags.Conditions.Value = 1;
                tags.Conditions.ListboxTop = 1;
            end
        end

        % create new contrasts list
        newcon = eye(nsubpreds);
        newcon = newcon(tags.Conditions.Value, :);
        newconnames = tags.Conditions.String;
        if ~iscell(newconnames)
            newconnames = cellstr(newconnames);
        end
        newconnames = newconnames(tags.Conditions.Value);

        % set contrasts and colors
        ne_gcfg.cc.(tstr).Config.contnames = newconnames;
        ne_gcfg.cc.(tstr).Config.contrasts = newcon;
        ne_gcfg.cc.(tstr).Config.barcolors = ...
            cc.Config.condcol(tags.Conditions.Value, :);

    % select conditions by regexp
    case {'condsel'}

        % enter a string
        vsel = inputdlg({'Please enter the string to regexp select conditions:'}, ...
            'NeuroElf - user input', 1, {'  .*'});
        if numel(vsel) ~= 1 || ~iscell(vsel)
            return;
        end
        connames = tags.Conditions.String;
        if ~iscell(connames)
            connames = cellstr(connames);
        end
        vsel = find(~cellfun('isempty', regexpi(connames, ddeblank(vsel{1}))));
        if ~isempty(vsel)
            tags.Conditions.Value = vsel(:);
            ne_glmplotbetasgui(0, 0, tstr, 'conditions');
        end
        return;

    % select contrast colors
    case {'contcolors'}

        % update colors and selection (if conditions selected)
        newcols = colorpicker(cc.Config.contcol, tags.Contrasts.String);
        ne_gcfg.cc.(tstr).Config.contcol = newcols;
        glm.RunTimeVars.ContrastColors = newcols;
        if ~isempty(tags.Contrasts.Value)
            ne_gcfg.cc.(tstr).Config.barcolors = newcols(tags.Contrasts.Value, :);
        else
            return;
        end

    % switch contrast(s)
    case {'contrasts'}

        % unselect conditions
        tags.Conditions.Value = [];

        % get new contrasts
        if isempty(tags.Contrasts.Value)
            tags.Contrasts.Value = 1;
        end
        newconnames = tags.Contrasts.String;
        if ~iscell(newconnames)
            newconnames = cellstr(newconnames);
        end
        newconnames = newconnames(tags.Contrasts.Value);

        % set contrasts and colors
        newcon = cat(2, rtv.Contrasts{tags.Contrasts.Value, 2})';
        newcon(:, end + 1) = 0;
        ne_gcfg.cc.(tstr).Config.contnames = newconnames;
        ne_gcfg.cc.(tstr).Config.contrasts = newcon;
        ne_gcfg.cc.(tstr).Config.barcolors = ...
            cc.Config.contcol(tags.Contrasts.Value, :);

    % copy data
    case {'copydata'}

        % assign fields of data in base workspace
        cdata = ne_gcfg.cc.(tstr).Data;
        df = fieldnames(cdata);
        for dfc = 1:numel(df)
            assignin('base', ['glmvb_' df{dfc}], cdata.(df{dfc}));
        end
        return;

    % copy plot
    case {'copyplot'}

        % first use saveaxes feature
        tfigname = [tempname '.fig'];
        ne_glmplotbetasgui(0, 0, tstr, 'saveaxes', tfigname);

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
        return;

    % switch covariate
    case {'covariates'}

        % set new covariate
        ne_gcfg.cc.(tstr).Config.covariate = ...
            rtv.CovariatesData(:, tags.Covariates.Value);

    % set covariate-based colors
    case {'covcolors'}

        % requires valid selection
        if tags.ColorByCov.Value < 2
            return;
        end
        covidx = tags.ColorByCov.Value - 1;
        if covidx > size(rtv.CovariatesColors, 1)
            return;
        end
        colrange = rtv.CovariatesColors{covidx, 1};

        % update colors and selection (if conditions selected)
        ncolrange = colorpicker(colrange);
        if ~isequal(size(ncolrange), size(colrange)) || isequal(colrange, ncolrange)
            return;
        end
        glm.RunTimeVars.CovariatesColors{covidx, 1} = ncolrange;
        
        % reselect covariate
        ne_glmplotbetasgui(0, 0, tstr, 'ColorByCov');
        return;

    % do groups (at all)
    case {'dogroups'}

        % get state
        dogroups = (tags.DoGroups.Value > 0);

        % set config
        if dogroups
            tags.Groups.Enable = 'on';
            ne_gcfg.cc.(tstr).Config.groups = ...
                tags.Groups.Value;
            if numel(tags.Groups.Value) < 2 || ...
               ~cc.Config.scgroups
                hFig.SetGroupEnabled('SGrps', 'on');
            else
                hFig.SetGroupEnabled('SGrps', 'off');
            end
        else
            tags.Groups.Enable = 'off';
            ne_gcfg.cc.(tstr).Config.groups = [];
            hFig.SetGroupEnabled('SGrps', 'on');
        end

    % keep updated
    case {'doupdate'}

        % switch flag
        ne_gcfg.cc.(tstr).Config.upplot = (tags.DoUpdate.Value > 0);

        % no further update?
        if tags.DoUpdate.Value < 1
            return;
        end

    % set error bar color(s)
    case {'ebarcolors'}

        % update what
        if numel(grp) ~= 1 || ...
           ~cc.Config.scgroups

            % set value
            ne_gcfg.cc.(tstr).Config.linecolor = ...
                colorpicker(cc.Config.linecolor, 'Error bar/line/scatter color');
            tags.EBarColors.CData = repmat(uint8(reshape( ...
                ne_gcfg.cc.(tstr).Config.linecolor, [1, 1, 3])), [12, 18]);
            cini.BetaPlot.ErrorBarColor = ne_gcfg.cc.(tstr).Config.linecolor;
        else
            ne_gcfg.cc.(tstr).Config.grpspecs{grp, 2} = ...
                colorpicker(cc.Config.grpspecs{grp, 2}, 'Line/scatter color');
        end

    % en/disable error bars
    case {'ebars'}

        % set value
        ne_gcfg.cc.(tstr).Config.drawse = (tags.EBars.Value > 0);

    % font size
    case {'fontsize'}

        % setting given
        if nargin < 5 || ...
           ~isa(varargin{5}, 'double') || ...
            numel(varargin{5}) ~= 1 || ...
            isinf(varargin{5}) || ...
            isnan(varargin{5}) || ...
            varargin{5} < 8 || ...
            varargin{5} > 60
            return;
        end

        % make setting
        ne_gcfg.cc.(tstr).Config.fontsize = round(varargin{5});
        cini.Satellites.FontSize = round(varargin{5});

        % disable checkmarks for controls
        tags.FontSize10pt.Checked = 'off';
        tags.FontSize12pt.Checked = 'off';
        tags.FontSize14pt.Checked = 'off';
        tags.FontSize16pt.Checked = 'off';
        tags.FontSize20pt.Checked = 'off';

        % check single control
        if any(varargin{5} == [10, 12, 14, 16, 20])
            fsf = sprintf('FontSize%dpt', varargin{5});
            tags.(fsf).Checked = 'on';
        end

        % call resize function
        ne_glmplotbetasrsz(0, 0, tstr);
        
    % set groups
    case {'groups'}

        % set groups
        ne_gcfg.cc.(tstr).Config.groups = ...
            tags.Groups.Value;
        if numel(tags.Groups.Value) < 2 || ...
           ~cc.Config.scgroups
            hFig.SetGroupEnabled('SGrps', 'on');
        else
            hFig.SetGroupEnabled('SGrps', 'off');
        end

    % toggle legend for bars
    case {'legendbars'}
        ne_gcfg.cc.(tstr).Config.legbars = ~cc.Config.legbars;
        if ne_gcfg.cc.(tstr).Config.legbars
            tags.LegendBars.Checked = 'on';
        else
            tags.LegendBars.Checked = 'off';
        end
        if ~strcmpi(cc.Config.type, 'bar')
            return;
        end

    % update position for legend
    case {'legendpos'}
        if nargin < 5 || ...
           ~ischar(varargin{5}) || ...
           ~any(strcmp(varargin{5}(:)', ...
                {'NorthEast', 'NorthWest', 'SouthEast', 'SouthWest'}))
            return;
        end
        ne_gcfg.cc.(tstr).Config.legpos = varargin{5}(:)';
        tags.LegPosNE.Checked = 'off';
        tags.LegPosNW.Checked = 'off';
        tags.LegPosSE.Checked = 'off';
        tags.LegPosSW.Checked = 'off';
        switch (varargin{5}(:)')
            case {'NorthEast'}
                tags.LegPosNE.Checked = 'on';
            case {'NorthWest'}
                tags.LegPosNW.Checked = 'on';
            case {'SouthEast'}
                tags.LegPosSE.Checked = 'on';
            case {'SouthWest'}
                tags.LegPosSW.Checked = 'on';
        end

    % toggle legend for scatters
    case {'legendscat'}
        ne_gcfg.cc.(tstr).Config.legscat = ~cc.Config.legscat;
        if ne_gcfg.cc.(tstr).Config.legscat
            tags.LegendScat.Checked = 'on';
        else
            tags.LegendScat.Checked = 'off';
        end
        if ~strcmpi(cc.Config.type, 'scatter')
            return;
        end

    % update options (menu)
    case {'options'}
        markertypes = {'d', 'o', 's', 'x', '*', '+', '.'};
        markertags = {'Diamond', 'Circle', 'Square', 'X', 'Asterisk', 'Plus', 'Period'};
        if cc.Config.remove0s
            tags.Remove0s.Checked = 'on';
        else
            tags.Remove0s.Checked = 'off';
        end
        cini.BetaPlot.Remove0s = cc.Config.remove0s;
        if cc.Config.robust
            tags.Robust.Checked = 'on';
        else
            tags.Robust.Checked = 'off';
        end
        cini.BetaPlot.Robust = cc.Config.robust;
        if cc.Config.scgroups && ...
            numel(cc.Config.groups) ~= 1
            hFig.SetGroupEnabled('SGrps', 'off');
            tags.ScLine.Checked = 'off';
            tags.ScQuadratic.Checked = 'off';
            tags.ScConfEllipse.Checked = 'off';
            tags.StatsOut.Checked = 'off';
        elseif cc.Config.scgroups
            hFig.SetGroupEnabled('SGrps', 'on');
            if cc.Config.grpspecs{tags.Groups.Value, 3}(2) && ...
               ~cc.Config.grpspecs{tags.Groups.Value, 3}(3)
                tags.ScLine.Checked = 'on';
            else
                tags.ScLine.Checked = 'off';
            end
            if all(cc.Config.grpspecs{tags.Groups.Value, 3}(2:3))
                tags.ScQuadratic.Checked = 'on';
            else
                tags.ScQuadratic.Checked = 'off';
            end
            if cc.Config.grpspecs{tags.Groups.Value, 3}(4)
                tags.ScConfEllipse.Checked = 'on';
            else
                tags.ScConfEllipse.Checked = 'off';
            end
            if cc.Config.grpspecs{tags.Groups.Value, 3}(1)
                tags.StatsOut.Checked = 'on';
            else
                tags.StatsOut.Checked = 'off';
            end
            for mt = 1:numel(markertags)
                if cc.Config.grpspecs{tags.Groups.Value, 1} ~= markertypes{mt}
                    mtc = 'off';
                else
                    mtc = 'on';
                end
                tags.(sprintf('ScMark%s', markertags{mt})).Checked = mtc;
            end
        else
            hFig.SetGroupEnabled('SGrps', 'on');
            if cc.Config.scline && ...
               ~cc.Config.quadratic
                tags.ScLine.Checked = 'on';
            else
                tags.ScLine.Checked = 'off';
            end
            if cc.Config.scline && ...
                cc.Config.quadratic
                tags.ScQuadratic.Checked = 'on';
            else
                tags.ScQuadratic.Checked = 'off';
            end
            if cc.Config.scellipse
                tags.ScConfEllipse.Checked = 'on';
            else
                tags.ScConfEllipse.Checked = 'off';
            end
            if cc.Config.scstats
                tags.StatsOut.Checked = 'on';
            else
                tags.StatsOut.Checked = 'off';
            end
            for mt = 1:numel(markertags)
                if cc.Config.scmarker ~= markertypes{mt}
                    mtc = 'off';
                else
                    mtc = 'on';
                end
                tags.(sprintf('ScMark%s', markertags{mt})).Checked = mtc;
            end
        end

        % return early
        return;

    % en/disable plotting of points
    case {'ppts'}

        % set value
        ppts = (tags.PPts.Value > 0);
        ne_gcfg.cc.(tstr).Config.plotpts = ppts;
        cini.BetaPlot.PlotPoints = ppts;

    % set new sampling radius
    case {'radius'}

        % get new radius value
        radval = tags.Radius.String;
        if isempty(radval) || ...
           ~all(radval >= '0' & radval <= '9')
            radval = '0';
        end
        radval = str2double(radval);
        if radval > 30
            radval = 30;
        end

        % compute voxels
        res = glm.Resolution;
        mrad = res * ceil(radval / res);
        mstp = -mrad:res:mrad;
        [xv, yv, zv] = ndgrid(mstp, mstp, mstp);
        xv = [xv(:), yv(:), zv(:)];
        xv(sqrt(sum(xv .* xv, 2)) > radval, :) = [];

        % put back into configuration
        ne_gcfg.cc.(tstr).Config.radius = radval;
        ne_gcfg.cc.(tstr).Config.radvox = xv;
        tags.Radius.String = sprintf('%d', radval);

    % rank-transform
    case {'ranktrans'}
        ne_gcfg.cc.(tstr).Config.ranktrans = ~cc.Config.ranktrans;
        if ne_gcfg.cc.(tstr).Config.ranktrans
            tags.RankTrans.Checked = 'on';
            if numel(grp) ~= 1 || ...
               ~cc.Config.scgroups
                ne_gcfg.cc.(tstr).Config.quadratic = false;
                ne_gcfg.cc.(tstr).Config.scellipse = false;
            elseif iscell(ne_gcfg.cc.(tstr).Config.grpspecs)
                for igc = 1:size(ne_gcfg.cc.(tstr).Config.grpspecs, 1)
                    ne_gcfg.cc.(tstr).Config.grpspecs{igc, 3}(3:4) = false;
                end
            end
        else
            tags.RankTrans.Checked = 'off';
        end

    % regress out global signal flag
    case {'reggs'}
        if isempty(ne_gcfg.cc.(tstr).Config.gsmap)
            ne_gcfg.cc.(tstr).Config.reggsx = false;
        else
            ne_gcfg.cc.(tstr).Config.reggsx = ~cc.Config.reggsx;
        end
        if ne_gcfg.cc.(tstr).Config.reggsx
            tags.RegGS.Checked = 'on';
        else
            tags.RegGS.Checked = 'off';
        end

    % switch remove0s flag
    case {'remove0s'}
        ne_gcfg.cc.(tstr).Config.remove0s = ~cc.Config.remove0s;
        cini.BetaPlot.Remove0s = ne_gcfg.cc.(tstr).Config.remove0s;
        if ne_gcfg.cc.(tstr).Config.remove0s
            tags.Remove0s.Checked = 'on';
        else
            tags.Remove0s.Checked = 'off';
        end

    % switch robust flag
    case {'robust'}
        ne_gcfg.cc.(tstr).Config.robust = ~cc.Config.robust;
        cini.BetaPlot.Robust = ne_gcfg.cc.(tstr).Config.robust;
        if ne_gcfg.cc.(tstr).Config.robust
            tags.Robust.Checked = 'on';
        else
            tags.Robust.Checked = 'off';
        end

    % save full axes view
    case {'saveaxes'}

        % get optvisible setting
        ovis = tags.OptVis.Value;

        % make setting
        if ovis > 0
            tags.OptVis.Value = 0;
            ne_glmplotbetasrsz(0, 0, tstr, 'OptVis');
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
            ne_glmplotbetasrsz(0, 0, tstr, 'OptVis');
        end
        drawnow;
        try
            delete(tags.Axes.Children);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

    % toggle conf ellipse
    case {'scatterconfellipse'}
        if numel(grp) ~= 1 || ...
           ~cc.Config.scgroups
            ne_gcfg.cc.(tstr).Config.scellipse = ~cc.Config.scellipse;
            if ne_gcfg.cc.(tstr).Config.scellipse
                if cc.Config.scline && ...
                    cc.Config.quadratic
                    ne_gcfg.cc.(tstr).Config.scline = false;
                    ne_gcfg.cc.(tstr).Config.quadratic = false;
                end
            else
            end
        else
            ne_gcfg.cc.(tstr).Config.grpspecs{grp, 3}(4) = ...
                ~cc.Config.grpspecs{grp, 3}(4);
            if all(ne_gcfg.cc.(tstr).Config.grpspecs{grp, 3}(2:4))
                ne_gcfg.cc.(tstr).Config.grpspecs{grp, 3}(2:3) = false;
            end
        end

    % toggle scatter filled markers
    case {'scatterfilledmarkers'}
        ne_gcfg.cc.(tstr).Config.scfmarker = ~cc.Config.scfmarker;
        if ne_gcfg.cc.(tstr).Config.scfmarker
            tags.ScFilledMarkers.Checked = 'on';
        else
            tags.ScFilledMarkers.Checked = 'off';
        end
        cini.BetaPlot.ScatterFilled = ne_gcfg.cc.(tstr).Config.scfmarker;

    % toggle scatter groups
    case {'scattergroups'}
        ne_gcfg.cc.(tstr).Config.scgroups = ~cc.Config.scgroups;
        if ne_gcfg.cc.(tstr).Config.scgroups
            tags.ScGroups.Checked = 'on';
            if numel(cc.Config.groups) < 2 || ...
                tags.DoGroups.Value == 0
                hFig.SetGroupEnabled('SGrps', 'on');
            else
                hFig.SetGroupEnabled('SGrps', 'off');
            end
        else
            tags.ScGroups.Checked = 'off';
            hFig.SetGroupEnabled('SGrps', 'on');
        end

    % toggle scatter reg line
    case {'scatterline'}
        if numel(grp) ~= 1 || ...
           ~cc.Config.scgroups
            if ~cc.Config.quadratic
                ne_gcfg.cc.(tstr).Config.scline = ~cc.Config.scline;
            else
                ne_gcfg.cc.(tstr).Config.scline = true;
            end
            if ne_gcfg.cc.(tstr).Config.scline
                ne_gcfg.cc.(tstr).Config.quadratic = false;
            end
        else
            if ~cc.Config.grpspecs{grp, 3}(3)
                ne_gcfg.cc.(tstr).Config.grpspecs{grp, 3}(2) = ...
                    ~cc.Config.grpspecs{grp, 3}(2);
            else
                ne_gcfg.cc.(tstr).Config.grpspecs{grp, 3}(2:3) = ...
                    [true, false];
            end
        end

    % set scatter marker
    case {'scattermarker'}

        % valid input
        markertypes = {'d', 'o', 's', 'x', '*', '+', '.'};
        if nargin < 5 || ...
           ~ischar(varargin{5}) || ...
            numel(varargin{5}) ~= 1 || ...
           ~any(strcmp(varargin{5}, markertypes))
            return;
        end

        % update what
        if numel(grp) ~= 1 || ...
           ~cc.Config.scgroups
            ne_gcfg.cc.(tstr).Config.scmarker = varargin{5};
        else
            ne_gcfg.cc.(tstr).Config.grpspecs{grp, 1} = varargin{5};
        end

    % scatter marker size
    case {'scattermarkersize'}
        
        % valid input
        if nargin < 5 || ...
           ~isa(varargin{5}, 'double') || ...
            numel(varargin{5}) ~= 1 || ...
            isinf(varargin{5}) || ...
            isnan(varargin{5}) || ...
            varargin{5} <= 0
            return;
        end

        % set size
        ne_gcfg.cc.(tstr).Config.scsize = varargin{5};
        
        % update menu
        tags.ScMarkerSize24.Checked = 'off';
        tags.ScMarkerSize36.Checked = 'off';
        tags.ScMarkerSize48.Checked = 'off';
        tags.ScMarkerSize60.Checked = 'off';
        tags.ScMarkerSize72.Checked = 'off';
        tags.ScMarkerSize96.Checked = 'off';
        tags.ScMarkerSize120.Checked = 'off';
        if any([24, 36, 48, 60, 72, 96, 120] == varargin{5})
            tags.(sprintf('ScMarkerSize%d', varargin{5})).Checked = 'on';
        end
        
    % toggle scatter reg quadratic
    case {'scatterquadratic'}
        if numel(grp) ~= 1 || ...
           ~cc.Config.scgroups
            if cc.Config.quadratic
                ne_gcfg.cc.(tstr).Config.scline = ~cc.Config.scline;
            else
                ne_gcfg.cc.(tstr).Config.scline = true;
            end
            if ne_gcfg.cc.(tstr).Config.scline
                tags.ScQuadratic.Checked = 'on';
                ne_gcfg.cc.(tstr).Config.scellipse = false;
                ne_gcfg.cc.(tstr).Config.quadratic = true;
                tags.ScLine.Checked = 'off';
                tags.ScConfEllipse.Checked = 'off';
            else
                tags.ScQuadratic.Checked = 'off';
            end
        else
            if cc.Config.grpspecs{grp, 3}(3)
                ne_gcfg.cc.(tstr).Config.grpspecs{grp, 3}(2) = ...
                    ~cc.Config.grpspecs{grp, 3}(2);
            else
                ne_gcfg.cc.(tstr).Config.grpspecs{grp, 3}(2:4) = ...
                    [true, true, false];
            end
        end

    % toggle show title and labels
    case {'showgeneralinfo'}
        ne_gcfg.cc.(tstr).Config.title = ~cc.Config.title;
        cini.BetaPlot.ShowTitle = ne_gcfg.cc.(tstr).Config.title;
        if ne_gcfg.cc.(tstr).Config.title
            tags.ShowGenInfo.Checked = 'on';
        else
            tags.ShowGenInfo.Checked = 'off';
        end
        ne_glmplotbetasrsz(0, 0, tstr);

    % toggle statsout
    case {'statsout'}
        if numel(grp) ~= 1 || ...
           ~cc.Config.scgroups
            ne_gcfg.cc.(tstr).Config.scstats = ~cc.Config.scstats;
            if ne_gcfg.cc.(tstr).Config.scstats
                tags.StatsOut.Checked = 'on';
            else
                tags.StatsOut.Checked = 'off';
            end
        else
            ne_gcfg.cc.(tstr).Config.grpspecs{grp, 3}(1) = ...
                ~cc.Config.grpspecs{grp, 3}(1);
        end

    % toggle sublabels
    case {'sublabels'}
        ne_gcfg.cc.(tstr).Config.sublabels = ~cc.Config.sublabels;
        if ne_gcfg.cc.(tstr).Config.sublabels
            tags.SubLabels.Checked = 'on';
        else
            tags.SubLabels.Checked = 'off';
        end

    % toggle sublines
    case {'sublines'}
        ne_gcfg.cc.(tstr).Config.plotlines = ~cc.Config.plotlines;
        if ne_gcfg.cc.(tstr).Config.plotlines
            tags.SubLines.Checked = 'on';
        else
            tags.SubLines.Checked = 'off';
        end

    % subject selection
    case {'subsel'}

        % get new subject selection
        [subsel, ssok] = listdlg( ...
            'PromptString', 'Please select subjects to use...', ...
            'ListString',   glm.Subjects, ...
            'InitialValue', cc.Config.subsel, ...
            'Name',         'NeuroElf - GLM beta plot request', ...
            'ListSize',     [min(200, max(640, 10 .* size(char(glm.Subjects), 2))), 300]);
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

    % switch type
    case {'type'}

        % set new type
        if cc.Tags.Type.Value == 1
            ne_gcfg.cc.(tstr).Config.type = 'bar';
            ne_gcfg.cc.(tstr).Config.scgroups = false;
            hFig.SetGroupEnabled('SGrps', 'on');
            tags.ScGroups.Checked = 'off';
        else
            ne_gcfg.cc.(tstr).Config.type = 'scatter';
            switch (cc.Tags.Type.Value)
                case {2}
                    ne_gcfg.cc.(tstr).Config.regress = 'none';
                case {3}
                    ne_gcfg.cc.(tstr).Config.regress = 'ols';
                case {4}
                    ne_gcfg.cc.(tstr).Config.regress = 'robust';
            end
            set(ax, 'XTick', [], 'XTickLabel', {});
            set(ax, 'XTickLabelMode', 'auto', 'XTickMode', 'auto');
        end

        % reset ranges
        ne_gcfg.cc.(tstr).Config.scrange = [-Inf, -Inf, Inf, Inf];
        ne_gcfg.cc.(tstr).Config.scrangea = [-1, -1, 1, 1];
        ne_gcfg.cc.(tstr).Config.xrange = [-Inf, Inf];
        tags.Axes.XLimMode = 'auto';
        tags.Axes.YLimMode = 'auto';

    % update range
    case {'xrange', 'yrange'}

        % get current and new range values
        cxlim = tags.Axes.XLim;
        cylim = tags.Axes.YLim;

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
            cxlim = tags.Axes.XLim;
            cylim = tags.Axes.YLim;
            if strcmpi(varargin{4}(:)', 'xrange')
                tags.Axes.XLim = cc.Config.scrangea([1, 3]);
                ne_gcfg.cc.(tstr).Config.scrange([1,3]) = [-Inf, Inf];
                ne_gcfg.cc.(tstr).Config.xrange = [-Inf, Inf];
                cxlim = tags.Axes.XLim;
            else
                ne_gcfg.cc.(tstr).Config.scrange([2,4]) = [-Inf, Inf];
                if strcmpi(cc.Config.type, 'bar')
                    tags.Axes.YLimMode = 'auto';
                else
                    tags.Axes.YLim = cc.Config.scrangea([2, 4]);
                end
                cylim = tags.Axes.YLim;
            end
        else
            ne_gcfg.cc.(tstr).Config.scrange = ...
                [cxlim(1), cylim(1), cxlim(2), cylim(2)];
            if strcmpi(varargin{4}(:)', 'xrange')
                tags.Axes.XLim = cxlim;
                ne_gcfg.cc.(tstr).Config.xrange = cxlim;
            else
                tags.Axes.YLim = cylim;
            end
        end

        % update text boxes
        tags.XFrom.String = sprintf('%.2f', cxlim(1));
        tags.XTo.String = sprintf('%.2f', cxlim(2));
        tags.YFrom.String = sprintf('%.2f', cylim(1));
        tags.YTo.String = sprintf('%.2f', cylim(2));
end

% reset X-range?
if ~any(strcmpi(varargin{4}(:)', ...
        {'backcolor'; 'condcolors'; 'contcolors'; 'ebarcolors'; ...
         'options'; 'ppts'; 'robust'; 'saveaxes'; ...
         'scattermarker'; 'showgeneralinfo'; 'xrange'; 'yrange'}))
    ne_gcfg.cc.(tstr).Config.xrange = [-Inf, Inf];
end

% return early
if nargin > 5 && ...
    numel(varargin{6}) == 1 && ...
    islogical(varargin{6}) && ...
   ~varargin{6}
    return;
end

% delete all handles
hf = fieldnames(ch);
for hc = 1:numel(hf)
    if ~any(strcmpi(hf{hc}, {'axes', 'barpos'})) && ...
       ~isempty(ch.(hf{hc}))
        try
            delete(ch.(hf{hc}));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        ch.(hf{hc}) = [];
    elseif strcmpi(hf{hc}, 'barpos')
        ch.(hf{hc}) = [];
    end
end
glmh.PlotHnd{pfi} = ch;

glm.SetHandle('PlotHnd', glmh.PlotHnd);

% title axes
if ne_gcfg.cc.(tstr).Config.title
    condval = tags.Conditions.Value;
    if strcmp(ne_gcfg.cc.(tstr).Config.type, 'bar')
        if isempty(condval)
            contval = tags.Contrasts.Value;
            if numel(contval) < 3
                connames = tags.Contrasts.String(contval);
            else
                connames = {sprintf('%d contrasts', numel(contval))};
            end
            ne_gcfg.cc.(tstr).Config.xlabel = 'Contrast';
        else
            if numel(condval) < 3
                connames = tags.Conditions.String(condval);
            else
                connames = {sprintf('%d conditions', numel(condval))};
            end
            ne_gcfg.cc.(tstr).Config.xlabel = 'Condition';
        end
        if ~isempty(connames)
            connames = strrep(sprintf('%s, ', connames{:}), '_', ' ');
            ne_gcfg.cc.(tstr).Config.titlestr = ...
                [ne_gcfg.cc.(tstr).Config.xlabel ' means: ' connames(1:end-2)];
        else
            ne_gcfg.cc.(tstr).Config.titlestr = ' ';
        end
    else
        if isempty(condval)
            contval = tags.Contrasts.Value;
            if ~isempty(contval)
                ne_gcfg.cc.(tstr).Config.titlestr = ...
                    strrep(sprintf('Scatter: %s ./. %s', ...
                    rtv.CovariatesNames{tags.Covariates.Value}, ...
                    tags.Contrasts.String{contval(1)}), '_', ' ');
            else
                ne_gcfg.cc.(tstr).Config.titlestr = ' ';
            end
        else
            ne_gcfg.cc.(tstr).Config.titlestr = ...
                strrep(sprintf('Scatter: %s ./. %s', ...
                rtv.CovariatesNames{tags.Covariates.Value}, ...
                tags.Conditions.String{condval(1)}), '_', ' ');
        end
        ne_gcfg.cc.(tstr).Config.xlabel = strrep( ...
            rtv.CovariatesNames{tags.Covariates.Value}, '_', '-');
    end
    switch glm.TransformationType
        case {1}
            ne_gcfg.cc.(tstr).Config.ylabel = 'z-scored betas';
        case {3}
            ne_gcfg.cc.(tstr).Config.ylabel = '% signal change';
        otherwise
            ne_gcfg.cc.(tstr).Config.ylabel = 'arbitraty units';
    end
else
    ne_gcfg.cc.(tstr).Config.titlestr = ' ';
    ne_gcfg.cc.(tstr).Config.xlabel = ' ';
    ne_gcfg.cc.(tstr).Config.ylabel = ' ';
end

% update
try
    ne_glmplotbetasup(0, 0, glm, 'fromdata', tstr);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
