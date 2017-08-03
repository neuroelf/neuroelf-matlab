% FUNCTION ne_mdmvoicondavgup: update MDM VOI condition average plot
function ne_mdmvoicondavgup(varargin)

% Version:  v1.1
% Build:    16021717
% Date:     Feb-17 2016, 5:13 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, 2015, 2016, Jochen Weber
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

% valid call?
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~isfield(ne_gcfg.cc, varargin{3}(:)')
    return;
end

% only once
tstr = varargin{3}(:)';
if any(strcmp(ne_gcfg.c.blockcb, [tstr '_up']))
    return;
end
ne_gcfg.c.blockcb{end+1} = [tstr '_up'];

% get configuration
cc = ne_gcfg.cc.(tstr).Config;
ax = cc.ax;
tags = ne_gcfg.cc.(tstr).Tags;
voiidx = tags.VOI.Value;
vn = tags.VOI.String(voiidx);
cval = tags.Conditions.Value;
clval = tags.Collapsed.Value;
ncl = numel(clval);
if (numel(cval) + ncl) == 2
    tags.Diff.Enable = 'on';
    condstr = tags.Conditions.String;
    collstr = tags.Collapsed.String;
    if ~iscell(condstr)
        condstr = cellstr(condstr);
    end
    if ~iscell(collstr)
        collstr = cellstr(collstr);
    end
    if ncl == 0
        tags.Diff1.Label = sprintf('%s - %s', condstr{cval(1)}, condstr{cval(2)});
        tags.Diff2.Label = sprintf('%s - %s', condstr{cval(2)}, condstr{cval(1)});
    elseif ncl == 1
        tags.Diff1.Label = sprintf('%s - %s', condstr{cval(1)}, collstr{clval(1)});
        tags.Diff2.Label = sprintf('%s - %s', collstr{clval(1)}, condstr{cval(1)});
    else
        tags.Diff1.Label = sprintf('%s - %s', collstr{clval(1)}, collstr{clval(2)});
        tags.Diff2.Label = sprintf('%s - %s', collstr{clval(2)}, collstr{clval(1)});
    end
    if tags.Type.Value < 2
        ne_gcfg.cc.(tstr).Config.diffmode = 0;
        cc.diffmode = 0;
    end
    tags.Diff1.Checked = 'off';
    tags.Diff2.Checked = 'off';
    if cc.diffmode == 1
        tags.Diff1.Checked = 'on';
    elseif cc.diffmode == 2
        tags.Diff2.Checked = 'on';
    end
else
    ne_gcfg.cc.(tstr).Config.diffmode = 0;
    cc.diffmode = 0;
    tags.Diff.Enable = 'off';
end
groups = tags.Groups.Value;
if ~isempty(groups)
    try
        mdmgrp = cc.mdm.RunTimeVars.Groups;
        groups = cat(1, mdmgrp{groups, 2})';
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        groups = [];
    end
end
robust = (tags.Robust.Value > 0);

% delete stuff
try
    delete(get(ax, 'Children'))
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
ne_gcfg.cc.(tstr).Config.bands = [];
ne_gcfg.cc.(tstr).Config.curves = [];
ne_gcfg.cc.(tstr).Config.legends = [];
ne_gcfg.cc.(tstr).Config.lines = [];
ne_gcfg.cc.(tstr).Config.texts = [];
curves = [];
ccurves = [];
leg = [];

% fix background color
set(ax, 'Color', (1 / 255) .* cc.axcol);

% leave already
if isempty(cval) && ...
    isempty(clval)
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, [tstr '_up'])) = [];
    return;
end

% construct arguments
ccols = (1 / 255) .* cc.condcol;
spa = cc.sealpha;

% no grouping?
if isempty(groups)
    groups = 1:size(cc.mtc, 4);
else
    groups = intersect(groups(:), (1:size(cc.mtc, 4))')';
end

% apply subject masking
groups = intersect(groups(:), cc.subsel(:));
gtype = tags.Type.Value;
if isempty(groups)
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, [tstr '_up'])) = [];
    return;
elseif numel(groups) == 1
    gtype = 1;
end

% smoothing kernel (for weighting)
smks = mean(cc.tr) / cc.avgopt.samptr;
smk = smoothkern(smks, 1e-6 * smks * smks);

% figure out which conditions are needed at all
cneeded = false(1, size(cc.mtc, 2));
cneeded(cval) = true;
for c = 1:ncl
    cneeded(cc.colls{clval(c), 2}) = true;
end
cneededb = zeros(1, numel(cneeded));
cneededb(cneeded) = 1:sum(cneeded);

% for now, trial-based baseline
nmtc = cc.mtc(:, cneeded, voiidx, groups, :) - ...
    cc.btc(ones(1, size(cc.mtc, 1)), cneeded, voiidx, groups, :);

% smoothing?
if cc.tsmooth > 0
    smtc = size(nmtc);
    nmtc = reshape(nmtc, smtc(1), prod(smtc(2:end)));
    str = 0.001 * cc.avgopt.samptr;
    nmtc = reshape(flexinterpn(nmtc, [Inf, Inf; 1, 1; 1, 1; size(nmtc)], ...
        {smoothkern(cc.tsmooth / str, 0, false, 'lanczos8'), [0; 1; 0]}, ...
        {1, 1}, 0), smtc);
end

% create weighting matrix
wnmtc = double(~isnan(nmtc));

% weighted (or robust)
if gtype == 3 || ...
    robust

    % compute variance
    vmtc = varc(nmtc, 1, true);
    vmtc(isinf(vmtc) | isnan(vmtc) | vmtc == 0) = 10;
    vmtc = min(1, (1 / sdist('normpdf', 0, 0, 1)) .* sdist('normpdf', vmtc, 0, 1));

    % down-weight where variance across time is too high
    if ~isempty(wnmtc)
        wnmtc = wnmtc .* repmat(vmtc, size(nmtc, 1), 1);
    end
end

% smooth with 1TR kernel
if ~isempty(wnmtc)
    smks = ceil(smks);
    smtc = [size(wnmtc), 1, 1, 1];
    wnmtc(end+1:end+smks, :, :, :, :) = repmat(wnmtc(end, :, :, :, :), smks, 1);
    wnmtc = reshape(flexinterpn(reshape(wnmtc, smtc(1) + smks, prod(smtc(2:end))), ...
        [Inf, Inf; 1 + smks, 1; 1, 1; smtc(1) + smks, prod(smtc(2:end))], ...
        {smk, [0;1;0]}, {1,1}), smtc(1:5));
end

% get required dimensions
dimt = size(nmtc, 1);
dims = size(nmtc, 4);
dimo = size(nmtc, 5);
avwin = cc.avgopt.avgbegin / 1000 + ((cc.avgopt.samptr / 1000) .* (0:(dimt-1))');

% no difference
if cc.diffmode == 0

    % if basic conditions
    if ~isempty(cval)

        % use function
        [mtc, mtcse, curves] = ne_mdmgrouping( ...
            nmtc(:, cneededb(cval), :, :, :), ...
            wnmtc(:, cneededb(cval), :, :, :), ...
            ax, avwin, ccols(cval, :), gtype, robust, ...
            tags.ECurves.Value > 0, tags.SDSE.Value > 0, cc.p05);
    end

    % loop over condidions
    lines = zeros(numel(cval) + numel(clval), 1);
    bands = lines;
    for c = 1:numel(cval)

        % use tcplot to plot into the same axes
        [lines(c + ncl), bands(c + ncl)] = ...
            tcplot(ax, avwin, mtc(:, c, 1), mtcse(:, c, 1), mtcse(:, c, 1), ...
            struct('color', ccols(cval(c), :), 'lwidth', 2.5, 'spalpha', spa, 'spline', false));
    end

    % loop over collapsings
    ccurves = repmat({zeros(0, 1)}, ncl, 1);
    for c = 1:ncl

        % get condition indices and color
        clidx = cc.colls{clval(c), 2};
        clcol = (1 / 255) .* cc.colls{clval(c), 3};

        % combine conditions in 5th dimension
        mtc = reshape(permute(nmtc(:, cneededb(clidx), :, :, :), ...
            [1, 3, 4, 5, 2]), [dimt, 1, 1, dims, dimo * numel(clidx)]);
        wmtc = reshape(permute(wnmtc(:, cneededb(clidx), :, :, :), ...
            [1, 3, 4, 5, 2]), size(mtc));

        % use function
        [mtc, mtcse, ccurves{c}] = ne_mdmgrouping(mtc, wmtc, ...
            ax, avwin, clcol, gtype, robust, ...
            tags.ECurves.Value > 0, tags.SDSE.Value > 0, cc.p05);

        % plot
        [lines(c), bands(c)] = ...
            tcplot(ax, avwin, mtc(:, 1), mtcse(:, 1), mtcse(:, 1), ...
            struct('color', clcol, 'lwidth', 2.5, 'spalpha', spa, 'spline', false));
    end
    if ~isempty(ccurves)
        ccurves = cat(1, ccurves{:});
    else
        ccurves = zeros(0, 1);
    end

    % legend texts
    condnames = strrep(cc.conds(cval), '_', ' ');
    collnames = strrep(cc.colls(clval, 1), '_', ' ');
    legnames = [collnames(:); condnames(:)];

% difference mode
else

    % two conditions
    if ncl == 0

        % get two conditions' data
        mtc1 = nmtc(:, cneededb(cval(1)), :, :, :);
        wmtc1 = wnmtc(:, cneededb(cval(1)), :, :, :);
        mtc2 = nmtc(:, cneededb(cval(2)), :, :, :);
        wmtc2 = wnmtc(:, cneededb(cval(2)), :, :, :);
        ccols = ccols(cval, :);
        if cc.diffmode == 1
            legnames = {sprintf('%s - %s', cc.conds{cval(1)}, cc.conds{cval(2)})};
        else
            legnames = {sprintf('%s - %s', cc.conds{cval(2)}, cc.conds{cval(1)})};
        end

    elseif ncl == 1

        % get condition's data
        mtc1 = nmtc(:, cneededb(cval), :, :, :);
        wmtc1 = wnmtc(:, cneededb(cval), :, :, :);

        % get collapsed data
        clidx = cc.colls{clval, 2};

        % combine conditions in 5th dimension
        mtc2 = reshape(permute(nmtc(:, cneededb(clidx), :, :, :), ...
            [1, 3, 4, 5, 2]), [dimt, 1, 1, dims, dimo * numel(clidx)]);
        wmtc2 = reshape(permute(wnmtc(:, cneededb(clidx), :, :, :), ...
            [1, 3, 4, 5, 2]), size(mtc2));

        % set colors
        ccols = [ccols(cval, :); (1 / 255) .* cc.colls{clval, 3}];

        % create name
        if cc.diffmode == 1
            legnames = {sprintf('%s - %s', cc.conds{cval}, cc.colls{clval, 1})};
        else
            legnames = {sprintf('%s - %s', cc.colls{clval, 1}, cc.conds{cval})};
        end

    else

        % get collapsed data
        clidx1 = cc.colls{clval(1), 2};
        clidx2 = cc.colls{clval(2), 2};
        mtc1 = reshape(permute(nmtc(:, cneededb(clidx1), :, :, :), ...
            [1, 3, 4, 5, 2]), [dimt, 1, 1, dims, dimo * numel(clidx1)]);
        wmtc1 = reshape(permute(wnmtc(:, cneededb(clidx1), :, :, :), ...
            [1, 3, 4, 5, 2]), size(mtc1));
        mtc2 = reshape(permute(nmtc(:, cneededb(clidx2), :, :, :), ...
            [1, 3, 4, 5, 2]), [dimt, 1, 1, dims, dimo * numel(clidx2)]);
        wmtc2 = reshape(permute(wnmtc(:, cneededb(clidx2), :, :, :), ...
            [1, 3, 4, 5, 2]), size(mtc2));
        ccols = (1 / 255) .* [cc.colls{clval(1), 3}; cc.colls{clval(2), 3}];
        if cc.diffmode == 1
            legnames = {sprintf('%s - %s', ...
                cc.colls{clval(1), 1}, cc.colls{clval(2), 1})};
        else
            legnames = {sprintf('%s - %s', ...
                cc.colls{clval(2), 1}, cc.colls{clval(1), 1})};
        end
    end

    % and use function
    [mtc, mtcse, curves, dcol] = ne_mdmgrouping(mtc1, wmtc1, ...
        ax, avwin, ccols, gtype, robust, ...
        tags.ECurves.Value > 0, tags.SDSE.Value > 0, cc.p05, ...
        mtc2, wmtc2, cc.diffmode);

    % use tcplot to plot into the same axes
    [lines, bands] = ...
        tcplot(ax, avwin, mtc, mtcse, mtcse, ...
        struct('color', dcol, 'lwidth', 2.5, 'spalpha', spa, 'spline', false));

    % good names
    legnames = strrep(legnames, '_', ' ');
end

% reorder children
if ~isempty(bands)
    if ne_gcfg.c.mlversion < 804
        set(bands, 'ZData', 0.5 .* ones(numel(get(bands(1), 'XData')), 1));
        set(lines, 'ZData', ones(numel(get(lines(1), 'XData')), 1));
    end
    set(ax, 'Children', [lines(:); ccurves; curves(:); bands(:)]);
    set(ax, 'ZLim', [-1, 2]);
end

% add legend
if ~isempty(lines)
    [leg, leghnd] = legend(ax, lines, legnames{:}, 'Location', cc.legloc);

    % adapt fontsize
    legtxt = leghnd(strcmpi(get(leghnd, 'Type'), 'text'));
    set(legtxt, 'Fontsize', cini.Satellites.FontSize);
    axcol = get(ax, 'Color');
    if isa(axcol, 'double') && ...
        numel(axcol) == 3 && ...
        sum(axcol) < 1.2
        set(legtxt, 'Color', [1, 1, 1]);
    end
end

% set title
txt = title(ax, strrep(vn{1}, '_', ' '));

% update handles
ne_gcfg.cc.(tstr).Config.bands = bands;
ne_gcfg.cc.(tstr).Config.curves = [ccurves; curves(:)];
ne_gcfg.cc.(tstr).Config.legends = leg;
ne_gcfg.cc.(tstr).Config.lines = lines;
ne_gcfg.cc.(tstr).Config.texts = txt;

% update range text boxes
axp = get(ax);
tags.XFrom.String = sprintf('%.2f', axp.XLim(1));
tags.XTo.String = sprintf('%.2f', axp.XLim(2));
tags.YFrom.String = sprintf('%.2f', axp.YLim(1));
tags.YTo.String = sprintf('%.2f', axp.YLim(2));

% clear block
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, [tstr '_up'])) = [];



% % % SUB FUNCTIONS % % %



function [mtc, mtcse, curves, ccols] = ne_mdmgrouping(mtc, wmtc, ax, avwin, ccols, gtype, robust, pcurves, sdse, p05, mtc2, wmtc2, d12)

% preset curves
curves = [];

% get size of input
smtc = [size(mtc), 1, 1, 1];

% FFX
if gtype == 1

    % remove subject as a factor
    mtc = reshape(mtc, [smtc(1:3), smtc(4) * smtc(5)]);
    wmtc = reshape(wmtc, size(mtc));

    % remove onsets without data
    owd = squeeze(all(all(all(isnan(mtc), 1), 2), 3));
    mtc(:, :, :, owd) = [];
    wmtc(:, :, :, owd) = [];

    % robust computation
    mtc = wmtc .* mtc;
    if robust
        mtc = replacerobmean(mtc, 4, true);
    end

    % compute SE measure
    mtcse = sqrt(varc(mtc, 4, true));

    % plot curves
    if pcurves
        curves = ne_mdmplotcurves(ax, avwin, mtc, ccols);
    end

    % build averages
    ges = sum(wmtc, 4);
    mtc(isnan(mtc)) = 0;
    mtc = sum(mtc, 4) ./ ges;

% RFX
else

    % build averages within subjects first
    mtc = mtc .* wmtc;
    wmtc = sum(wmtc, 5);
    mtc(isnan(mtc)) = 0;
    mtc = sum(mtc, 5) ./ wmtc;

    % for further computation all-0 responses are NaN!
    mtc(:, all(mtc(:, :) == 0, 1)) = NaN;

    % robust computation
    if robust
        mtc = replacerobmean(mtc, 4, true);
    end

    % differences mode
    if nargin > 12

        % also perform computation on second set
        mtc2 = mtc2 .* wmtc2;
        wmtc2 = sum(wmtc2, 5);
        mtc2(isnan(mtc2)) = 0;
        mtc2 = sum(mtc2, 5) ./ wmtc2;

        % for further computation all-0 responses are NaN!
        mtc2(:, all(mtc2(:, :) == 0, 1)) = NaN;

        % robust computation
        if robust
            mtc2 = replacerobmean(mtc2, 4, true);
        end

        % depending on direction
        if d12 == 1
            mtc = mtc - mtc2;
            ccols = limitrangec(ccols(1, :) - ccols(2, :), 0, 1, 0);
        else
            mtc = mtc2 - mtc;
            ccols = limitrangec(ccols(2, :) - ccols(1, :), 0, 1, 0);
        end

        % for weighted RFX
        if gtype == 3
            wmtc = 1 ./ (1 ./ sqrt(wmtc) + 1 ./ sqrt(wmtc2));
            wmtc(isinf(wmtc) | isnan(wmtc)) = 0;
        end
    end

    % compute SE measures
    mtcse = sqrt(varc(mtc, 4, true));

    % plot curves
    if pcurves > 0
        curves = ne_mdmplotcurves(ax, avwin, mtc, ccols);
    end

    % build averages along subject dimension
    be = (isnan(mtc) | mtc == 0);
    mtc(be) = 0;
    if gtype < 3
        ges = sum(~be, 4);
        mtc = sum(mtc, 4) ./ ges;
    else
        wmtc = sqrt(wmtc);
        ges = sum(wmtc, 4) ./ max(wmtc(:));
        mtc = sum(wmtc .* mtc, 4) ./ sum(wmtc, 4);
    end
end

% SE instead of SD
if sdse

    % actual CI
    if p05
        mtcse = (1.96 ./ sqrt(ges)) .* mtcse;

    % just SE
    else
        mtcse = (1 ./ sqrt(ges)) .* mtcse;
    end
end

function curves = ne_mdmplotcurves(ax, avwin, mtc, ccols)
curves = cell(1, size(mtc, 2));
axcol = get(ax, 'Color');
if isa(axcol, 'double') && ...
    numel(axcol) == 3
    ccols = 0.5 .* (ccols + axcol(ones(1, size(ccols, 1)), :));
end
for c = 1:numel(curves)
    curves{c} = plot(ax, avwin, squeeze(mtc(:, c, 1, ~isnan(mtc(1, c, 1, :)))));
    set(curves{c}, 'Color', ccols(c, :), 'LineWidth', 0.25);
    curves{c} = curves{c}(:);
end
curves = cat(1, curves{:});
