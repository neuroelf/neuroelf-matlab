function varargout = ne_setcstatmap(varargin)
% ne_setcstatmap  - set current StatsVarIdx map(s)
%
% FORMAT:       ne_setcstatmap(SRC, EVT [, mapidx])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       mapidx      optional indices to set for current stats object
%
% No output fields.
%
% Example:
%
%     ne_setcstatmap(0, 0, 4);
%
%     sets the 4th map in the currently selected stats object as current
%     map.
%
% Note: this function issues ne_setslicepos() after settings the map; so
%       unless ne_gcfg.fcfg.noupdate == true, the GUI will update as well.

% Version:  v1.1
% Build:    16052922
% Date:     May-29 2016, 10:00 PM EST
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% clear index (in case we fail)
ne_gcfg.fcfg.StatsVarIdx = [];

% get handles
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;
cinist = ne_gcfg.c.ini.Statistics;
stvar = cc.StatsVar;

% hide stats texts
set(ch.CorStatsText, 'Visible', 'off', 'String', ' ');
set(ch.ZoomStatsText, 'Visible', 'off', 'String', ' ');

% disable SVC menu
set(ch.SVCEntries, 'Enable', 'off');

% and if we're supposed to be disabled
if strcmpi(ch.StatsVar.Enable, 'off') || ~isxff(stvar, true)

    % then set other controls also disabled and return
    ch.MainFig.SetGroupEnabled('SLoaded', 'off');
    ch.MainFig.SetGroupEnabled('SLdNVMP', 'off');

    % update render UI
    if isstruct(ne_gcfg.fcfg.Render) && isfield(ne_gcfg.fcfg.Render, 'hFig') && ...
        isxfigure(ne_gcfg.fcfg.Render.hFig, true)
        ne_gcfg.fcfg.Render.stvar = struct('FileType', 'NONE');
        ne_gcfg.fcfg.Render.stvix = [];
        ne_gcfg.fcfg.Render.stvixo = [];
        ne_gcfg.fcfg.Render.hFig.SetGroupEnabled('StVar', 'off');
        ne_gcfg.fcfg.Render.hFig.SetGroupEnabled('StVMP', 'off');
        ne_gcfg.fcfg.Render.hTag.ED_render_stvarfile.String = '<none>';
    end
    return;
end

% otherwise begin with default "loaded" configuration
ch.MainFig.SetGroupEnabled('SLoaded', 'on');
ch.MainFig.SetGroupEnabled('SLdNVMP', 'off');
ch.MainFig.SetGroupEnabled('SngVMP', 'off');

% input
stvix = [];
if nargin > 2 && isa(varargin{3}, 'double') && ...
   ~any(isinf(varargin{3}(:)) | isnan(varargin{3}(:)) | varargin{3}(:) < 1)
    stvix = round(varargin{3}(:));
end

% get index of dropdown
if isempty(stvix)
    stvix = ch.StatsVarMaps.Value;
end

% and set to configuration
stvarmnames = stvar.MapNames;
stvix(stvix < 1 | stvix > numel(stvarmnames)) = [];
cc.StatsVar.RunTimeVars.MapSelection = {stvarmnames(stvix), stvix(:)};
ne_gcfg.fcfg.StatsVarIdx = stvix;
ch.StatsVarMaps.Value = stvix;
nstvix = numel(stvix);

% stats text
stms = stvar.Map(stvix);
stlt = zeros(1, nstvix);
stut = zeros(1, nstvix);
stti = zeros(3, nstvix);
for umc = 1:nstvix
    stlt(umc) = stms(umc).LowerThreshold;
    stut(umc) = stms(umc).UpperThreshold;
    stti(:, umc) = [stms(umc).Type; stms(umc).DF1; stms(umc).DF2];
end
if ~isempty(stvix) && all(stlt == stlt(1)) && all(stut == stut(1)) && ...
    isequal(stti, repmat(stti(:, 1), 1, nstvix))
    sttyp = stti(1);
    stdf1 = stti(2);
    stdf2 = stti(3);
    sdtxt = sprintf('%%.%df', cinist.ShowThreshDecimals);
    switch (sttyp)
        case 1
            sttl = sprintf(['t_{[%d]}>=' sdtxt], stdf1, stlt(1));
            sttu = sprintf(['t_{[%d]}<' sdtxt], stdf1, stut(1));
        case 2
            sttl = sprintf(['r_{df=%d}>=' sdtxt], stdf1, stlt(1));
            sttu = sprintf(['r_{df=%d}<' sdtxt], stdf1, stut(1));
        case 4
            sttl = sprintf(['F_{[%d,%d]}>=' sdtxt], stdf1, stdf2, stlt(1));
            sttu = sprintf(['F_{[%d,%d]}<' sdtxt], stdf1, stdf2, stut(1));
        otherwise
            sttl = '';
    end
    if ~isempty(sttl)
        set(ch.CorStatsText(1), 'String', sttl);
        set(ch.CorStatsText(2), 'String', sttu);
        set(ch.ZoomStatsText(1), 'String', sttl);
        set(ch.ZoomStatsText(2), 'String', sttu);
        if cc.page == 1 && cinist.ShowThreshBars && cinist.ShowThreshText
            set(ch.CorStatsText, 'Visible', 'on');
        elseif cc.page == 2 && cinist.ShowThreshBars && cinist.ShowThreshText
            set(ch.ZoomStatsText, 'Visible', 'on');
        end
    end
end

% if multiple are selected
ne_gcfg.c.title{1, 2} = ch.StatsVar.String{ch.StatsVar.Value};
if numel(stvix) ~= 1

    % udpate title
    ne_gcfg.c.title{1, 3} = sprintf('(%d maps)', numel(stvix));
    ne_updatename;

    % then we won't update any controls...
    ne_setslicepos;
    return;
end

% we definitely have one map only!
ch.MainFig.SetGroupEnabled('SngVMP', 'on');
ch.MainFig.SetGroupEnabled('SLdNVMP', 'off');

% re-read the configuration
cc = ne_gcfg.fcfg;
sttyp = lower(stvar.Filetype);

% also enable SVC?
if strcmp(sttyp, 'vmp') && isfield(stvar.Map(stvix), 'RunTimeVars') && ...
    numel(stvar.Map(stvix).RunTimeVars) == 1 && ...
    isstruct(stvar.Map(stvix).RunTimeVars) && ...
    isfield(stvar.Map(stvix).RunTimeVars, 'FWHMResEst') && ...
   ~isempty(stvar.Map(stvix).RunTimeVars.FWHMResEst) && ...
    isfield(stvar.Map(stvix).RunTimeVars, 'FWHMResImg') && ...
    isequal(size(stvar.Map(stvix).VMPData), size(stvar.Map(stvix).RunTimeVars.FWHMResImg))
    set(ch.SVCEntries, 'Enable', 'on');
end

% update figure title
ne_gcfg.c.title{1, 3} = ch.StatsVarMaps.String{ch.StatsVarMaps.Value};
ne_updatename(0, 0, false);

% as default, we use LUT coloring scheme
ch.Stats.UseLUT.RadioGroupSetOne;

% the rest depends on filetype
switch (sttyp)

    % for VMPs
    case {'ava', 'cmp', 'glm', 'hdr', 'head', 'vmp', 'vtc'}

        % do we need clustering ?
        if strcmp(sttyp, 'vmp') && isempty(stvar.Map(stvix).VMPDataCT) && ...
            stvar.Map(stvix).EnableClusterCheck > 0
            stvar.ClusterTable(stvix, []);
        end

        % get Map shortcut handle
        stmap = stvar.Map(stvix);

        % if invalid
        if isempty(stmap.LowerThreshold) || isempty(stmap.UpperThreshold)

            % slice data (forces thresolds)
            ntio = transimg(16, 16);
            stvar.SliceToTransimg([0, 0, 0], ntio, struct( ...
                'frame', [-8, -8, -8; 7.99, 7.99, 7.99], 'dir', 'sag', ...
                'mapvol', stvix, 'type', 'rgb'));

            % then delete transimg
            delete(ntio);

            % and re-get map
            stmap = stvar.Map(stvix);
        end

        % and make sure to update the enabled flags
        ch.MainFig.SetGroupEnabled('SLdNVMP', 'on');

        % set the thresholds in the figure config
        ne_gcfg.fcfg.StatsVarkThr = stmap.ClusterSize;
        ch.Stats.PosTail.Value = double(mod(stmap.ShowPositiveNegativeFlag, 2) > 0);
        ch.Stats.NegTail.Value = double(stmap.ShowPositiveNegativeFlag > 1);
        spvals = cellfun(@str2double, ch.Stats.PThresh.String);

        % and also set the parameters of the stats
        tvals = [];
        switch (stmap.Type)
            case 1
                ne_gcfg.fcfg.StatsVarPar = {'t', stmap.DF1, 0};
                if stmap.ShowPositiveNegativeFlag == 3
                    tvals = -sdist('tinv', 0.5 .* spvals, stmap.DF1);
                else
                    tvals = -sdist('tinv', spvals, stmap.DF1);
                end
            case 4
                ne_gcfg.fcfg.StatsVarPar = {'F', stmap.DF1, stmap.DF2};
                tvals = sdist('finv', spvals, stmap.DF1, stmap.DF2, true);
            case 2
                ne_gcfg.fcfg.StatsVarPar = {'r', stmap.DF1, 0};
                if stmap.ShowPositiveNegativeFlag == 3
                    tvals = correlinvtstat(-sdist('tinv', 0.5 .* spvals, stmap.DF1), stmap.DF1);
                else
                    tvals = correlinvtstat(-sdist('tinv', spvals, stmap.DF1), stmap.DF1);
                end
            case 9
                ne_gcfg.fcfg.StatsVarPar = {'m', stmap.DF1, 0};
            case 30
                ne_gcfg.fcfg.StatsVarPar = {'a', stmap.DF1, 0};
            otherwise
                ne_gcfg.fcfg.StatsVarPar{1} = '!';
        end
        if ~isempty(tvals)
            tvali = findfirst(abs(tvals - stmap.LowerThreshold) <= 0.0005);
            if ~isempty(tvali)
                ch.Stats.PThresh.Value = tvali;
            end
        end

        % update controls on figure -> cluster threshold (size and status)
        set(ch.Stats.kThresh, 'String', sprintf('%d', ...
            max(1, floor(stmap.ClusterSize))));
        ch.Stats.UsekThr.Value = double(stmap.EnableClusterCheck ~= 0);

        % LUT/RGB choice
        if stmap.UseRGBColor ~= 0
            ch.Stats.UseRGB.RadioGroupSetOne;
        end

        % colors (buttons for colorpicker)
        bcolor = min(255, max(0, stmap.RGBLowerThreshPos(:)));
        ch.Stats.RGBLPos.BackgroundColor = (1 / 255) .* bcolor';
        bcolor = min(255, max(0, stmap.RGBUpperThreshPos(:)));
        ch.Stats.RGBUPos.BackgroundColor = (1 / 255) .* bcolor';
        bcolor = min(255, max(0, stmap.RGBLowerThreshNeg(:)));
        ch.Stats.RGBLNeg.BackgroundColor = (1 / 255) .* bcolor';
        bcolor = min(255, max(0, stmap.RGBUpperThreshNeg(:)));
        ch.Stats.RGBUNeg.BackgroundColor = (1 / 255) .* bcolor';

    % for AVA files, etc.
    otherwise

        % get thresholds for current map
        try
            thr = stvar.RunTimeVars.Thresholds;
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            thr = nan(1, 5);
        end
        thrl = min(size(thr, 1), stvix);
        thrv = thr(thrl, :);

        % if invalid
        if any(isnan(thrv))

            % slice data (forces thresolds)
            ntio = transimg(16, 16);
            stvar.SliceToTransimg([0, 0, 0], ntio, struct( ...
                'frame', [-8, -8, -8; 7.99, 7.99, 7.99], 'dir', 'sag', ...
                'mapvol', stvix, 'type', 'rgb'));

            % then reget thr
            delete(ntio);
        end
end

% re-read configuration (after changes from VMP)
cc = ne_gcfg.fcfg;

% but if not GLM
if ~strcmpi(stvar.Filetype, 'glm') || stvar.ProjectTypeRFX == 0
    ch.StatsVarProject.Enable = 'off';
else
    ch.StatsVarProject.Enable = 'on';
end

% set thresholds to text controls
ch.Stats.LThresh.String = sprintf('%.4f', stmap.LowerThreshold);
ch.Stats.UThresh.String = sprintf('%.4f', stmap.UpperThreshold);

% update render page
if isstruct(ne_gcfg.fcfg.Render) && ...
    isfield(ne_gcfg.fcfg.Render, 'hFig') && ...
    isxfigure(ne_gcfg.fcfg.Render.hFig, true)

    % get name
    stvarfile = stvar.FilenameOnDisk(2);
    if isempty(stvarfile)
        stvarfile = sprintf('<untitled.%d>', lower(stvar.Filetype));
    elseif numel(stvarfile) > 51
        stvarfile = [stvarfile(1:24) '...' stvarfile(end-23:end)];
    end
    hRendFig = ne_gcfg.fcfg.Render.hFig;
    hRendFig.SetGroupEnabled('StVMP', 'off');
    hRendFig.SetGroupEnabled('StVar', 'on');
    if isxff(stvar, {'cmp', 'hdr', 'head', 'vmp'})
        hRendFig.SetGroupEnabled('StVMP', 'on');
    end
    ne_gcfg.fcfg.Render.hTag.ED_render_stvarfile.String = ['  ' stvarfile];

    % make settings
    ne_gcfg.fcfg.Render.stvar = stvar;
    ne_gcfg.fcfg.Render.stvix = cc.StatsVarIdx;
    ne_gcfg.fcfg.Render.stvixo = cc.StatsVarIdx;

    % update UI
    if ne_gcfg.fcfg.page == 4
        ne_render_setview;
        return;
    end
end

% udpate screen
ne_setslicepos;
