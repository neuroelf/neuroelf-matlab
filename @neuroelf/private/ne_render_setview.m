function varargout = ne_render_setview(varargin)
% ne_render_setview  - update render page view
%
% FORMAT:       ne_render_setview([SRC, EVT, window, viewpt])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       window      window specifier (used to target satellite windows)
%       viewpt      1x5 double with
%        (1)        1x1 double zenith angle (between -90 and 90)
%        (2)        1x1 double azimuth (around Z) angle (180 = left hemi)
%        (3:4)      1x2 translation
%        (5)        1x1 zoom
%
% No output fields. (will be set to [])
%
% Example:
%
%     neuroelf_gui('render_setview', 'main', [0, 180, -20, 0, 1.2]);
%
%     this sets the render (camera) viewpoint in the main GUI window to
%     a left-hemisphere-view with a [-20, 0] translation and a zoom factor
%     of 1.2

% Version:  v1.1
% Build:    16052611
% Date:     May-26 2016, 11:16 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2015, 2016, Jochen Weber
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

% render not open
if ~isfield(ne_gcfg.fcfg, 'Render') || ...
   ~isstruct(ne_gcfg.fcfg.Render) || ...
   ~isfield(ne_gcfg.fcfg.Render, 'roty')
    try
        ne_render;
    catch ne_eo;
        rethrow(ne_eo);
    end
end

% get handles
if nargin < 3 || ~ischar(varargin{3}) || ~isfield(ne_gcfg.cc, varargin{3})
    fromroot = true;
    ch = ne_gcfg.h;
    chfig = ch.MainFig;
    pc = ne_gcfg.fcfg;
    fPos = [ne_gcfg.fcfg.surfpos; ne_gcfg.fcfg.histpos; ne_gcfg.fcfg.tcpos];
    tio = ne_gcfg.tio.imRnd;
    xPos = get(ne_gcfg.h.Surface, 'Position');

    % TC undocked?
    ccf = fieldnames(ne_gcfg.cc);
    if ~isempty(ccf)
        ccf(cellfun('isempty', regexp(ccf, '^TC'))) = [];
    end
    if ~isempty(ccf) && numel(ne_gcfg.cc.(ccf{1})) == 1 && isstruct(ne_gcfg.cc.(ccf{1})) && ...
        isfield(ne_gcfg.cc.(ccf{1}), 'Satellite') && isxfigure(ne_gcfg.cc.(ccf{1}).Satellite, true)
        tcpct = ccf{1};
        tcpch = ne_gcfg.cc.(ccf{1});
        ch.TCPlot = tcpch.TCPlot;
        ch.TCPlotChild = tcpch.TCPlotChild;
        ch.TCPlotChildren = tcpch.TCPlotChildren;
        ch.TCPlotDiscards = tcpch.TCPlotDiscards;
    else
        tcpct = '';
    end

    % hide timecourse plot (unless, see below)
    ch.TCPlot.Visible = 'off';
    ch.MainFig.SetGroupEnabled('TCPlot', 'off');
    if ~isempty(tcpct)
        tcpch.Satellite.Visible = 'off';
    end
else
    fromroot = false;
    ch = ne_gcfg.cc.(varargin{3});
    chfig = ch.Satellite;
    pc = ch.Config;
    fPos = ch.Config.surfpos;
    tio = pc.tio.imRnd;
    xPos = get(ne_gcfg.cc.(varargin{3}).Surface, 'Position');
    tcpct = '';
end
xPoshsz = floor(0.5 .* xPos(3:4));
xPosmid = xPos(1:2) + xPoshsz;

% get mouse position
nPos = chfig.CurrentPoint;

% get ini settings
ci = ne_gcfg.c.ini.Render;

% position of controls on which clicks allows updating of position
cc = pc.srfcfg;

% render config
rc = ne_gcfg.fcfg.Render;

% currently configured position
cpos = [cc.angley, cc.anglex, cc.trans(2:3), cc.zoom];

% update with input
if nargin > 3 && isa(varargin{4}, 'double') && numel(varargin{4}) == 5 && ...
   ~any(isinf(varargin{4}) | isnan(varargin{4}))
    cpos = varargin{4}(:)';
    if fromroot
        ne_gcfg.fcfg.srfcfg.angley = cpos(1);
        ne_gcfg.fcfg.srfcfg.anglex = cpos(2);
        ne_gcfg.fcfg.srfcfg.trans = [0, cpos(3:4)];
        ne_gcfg.fcfg.srfcfg.zoom = cpos(5);
        pc = ne_gcfg.fcfg;
    else
        ne_gcfg.cc.(varargin{3}).Config.srfcfg.angley = cpos(1);
        ne_gcfg.cc.(varargin{3}).Config.srfcfg.anglex = cpos(2);
        ne_gcfg.cc.(varargin{3}).Config.srfcfg.trans = [0, cpos(3:4)];
        ne_gcfg.cc.(varargin{3}).Config.srfcfg.zoom = cpos(5);
        ch = ne_gcfg.cc.(varargin{3});
        pc = ch.Config;
    end
    cc = pc.srfcfg;
end

% get variable that is to be sliced
svar = pc.SliceVar;
if ~isxff(svar, true)
    svar = struct('FileType', '', 'RunTimeVars', struct, 'Handles', struct, ...
        'BoundingBox', struct('Dim', [256, 256, 256]));
else
    if ~isfield(svar.RunTimeVars, 'AlphaTable')
        atable = (0:(1/255):1)';
        atable = 1 - (1 - atable) .^ (2 .* atable);
        svar.RunTimeVars.AlphaTable = atable;
    end
    if ~isfield(svar.RunTimeVars, 'SliceRanges')
        svar.RunTimeVars.SliceRanges = sliceranges(svar.GetVolume(1));
    end
end

% get some more shorthands
svartyp = svar.FileType;
rtv = svar.RunTimeVars;
stvar = pc.StatsVar;
svi = pc.StatsVarIdx;
if fromroot && isxff(stvar, 'vtc') && ~isempty(svi)
    ch.TCPlot.Visible = 'on';
    ch.MainFig.SetGroupEnabled('TCPlot', 'on');
    if ~isempty(tcpct)
        tcpch.Satellite.Visible = 'on';
    end
end

% adapt bounding box?
if nargin > 4 && isa(varargin{5}, 'double') && isequal(size(varargin{5}), [3, 2]) && ...
   ~any(isinf(varargin{5}(:)) | isnan(varargin{5}(:))) && ...
   ~any(varargin{5}(:, 1) < 1 | varargin{5}(:, 2) > rtv.RenderBBoxFull(:, 2) | varargin{5}(:, 1) >= varargin{5}(:, 2))
    svar.RunTimeVars.RenderBBox = round(varargin{5});
    rtv = svar.RunTimeVars;
end

% do we need a hit-test
allowpreview = true;
warpip = '';
if any(ne_gcfg.c.btdown == gcbf) && isempty(pc.mpos.ddat)

    % make hit-test
    cobj = findfirst( ...
        fPos(:, 1) <= nPos(1) & fPos(:, 2) <= nPos(2) & ...
        fPos(:, 3) >  nPos(1) & fPos(:, 4) >  nPos(2));

    % object not hit?
    if isempty(cobj) || cobj == 2

        % pretend it's not down (no further updates)
        ne_gcfg.c.btdown = [];

        % for the histogram (min/max) sliders
        if ~isempty(cobj) && ...
            any(strcmpi(svartyp, {'fmr', 'hdr', 'head', 'vmr', 'vtc'}))

            % y axes and RunTimeVars
            cp = nPos - fPos(cobj, 1:2);
            hpos = cp(2) / 255;
            if strcmp(svartyp, 'vmr') && ...
                rtv.ShowV16
                swl = rtv.ScalingWindowLim16;
            else
                swl = rtv.ScalingWindowLim;
            end

            % determine which boundary
            if pc.histset == 0

                % the closer one
                if abs(hpos - pc.histval(1)) < abs(hpos - pc.histval(2))
                    ne_gcfg.fcfg.histset = 1;
                else
                    ne_gcfg.fcfg.histset = 2;
                end
            end
            hset = ne_gcfg.fcfg.histset;

            % which histogram line and limit to update?
            if hset == 1

                % compute new value
                hpos = min(hpos, pc.histval(2) - 1/256);
                ne_gcfg.fcfg.histval(1) = hpos;

                % and set
                set(ch.HistLine1, 'YData', 0.002 + 0.996 * [hpos; hpos]);
            else
                hpos = max(hpos, pc.histval(1) + 1/256);
                ne_gcfg.fcfg.histval(2) = hpos;
                set(ch.HistLine2, 'YData', 0.002 + 0.996 * [hpos; hpos]);
            end

            % compute actual value
            hval = hpos * swl(2) + (1 - hpos) * swl(1);

            % and update
            if strcmp(svartyp, 'vmr') && ...
                rtv.ShowV16
                svar.RunTimeVars.ScalingWindow16(hset) = hval;
                swn = svar.RunTimeVars.ScalingWindow16;
            else
                svar.RunTimeVars.ScalingWindow(hset) = hval;
                swn = svar.RunTimeVars.ScalingWindow;
            end

            % as well as image
            hid = (0:255)';
            swld = abs(swl(2) - swl(1));
            swnd = abs(swn(2) - swn(1));
            hid = uint8(min(255, max(0, round( ...
                (hid - (256 / swld) * (swn(1) - swl(1))) * (swld / swnd)))));
            set(ch.HistImage, 'CData', repmat(hid, [1, 16, 3]));

        % otherwise
        else

            % return early
            return;
        end

    % first object -> update movement
    elseif cobj == 1

        % store config in global variable
        if fromroot
            ne_gcfg.fcfg.mpos.ddat = {3, nPos, cc};
        else
            ne_gcfg.cc.(varargin{3}).Config.mpos.ddat = {3, nPos, cc};
        end

    % time-course update
    elseif fromroot && isxff(stvar, 'vtc') && ~isempty(svi)

        % remove statsvar handle content
        svar.SetHandle('RenderSVol', cell(0, 5));

        % update SubVol
        cp = nPos(1) - fPos(cobj, 1);
        nvol = stvar.RunTimeVars.NrOfVolumesPerTC;
        tsvalue = min(nvol, max(1, ((nvol + 2) * cp / 538) + 1));
        stvar.RunTimeVars.SubMapVol = tsvalue;
        tfrom = stvar.RunTimeVars.AvgWindowFrom;
        tstep = stvar.RunTimeVars.AvgWindowStep;
        xtpos = 0.001 * (tstep * (stvar.RunTimeVars.SubMapVol - 1) - tfrom);
        if numel(ch.TCPlotChildren) > 2
            set(ch.TCPlotChildren(end), 'XData', [xtpos, xtpos]);
        end

        % force correct update (but only once)
        allowpreview = false;

        % pretend it's not down (no further updates)
        ne_gcfg.c.btdown = [];
    end

% we still need to update the position
elseif any(ne_gcfg.c.btdown == gcbf)

    % get original position and configuration
    oPos = pc.mpos.ddat{2};
    occ = pc.mpos.ddat{3};

    % depending on modifiers hit
    if isempty(pc.mpos.mods)

        % compute new angles
        cc.anglex = mod(occ.anglex + oPos(1) - nPos(1), 360);
        cc.angley = min(90, max(-90, occ.angley + ...
            0.5 * (oPos(2) - nPos(2))));

    % shifting (translation)
    elseif numel(pc.mpos.mods) == 1 && strcmpi(pc.mpos.mods{1}, 'shift')

        % compute new translation
        cc.trans = min(256, max(-256, occ.trans + [0, nPos - oPos]));

    % zooming
    elseif numel(pc.mpos.mods) == 1 && strcmpi(pc.mpos.mods{1}, 'alt')

        % compute new translation
        cc.zoom = min(5, max(0.2, occ.zoom * (1.01 ^ (round(oPos(2) - nPos(2))))));

    % slicing range x/y
    elseif numel(pc.mpos.mods) == 1 && strcmpi(pc.mpos.mods{1}, 'control')

        % change slicing range
        allowpreview = false;
        warpip = 'linear';

        % compute difference between midpoints
        xdiff = rtv.RenderBBoxFull(1, 2) * min(1, max(-1, (nPos(1) - xPosmid(1)) / xPoshsz(1)));
        ydiff = rtv.RenderBBoxFull(2, 2) * min(1, max(-1, (nPos(2) - xPosmid(2)) / xPoshsz(2)));

        % apply to bounding box
        if xdiff < 0
            svar.RunTimeVars.RenderBBox(1, :) = round( ...
                [min(rtv.RenderBBoxFull(1, 2), -xdiff), rtv.RenderBBoxFull(1, 2)]);
        elseif xdiff > 0
            svar.RunTimeVars.RenderBBox(1, :) = round( ...
                [1, max(1, rtv.RenderBBoxFull(1, 2) - xdiff)]);
        end
        if ydiff < 0
            svar.RunTimeVars.RenderBBox(2, :) = round( ...
                [min(rtv.RenderBBoxFull(2, 2), -ydiff), rtv.RenderBBoxFull(2, 2)]);
        elseif ydiff > 0
            svar.RunTimeVars.RenderBBox(2, :) = round( ...
                [1, max(1, rtv.RenderBBoxFull(2, 2) - ydiff)]);
        end

    % slicing range z
    elseif numel(pc.mpos.mods) == 2 && ...
        any(strcmpi(pc.mpos.mods, 'control')) && any(strcmpi(pc.mpos.mods, 'alt'))

        % change slicing range
        allowpreview = false;
        warpip = 'linear';

        % compute difference between midpoints
        xdiff = rtv.RenderBBoxFull(1, 2) * min(1, max(-1, (nPos(1) - xPosmid(1)) / xPoshsz(1)));
        ydiff = max(0, xPos(4) - nPos(2)) ./ xPos(4);

        % apply to bounding box
        if xdiff < 0
            svar.RunTimeVars.RenderBBox(3, :) = round( ...
                [min(rtv.RenderBBoxFull(3, 2), -xdiff), rtv.RenderBBoxFull(3, 2)]);
        elseif xdiff > 0
            svar.RunTimeVars.RenderBBox(3, :) = round( ...
                [1, max(1, rtv.RenderBBoxFull(3, 2) - xdiff)]);
        end
        if fromroot
            ne_gcfg.fcfg.stalphared = 4 * ydiff;
            pc = ne_gcfg.fcfg;
        else
            ne_gcfg.cc.(varargin{3}).Config.stalphared = 4 * ydiff;
            pc = ne_gcfg.cc.(varargin{3}).Config;
        end
    end

    % re-store config
    if fromroot
        ne_gcfg.fcfg.srfcfg = cc;
    else
        ne_gcfg.cc.(varargin{3}).Config.srfcfg = cc;
    end
end

% update Render UI
ne_gcfg.fcfg.Render.roty = cc.angley;
ne_gcfg.fcfg.Render.rotz = cc.anglex;
ne_gcfg.fcfg.Render.hTag.ED_render_yzenith.String = sprintf('%.0f', cc.angley);
ne_gcfg.fcfg.Render.hTag.ED_render_zazimuth.String = sprintf('%.0f', cc.anglex);

% preview mode
if nargin > 3 && ischar(varargin{4}) && strcmpi(varargin{4}(:)', 'preview')
    preview = true & allowpreview;
else
    preview = false;
end

% get handle of StatsVar and current selection
svarh = svar.Handles;
if isfield(svarh, 'Underlay')
    uvar = svarh.Underlay;
else
    uvar = struct('FileType', '', 'RunTimeVars', struct, 'Handles', struct);
end
if numel(uvar) == 1 && isxff(uvar, {'hdr', 'head', 'vmr'})
    urtv = uvar.RunTimeVars;
    svol = {uvar, 1, struct( ...
        'alpha',  (0:0.25:1)', ...
        'collut', (0:0.25:1)' * ones(1, 3), ...
        'mapvol', 1, ...
        'max',    urtv.ScalingWindow(2), ...
        'min',    urtv.ScalingWindow(1))};
else
    svol = cell(1, 3);
end

% generate svol
if numel(stvar) == 1 && isxff(stvar, {'glm', 'hdr', 'head', 'vmp', 'vtc'}) && ~isempty(svi)
    svol = [svol; cell(2 * numel(svi), 3)];
    stmap = stvar.Map(svi);
    for svc = numel(svi):-1:1
        staval = stmap(svc).TransColorFactor;
        if ci.StatsForceTransparent
            staval = -abs(staval);
        end
        if abs(staval) > ci.StatsMaxAlpha
            staval = sign(staval) * ci.StatsMaxAlpha;
        end
        stclut = stmap(svc).OverlayColors;
        stnlut = size(stclut, 1);
        if staval >= 0
            stalpha = [0; repmat(min(1, staval), 256, 1)];
        else
            staval = -staval;
            stalpha = (0:(staval/4):staval)';
        end
        if stmap(svc).ShowPositiveNegativeFlag > 1
            svol(2*svc+1, :) = {stvar, svi(svc), struct( ...
                'alpha',  stalpha(end:-1:1, 1), ...
                'collut', stclut(end:-1:(round(0.5 * stnlut)+1), :), ...
                'mapvol', svi(svc), ...
                'max',    -stmap(svc).LowerThreshold, ...
                'min',    -stmap(svc).UpperThreshold)};
        end
        if mod(stmap(svc).ShowPositiveNegativeFlag, 2) == 1
            svol(2*svc, :) = {stvar, svi(svc), struct( ...
                'alpha',  stalpha, ...
                'collut', stclut(1:round(0.5 * stnlut), :), ...
                'mapvol', svi(svc), ...
                'max',    stmap(svc).UpperThreshold, ...
                'min',    stmap(svc).LowerThreshold)};
        end
    end
end
svol(cellfun('isempty', svol(:, 1)), :) = [];

% gray-scale lut
if ~isempty(pc.graylut)
    collut = pc.graylut;
else
    collut = [rc.locol; rc.hicol];
end

% and do the work if either page 1 or 2 is shown and SliceVar is valid
if numel(svar) == 1 && isxff(svar, true)

    % and (for temporal vars) the volume number
    svft = lower(svar.Filetype);
    if strcmp(svft, 'vtc')
        volnum = round(ch.Coord.TempSlider.Value);
    elseif strcmp(svft, 'hdr') && svar.NrOfVolumes > 1
        volnum = round(ch.Coord.TempSlider.Value);
    else
        volnum = 1;
    end

    % compute mview
    mview = spmtrf([0, -2 * cpos(3), 2 * cpos(4)]) * ...
        spmtrf([0, 0, 0], (pi / 180) .* [0, -cpos(1), 0]) * ...
        spmtrf([0, 0, 0], (pi / 180) .* [0, 0, cpos(2) - 180]) * ...
        spmtrf([0, 0, 0], [0, 0, 0], cpos(5) .* [2, 2, 2]);

    % sample requested directions and get sampling voxel value and coord
    ropts = struct( ...
        'alpha', (1 / (1 + pc.stalphared * double(~isempty(svi)))) .* rtv.AlphaTable, ...
        'backcolor', rc.bgcol, 'bbox', rtv.RenderBBox, ...
        'collut', collut, 'layer', 1, 'mapvol', volnum, 'preview', preview, ...
        'sranges', rtv.SliceRanges, 'svinterp', rc.imetha, 'svol', {svol}, ...
        'svolafac', 1 / (1 + pc.stalphared * numel(svi)), 'warpip', warpip);
    if isfield(rtv, 'AlphaVolume') && ischar(rtv.AlphaVolume) && ...
        numel(rtv.AlphaVolume) == 24 && ~isempty(regexpi(rtv.AlphaVolume(:)', '^[a-f0-9]+$'))
        try
            avol = xff(rtv.AlphaVolume(:)');
            ropts.avol = avol;
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
    end
    if isfield(rtv, 'ShowV16') && islogical(rtv.ShowV16) && numel(rtv.ShowV16) == 1
        ropts.v16 = rtv.ShowV16;
    end
    if ~preview
        chptr = chfig.Pointer;
        chfig.Pointer = 'watch';
        drawnow;
    end
    svar.RenderToTransimg(mview, tio, ropts);
    if ~preview
        chfig.Pointer = chptr;
    end

% no good variable
else
    setlayer(tio, 1, zeros(512, 512), 1);
end

% update
display(render(tio));
ne_gcfg.c.rpreview = preview;

% remote
if ne_gcfg.c.remote

    % grab screenshot and write to images folder
    simg = tio.Rendered;
    ifmt = lower(ne_gcfg.c.ini.Remote.ImageFormat);
    if strcmp(ifmt, 'jpg')
        iqual = {'Quality', ne_gcfg.c.ini.Remote.ImageJPGQuality};
    else
        iqual = {};
    end
    ipath = [neuroelf_path('remote') '/images'];
    try
        imwrite(simg, sprintf('%s/render.%s', ipath, ifmt), iqual{:});
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
