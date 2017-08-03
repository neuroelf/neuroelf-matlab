function ne_keypress(src, ke, varargin)
% ne_keypress  - handle keyboard (if no control has focus)
%
% FORMAT:       ne_keypress(SRC, EVT, windowid)
%
% Input fields:
%
%       SRC, EVT    Matlab-internal handle and keyboard event information
%       windowid    either empty or a valid satellite window ID
%
% No output fields:
%
% Notes:
%
%   this function should not be called manually but is triggered by
%   keyboard presses while no control has the focus

% Version:  v1.1
% Build:    16052921
% Date:     May-29 2016, 9:52 PM EST
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

% get Key and Modifier from keyboard event (see Matlab docu!)
kk = ke.Key;
mn = ke.Modifier;

% not further allowed while doing modal stuff
if ne_gcfg.c.incb && (~strcmpi(kk, 'r') || ~isempty(mn))
    return;
end
ne_gcfg.c.incb = true;
ne_gcfg.fcfg.mods = mn;

% get configuration and handles
if nargin < 3 || ~ischar(varargin{1}) || ~isfield(ne_gcfg.cc, varargin{1})
    cc = ne_gcfg.fcfg;
    ch = ne_gcfg.h;
else
    tsat = varargin{1};
    ch = ne_gcfg.cc.(tsat);
    cc = ch.Config;
end

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

% get current position
cpos = cc.cpos;
cpag = cc.page;
cstp = cc.cstep;
cori = 2 * double(lower(cc.orient(1)) == 'r') - 1;

% for the main fig
if ne_gcfg.h.MainFigMLH == src

    % handle events without modifiers
    if ~any(km)
        switch (lower(kk))

            % cursor movements -> update cpos/axes props and update window
            case 'downarrow'
                if cpag < 3
                    ne_gcfg.fcfg.cpos(3) = min(127.9999, max(-128, cpos(3) - cstp));
                    ne_setslicepos(0, 0, [], 'OnCursor');
                elseif any([3, 4] == cpag)
                    ne_gcfg.fcfg.srfcfg.angley = min(90, max(-90, ...
                        5 * round(0.2 * ne_gcfg.fcfg.srfcfg.angley) + 5));
                    if cpag == 3
                        ne_setsurfpos;
                    else
                        ne_render_setview;
                    end
                end
            case 'leftarrow'
                if cpag < 3
                    ne_gcfg.fcfg.cpos(1) = min(127.9999, max(-128, cpos(1) + cori * cstp));
                    ne_setslicepos(0, 0, [], 'OnCursor');
                elseif any([3, 4] == cpag)
                    ne_gcfg.fcfg.srfcfg.anglex = mod(...
                        5 * round(0.2 * ne_gcfg.fcfg.srfcfg.anglex) + 5, 360);
                    if cpag == 3
                        ne_setsurfpos;
                    else
                        ne_render_setview;
                    end
                end
            case 'rightarrow'
                if cpag < 3
                    ne_gcfg.fcfg.cpos(1) = min(127.9999, max(-128, cpos(1) - cori * cstp));
                    ne_setslicepos(0, 0, [], 'OnCursor');
                elseif any([3, 4] == cpag)
                    ne_gcfg.fcfg.srfcfg.anglex = mod(...
                        5 * round(0.2 * ne_gcfg.fcfg.srfcfg.anglex) - 5, 360);
                    if cpag == 3
                        ne_setsurfpos;
                    else
                        ne_render_setview;
                    end
                end
            case 'uparrow'
                if cpag < 3
                    ne_gcfg.fcfg.cpos(3) = min(127.9999, max(-128, cpos(3) + cstp));
                    ne_setslicepos(0, 0, [], 'OnCursor');
                elseif any([3, 4] == cpag)
                    ne_gcfg.fcfg.srfcfg.angley = min(90, max(-90, ...
                        5 * round(0.2 * ne_gcfg.fcfg.srfcfg.angley) - 5));
                    if cpag == 3
                        ne_setsurfpos;
                    else
                        ne_render_setview;
                    end
                end

            % toggle lines on/off for VMR Browser
            case 'a'
                sl = ~cc.chair;
                ne_gcfg.fcfg.chair = sl;
                if sl
                    sl = 'on';
                else
                    sl = 'off';
                end

                % depending on currently shown page
                if cpag == 1
                    set([ch.CorLineX, ch.CorLineY, ...
                         ch.SagLineX, ch.SagLineY, ...
                         ch.TraLineX, ch.TraLineY], 'Visible', sl);
                    set([ch.ZoomLineX, ch.ZoomLineY], 'Visible', 'off');
                elseif cpag == 2
                    set([ch.CorLineX, ch.CorLineY, ...
                         ch.SagLineX, ch.SagLineY, ...
                         ch.TraLineX, ch.TraLineY], 'Visible', 'off');
                    set([ch.ZoomLineX, ch.ZoomLineY], 'Visible', sl);
                else
                    set([ch.CorLineX, ch.CorLineY, ...
                         ch.SagLineX, ch.SagLineY, ...
                         ch.TraLineX, ch.TraLineY, ...
                         ch.ZoomLineX, ch.ZoomLineY], 'Visible', 'off');
                end

            % browse mode
            case 'b'
                ne_gcfg.fcfg.paint.mode = sign(ne_gcfg.fcfg.paint.mode) * 1;

            % toggle cluster check
            case 'c'
                if strcmpi(ch.Stats.UsekThr.Enable, 'on')
                    ch.Stats.UsekThr.Value = double(ch.Stats.UsekThr.Value == 0);
                    ne_setstatthrccheck;
                end

            % drawing mode
            case 'd'
                ne_gcfg.fcfg.paint.mode = sign(ne_gcfg.fcfg.paint.mode) * 2;

            % surface display mode
            case {'f', 'w'}

                % only valid for page 3
                if cpag ~= 3 || isempty(ch.Scenery.UserData) || ...
                    numel(ch.Scenery.Value) ~= 1
                    ne_gcfg.c.incb = false;
                    return;
                end

                % get surface
                svars = ch.Scenery.UserData;
                srf = svars{ch.Scenery.Value, 4};
                srfh = handles(srf);
                srfp = srfh.SurfProps;

                % update mode
                if ~strcmpi(srfp{5}, kk)
                    srfp{5} = lower(kk);
                    srf.SetHandle('SurfProps', srfp);

                    % update
                    ne_setsurfpos;
                end

            % gradient mode
            case 'g'
                if cpag == 3 || ~isxff(cc.SliceVar, true)
                    return;
                end
                if cc.SliceVar.RunTimeVars.DisplayType(1) == 'g'
                    cc.SliceVar.RunTimeVars.DisplayType = ...
                        ne_gcfg.c.ini.MainFig.GradientType;
                else
                    cc.SliceVar.RunTimeVars.DisplayType = 'gray';
                end
                if cpag < 3
                    ne_setslicepos;
                end

            % interpolation
            case 'i'
                ch.Interpolate.Value = double(ch.Interpolate.Value == 0);
                if cpag < 3
                    ne_setslicepos;
                end

            % map joining
            case 'j'
                ne_gcfg.fcfg.join = ~ne_gcfg.fcfg.join;
                if cpag < 3
                    ne_setslicepos;
                elseif cpag == 3
                    if isxff(ne_gcfg.fcfg.SurfStatsVar, true)
                        ne_setcsrfstatmap;
                    end
                    ne_setsurfpos(0, 0, true);
                end

            % k-thresholding
            case 'k'
                ch.Stats.UsekThr.Value = double(ch.Stats.UsekThr.Value == 0);
                ne_setstatthrccheck;

            % split into local max./min. for clustertable
            case 'l'
                ne_gcfg.fcfg.localmax = ~ne_gcfg.fcfg.localmax;
                ch.SplitLocalMax.Value = double(ne_gcfg.fcfg.localmax);

            % interpolation method
            case 'm'
                if strcmpi(cc.imethod, 'linear')
                    ne_gcfg.fcfg.imethod = 'cubic';
                else
                    ne_gcfg.fcfg.imethod = 'linear';
                end
                if cpag < 3
                    ne_setslicepos;
                end

            % negate alpha thresholding setting
            case 'n'

                % which page
                if cpag < 3 || cpag == 4

                    % requires valid stats
                    if ~isxff(cc.StatsVar)
                        return;
                    end

                    % for single maps
                    if numel(cc.StatsVarIdx) == 1

                        % negate alpha
                        cc.StatsVar.Map(cc.StatsVarIdx).TransColorFactor = ...
                                min(1, -cc.StatsVar.Map(cc.StatsVarIdx).TransColorFactor);

                    % for multiple maps
                    elseif numel(cc.StatsVarIdx) > 1
                        newalpha = -2 / numel(cc.StatsVarIdx);
                        for mcc = 1:numel(cc.StatsVarIdx)
                            cc.StatsVar.Map(cc.StatsVarIdx(mcc)).TransColorFactor = newalpha;
                        end
                    end
                    if cpag < 3
                        ne_setslicepos;
                    else
                        if ~isempty(cc.StatsVarIdx)
                            ne_setcstatmap;
                        end
                    end
                elseif cpag == 3 && numel(cc.SurfStatsVar) == 1 && ...
                    isxff(cc.SurfStatsVar, {'fsmf', 'glm', 'smp'}) && ...
                    ~isempty(cc.SurfStatsVarIdx)
                    ustvar = cc.SurfStatsVar;
                    ustvix = cc.SurfStatsVarIdx;
                    fullv = (ustvar.Map(ustvix(1)).TransColorFactor < 0 || ...
                        numel(ustvix) < 3);
                    for uc = 1:numel(ustvix)
                        if fullv
                            newv = min(1, -ustvar.Map(ustvix(uc)).TransColorFactor);
                        else
                            newv = -2 / numel(ustvix);
                        end
                        ustvar.Map(ustvix(uc)).TransColorFactor = newv;
                    end
                    ne_setcsrfstatmap;
                end

            % orientation
            case 'o'
                if cc.orient == 'r'
                    ne_setoption(0, 0, 'orientation', 'n');
                else
                    ne_setoption(0, 0, 'orientation', 'r');
                end

            % show "points" (markeredgecolor) for surfaces
            case 'p'

                % only valid for page 3
                if cpag ~= 3 || isempty(ch.Scenery.UserData) || ...
                    numel(ch.Scenery.Value) ~= 1
                    ne_gcfg.c.incb = false;
                    return;
                end

                % get surface
                svars = ch.Scenery.UserData;
                srf = svars{ch.Scenery.Value, 4};
                srfh = handles(srf);
                srfp = srfh.SurfProps;

                % update mode
                if lower(srfp{7}(1)) == 'n'
                    srfp{7} = 'flat';
                else
                    srfp{7} = 'none';
                end
                srf.SetHandle('SurfProps', srfp);

                % update
                ne_setsurfpos(0, 0, true);

            % reset view
            case 'r'

                % try to restore as many settings as possible/necessary
                ne_gcfg.c.btdown = [];
                ne_gcfg.c.btdoup = false;
                ne_gcfg.c.incb = false;
                ne_gcfg.fcfg.chair = 1;
                ne_gcfg.fcfg.cpos = [0, 0, 0];
                ne_gcfg.fcfg.graylut = [];
                ne_gcfg.fcfg.imethod = 'cubic';
                ne_gcfg.fcfg.join = true;
                ne_gcfg.fcfg.mpos.mods = {};
                ne_gcfg.c.ini.MainFig.Orientation = 'N';
                ne_gcfg.fcfg.orient = 'n';
                set(ne_gcfg.h.MainFig.TagStruct.UIM_NeuroElf_OrientNeuro.MLHandle, ...
                    'Checked', 'on');
                set(ne_gcfg.h.MainFig.TagStruct.UIM_NeuroElf_OrientRadio.MLHandle, ...
                    'Checked', 'off');
                ne_gcfg.fcfg.paint.mode = 1;
                ne_gcfg.fcfg.paint.over = [0, 32767];
                ne_gcfg.fcfg.srfcfg.anglex = 180;
                ne_gcfg.fcfg.srfcfg.angley = 0;
                ne_gcfg.fcfg.srfcfg.trans = [0, 0, 0];
                ne_gcfg.fcfg.srfcfg.time = 0;
                ne_gcfg.fcfg.srfcfg.zoom = 1;
                ne_gcfg.fcfg.stalphared = 2;
                ne_gcfg.fcfg.strans = eye(4);
                ne_gcfg.fcfg.strrot = [0, 0, 0];
                ne_gcfg.fcfg.strscl = [1, 1, 1];
                ne_gcfg.fcfg.strtra = [0, 0, 0];
                ne_gcfg.fcfg.strzoom = false;
                ne_gcfg.fcfg.szoom = false;
                ne_gcfg.fcfg.zoom = 0;
                if isxff(cc.StatsVar, true) && numel(cc.StatsVarIdx) == 1
                    cc.StatsVar.Map(cc.StatsVarIdx).TransColorFactor = 1;
                end
                ch.Interpolate.Value = 1;
                ne_setsceneproj(0, 0, 'orthographic');
                ne_setoption(0, 0, 'surfbgcol', [0, 0, 0]);

                % and also reset surface properties to normal
                scu = ne_gcfg.h.Scenery.UserData;
                for sco = 1:size(scu, 1)
                    if isxff(scu{sco, 4}, 'srf')
                        if isempty(scu{sco, 4}.TriangleVertex)
                            scu{sco, 4}.SetHandle('SurfProps', ...
                                {[0, 0, 0], [0, 0, 0], [1, 1, 1], 1, 'w', [], 'none'});
                        else
                            scu{sco, 4}.SetHandle('SurfProps', ...
                                {[0, 0, 0], [0, 0, 0], [1, 1, 1], 1, 'f', [], 'none'});
                        end
                        try
                            btc_meshcolor(scu{sco, 4});
                        catch ne_eo;
                            ne_gcfg.c.lasterr = ne_eo;
                        end
                    end
                end

                % then show page (which also updates the window)
                if cpag < 3

                    % and update the scaling window of the current dataset
                    if isxff(cc.SliceVar, true) && ...
                        isfield(cc.SliceVar.RunTimeVars, 'ScalingWindow') && ...
                        isfield(cc.SliceVar.RunTimeVars, 'ScalingWindowLim')
                        cc.SliceVar.RunTimeVars.ScalingWindow = ...
                            cc.SliceVar.RunTimeVars.ScalingWindowLim;
                        ne_setcvar;
                    end
                    ne_showpage(0, 0, 1);
                elseif cpag == 3
                    ne_setsurfpos(0, 0, 1);
                elseif cpag == 4
                    if isxff(cc.SliceVar, true)
                        cc.SliceVar.RunTimeVars.RenderBBox = ...
                            cc.SliceVar.RunTimeVars.RenderBBoxFull;
                        cc.SliceVar.RunTimeVars.ScalingWindow = ...
                            cc.SliceVar.RunTimeVars.ScalingWindowLim;
                        ne_setcvar;
                    end
                    ne_render_setview;
                end

                % reset size
                if ne_gcfg.c.mlversion >= 804
                    ne_swapfullsize(0, 0, [1002, 696]);
                end
                ne_swapfullsize(0, 0, [1000, 696]);

                % reset pointer
                ne_gcfg.h.MainFig.Pointer = 'arrow';
                drawnow;

            % size adaptation
            case 's'

                % make call
                ne_swapfullsize(0, 0, 'swap');

            % zoomed view
            case 't'

                % if either page 1 or 2 is shown
                if cpag > 0 && cpag < 3

                    % switch through the 4 zoom modes
                    ne_gcfg.fcfg.zoom = mod(cc.zoom + 1, 4);

                    % then either show 3-slice (page 1) or zoom (page 2)
                    if ne_gcfg.fcfg.zoom == 0
                        ne_showpage(0, 0, 1);
                    else
                        ne_showpage(0, 0, 2);
                    end

                % for surface
                elseif cpag == 3 && ~isempty(ch.Scenery.UserData) && ...
                    numel(ch.Scenery.Value) == 1

                    % get surface
                    svars = ch.Scenery.UserData;
                    srf = svars{ch.Scenery.Value, 4};
                    srfh = handles(srf);
                    srfp = srfh.SurfProps;
                    srfi = ne_gcfg.c.ini.Surface;

                    % update alpha value
                    if srfp{4} < srfi.Alpha
                        nalpha = srfi.Alpha;
                    else
                        nalpha = srfi.TransparentAlpha;
                    end
                    srfp{4} = nalpha;
                    srf.SetHandle('SurfProps', srfp);

                    % update
                    ne_setsurfpos;
                end

            % toggle draw/undo drawing mode
            case 'u'
                ch.DrawUndo.Value = 1 - double(ch.DrawUndo.Value > 0);
                drawnow;
                ne_setdrawmode(0, 0, -1);

            % auto extract cluster betas
            case 'x'
                switch lower(ne_gcfg.c.ini.Statistics.ExtractOnSelect)
                    case 'manual'
                        ne_setoption(0, 0, 'extonselect', 'single');
                    case 'multi'
                        ne_setoption(0, 0, 'extonselect', 'manual');
                    case 'single'
                        ne_setoption(0, 0, 'extonselect', 'multi');
                end

            % toggle brain zooming flag
            case 'z'
                if cpag < 3
                    ne_gcfg.fcfg.szoom = ~ne_gcfg.fcfg.szoom;
                    ne_setslicepos;
                end

            % positive tail
            case '1'
                ch.Stats.PosTail.Value = double(ch.Stats.PosTail.Value == 0);
                ne_setstatthrtails;

            % negative tail
            case '2'
                ch.Stats.NegTail.Value = double(ch.Stats.NegTail.Value == 0);
                ne_setstatthrtails;

            % 3D drawing mode
            case '3'
                ne_gcfg.fcfg.paint.mode = sign(ne_gcfg.fcfg.paint.mode) * 3;
        end

    % handle events with ALT
    elseif km(1)
        switch (lower(kk))

            % alpha blending
            case 'downarrow'
                if cpag < 3
                    if isxff(cc.StatsVar, {'ava', 'cmp', 'glm', 'hdr', 'head', 'vmp'})
                        if numel(cc.StatsVarIdx) == 1
                            cc.StatsVar.Map(cc.StatsVarIdx).TransColorFactor = ...
                                max(-2, cc.StatsVar.Map(cc.StatsVarIdx).TransColorFactor - 0.1);
                        else
                            alphadiff = -0.2 / numel(cc.StatsVarIdx);
                            for mcc = 1:numel(cc.StatsVarIdx)
                                cc.StatsVar.Map(cc.StatsVarIdx(mcc)).TransColorFactor = ...
                                    max(-2, cc.StatsVar.Map(cc.StatsVarIdx(mcc)).TransColorFactor + alphadiff);
                            end
                        end
                        ne_setslicepos;
                    end
                elseif cpag == 3 && ...
                    numel(cc.SurfStatsVar) == 1 && ...
                    isxff(cc.SurfStatsVar, {'fsmf', 'glm', 'smp'}) && ...
                    numel(cc.SurfStatsVarIdx) == 1
                    cc.SurfStatsVar.Map(cc.SurfStatsVarIdx).TransColorFactor = ...
                        max(-2, cc.SurfStatsVar.Map(cc.SurfStatsVarIdx).TransColorFactor - 0.1);
                    ne_setcsrfstatmap;
                end
            case 'uparrow'
                if cpag < 3
                    if isxff(cc.StatsVar, {'ava', 'cmp', 'glm', 'hdr', 'head', 'vmp'})
                        if numel(cc.StatsVarIdx) == 1
                            cc.StatsVar.Map(cc.StatsVarIdx).TransColorFactor = ...
                                min(1, cc.StatsVar.Map(cc.StatsVarIdx).TransColorFactor + 0.1);
                        else
                            alphadiff = 0.2 / numel(cc.StatsVarIdx);
                            for mcc = 1:numel(cc.StatsVarIdx)
                                cc.StatsVar.Map(cc.StatsVarIdx(mcc)).TransColorFactor = ...
                                    min(1, cc.StatsVar.Map(cc.StatsVarIdx(mcc)).TransColorFactor + alphadiff);
                            end
                        end
                        ne_setslicepos;
                    end
                elseif cpag == 3 && ...
                    numel(cc.SurfStatsVar) == 1 && ...
                    isxff(cc.SurfStatsVar, {'fsmf', 'glm', 'smp'}) && ...
                    numel(cc.SurfStatsVarIdx) == 1
                    cc.SurfStatsVar.Map(cc.SurfStatsVarIdx).TransColorFactor = ...
                        min(1, cc.SurfStatsVar.Map(cc.SurfStatsVarIdx).TransColorFactor + 0.1);
                    ne_setcsrfstatmap;
                end

            % view selection
            case '1'
                ne_setview(0, 0, 1, 0);
            case '2'
                ne_setview(0, 0, 2, 1);
            case '3'
                ne_setview(0, 0, 2, 2);
            case '4'
                ne_setview(0, 0, 2, 3);
            case '5'
                ne_setview(0, 0, 3, 0);
            case '6'
                if ~isstruct(ne_gcfg.h.Render) || numel(ne_gcfg.h.Render) ~= 1 || ...
                   ~isfield(ne_gcfg.h.Render, 'RendFig') || ...
                   ~isxfigure(ne_gcfg.h.Render.RendFig, true)
                    ne_render;
                else
                    ne_render_setview;
                    ne_showpage(0, 0, 4);
                end

            % max-size adaptation
            case 's'

                % get current GUI size
                gsize = ne_gcfg.h.MainFig.Position(3:4);

                % increase to next 512x512 fuller size
                gsize = cc.fullsize + 512 .* ceil((1 + gsize - cc.fullsize) / 512);

                % make call
                ne_swapfullsize(0, 0, gsize);
        end

    % handle events with CTRL
    elseif km(2)
        switch (lower(kk))

            % switch linking docked slicevars
            case 'l'
                ne_togglelinked;

            % reset TrfPlus field
            case 'r'
                ne_gcfg.fcfg.SliceVar.RunTimeVars.TrfPlus = eye(4);
                if cpag < 3
                    ne_setslicepos;
                end

            % switch dataset
            case 's'
                if cpag < 3
                    svcsel = ch.SliceVar.Value;
                    svavail = ch.SliceVar.String;
                    if ~iscell(svavail)
                        svavail = cellstr(svavail);
                    end
                    if numel(svavail) > 1
                        ch.SliceVar.Value = 1 + mod(svcsel, numel(svavail));
                        ne_setcvar;
                    end
                end

            % close currect SliceVar/SurfVar
            case 'w'
                if cpag ~= 3
                    if numel(cc.SliceVar) == 1 && isxff(cc.SliceVar, true)
                        ne_closefile(0, 0, 'SliceVar');
                    end
                else
                    if numel(cc.SurfVar) == 1 && isxff(cc.SurfVar, 'srf')
                        ne_closefile(0, 0, 'SurfVar');
                    end
                end

            % surface zoom
            case 'downarrow'
                if cpag == 3
                    ne_gcfg.fcfg.srfcfg.zoom = min(20, max(0.05, ...
                        ne_gcfg.fcfg.srfcfg.zoom * 0.99));
                    ne_setsurfpos;
                end
            case 'uparrow'
                if cpag == 3
                    ne_gcfg.fcfg.srfcfg.zoom = min(20, max(0.05, ...
                        ne_gcfg.fcfg.srfcfg.zoom / 0.99));
                    ne_setsurfpos;
                end
        end


    % handle events with SHIFT
    elseif km(3)
        switch (lower(kk))

            % cursor movements
            case 'downarrow'
                if cpag < 3
                    ne_gcfg.fcfg.cpos(2) = min(127.9999, max(-128, cpos(2) - cstp));
                    ne_setslicepos(0, 0, [], 'OnCursor');
                elseif any([3, 4] == cpag)
                    ne_gcfg.fcfg.srfcfg.trans(3) = ...
                        5 * round(0.2 * ne_gcfg.fcfg.srfcfg.trans(3)) - 5;
                    if cpag == 3
                        ne_setsurfpos;
                    else
                        ne_render_setview;
                    end
                end
            case 'leftarrow'
                if cpag < 3
                    tsval = floor(ch.Coord.TempSlider.Value);
                    if cc.tcplot && tsval > 1
                        tsval = tsval - 1;
                        ne_gcfg.h.Coord.TempSlider.Value = tsval;
                        ne_gcfg.h.Coord.Temp.String = sprintf('%d', tsval);
                        ne_setslicepos;
                    end
                elseif any([3, 4] == cpag)
                    ne_gcfg.fcfg.srfcfg.trans(2) = ...
                        5 * round(0.2 * ne_gcfg.fcfg.srfcfg.trans(2)) - 5;
                    if cpag == 3
                        ne_setsurfpos;
                    else
                        ne_render_setview;
                    end
                end
            case 'rightarrow'
                if cpag < 3
                    tsval = floor(ch.Coord.TempSlider.Value);
                    if cc.tcplot && tsval < ch.Coord.TempSlider.Max
                        tsval = tsval + 1;
                        ne_gcfg.h.Coord.TempSlider.Value = tsval;
                        ne_gcfg.h.Coord.Temp.String = sprintf('%d', tsval);
                        ne_setslicepos;
                    end
                elseif any([3, 4] == cpag)
                    ne_gcfg.fcfg.srfcfg.trans(2) = ...
                        5 * round(0.2 * ne_gcfg.fcfg.srfcfg.trans(2)) + 5;
                    if cpag == 3
                        ne_setsurfpos;
                    else
                        ne_render_setview;
                    end
                end
            case 'uparrow'
                if cpag < 3
                    ne_gcfg.fcfg.cpos(2) = min(127.9999, max(-128, cpos(2) + cstp));
                    ne_setslicepos(0, 0, [], 'OnCursor');
                elseif any([3, 4] == cpag)
                    ne_gcfg.fcfg.srfcfg.trans(3) = ...
                        5 * round(0.2 * ne_gcfg.fcfg.srfcfg.trans(3)) + 5;
                    if cpag == 3
                        ne_setsurfpos;
                    else
                        ne_render_setview;
                    end
                end

            % SHIFT-keyboard commands

            % ShowThreshBars toggle
            case 'b'
                ne_setoption(0, 0, 'showthreshbars');

            % copy current map (only for single VMP map)
            case 'c'
                if numel(ne_gcfg.fcfg.StatsVar) == 1 && ...
                    isxff(ne_gcfg.fcfg.StatsVar, 'vmp') && ...
                    numel(ne_gcfg.fcfg.StatsVarIdx) == 1

                    % copy
                    ne_gcfg.fcfg.StatsVar.Map = [ ...
                        lsqueeze(ne_gcfg.fcfg.StatsVar.Map(1:ne_gcfg.fcfg.StatsVarIdx)); ...
                        lsqueeze(ne_gcfg.fcfg.StatsVar.Map(ne_gcfg.fcfg.StatsVarIdx:end))];

                    % display
                    ne_setcstats;
                    ne_setcstatmap(0, 0, ne_gcfg.fcfg.StatsVarIdx + 1);
                end

            % instant seed-corr mode
            case 'i'
                ne_setoption(0, 0, 'instseedcorr');

            % join-MD mode
            case 'j'
                ne_setoption(0, 0, 'joinmd2');
                if cpag < 3
                    ne_setslicepos;
                end
                ne_setsurfpos(0, 0, true);

            % cycle through maps
            case 'm'

                % StatsVarMaps
                if cpag < 3

                    % get current map selection
                    svms = ch.StatsVarMaps.String;
                    svmi = ch.StatsVarMaps.Value;

                    % if only one map
                    if numel(svmi) == 1
                        if ~iscell(svms)
                            svms = cellstr(svms);
                        end

                        % cycle + 1
                        ch.StatsVarMaps.Value = 1 + mod(svmi, numel(svms));
                        ne_setcstatmap;
                    end

                % SurfStatsVarMaps
                elseif cpag == 3
                    svms = ch.SurfStatsVarMaps.String;
                    svmi = ch.SurfStatsVarMaps.Value;
                    if numel(svmi) == 1
                        if ~iscell(svms)
                            svms = cellstr(svms);
                        end
                        ch.SurfStatsVarMaps.Value = 1 + mod(svmi, numel(svms));
                        ne_setcsrfstatmap;
                    end
                end

            % screenshot
            case 's'

                % create screenshot
                if cpag == 3 || ne_gcfg.c.ostype.ispc
                    ne_screenshot(0, 0, '', '', 'high-q');
                else
                    ne_screenshot(0, 0, '');
                end

            % ShowThreshText toggle
            case 't'
                ne_setoption(0, 0, 'showthreshtext');

            % dataset selection
            case {'1', '2', '3', '4', '5', '6', '7', '8', '9'}
                slvarsel = str2double(kk);
                if size(ch.SliceVar.UserData, 1) >= slvarsel
                    ch.SliceVar.Value = slvarsel;
                    ne_setcvar;
                    if cpag == 3
                        ne_showpage(0, 0, 3);
                    end
                else
                    stvarsel = slvarsel - size(ch.SliceVar.UserData, 1);
                    if size(ch.StatsVar.UserData, 1) >= stvalsel
                        ch.StatsVar.Value = stvarsel;
                        ne_setcstats;
                        if cpag == 3
                            ne_showpage(0, 0, 3);
                        end
                    end
                end
        end
    end

% for satellite windows (other than TCPlot)
elseif ~strcmpi(cc.sattype, 'tcplot')

    % make sure cstep has 3 values!
    if numel(cstp) == 1
        cstp = cstp([1, 1, 1]);
    end

    % handle events without modifiers
    if ~any(km)
        switch (lower(kk))

            % cursor movements -> update cpos/axes props and update window
            case 'downarrow'
                if cpag < 3
                    ne_gcfg.cc.(tsat).Config.cpos(3) = ...
                        min(127.9999, max(-128, cpos(3) - cstp(3)));
                    ne_setsatslicepos(0, 0, tsat);
                elseif any([3, 4] == cpag)
                    ne_gcfg.cc.(tsat).Config.srfcfg.angley = ...
                        min(90, max(-90, 5 * round(0.2 * cc.srfcfg.angley) + 5));
                    if cpag == 3
                        ne_setsurfpos(0, 0, tsat);
                    else
                        ne_render_setview(0, 0, tsat);
                    end
                end
            case 'leftarrow'
                if cpag < 3
                    ne_gcfg.cc.(tsat).Config.cpos(1) = ...
                        min(127.9999, max(-128, cpos(1) + cori * cstp(1)));
                    ne_setsatslicepos(0, 0, tsat);
                elseif any([3, 4] == cpag)
                    ne_gcfg.cc.(tsat).Config.srfcfg.anglex = ...
                        mod(5 * round(0.2 * cc.srfcfg.anglex) + 5, 360);
                    if cpag == 3
                        ne_setsurfpos(0, 0, tsat);
                    else
                        ne_render_setview(0, 0, tsat);
                    end
                end
            case 'rightarrow'
                if cpag < 3
                    ne_gcfg.cc.(tsat).Config.cpos(1) = ...
                        min(127.9999, max(-128, cpos(1) - cori * cstp(1)));
                    ne_setsatslicepos(0, 0, tsat);
                elseif any([3, 4] == cpag)
                    ne_gcfg.cc.(tsat).Config.srfcfg.anglex = ...
                        mod(5 * round(0.2 * cc.srfcfg.anglex) - 5, 360);
                    if cpag == 3
                        ne_setsurfpos(0, 0, tsat);
                    else
                        ne_render_setview(0, 0, tsat);
                    end
                end
            case 'uparrow'
                if cpag < 3
                    ne_gcfg.cc.(tsat).Config.cpos(3) = ...
                        min(127.9999, max(-128, cpos(3) + cstp(3)));
                    ne_setsatslicepos(0, 0, tsat);
                elseif any([3, 4] == cpag)
                    ne_gcfg.cc.(tsat).Config.srfcfg.angley = ...
                        min(90, max(-90, 5 * round(0.2 * cc.srfcfg.angley) - 5));
                    if cpag == 3
                        ne_setsurfpos(0, 0, tsat);
                    else
                        ne_render_setview(0, 0, tsat);
                    end
                end

            % toggle lines on/off for VMR Browser
            case 'a'
                sl = ~cc.chair;
                ne_gcfg.cc.(tsat).Config.chair = sl;
                if sl
                    sl = 'on';
                else
                    sl = 'off';
                end

                % depending on currently shown page
                if cpag == 1
                    set([ch.CorLineX, ch.CorLineY, ...
                         ch.SagLineX, ch.SagLineY, ...
                         ch.TraLineX, ch.TraLineY], 'Visible', sl);
                    set([ch.ZoomLineX, ch.ZoomLineY], 'Visible', 'off');
                elseif cpag == 2
                    set([ch.CorLineX, ch.CorLineY, ...
                         ch.SagLineX, ch.SagLineY, ...
                         ch.TraLineX, ch.TraLineY], 'Visible', 'off');
                    set([ch.ZoomLineX, ch.ZoomLineY], 'Visible', sl);
                else
                    set([ch.CorLineX, ch.CorLineY, ...
                         ch.SagLineX, ch.SagLineY, ...
                         ch.TraLineX, ch.TraLineY, ...
                         ch.ZoomLineX, ch.ZoomLineY], 'Visible', 'off');
                end

            % background color (surfaces)
            case 'b'
                if cpag == 3
                    oldcolor = round(255 .* get(ch.Surface, 'Color'));
                    newcolor = colorpicker(oldcolor, {'Axes background'});
                    if ~isequal(newcolor, oldcolor)
                        ne_satsetcolor(0, 0, tsat, newcolor);
                    end
                end

            % gradient mode
            case 'g'
                if cpag == 3 || ~isxff(cc.SliceVar, true)
                    return;
                end
                if cc.SliceVar.RunTimeVars.DisplayType(1) == 'g'
                    cc.SliceVar.RunTimeVars.DisplayType = ...
                        ne_gcfg.c.ini.MainFig.GradientType;
                else
                    cc.SliceVar.RunTimeVars.DisplayType = 'gray';
                end
                if cpag < 3
                    ne_setsatslicepos(0, 0, tsat);
                end

            % interpolation
            case 'i'
                ne_gcfg.cc.(tsat).Interpolate.Value = double(ch.Interpolate.Value == 0);
                ne_setsatslicepos(0, 0, tsat);

            % interpolation method
            case 'm'
                if strcmpi(cc.imethod, 'linear')
                    ne_gcfg.cc.(tsat).Config.imethod = 'cubic';
                else
                    ne_gcfg.cc.(tsat).Config.imethod = 'linear';
                end
                ne_setsatslicepos(0, 0, tsat);

            % orientation
            case 'o'
                if cc.orient == 'r'
                    ne_setoption(0, 0, 'orientation', 'n', tsat);
                else
                    ne_setoption(0, 0, 'orientation', 'r', tsat);
                end

            % show "points" (markeredgecolor) for surfaces
            case 'p'
                
            % zoomed view
            case 't'

                % if either page 1 or 2 is shown
                if cpag == 1 || cpag == 2

                    % switch through the 4 zoom modes
                    ne_gcfg.cc.(tsat).Config.zoom = mod(cc.zoom + 1, 4);

                    % then either show 3-slice (page 1) or zoom (page 2)
                    if ne_gcfg.cc.(tsat).Config.zoom == 0
                        ne_showsatpage(0, 0, tsat, 1);
                    else
                        ne_showsatpage(0, 0, tsat, 2);
                    end
                end

            % drawing also changeable for satellites (if linked)
            case 'u'
                if (cpag == 1 || cpag == 2) && ne_gcfg.c.linked
                    ne_gcfg.h.DrawUndo.Value = 1 - double(ne_gcfg.h.DrawUndo.Value > 0);
                    drawnow;
                    ne_setdrawmode(0, 0, -1);
                end
        end

    % handle events with ALT
    %elseif km(1)

    % handle events with CTRL
    elseif km(2)
        switch (lower(kk))

            % alpha blending
            case 'downarrow'
                if any([3, 4] == cpag)
                    ne_gcfg.cc.(tsat).Config.srfcfg.zoom = ...
                        min(5, max(0.2, cc.srfcfg.zoom * 0.99));
                    if cpag == 3
                        ne_setsurfpos(0, 0, tsat);
                    else
                        ne_render_setview(0, 0, tsat);
                    end
                end
            case 'uparrow'
                if any([3, 4] == cpag)
                    ne_gcfg.cc.(tsat).Config.srfcfg.zoom = ...
                        min(5, max(0.2, cc.srfcfg.zoom / 0.99));
                    if cpag == 3
                        ne_setsurfpos(0, 0, tsat);
                    else
                        ne_render_setview(0, 0, tsat);
                    end
                end

            % KEYBOARD shortcuts

            % switch linking docked slicevars
            case 'l'
                ne_togglelinked;

        end


    % handle events with SHIFT
    elseif km(3)
        switch (lower(kk))

            % cursor movements
            case 'downarrow'
                if cpag < 3
                    ne_gcfg.cc.(tsat).Config.cpos(2) = ...
                        min(127.9999, max(-128, cpos(2) - cstp(2)));
                    ne_setsatslicepos(0, 0, tsat);
                elseif any([3, 4] == cpag)
                    ne_gcfg.cc.(tsat).Config.srfcfg.trans(3) = ...
                        5 * round(0.2 * cc.srfcfg.trans(3)) - 5;
                    if cpag == 3
                        ne_setsurfpos(0, 0, tsat);
                    else
                        ne_render_setview(0, 0, tsat);
                    end
                end
            case 'leftarrow'
                if cpag < 3
                    tsval = floor(ch.Coord.TempSlider.Value);
                    if cc.tcplot && tsval > 1
                        tsval = tsval - 1;
                        ne_gcfg.cc.(tsat).Coord.TempSlider.Value = tsval;
                        ne_setsatslicepos(0, 0, tsat);
                    end
                elseif any([3, 4] == cpag)
                    ne_gcfg.cc.(tsat).Config.srfcfg.trans(2) = ...
                        5 * round(0.2 * cc.srfcfg.trans(2)) - 5;
                    if cpag == 3
                        ne_setsurfpos(0, 0, tsat);
                    else
                        ne_render_setview(0, 0, tsat);
                    end
                end
            case 'rightarrow'
                if cpag < 3
                    tsval = floor(ch.Coord.TempSlider.Value);
                    if cc.tcplot && tsval < ch.Coord.TempSlider.Max
                        tsval = tsval + 1;
                        ne_gcfg.cc.(tsat).Coord.TempSlider.Value = tsval;
                        ne_setsatslicepos(0, 0, tsat);
                    end
                elseif any([3, 4] == cpag)
                    ne_gcfg.cc.(tsat).Config.srfcfg.trans(2) = ...
                        5 * round(0.2 * cc.srfcfg.trans(2)) + 5;
                    if cpag == 3
                        ne_setsurfpos(0, 0, tsat);
                    else
                        ne_render_setview(0, 0, tsat);
                    end
                end
            case 'uparrow'
                if cpag < 3
                    ne_gcfg.cc.(tsat).Config.cpos(2) = ...
                        min(127.9999, max(-128, cpos(2) + cstp(2)));
                    ne_setsatslicepos(0, 0, tsat);
                elseif any([3, 4] == cpag)
                    ne_gcfg.cc.(tsat).Config.srfcfg.trans(3) = ...
                        5 * round(0.2 * cc.srfcfg.trans(3)) + 5;
                    if cpag == 3
                        ne_setsurfpos(0, 0, tsat);
                    else
                        ne_render_setview(0, 0, tsat);
                    end
                end

            % ShowThreshBars toggle
            case 'b'
                ne_setoption(0, 0, 'showthreshbars');

            % gradient mode for underlay
            case 'g'
                ne_gcfg.cc.(tsat).Config.gradientu = ~cc.gradientu;
                if cpag < 3
                    ne_setsatslicepos(0, 0, tsat);
                end

            % cycle through maps
            case {'m', 'n'}

                % StatsVarMaps
                if cpag < 3

                    % get current index
                    svmi = cc.StatsVarIdx;

                    % if only one map
                    if numel(svmi) == 1 && numel(cc.StatsVar) == 1 && isxff(cc.StatsVar)

                        % for VMPs
                        if isxff(cc.StatsVar, 'vmp')
                            svms = numel(cc.StatsVar.Map);

                        % for HDRs
                        elseif isxff(cc.StatsVar, 'hdr')
                            svms = size(cc.StatsVar.VoxelData, 4);

                        % for HEADs
                        elseif isxff(cc.StatsVar, 'head')
                            svms = numel(cc.StatsVar.Brick);

                        % otherwise nothing
                        else
                            svms = 1;
                        end

                        % cycle +/- 1
                        if lower(kk) == 'm'
                            ne_gcfg.cc.(tsat).Config.StatsVarIdx = 1 + mod(svmi, svms);
                        else
                            ne_gcfg.cc.(tsat).Config.StatsVarIdx = 1 + mod(svmi + svms - 2, svms);
                        end

                        % update title
                        mapnames = cc.StatsVar.MapNames(ne_gcfg.c.extmapnames);
                        stmap = sprintf(' (Map %d: %s)', ne_gcfg.cc.(tsat).Config.StatsVarIdx, ...
                            mapnames{ne_gcfg.cc.(tsat).Config.StatsVarIdx});
                        [stp, stf] = fileparts(cc.StatsVar.FilenameOnDisk);
                        if isempty(stf)
                            stf = 'Unsaved';
                        end
                        stname = sprintf(' - %s%s', stf, stmap);
                        if isxff(cc.SliceVar)
                            [slp, slf] = fileparts(cc.SliceVar.FilenameOnDisk);
                        else
                            slf = 'Empty';
                        end
                        ne_gcfg.cc.(tsat).Satellite.Name = sprintf('NeuroElf GUI - %s%s', slf, stname);

                        % and update
                        set(ch.CorStatsText, 'String', ' ');
                        set(ch.ZoomStatsText, 'String', ' ');
                        ne_setsatslicepos(0, 0, tsat);
                    end

                % SurfStatsVarMaps
                else

                    % get list of shown surfaces (with a Matlab UI handle)
                    srflist = find(~cellfun('isempty', ne_gcfg.cc.(tsat).Scenery.UserData(:, 5)));

                    % only continue if one single surface
                    if numel(srflist) == 1

                        % get this surface and stats handles
                        srfsh = ne_gcfg.cc.(tsat).Scenery.UserData(srflist, :);
                        srf = srfsh{4};
                        srfh = srfsh{5};
                        srfsh = srfsh{6};

                        % only continue with one map
                        if iscell(srfsh) && numel(srfsh) == 2 && ...
                            isxff(srfsh{1}, true) && numel(srfsh{2}) == 1
                            ne_gcfg.cc.(tsat).Config.SurfStatsVar = srfsh{1};
                            svmi = srfsh{2};
                            switch (lower(srfsh{1}.Filetype))
                                case {'fsmf', 'glm', 'smp'}
                                    svms = numel(srfsh{1}.Map);
                                otherwise
                                    svms = 1;
                            end
                            if lower(kk) == 'n'
                                svmi = 1 + mod(svmi, svms);
                            else
                                svmi = 1 + mod(svmi + svms - 2, svms);
                            end
                            ne_gcfg.cc.(tsat).Config.SurfStatsVarIdx = svmi;
                            ne_gcfg.cc.(tsat).Scenery.UserData{srflist, 6}{2} = svmi;
                            btc_meshcolor(srf, false, tsat, srfh);
                        end
                    end
                end

            % renderer
            case 'r'

                % switch renderer
                if strcmpi(ne_gcfg.cc.(tsat).Satellite.Renderer(1), 'z')
                    ne_gcfg.cc.(tsat).Config.renderer = 'OpenGL';
                    ne_gcfg.cc.(tsat).Satellite.Renderer = 'OpenGL';
                else
                    ne_gcfg.cc.(tsat).Config.renderer = 'zbuffer';
                    ne_gcfg.cc.(tsat).Satellite.Renderer = 'zbuffer';
                end
                drawnow;

            % screenshot
            case 's'

                % create screenshot
                if cpag == 3 || ne_gcfg.c.ostype.ispc
                    ne_screenshot(0, 0, tsat, '', 'high-q');
                else
                    ne_screenshot(0, 0, tsat);
                end

            % ShowThreshText toggle
            case 't'
                ne_setoption(0, 0, 'showthreshtext');
        end
    end

% for TCPlot satellites
else

    % no modifier keys
    if ~any(km)
        switch (lower(kk))
            case 's'
                scol = ch.Satellite.Color;
                ch.Satellite.Color = [1, 1, 1];
                ne_screenshot(0, 0, tsat, '');
                ch.Satellite.Color = scol;
        end

    % with shift
    elseif km(3)
        switch (lower(kk))
            case 's'
                scol = ch.Satellite.Color;
                ch.Satellite.Color = [1, 1, 1];
                ne_screenshot(0, 0, tsat, '', 'high-q');
                ch.Satellite.Color = scol;
        end
    end
end

% open up keypress again
ne_gcfg.c.incb = false;
