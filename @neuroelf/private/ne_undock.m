function varargout = ne_undock(varargin)
% ne_undock  - undock the current document in satellite window
%
% FORMAT:       [hSat, tags, iSat] = ne_undock(SRC, EVT)
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%
% Output fields:
%
%       hSat        xfigure handle of satellite window
%       tags        1x1 structure with further handles and settings
%       iSat        string representing the satellite ID
%
% Example:
%
%       [surface_viewer, ~, surface_id] = ne_undock(0, 0);
%       neuroelf_gui('satresize', surface_id, [640, 480]);
%       neuroelf_gui('setsurfpos', surface_id, {150, 15, [0, -30, 0], 1.25, 0});

% Version:  v1.1
% Build:    16052822
% Date:     May-28 2016, 10:03 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;
cini = ne_gcfg.c.ini;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% only for pages 1 through 4
if ~any(cc.page == (1:4))
    return;
end

% echo
if ne_gcfg.c.echo
    ne_echo('neuroelf_gui(''undock'');');
end

% main figure size
mfsz = ch.MainFig.Position(3:4);
bsc = fieldnames(ne_gcfg.cc);

% open satellite window
hSat = xfigure([neuroelf_path('tfg') '/ne_satellite.tfg']);
tSat = hSat.Tag;
iSat = tSat(1:8);

% create tags
tSat = hSat.TagStruct;
htag = struct;
htag.Satellite = hSat;
htag.SatelliteMLH = hSat.MLHandle;

% add crosshair lines to images axes objects -> SAG, COR, TRA, Zoom
chax = tSat.(sprintf('AX_%s_Slice_SAG', iSat)).MLHandle;
set(chax, 'Units', 'pixels');
htag.SagAxes = chax;
htag.SagLineX = line([0; 0.999], [0.5; 0.5], 'Color', cc.chcol, 'Parent', chax);
htag.SagLineY = line([0.5; 0.5], [0.001; 0.999], 'Color', cc.chcol, 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');
chax = tSat.(sprintf('AX_%s_Slice_COR', iSat)).MLHandle;
set(chax, 'Units', 'pixels');
htag.CorAxes = chax;
htag.CorLineX = line([0; 0.999], [0.5; 0.5], 'Color', cc.chcol, 'Parent', chax);
htag.CorLineY = line([0.5; 0.5], [0.001; 0.999], 'Color', cc.chcol, 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');
htag.CorStatsText = [text(0.99, 0.22, get(ch.CorStatsText(1), 'String'), ...
    'Parent', chax, 'Color', [1, 1, 1], 'FontSize', 10, 'HorizontalAlignment', 'right'), ...
    text(0.99, 0.77, get(ch.CorStatsText(2), 'String'), ...
    'Parent', chax, 'Color', [1, 1, 1], 'FontSize', 10, 'HorizontalAlignment', 'right')];
chax = tSat.(sprintf('AX_%s_Slice_TRA', iSat)).MLHandle;
set(chax, 'Units', 'pixels');
htag.TraAxes = chax;
htag.TraLineX = line([0; 0.999], [0.5; 0.5], 'Color', cc.chcol, 'Parent', chax);
htag.TraLineY = line([0.5; 0.5], [0.001; 0.999], 'Color', cc.chcol, 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');
chax = tSat.(sprintf('AX_%s_Slice_Zoom', iSat)).MLHandle;
set(chax, 'Units', 'pixels');
htag.ZoomAxes = chax;
htag.ZoomLineX = line([0; 0.999], [0.5; 0.5], 'Color', cc.chcol, 'Parent', chax);
htag.ZoomLineY = line([0.5; 0.5], [0.001; 0.999], 'Color', cc.chcol, 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');
htag.ZoomStatsText = [text(0.99, 0.22, get(ch.ZoomStatsText(1), 'String'), 'Parent', chax, 'Color', [1, 1, 1], ...
    'FontSize', 12, 'HorizontalAlignment', 'right'), ...
    text(0.99, 0.78, get(ch.ZoomStatsText(2), 'String'), 'Parent', chax, 'Color', [1, 1, 1], ...
    'FontSize', 12, 'HorizontalAlignment', 'right')];

% time-dim and interpolation controls replacement
htag.Coord = struct('TempSlider', struct( ...
    'Max',   ch.Coord.TempSlider.Max, ...
    'Value', ch.Coord.TempSlider.Value));
htag.Interpolate = struct('Value', ch.Interpolate.Value);
htag.Surface = tSat.(sprintf('AX_%s_Slice_Zoom', iSat)).MLHandle;
htag.SurfaceStatsBar = [];

% callbacks
hSat.CloseRequestFcn = {@ne_closesatwindow, iSat};
hSat.KeyPressFcn = {@ne_keypress, iSat};
hSat.KeyReleaseFcn = @ne_keyrelease;
hSat.WindowButtonDownFcn = @ne_btdown;
hSat.WindowButtonMotionFcn = @ne_btmove;
hSat.WindowButtonUpFcn = {@ne_btup, iSat};
hSat.Renderer = cc.renderer;

% set configuration
htag.Config = cc;
htag.Config.plp = [];
htag.Config.plpcfg = struct( ...
    'bcolor',  [], ...
    'color',   cini.PLPPlot.Color, ...
    'cond',    '', ...
    'range',   cini.PLPPlot.Distance, ...
    'sel',     [], ...
    'symsize', cini.PLPPlot.SymbolSize);
htag.Config.plpph = [];
htag.Config.sattag = iSat;

% for slice-based view
if cc.page == 1 || ...
    cc.page == 2

    % get selected SliceVar
    slvar = cc.SliceVar;
    
    % update cstep?
    if numel(slvar) == 1 && ...
        isxff(slvar, 'vmr')
        htag.Config.cstep = slvar.BoundingBox.ResXYZ([3, 1, 2]);
    end

    % record type, create transimg objects, and keep required settings
    htag.Config.sattype = 'slice';
    htag.Config.tio = struct( ...
        'imSag', transimg(256, 256), ...
        'imCor', transimg(256, 256), ...
        'imTra', transimg(256, 256), ...
        'imSlZ', transimg(512, 512));
    ne_gcfg.tio.satSag = htag.Config.tio.imSag;
    ne_gcfg.tio.satCor = htag.Config.tio.imCor;
    ne_gcfg.tio.satTra = htag.Config.tio.imTra;
    ne_gcfg.tio.satSlz = htag.Config.tio.imSlZ;
    tio = htag.Config.tio;
    sethandle(tio.imSag, ...
        get(tSat.(sprintf('IM_%s_Slice_SAG', iSat)).MLHandle, 'Children'));
    sethandle(tio.imCor, ...
        get(tSat.(sprintf('IM_%s_Slice_COR', iSat)).MLHandle, 'Children'));
    sethandle(tio.imTra, ...
        get(tSat.(sprintf('IM_%s_Slice_TRA', iSat)).MLHandle, 'Children'));
    sethandle(tio.imSlZ, ...
        get(tSat.(sprintf('IM_%s_Slice_Zoom', iSat)).MLHandle, 'Children'));
    htag.Config.slicepos = [ ...
        tSat.(sprintf('IM_%s_Slice_SAG', iSat)).Position; ...
        tSat.(sprintf('IM_%s_Slice_COR', iSat)).Position; ...
        tSat.(sprintf('IM_%s_Slice_TRA', iSat)).Position];
    htag.Config.slicepos(:, 3:4) = ...
        htag.Config.slicepos(:, 1:2) + htag.Config.slicepos(:, 3:4);
    htag.Config.zslicepos = ...
        tSat.(sprintf('IM_%s_Slice_Zoom', iSat)).Position;
    htag.Config.zslicepos(3:4) = ...
        htag.Config.zslicepos(1:2) + htag.Config.zslicepos(3:4);
    htag.Config.zslicepos = htag.Config.zslicepos([1, 1, 1], :);

    % text tags
    htag.Text.BVSX = tSat.(sprintf('ED_%s_BVSX', iSat));
    htag.Text.BVSY = tSat.(sprintf('ED_%s_BVSY', iSat));
    htag.Text.BVSZ = tSat.(sprintf('ED_%s_BVSZ', iSat));
    htag.Text.TALX = tSat.(sprintf('ED_%s_TALX', iSat));
    htag.Text.TALY = tSat.(sprintf('ED_%s_TALY', iSat));
    htag.Text.TALZ = tSat.(sprintf('ED_%s_TALZ', iSat));
    htag.Text.Values = tSat.(sprintf('TX_%s_SValues', iSat));
    htag.Coord.TEdX = mlhandle(htag.Text.TALX);
    htag.Coord.TEdY = mlhandle(htag.Text.TALY);
    htag.Coord.TEdZ = mlhandle(htag.Text.TALZ);
    htag.Coord.VEdX = mlhandle(htag.Text.BVSX);
    htag.Coord.VEdY = mlhandle(htag.Text.BVSY);
    htag.Coord.VEdZ = mlhandle(htag.Text.BVSZ);
    htag.Config.txtpos = [ ...
        htag.Coord.TEdX, ...
        htag.Coord.TEdY, ...
        htag.Coord.TEdZ, ...
        htag.Coord.VEdX, ...
        htag.Coord.VEdY, ...
        htag.Coord.VEdZ];
    htag.Text.TALX.Callback = {@ne_setwindowpos, iSat};
    htag.Text.TALY.Callback = {@ne_setwindowpos, iSat};
    htag.Text.TALZ.Callback = {@ne_setwindowpos, iSat};
    htag.Text.BVSX.Callback = {@ne_setwindowpos, iSat};
    htag.Text.BVSY.Callback = {@ne_setwindowpos, iSat};
    htag.Text.BVSZ.Callback = {@ne_setwindowpos, iSat};

    % copy over some more config
    htag.Stats = struct( ...
        'NegTail', struct('Value', ch.Stats.NegTail.Value), ...
        'PosTail', struct('Value', ch.Stats.PosTail.Value));

    % show correct page
    ne_gcfg.cc.(iSat) = htag;
    ne_showsatpage(0, 0, iSat, cc.page);

    % if slice and stats are valid xff
    if isxff(htag.Config.StatsVar, true)
        if numel(htag.Config.StatsVarIdx) == 1
            mapnames = htag.Config.StatsVar.MapNames(ne_gcfg.c.extmapnames);
            stmap = sprintf(' (Map %d: %s)', htag.Config.StatsVarIdx, ...
                mapnames{htag.Config.StatsVarIdx});
        else
            stmap = sprintf(' (%d maps%s)', numel(htag.Config.StatsVarIdx), ...
                sprintf(', %d', htag.Config.StatsVarIdx));
        end
        [stp, stf] = fileparts(htag.Config.StatsVar.FilenameOnDisk(true));
        if isempty(stf)
            stf = 'Unsaved';
        end
        stname = sprintf(' - %s%s', stf, stmap);
    else
        stname = '';
    end
    if isxff(slvar, true)
        [slp, slf] = fileparts(slvar.FilenameOnDisk(true));
        if isempty(slf)
            slf = 'Unsaved';
        end
    else
        slf = 'Empty';
    end
    hSat.Name = sprintf('NeuroElf GUI - %s%s', slf, stname);
    
    % force renderer to 'zbuffer'
    hSat.Renderer = 'zbuffer';

% for surface views
elseif cc.page == 3

    % make settings
    htag.Config.sattype = 'surf';
    htag.Scenery = get(ne_gcfg.h.Scenery.MLHandle);
    srf = htag.Surface;
    srfcfg = cini.Surface;
    srfbcl = (1 / 255) .* srfcfg.BackgroundColor(:)';
    set(srf, 'Color', srfbcl);
    set(srf, 'View', [90, 0]);
    slim = [-128, 128];
    set(srf, 'XLim', 4 * slim, 'YLim', slim, 'ZLim', slim);
    set(srf, 'XTick', [], 'YTick', [], 'ZTick', []);
    set(srf, 'Units', 'normalized');
    for lc = 1:numel(srfcfg.Lights)
        light('Parent', srf, 'Position', ...
            srfcfg.Lights{lc}, 'Color', (1 / 255) .* srfcfg.LightColors{lc});
    end
    set(srf, 'XColor', [0, 0, 0], 'YColor', [0, 0, 0], 'ZColor', [0, 0, 0]);
    sbarp = get(ne_gcfg.h.SurfaceStatsBar);
    ssbv = sbarp.Vertices;
    ssbf = sbarp.Faces;
    htag.SurfaceStatsBar = ...
        patch(ssbv(:, 1), ssbv(:, 2), ssbv(:, 3), zeros(size(ssbv, 1), 1), ...
        'FaceColor', 'none', 'EdgeColor', 'none', 'Parent', srf, 'Visible', 'off');
    set(htag.SurfaceStatsBar, 'Faces', ssbf, 'FaceVertexCData', sbarp.FaceVertexCData, ...
        'FaceColor', 'flat');
    set(htag.SurfaceStatsBar, 'Visible', sbarp.Visible);
    stextxpos = 0.015 * slim(1) + 0.985 * slim(2);
    stextypos1 = 0.78 * slim(1) + 0.22 * slim(2);
    stextypos2 = 0.22 * slim(1) + 0.78 * slim(2);
    htag.SurfaceStatsText = [text(0, stextxpos, stextypos1, ' ', 'Parent', srf, ...
        'Color', [1, 1, 1], 'FontSize', 12, 'HorizontalAlignment', 'right'), ...
        text(0, stextxpos, stextypos2, ' ', 'Parent', srf, 'Color', [1, 1, 1], ...
        'FontSize', 12, 'HorizontalAlignment', 'right')];
    htag.SurfaceTransform = hgtransform('Parent', srf);
    htag.Config.surfpos = hSat.Position;
    htag.Config.surfpos(1:2) = 0;
    hSat.Units = 'normalized';
    set(srf, 'Position', [0, 0, 1, 1]);
    hSat.Units = 'pixels';

    % add all visible surfaces
    scu = htag.Scenery.UserData;
    sci = htag.Scenery.Value;
    if isempty(scu)
        sci = [];
    else
        scu{1, 6} = [];
    end
    for scc = sci(:)'

        % get xff (surface), current surface handle (mesh) props
        f = scu{scc, 4};
        fh = handles(f);
        fmhp = get(fh.Surface);
        if isempty(fmhp.FaceVertexAlphaData)
            fmhp.FaceVertexAlphaData = 1;
        end

        % get coordinates and normals (with current main window config)
        [p, pn] = btc_meshcn(f, cc.srfcfgload, ...
            ~strcmpi(fmhp.FaceColor, 'none') || ...
            ~strcmpi(cc.renderer, 'opengl'));

        % create sub-group
        hgtrf = hgtransform('Parent', htag.SurfaceTransform);

        % create surface
        hold(srf, 'on');
        if ~isempty(f.TriangleVertex)
            hsrf = patch( ...
                'Faces', f.TriangleVertex(:, [1, 3, 2]), ...
                'Vertices', p, ...
                'VertexNormals', pn, ...
                'FaceVertexAlphaData', fmhp.FaceVertexAlphaData, ...
                'FaceVertexCData', fmhp.FaceVertexCData, ...
                'FaceColor', fmhp.FaceColor, ...
                'EdgeColor', fmhp.EdgeColor, ...
                'Parent', hgtrf);

            % and make some more initial settings
            set(hsrf, 'AmbientStrength', fmhp.AmbientStrength, ...
                      'BackFaceLighting', fmhp.BackFaceLighting, ...
                      'ButtonDownFcn', @ne_btdown, ...
                      'DiffuseStrength', fmhp.DiffuseStrength, ...
                      'EdgeLighting', fmhp.EdgeLighting, ...
                      'FaceAlpha', fmhp.FaceAlpha, ...
                      'AlphaDataMapping', fmhp.AlphaDataMapping, ...
                      'FaceLighting', fmhp.FaceLighting, ...
                      'LineStyle', fmhp.LineStyle, ...
                      'Marker', fmhp.Marker, ...
                      'MarkerSize', fmhp.MarkerSize, ...
                      'MarkerEdgeColor', fmhp.MarkerEdgeColor, ...
                      'MarkerFaceColor', fmhp.MarkerFaceColor, ...
                      'SpecularStrength', fmhp.SpecularStrength, ...
                      'SpecularExponent', fmhp.SpecularExponent, ...
                      'SpecularColorReflectance', fmhp.SpecularColorReflectance, ...
                      'UserData', struct('SRF', f), 'Visible', 'on');
        else
            hsrf = patch(p(:, 1), p(:, 2), p(:, 3), 0, 'Parent', hgtrf);
            set(hsrf, ...
                'CData', fmhp.CData, ...
                'EdgeColor', fmhp.EdgeColor, ...
                'EdgeAlpha', fmhp.EdgeAlpha, ...
                'FaceAlpha', fmhp.FaceAlpha, ...
                'FaceColor', fmhp.FaceColor, ...
                'AmbientStrength', fmhp.AmbientStrength, ...
                'DiffuseStrength', fmhp.DiffuseStrength, ...
                'SpecularStrength', fmhp.SpecularStrength, ...
                'SpecularExponent', fmhp.SpecularExponent, ...
                'SpecularColorReflectance', fmhp.SpecularColorReflectance);
            set(hsrf, 'Visible', 'on');
        end

        % add to userdata
        scu(scc, 5:6) = {hsrf, fh.Stats};
    end

    % update UserData field (not in handle; now a struct!)
    htag.Scenery.UserData = scu;

    % show correct page, set slice pos
    ne_gcfg.cc.(iSat) = htag;
    ne_showsatpage(0, 0, iSat, cc.page);
    ne_setcsrfstatbars(0, 0, iSat);

    % delete zoomline axes
    delete(htag.ZoomLineX);
    delete(htag.ZoomLineY);

% for render view
elseif cc.page == 4

    % get selected SliceVar
    slvar = cc.SliceVar;

    % record type, create transimg objects, and keep required settings
    htag.Config.sattype = 'render';
    htag.Config.tio = struct('imRnd', transimg(528, 528));
    tio = htag.Config.tio;
    ne_gcfg.tio.satRnd = tio.imRnd;
    sethandle(tio.imRnd, get(tSat.(sprintf('IM_%s_Slice_Rend', iSat)).MLHandle, 'Children'));
    htag.Config.surfpos = hSat.Position;
    htag.Config.surfpos(1:2) = 0;
    set(tSat.(sprintf('AX_%s_Slice_Rend', iSat)).MLHandle, 'XTick', [], 'YTick', [], 'ZTick', []);

    % show correct page
    ne_gcfg.cc.(iSat) = htag;
    ne_showsatpage(0, 0, iSat, cc.page);

    % if slice and stats are valid xff
    if isxff(htag.Config.StatsVar, true)
        if numel(htag.Config.StatsVarIdx) == 1
            mapnames = htag.Config.StatsVar.MapNames(ne_gcfg.c.extmapnames);
            stmap = sprintf(' (Map %d: %s)', htag.Config.StatsVarIdx, ...
                mapnames{htag.Config.StatsVarIdx});
        else
            stmap = sprintf(' (%d maps%s)', numel(htag.Config.StatsVarIdx), ...
                sprintf(', %d', htag.Config.StatsVarIdx));
        end
        [stp, stf] = fileparts(htag.Config.StatsVar.FilenameOnDisk(true));
        if isempty(stf)
            stf = 'Unsaved';
        end
        stname = sprintf(' - %s%s', stf, stmap);
    else
        stname = '';
    end
    if isxff(slvar, true)
        [slp, slf] = fileparts(slvar.FilenameOnDisk(true));
        if isempty(slf)
            slf = 'Unsaved';
        end
    else
        slf = 'Empty';
    end
    hSat.Name = sprintf('NeuroElf GUI - %s%s', slf, stname);

    % force renderer to 'zbuffer'
    hSat.Renderer = 'zbuffer';

    % delete zoomline axes
    delete(htag.ZoomLineX);
    delete(htag.ZoomLineY);
end

% position
lastpos = cini.Satellites.Position;
if all(lastpos == -1)
    lastpos = ne_gcfg.h.MainFig.Position;
    lastpos = lastpos(1:2) + [lastpos(3) + 10, lastpos(4) - 508];
end
for fc = numel(bsc):-1:1
    if ~isempty(regexpi(bsc{fc}, '^bs')) && isfield(ne_gcfg.cc.(bsc{fc}), 'SatelliteMLH') && ...
        numel(ne_gcfg.cc.(bsc{fc}).SatelliteMLH) == 1 && ishandle(ne_gcfg.cc.(bsc{fc}).SatelliteMLH) && ...
        (isnumeric(ne_gcfg.cc.(bsc{fc}).SatelliteMLH) || isvalid(ne_gcfg.cc.(bsc{fc}).SatelliteMLH))
        testpos = get(ne_gcfg.cc.(bsc{fc}).SatelliteMLH, 'Position');
        if testpos(1:2) == lastpos
            lastpos = lastpos + [40, -40];
            break;
        end
    end
end
hSat.Position(1:2) = lastpos;

% resize
if any(mfsz > cc.fullsize)
    ne_satresize(0, 0, iSat, hSat.Position(3:4) + mfsz - cc.fullsize);
end

% finally, make satellite visible
hSat.HandleVisibility = 'callback';
hSat.Visible = 'on';

% allow resize
set(hSat.MLHandle, 'Resize', 'on');
hSat.ResizeFcn = {@ne_satresize, iSat};

% for output
if nargout > 0
    varargout{1} = hSat;
    if nargout > 1
        varargout{2} = htag;
        if nargout > 2
            varargout{3} = iSat;
        end
    end
end
