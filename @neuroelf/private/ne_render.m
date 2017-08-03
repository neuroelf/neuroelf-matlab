% FUNCTION ne_render: create 3D rendering
function varargout = ne_render(varargin)

% Version:  v1.0
% Build:    14103015
% Date:     Oct-30 2014, 3:17 PM EST
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% show render page
ne_showpage(0, 0, 4);

% only open if not open
ch = ne_gcfg.h;
if isstruct(ch.Render) && ...
    isfield(ch.Render, 'RendFig') && ...
    numel(ch.Render.RendFig) == 1 && ...
    isxfigure(ch.Render.RendFig, true)
    figure(ch.Render.RendFig.MLHandle);
    ne_render_setview;
    drawnow;
    return;
end

% if slvar invalid, do nothing
cc = ne_gcfg.fcfg;
if numel(cc.SliceVar) ~= 1 || ...
   ~isxff(cc.SliceVar)
    return;
end
slvar = cc.SliceVar;
slvarfile = slvar.FilenameOnDisk(2);
slvarrtv = slvar.RunTimeVars;
if isempty(slvarfile)
    slvarfile = sprintf('<untitled.%d>', lower(slvar.Filetype));
elseif numel(slvarfile) > 51
    slvarfile = [slvarfile(1:24) '...' slvarfile(end-23:end)];
end

% check other slicing vars
slvars = ne_gcfg.h.SliceVar.UserData;
slvarsn = {'<same>'};
slvarso = {[]};
for slvc = 1:size(slvars, 1)
    if isxff(slvars{slvc, 4}, true) && ...
        slvars{slvc, 4} ~= slvar
        slvarsn{end+1} = slvars{slvc, 1};
        slvarso{end+1} = slvars{slvc, 4};
    end
end

% echo
if ne_gcfg.c.echo
    ne_echo('neuroelf_gui(''render_ex'')');
end

% open configuration window
try
    hRendFig = xfigure([neuroelf_path('tfg') '/ne_render.tfg']);
    hTag = hRendFig.TagStruct;
    ne_gcfg.h.Render = struct( ...
        'RendFig', hRendFig, ...
        'h', hTag);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Error opening render configuration.', ...
        'NeuroElf GUI - warning', 'modal'));
    return;
end

% put figure into struct
ci = ne_gcfg.c.ini.Render;
ne_gcfg.fcfg.Render = struct( ...
    'actdepth',  ci.ActivationDepth, ...
    'agalpha',   ci.GrayAlphaFactor, ...
    'bgcol',     ci.BackColor, ...
    'colblend',  ci.ColorBlend, ...
    'dist',      ci.Distance, ...
    'filename',  ci.Filename, ...
    'filenmov',  ci.FilenameMovie, ...
    'filewmov',  ci.FileWriteMovie, ...
    'filewrt',   ci.FileWrite, ...
    'frame',     ci.Frame, ...
    'galpha',    ci.GrayAlpha, ...
    'hFig',      hRendFig, ...
    'hTag',      hTag, ...
    'hicol',     ci.HighColor, ...
    'imeth',     ci.StatsInterp, ...
    'imetha',    ci.AnatInterp, ...
    'join',      cc.join, ...
    'joinulay',  cc.joinulay, ...
    'locol',     ci.LowColor, ...
    'movhdfrom', ci.MovieHDFrom, ...
    'proty',     0, ...
    'protz',     0, ...
    'res',       ci.Resolution, ...
    'roty',      round(cc.srfcfg.angley), ...
    'rotz',      round(cc.srfcfg.anglex), ...
    'showinfig', ci.ShowWhileRendering, ...
    'slavar',    [], ...
    'slfrom',    ci.SliceFrom, ...
    'slstep',    ci.SliceStep, ...
    'slto',      ci.SliceTo, ...
    'slvar',     slvar, ...
    'smstat',    ci.SmoothStatBorder, ...
    'smstatk',   ci.SmoothStatBorderKeep, ...
    'stalp',     ci.StatsAlpha, ...
    'stminaalp', ci.StatsMinAnaAlpha, ...
    'stmulaalp', ci.StatsMultWithAnaAlpha, ...
    'stsmooth',  ci.SmoothStats, ...
    'stvar',     cc.StatsVar, ...
    'stvix',     cc.StatsVarIdx, ...
    'stvixo',    cc.StatsVarIdx);
rcfg = ne_gcfg.fcfg.Render;

% set color buttons
hTag.BT_render_btcol.CData = ...
    repmat(reshape(uint8(ci.BackColor), [1, 1, 3]), 16, 22);
hTag.BT_render_hicol.CData = ...
    repmat(reshape(uint8(ci.HighColor), [1, 1, 3]), 16, 22);
hTag.BT_render_locol.CData = ...
    repmat(reshape(uint8(ci.LowColor), [1, 1, 3]), 16, 22);
if strcmpi(ci.ColorBlend, 'hsv')
    hTag.DD_render_colblend.Value = 2;
else
    hTag.DD_render_colblend.Value = 1;
end

% LUT rendering?
if isfield(slvarrtv, 'GrayScaleLUT') && ...
    isequal(size(slvarrtv.GrayScaleLUT), [256, 3]) && ...
    isa(slvarrtv.GrayScaleLUT, 'double')
    hTag.DD_render_colblend.String = {'RGB blend'; 'HSV blend'; 'LUT coloring'};
    hTag.DD_render_colblend.Value = 3;
    ne_gcfg.fcfg.Render.colblend = 'lut';
end

% set initial values -> slicing object(s)
hTag.ED_render_slvarfile.String = ['  ' slvarfile];
hTag.DD_render_slavarfile.String = slvarsn;
hTag.DD_render_slavarfile.UserData = slvarso;
if numel(slvarsn) > 1
    hRendFig.SetGroupEnabled('SlAlpha', 'on');
end
hTag.ED_render_agalpha.String = sprintf('%.3g', rcfg.agalpha);

% -> slicing range, rotation, perspective, resolution
hTag.ED_render_xfrom.String = sprintf('%.6g', rcfg.slfrom);
hTag.ED_render_xto.String = sprintf('%.6g', rcfg.slto);
hTag.ED_render_xstep.String = sprintf('%.4g', rcfg.slstep);
hTag.ED_render_yzfrom.String = sprintf('%d', rcfg.frame(1));
hTag.ED_render_yzto.String = sprintf('%d', rcfg.frame(2));
hTag.ED_render_yzenith.String = sprintf('%d', rcfg.roty);
hTag.ED_render_zazimuth.String = sprintf('%d', rcfg.rotz);
hTag.ED_render_dist.String = sprintf('%d', rcfg.dist);
hTag.ED_render_res.String = sprintf('%d', rcfg.res);

% -> stats object
if numel(rcfg.stvar) ~= 1 || ...
   ~isxff(rcfg.stvar)
    stvarfile = '';
else
    hRendFig.SetGroupEnabled('StVar', 'on');
    stvarfile = rcfg.stvar.FilenameOnDisk(2);
    if isempty(stvarfile)
        stvarfile = sprintf('<untitled.%s>', lower(rcfg.stvar.Filetype));
    elseif numel(stvarfile) > 51
        stvarfile = [stvarfile(1:24) '...' stvarfile(end-23:end)];
    end
    if isxff(rcfg.stvar, {'cmp', 'hdr', 'head', 'vmp'})
        hRendFig.SetGroupEnabled('StVMP', 'on');
    end
end
if ~isempty(stvarfile)
    hTag.ED_render_stvarfile.String = ['  ' stvarfile];
end
hTag.ED_render_smstat.String = sprintf('%.3g', rcfg.smstat);
hTag.CB_render_smstatk.Value = double(rcfg.smstatk);

% -> alpha values map visibility/smoothing controls
hTag.ED_render_galpha.String = sprintf('%.3g', rcfg.galpha);
hTag.ED_render_stalp.String = sprintf('%.3g', rcfg.stalp);
hTag.CB_render_join.Value = double(rcfg.join);
hTag.CB_render_stmulaalp.Value = double(rcfg.stmulaalp);
hTag.ED_render_stminaalp.String = sprintf('%.2g', rcfg.stminaalp);
hTag.ED_render_stsmooth.String = sprintf('%.1g', rcfg.stsmooth);
hTag.ED_render_actdepth.String = sprintf('%.1g', rcfg.actdepth);

% -> interpolation types
interpm = {'linear', 'cubic', 'lanczos3'};
interpi = find(strcmpi(rcfg.imetha, interpm));
if isempty(interpi)
    interpi = 1;
    rcfg.imetha = interpm{interpi};
end
hTag.DD_render_avarinterp.Value = interpi;
interpi = find(strcmpi(rcfg.imeth, interpm));
if isempty(interpi)
    interpi = 1;
    rcfg.imeth = interpm{interpi};
end
hTag.DD_render_stvarinterp.Value = interpi;

% -> processing controls
hTag.CB_render_showinfig.Value = double(rcfg.showinfig);
hTag.ED_render_filename.String = ['  ' rcfg.filename];
hTag.CB_render_filewrt.Value = double(rcfg.filewrt);
hTag.ED_render_movfilename.String = ['  ' rcfg.filenmov];
hTag.CB_render_movfilewrt.Value = double(rcfg.filewmov);

% set callbacks
hTag.BT_render_bgcol.Callback = {@ne_render_updateui, 'setcolor', 'bgcol'};
hTag.BT_render_locol.Callback = {@ne_render_updateui, 'setcolor', 'locol'};
hTag.BT_render_hicol.Callback = {@ne_render_updateui, 'setcolor', 'hicol'};
hTag.DD_render_colblend.Callback = {@ne_render_updateui, 'check'};
hTag.DD_render_slavarfile.Callback = {@ne_render_updateui, 'slavar'};
hTag.ED_render_agalpha.Callback = {@ne_render_updateui, 'agalpha'};
hTag.ED_render_xfrom.Callback = {@ne_render_updateui, 'xfrom'};
hTag.ED_render_xto.Callback = {@ne_render_updateui, 'xto'};
hTag.ED_render_xstep.Callback = {@ne_render_updateui, 'xstep'};
hTag.ED_render_yzfrom.Callback = {@ne_render_updateui, 'yzfrom'};
hTag.ED_render_yzto.Callback = {@ne_render_updateui, 'yzto'};
hTag.ED_render_yzenith.Callback = {@ne_render_updateui, 'yrot'};
hTag.ED_render_zazimuth.Callback = {@ne_render_updateui, 'zrot'};
hTag.ED_render_dist.Callback = {@ne_render_updateui, 'dist'};
hTag.ED_render_res.Callback = {@ne_render_updateui, 'res'};
hTag.BT_render_smoothstat.Callback = @ne_render_smoothstats;
hTag.ED_render_smstat.Callback = {@ne_render_updateui, 'smstat'};
hTag.CB_render_smstatk.Callback = {@ne_render_updateui, 'check'};
hTag.ED_render_galpha.Callback = {@ne_render_updateui, 'galpha'};
hTag.ED_render_stalp.Callback = {@ne_render_updateui, 'stalp'};
hTag.CB_render_join.Callback = {@ne_render_updateui, 'check'};
hTag.CB_render_stmulaalp.Callback = {@ne_render_updateui, 'check'};
hTag.ED_render_stminaalp.Callback = {@ne_render_updateui, 'stminaalp'};
hTag.ED_render_stsmooth.Callback = {@ne_render_updateui, 'stsmooth'};
hTag.ED_render_actdepth.Callback = {@ne_render_updateui, 'actdepth'};
hTag.DD_render_avarinterp.Callback = {@ne_render_updateui, 'check'};
hTag.DD_render_stvarinterp.Callback = {@ne_render_updateui, 'check'};
hTag.CB_render_showinfig.Callback = {@ne_render_updateui, 'check'};
hTag.ED_render_filename.Callback = {@ne_render_updateui, 'filename'};
hTag.CB_render_filewrt.Callback = {@ne_render_updateui, 'check'};
hTag.ED_render_movfilename.Callback = {@ne_render_updateui, 'filenmov'};
hTag.CB_render_movfilewrt.Callback = {@ne_render_updateui, 'check'};
hTag.BT_render_cancel.Callback = @ne_render_closeui;
hTag.BT_render_ex.Callback = @ne_render_ex;
hRendFig.CloseRequestFcn = @ne_render_closeui;

% put to last position (if available)
try
    lastpos = ne_gcfg.c.ini.Children.RenderPosition;
    if any(lastpos ~= -1)
        hRendFig.Position(1:2) = lastpos;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% set figure visible and let the user interact...
hRendFig.HandleVisibility = 'callback';
hRendFig.Visible = 'on';
hRendFig.Pointer = 'watch';
drawnow;

% and update view
try
    ne_render_setview;
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% restore pointer
hRendFig.Pointer = 'arrow';
drawnow;
