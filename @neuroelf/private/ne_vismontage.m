% FUNCTION ne_vismontage: create screenshot (in new figure) of slicing
function varargout = ne_vismontage(varargin)

% Version:  v1.0
% Build:    15122814
% Date:     Dec-28 2015, 2:16 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2015, Jochen Weber
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

% only open if not open
ch = ne_gcfg.h;
if isstruct(ch.VisMontage) && ...
    isfield(ch.VisMontage, 'VMFig') && ...
    numel(ch.VisMontage.VMFig) == 1 && ...
    isxfigure(ch.VisMontage.VMFig)
    figure(ch.VisMontage.VMFig.MLHandle);
    drawnow;
    return;
end

% echo
if ne_gcfg.c.echo
    ne_echo('neuroelf_gui(''vismontage'');');
end

% open configuration window
try
    hVisFig = xfigure([neuroelf_path('tfg') '/ne_vismontage.tfg']);
    hTag = hVisFig.TagStruct;
    ne_gcfg.h.VisMontage = struct( ...
        'VMFig', hVisFig, ...
        'h', hTag);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Error opening montage configuration.', ...
        'NeuroElf GUI - warning', 'modal'));
    return;
end

% put figure into struct
ne_gcfg.fcfg.VisMontage = struct( ...
    'dir',    'sag', ...
    'frame',  [128, 128, 128; -127.9999, -127.9999, -127.9999], ...
    'hFig',   hVisFig, ...
    'hTag',   hTag, ...
    'layout', [1, 1]);

% try to get resolution of stats var
stvar = ne_gcfg.fcfg.StatsVar;
stres = 3;
if ~isempty(stvar) && ...
    isxff(stvar, true)
    try
        stres = round(mean(stvar.BoundingBox.ResXYZ));
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        stres = 3;
    end
end

% set initial values
steps = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; 'S'};
hTag.DD_vismontage_dir.String = {'sagittal'; 'coronal'; 'transversal'};
hTag.DD_vismontage_xstep.String = steps;
hTag.DD_vismontage_xstep.Value = stres;
hTag.DD_vismontage_xstep.Enable = 'on';
hTag.DD_vismontage_ystep.String = steps;
hTag.DD_vismontage_ystep.Value = stres;
hTag.DD_vismontage_zstep.String = steps;
hTag.DD_vismontage_zstep.Value = stres;
hTag.DD_vismontage_pixpvox.String = {'0.5'; '1'; '2'; '3'; '4'; '6'; '8'};
hTag.DD_vismontage_pixpvox.Value = 2;
hTag.DD_vismontage_imgborder.String = [{'0'}; steps(1:10, 1)];
hTag.ED_vismontage_filename.Enable = 'off';
hTag.DD_vismontage_brainbox.String = {'Full box', 'Stats box', 'MNI brain', 'AFNI brain'};

% set callbacks
hTag.DD_vismontage_dir.Callback = {@ne_vismontage_updateui, 'dir'};
hTag.ED_vismontage_xfrom.Callback = {@ne_vismontage_updateui, 'xfrom'};
hTag.ED_vismontage_xto.Callback = {@ne_vismontage_updateui, 'xto'};
hTag.DD_vismontage_xstep.Callback = {@ne_vismontage_updateui, 'xstep'};
hTag.ED_vismontage_yfrom.Callback = {@ne_vismontage_updateui, 'yfrom'};
hTag.ED_vismontage_yto.Callback = {@ne_vismontage_updateui, 'yto'};
hTag.DD_vismontage_ystep.Callback = {@ne_vismontage_updateui, 'ystep'};
hTag.ED_vismontage_zfrom.Callback = {@ne_vismontage_updateui, 'zfrom'};
hTag.ED_vismontage_zto.Callback = {@ne_vismontage_updateui, 'zto'};
hTag.DD_vismontage_zstep.Callback = {@ne_vismontage_updateui, 'zstep'};
hTag.CB_vismontage_anatransp.Callback = {@ne_vismontage_bgcol, 'chbox'};
hTag.BT_vismontage_anabackgc.Callback = {@ne_vismontage_bgcol, 'pick'};
hTag.RB_vismontage_showinfig.Callback = {@ne_vismontage_updateui, 'show'};
hTag.BT_vismontage_fontcolor.Callback = {@ne_vismontage_updateui, 'fontcolor'};
hTag.RB_vismontage_writefile.Callback = {@ne_vismontage_updateui, 'write'};
hTag.DD_vismontage_brainbox.Callback = {@ne_vismontage_updateui, 'brain'};
hTag.BT_vismontage_cancel.Callback = @ne_vismontage_closeui;
hTag.BT_vismontage_create3s.Callback = @ne_vismontage_create3s;
hTag.BT_vismontage_create.Callback = @ne_vismontage_create;
hVisFig.CloseRequestFcn = @ne_vismontage_closeui;

% pretend we're updating xfrom
hTag.DD_vismontage_brainbox.Value = 3;
ne_vismontage_updateui(0, 0, 'brain');
ne_vismontage_updateui(0, 0, 'xfrom');

% put to last position (if available)
try
    lastpos = ne_gcfg.c.ini.Children.VisMontagePosition;
    if any(lastpos ~= -1)
        hVisFig.Position(1:2) = lastpos;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% set figure visible and let the user interact...
hVisFig.HandleVisibility = 'callback';
hVisFig.Visible = 'on';
