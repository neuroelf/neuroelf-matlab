function varargout = ne_tcundock(varargin)
% ne_tcundock  - undock the time course plot into (resizable) window UI
%
% FORMAT:       [hSat, tags, iSat] = ne_tcundock(SRC, EVT)
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
% Build:    16052611
% Date:     May-26 2016, 11:08 AM EST
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% if already undocked, make visible and bring forward
cf = fieldnames(ne_gcfg.cc);
if ~isempty(cf)
    cf(cellfun('isempty', regexp(cf, '^TC'))) = [];
end
for fc = 1:numel(cf)
    ccc = ne_gcfg.cc.(cf{fc});
    if numel(ccc) == 1 && isstruct(ccc) && isfield(ccc, 'Satellite') && ...
        numel(ccc.Satellite) == 1 && isxfigure(ccc.Satellite, true)
        ccc.Satellite.Visible = 'on';
        figure(ccc.Satellite.MLHandle);
        drawnow;
        return;
    else
        ne_gcfg.cc = rmfield(ne_gcfg.cc, cf{fc});
    end
end

% only for pages 1, 2, and 4
if ~any(cc.page == [1, 2, 4])
    return;
end

% if slice and stats are valid xff
if isxff(cc.SliceVar, {'fmr', 'hdr', 'head', 'vtc'}) && cc.SliceVar.NrOfVolumes > 1
    tcvar = cc.SliceVar;
elseif isxff(cc.StatsVar, {'vtc'})
    tcvar = cc.StatsVar;
else
    return;
end

% echo
if ne_gcfg.c.echo
    ne_echo('neuroelf_gui(''tcundock'');');
end

% open satellite window
hSat = xfigure([neuroelf_path('tfg') '/ne_tcplot.tfg']);
tSat = hSat.Tag;
iSat = tSat(1:8);

% create tags
tSat = hSat.TagStruct;
htag = struct;
htag.Satellite = hSat;
htag.SatelliteMLH = hSat.MLHandle;
htag.TCPlot = tSat.(sprintf('AX_%s_TCPlot', iSat));
tcplot = htag.TCPlot;
tcplotmlh = tcplot.MLHandle;
htag.TCPlotChild = plot(tcplotmlh, (1:120)', zeros(120, 1));
set(tcplotmlh, 'HandleVisibility', 'off');
hold(tcplotmlh, 'on');
htag.TCPlotChildren = [];
htag.TCPlotDiscards = image([1, 120], [0, 0], ...
    uint8(repmat(reshape([255, 0, 0], [1, 1, 3]), 1, 120)), 'Parent', tcplotmlh);
set(htag.TCPlotDiscards, 'AlphaData', zeros(1, 120));

% callbacks
hSat.CloseRequestFcn = {@ne_closesatwindow, iSat};
hSat.KeyPressFcn = {@ne_keypress, iSat};
hSat.KeyReleaseFcn = @ne_keyrelease;
hSat.ResizeFcn = {@ne_satresize, iSat};
hSat.WindowButtonDownFcn = {@ne_tcbtdown, iSat};
hSat.WindowButtonMotionFcn = {@ne_tcbtmove, iSat};
hSat.WindowButtonUpFcn = {@ne_btup, iSat};
hSat.Renderer = cc.renderer;

% set configuration
htag.Config = cc;
htag.Config.mods = {};
htag.Config.mpos = struct('cur', [0, 0], 'ddat', {{}}, 'down', [-1, -1], 'last', [0, 0], 'mods', {{}});
htag.Config.sattag = iSat;
htag.Config.tcvar = tcvar;

% record type, create transimg objects, and keep required settings
htag.Config.sattype = 'tcplot';

% show correct page
ne_gcfg.cc.(iSat) = htag;

% filename for window title
[slp, slf, slfe] = fileparts(tcvar.FilenameOnDisk);
if isempty(slf)
    slf = 'unsaved';
    slfe = ['.' lower(tcvar.Filetype)];
end
slf = [slf slfe];
hSat.Name = sprintf('NeuroElf - TC plot: %s', slf);

% force renderer to 'zbuffer'
hSat.Renderer = 'zbuffer';

% allow resize
set(hSat.MLHandle, 'Resize', 'on');

% call to setslicepos (to re-draw TCPlot)
ne_setslicepos;
ch.TCPlot.Visible = 'off';

% finally, make satellite visible
hSat.HandleVisibility = 'callback';
hSat.Visible = 'on';

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
