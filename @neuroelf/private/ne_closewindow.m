function varargout = ne_closewindow(varargin)
% ne_closewindow  - request and then close the main UI window
%
% FORMAT:       ne_closewindow(SRC, EVT [, clearflag])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       clearflag   if given, do not show request dialog; additionally,
%                   if set to 'Yes & clear loaded', clear objects
%
% No output fields.
%
% Example:
%
%     ne_closewindow(0, 0, 'yes')

% Version:  v1.1
% Build:    16060912
% Date:     Jun-09 2016, 12:08 PM EST
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% initial number of arguments
nnargin = nargin;

% reject if in callback or procedure
if ~isempty(ne_gcfg.c.blockcb)
    asw = questdlg(['A component depending on this window has not been closed!', ...
        char([10, 10]), 'Are you sure you want to close NeuroElf GUI?'], ...
        'NeuroElf GUI close request...', 'No', 'Yes', 'Yes & clear loaded', 'No');
    if strcmpi(asw, 'no')
        return;
    end
    nnargin = 3;
    varargin{3} = asw;
end

% without argument, ask: are you sure...?
if nnargin < 3
    asw = questdlg('Are you sure you want to close NeuroElf GUI?', ...
        'NeuroElf GUI close request...', 'No', 'Yes', 'Yes & clear loaded', 'Yes');
    if strcmpi(asw, 'no')
        return;
    end

% if a char argument is given
elseif ischar(varargin{3})
    asw = varargin{3}(:)';

% otherwise assume only to close but not clear vars in workspace
else
    asw = 'Yes';
end

% since we're closing shop, hide all menus, disable controls
ch = ne_gcfg.h;
fch = ch.MainFig.Children;
fcht = get(fch, 'Type');
set(fch(strcmpi(fcht, 'uimenu')), 'Visible', 'off');
set(fch(strcmpi(fcht, 'uicontrol')), 'Enable', 'off');

% disallow further close/resize requests
ch.MainFig.CloseRequestFcn = '';
ch.MainFig.ResizeFcn = '';
ch.MainFig.Resize = 'off';

% and show only progress bar
mainpos = ch.MainFig.Position;
progpos = ch.Progress.Position;
ch.MainFig.Position = [mainpos(1:2) + 0.5 .* (mainpos(3:4) - progpos(3:4)), ...
    progpos(3:4) + 12];
ch.Progress.Position = progpos - [6, 6, 0, 0];
ch.MainFig.Name = 'NeuroElf GUI - cleaning up...';
drawnow;

% close any referenced files
srvaru = ch.StatsVarRefs.UserData;
if ~isempty(srvaru) && ...
    numel(srvaru{1}) == 1 && ...
    isxff(srvaru{1})
    try
        clearxffobjects(srvaru(1));
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
srvaru = ch.SurfStatsVarRefs.UserData;
if ~isempty(srvaru) && ...
    numel(srvaru{1}) == 1 && ...
    isxff(srvaru{1})
    try
        clearxffobjects(srvaru(1));
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end

% close the matlab handle children (if any) and the main figure
try
    delete(ch.fChild);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% also close satellites
bfc = fieldnames(ne_gcfg.cc);
for bfn = bfc(:)'
    hc = ne_gcfg.cc.(bfn{1});
    if isstruct(hc) && ...
        numel(hc) == 1 && ...
        isfield(hc, 'Satellite') && ...
        numel(hc.Satellite) == 1 && ...
    	isxfigure(hc.Satellite, true)
        ne_closesatwindow(0, 0, bfn{1});
    end
end

% then iterate over any xfigure children
bfc = fieldnames(ch.Children);
for bfn = bfc(:)'
    if ~isempty(ch.Children.(bfn{1}))
        hCFig = ch.Children.(bfn{1});
        clrq = hCFig.CloseRequestFcn;
        figure(hCFig.MLHandle);
        if ischar(clrq) && ...
           ~isempty(clrq)
            evalin('base', clrq);
        elseif isa(clrq, 'function_handle')
            feval(clrq, hCFig.MLHandle, 0);
        elseif iscell(clrq) && ...
           ~isempty(clrq) && ...
            isa(clrq{1}, 'function_handle')
            feval(clrq{1}, hCFig.MLHandle, 0, clrq{2:end});
        else
            hCFig.Delete;
        end
    end
end

% then close special sub-UIs
if isstruct(ch.AC) && ...
    isfield(ch.AC, 'ACFig') && ...
    numel(ch.AC.ACFig) == 1 && ...
    isxfigure(ch.AC.ACFig, true)
    try
        ne_ancova_ui(0, 0, 'closeui');
        drawnow;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
if isstruct(ch.AP) && ...
    isfield(ch.AP, 'APFig') && ...
    numel(ch.AP.APFig) == 1 && ...
    isxfigure(ch.AP.APFig, true)
    try
        ne_audioplayer(0, 0, 'close');
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
if isstruct(ch.RM) && ...
    isfield(ch.RM, 'RMFig') && ...
    numel(ch.RM.RMFig) == 1 && ...
    isxfigure(ch.RM.RMFig, true)
    try
        ne_rm_closeui;
        drawnow;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
if isstruct(ch.CM) && ...
    isfield(ch.CM, 'CMFig') && ...
    numel(ch.CM.CMFig) == 1 && ...
    isxfigure(ch.CM.CMFig, true)
    try
        ne_cm_closeui;
        drawnow;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
if isstruct(ch.MDM) && ...
    isfield(ch.MDM, 'MDMFig') && ...
    numel(ch.MDM.MDMFig) == 1 && ...
    isxfigure(ch.MDM.MDMFig, true)
    try
        ne_mdm_closeui;
        drawnow;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
if isstruct(ch.MKDA) && ...
    isfield(ch.MKDA, 'MKDAFig') && ...
    numel(ch.MKDA.MKDAFig) == 1 && ...
    isxfigure(ch.MKDA.MKDAFig, true)
    try
        ne_mkda_closeui;
        drawnow;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
if isstruct(ch.Render) && ...
    isfield(ch.Render, 'RendFig') && ...
    numel(ch.Render.RendFig) == 1 && ...
    isxfigure(ch.Render.RendFig, true)
    try
        ne_render_closeui;
        drawnow;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
if isstruct(ch.SurfMontage) && ...
    isfield(ch.SurfMontage, 'SMFig') && ...
    numel(ch.SurfMontage.SMFig) == 1 && ...
    isxfigure(ch.SurfMontage.SMFig, true)
    try
        ne_surfmontage_closeui;
        drawnow;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
if isstruct(ch.VisMontage) && ...
    isfield(ch.VisMontage, 'VMFig') && ...
    numel(ch.VisMontage.VMFig) == 1 && ...
    isxfigure(ch.VisMontage.VMFig, true)
    try
        ne_vismontage_closeui;
        drawnow;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end

% record the last known position
ne_gcfg.c.ini.MainFig.Position = mainpos(1:2);
ne_gcfg.c.ini.MainFig.Size = mainpos(3:4);

% clear storage
if strcmpi(asw, 'yes & clear loaded')

    % go through variables in workspace
    sf = fieldnames(ne_gcfg.wc);
    for sfc = 1:numel(sf)

        % if loaded by tool, and still a valid object
        if ne_gcfg.wc.(sf{sfc}) && isfield(ne_gcfg.w, sf{sfc}) && ...
            numel(ne_gcfg.w.(sf{sfc})) == 1 && isxff(ne_gcfg.w.(sf{sfc}), true)
            try
                sfsfc = ne_gcfg.w.(sf{sfc});
                if ~ne_closefile(0, 0, sfsfc, 'final')
                    sfsfc.ClearObject;
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        elseif isfield(ne_gcfg.w, sf{sfc}) && numel(ne_gcfg.w.(sf{sfc})) == 1 && ...
            isxff(ne_gcfg.w.(sf{sfc}), true)
            try
                ne_closefile(0, 0, ne_gcfg.w.(sf{sfc}), 'final');
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end
    end

% only closing in GUI
elseif strcmpi(asw, 'yes')
    sf = fieldnames(ne_gcfg.w);
    for sfc = 1:numel(sf)
        if numel(ne_gcfg.w.(sf{sfc})) == 1 && isxff(ne_gcfg.w.(sf{sfc}), true)
            ne_closefile(0, 0, ne_gcfg.w.(sf{sfc}), 'final');
        end
    end
end

% shut down remote
if ne_gcfg.c.remote
    ne_remote(0, 0, 'unlisten');
    drawnow;
end

% then the main figure
try
    ch.MainFig.Delete;
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% clear LUT
try
    if isxff(ne_gcfg.lut)
        ne_gcfg.lut.ClearObject;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% clear POI
try
    if isxff(ne_gcfg.poi)
        ne_gcfg.poi.ClearObject;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% clear VOI
try
    if isxff(ne_gcfg.voi)
        voih = handles(ne_gcfg.voi);
        if isfield(voih, 'RGBImage') && ...
            numel(voih.RGBImage) == 1 && ...
            isxff(voih.RGBImage, 'hdr')
            voih.RGBImage.ClearObject;
        end
        ne_gcfg.voi.ClearObject;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% clear transimg objects
delete(ne_gcfg.tio.imSag);
delete(ne_gcfg.tio.imCor);
delete(ne_gcfg.tio.imTra);
delete(ne_gcfg.tio.imSlZ);
delete(ne_gcfg.tio.imRnd);

% save and release ini file
try
    ne_gcfg.c.ini.Save;
    ne_gcfg.c.ini.Release('force');
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% if working in remote mode
if ne_gcfg.c.remote

    % just keep stub
    ne_gcfg = struct( ...
        'c',    ne_gcfg.c, ...
        'cc',   [], ...
        'fcfg', [], ...
        'h',    ne_gcfg.h, ...
        'lut',  [], ...
        'poi',  [], ...
        'tio',  [], ...
        'voi',  [], ...
        'w',    ne_gcfg.w, ...
        'wc',   []);

    % but END remote after next turn!
    ne_gcfg.c.remote = false;

% not in remote mode
else

    % clear global variable completely
    ne_gcfg(:) = [];
end
