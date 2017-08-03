% FUNCTION ne_surfmontage: create surface screenshot (in new figure)
function varargout = ne_surfmontage(varargin)

% Version:  v1.1
% Build:    16053009
% Date:     May-30 2016, 9:40 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015, 2016, Jochen Weber
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
if isstruct(ch.SurfMontage) && isfield(ch.SurfMontage, 'SMFig') && ...
    numel(ch.SurfMontage.SMFig) == 1 && isxfigure(ch.SurfMontage.SMFig)
    figure(ch.SurfMontage.SMFig.MLHandle);
    drawnow;
    return;
end
ci = ne_gcfg.c.ini;

% echo
if ne_gcfg.c.echo
    ne_echo('neuroelf_gui(''surfmontage'');');
end

% open configuration window
try
    hSurfFig = xfigure([neuroelf_path('tfg') '/ne_surfmontage.tfg']);
    hTag = hSurfFig.TagStruct;
    ne_gcfg.h.SurfMontage = struct('SMFig', hSurfFig, 'h', hTag);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Error opening montage configuration.', ...
        'NeuroElf GUI - warning', 'modal'));
    return;
end

% load data from ini
mc = ci.SurfMontage;
cc = ci.SurfMontageConfigs;
ec = ci.SurfMontageElems;

% put figure into struct
ne_gcfg.fcfg.SurfMontage = struct('cc', cc, 'ec', ec, 'hFig', hSurfFig, 'hTag', hTag, 'mc', mc);

% set callbacks
hTag.DD_surfmontage_configs.Callback = {@ne_surfmontage_ui, 'setconfig'};
hTag.BT_surfmontage_addcfg.Callback = {@ne_surfmontage_ui, 'addconfig'};
hTag.BT_surfmontage_delcfg.Callback = {@ne_surfmontage_ui, 'deleteconfig'};
hTag.BT_surfmontage_rencfg.Callback = {@ne_surfmontage_ui, 'renameconfig'};
hTag.BT_surfmontage_cpycfg.Callback = {@ne_surfmontage_ui, 'copyconfig'};
hTag.ED_surfmontage_imagex.Callback = {@ne_surfmontage_ui, 'checkfield', 'imagex'};
hTag.ED_surfmontage_imagey.Callback = {@ne_surfmontage_ui, 'checkfield', 'imagey'};
hTag.BT_surfmontage_bgcol.Callback = {@ne_surfmontage_ui, 'bgcolor'};
hTag.BT_surfmontage_bgcolc.Callback = {@ne_surfmontage_ui, 'bgcolor'};
hTag.CB_surfmontage_stbars.Callback = {@ne_surfmontage_ui, 'checkfield', 'stbars'};
hTag.LB_surfmontage_celems.Callback = {@ne_surfmontage_ui, 'selectelement'};
hTag.BT_surfmontage_delfelms.Callback = {@ne_surfmontage_ui, 'deletefromelements'};
hTag.BT_surfmontage_addtoelms.Callback = {@ne_surfmontage_ui, 'addtoelements'};
hTag.BT_surfmontage_downelm.Callback = {@ne_surfmontage_ui, 'moveelementdown'};
hTag.BT_surfmontage_upelm.Callback = {@ne_surfmontage_ui, 'moveelementup'};
hTag.DD_surfmontage_element.Callback = {@ne_surfmontage_ui, 'setelement'};
hTag.BT_surfmontage_addelem.Callback = {@ne_surfmontage_ui, 'addelement'};
hTag.BT_surfmontage_delelem.Callback = {@ne_surfmontage_ui, 'deleteelement'};
hTag.ED_surfmontage_smniter.Callback = {@ne_surfmontage_ui, 'checkfield', 'smniter'};
hTag.ED_surfmontage_smforce.Callback = {@ne_surfmontage_ui, 'checkfield', 'smforce'};
hTag.ED_surfmontage_infniter.Callback = {@ne_surfmontage_ui, 'checkfield', 'infniter'};
hTag.ED_surfmontage_infforce.Callback = {@ne_surfmontage_ui, 'checkfield', 'infforce'};
hTag.ED_surfmontage_transx.Callback = {@ne_surfmontage_ui, 'checkfield', 'transx'};
hTag.ED_surfmontage_transy.Callback = {@ne_surfmontage_ui, 'checkfield', 'transy'};
hTag.ED_surfmontage_azimuth.Callback = {@ne_surfmontage_ui, 'checkfield', 'azimuth'};
hTag.ED_surfmontage_zenith.Callback = {@ne_surfmontage_ui, 'checkfield', 'zenith'};
hTag.ED_surfmontage_zoom.Callback = {@ne_surfmontage_ui, 'checkfield', 'zoom'};
hTag.ED_surfmontage_time.Callback = {@ne_surfmontage_ui, 'checkfield', 'time'};
hTag.ED_surfmontage_elemx.Callback = {@ne_surfmontage_ui, 'checkfield', 'elemx'};
hTag.ED_surfmontage_elemy.Callback = {@ne_surfmontage_ui, 'checkfield', 'elemy'};
hTag.ED_surfmontage_offsetx.Callback = {@ne_surfmontage_ui, 'checkfield', 'offsetx'};
hTag.ED_surfmontage_offsety.Callback = {@ne_surfmontage_ui, 'checkfield', 'offsety'};
hTag.RB_surfmontage_showinfig.Callback = {@ne_surfmontage_ui, 'targetfigure'};
hTag.RB_surfmontage_writefile.Callback = {@ne_surfmontage_ui, 'targetfile'};
hTag.ED_surfmontage_filename.Callback = {@ne_surfmontage_ui, 'checkfield', 'filename'};
hTag.CB_surfmontage_samplevmp.Callback = {@ne_surfmontage_ui, 'smpvmp'};
hTag.ED_surfmontage_lowthresh.Callback = {@ne_surfmontage_ui, 'checkfield', 'lowthresh'};
hTag.ED_surfmontage_uppthresh.Callback = {@ne_surfmontage_ui, 'checkfield', 'uppthresh'};
hTag.BT_surfmontage_cancel.Callback = @ne_surfmontage_closeui;
hTag.BT_surfmontage_create.Callback = @ne_surfmontage_create;
hSurfFig.CloseRequestFcn = @ne_surfmontage_closeui;

% generate list of configs
ccf = fieldnames(cc);
for cci = 1:numel(ccf)
    ccf{cci} = cc.(ccf{cci}){1};
end
if isempty(ccf)
    ccf = {'<no configurations available>'};
end
hTag.DD_surfmontage_configs.String = ccf;
hTag.DD_surfmontage_configs.Value = 1;
if strcmp(ccf{1}, '<no configurations available>')
    hSurfFig.SetGroupEnabled('Configs', 'off');
else
    hSurfFig.SetGroupEnabled('Configs', 'on');    
end
ne_surfmontage_ui(0, 0, 'genelemlist');
hTag.DD_surfmontage_element.Value = 1;
if strcmp(ccf{1}, '<no configurations available>')
    hSurfFig.SetGroupEnabled('Configs', 'off');
else
    hSurfFig.SetGroupEnabled('Configs', 'on');    
end
ecf = hTag.DD_surfmontage_element.String;
if strcmp(ecf{1}, '<no surface elements available>')
    hSurfFig.SetGroupEnabled('Elem', 'off');
else
    hSurfFig.SetGroupEnabled('Elem', 'on');
end

% set global fields content
if mc.WriteToFile
    hSurfFig.RadioGroupSetOne('VisMOut', 2);
    hTag.ED_surfmontage_filename.Enable = 'on';
else
    hSurfFig.RadioGroupSetOne('VisMOut', 1);
    hTag.ED_surfmontage_filename.Enable = 'off';
end
hTag.ED_surfmontage_filename.String = [' ' mc.WriteFilename];
hTag.CB_surfmontage_samplevmp.Value = double(mc.SampleVMP);
hTag.CB_surfmontage_voisonly.Value = double(mc.RestrictVMPsToVOIs);
hTag.CB_surfmontage_lowthresh.Value = double(mc.MultiplyThresholds);
hTag.ED_surfmontage_lowthresh.String = sprintf('%g', mc.ThresholdFactors(1));
hTag.ED_surfmontage_uppthresh.String = sprintf('%g', mc.ThresholdFactors(2));

% update configs and everything
ne_surfmontage_ui(0, 0, 'setconfig');

% put to last position (if available)
try
    lastpos = ne_gcfg.c.ini.Children.SurfMontagePosition;
    if any(lastpos ~= -1)
        hSurfFig.Position(1:2) = lastpos;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% set figure visible and let the user interact...
hSurfFig.HandleVisibility = 'callback';
hSurfFig.Visible = 'on';
