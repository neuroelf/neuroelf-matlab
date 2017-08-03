% FUNCTION ne_rm_open: loads the RFX mediation UI
function varargout = ne_rm_open(varargin)

% Version:  v1.1
% Build:    16012618
% Date:     Jan-26 2016, 6:41 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% already open?
if isfield(ne_gcfg.h, 'RM') && ...
    isstruct(ne_gcfg.h.RM) && ...
    isfield(ne_gcfg.h.RM, 'RMFig') && ...
    numel(ne_gcfg.h.RM.RMFig) == 1 && ...
    isxfigure(ne_gcfg.h.RM.RMFig, true)
    ne_gcfg.h.RM.RMFig.Visible = 'on';
    figure(ne_gcfg.h.RM.RMFig.MLHandle);
    return;
end

% blocked?
if any(strcmp(ne_gcfg.c.blockcb, 'rm_open'))
    return;
end

% echo
if ne_gcfg.c.echo
    ne_echo('neuroelf_gui(''rm_open'');');
end

% get list of loaded GLMs
nelfo = xff(0, 'objects');
isglm = false(size(nelfo));
for oc = 1:numel(nelfo)
    if strcmpi(nelfo(oc).Filetype, 'glm') && ...
        isfield(nelfo(oc).RunTimeVars, 'Contrasts') && ...
        isfield(nelfo(oc).RunTimeVars, 'CovariatesNames') && ...
       ~isempty(nelfo(oc).RunTimeVars.Contrasts) && ...
       ~isempty(nelfo(oc).RunTimeVars.CovariatesNames)
        isglm(oc) = true;
    end
end
nelfo(~isglm) = [];

% return if no GLMs found
if isempty(nelfo)
    uiwait(msgbox('A GLM file with configured contrasts and covariates must be loaded.', ...
        'NeuroElf GUI - info', 'modal'));
    return;
end

% load contrast manager
try
    hFig = xfigure([neuroelf_path('tfg') '/ne_rfxmediation.tfg']);
    ne_gcfg.h.RM.RMFig = hFig;
    tags = hFig.TagStruct;

    % block further callbacks
    ne_gcfg.c.blockcb{end+1} = 'rm_open';

    % last selected GLM is currently empty!
    hFig.UserData = struct('lastglm', []);

    % get required controls
    ch = struct;
    ch.GLMs = tags.DD_NeuroElf_RM_GLMs;
    ch.Subjects = tags.LB_NeuroElf_RM_subjects;
    ch.StrBoot = tags.RB_NeuroElf_RM_abboot;
    ch.StrMCMAM = tags.RB_NeuroElf_RM_abmcmam;
    ch.StrSobel = tags.RB_NeuroElf_RM_absobel;
    ch.Robust = tags.CB_NeuroElf_RM_robust;
    ch.Ztrans = tags.CB_NeuroElf_RM_ztrans;
    ch.BootNumSmp = tags.ED_NeuroElf_RM_numsmp;
    ch.BootReuse = tags.CB_NeuroElf_RM_bootre;
    ch.BootSEPerc = tags.RB_NeuroElf_RM_bsperc;
    ch.BootSEBCa = tags.RB_NeuroElf_RM_bsbca;
    ch.BootSEVar = tags.RB_NeuroElf_RM_bsvar;
    ch.XType = tags.DD_NeuroElf_RM_Xtype;
    ch.XList = tags.LB_NeuroElf_RM_Xlist;
    ch.YType = tags.DD_NeuroElf_RM_Ytype;
    ch.YList = tags.LB_NeuroElf_RM_Ylist;
    ch.MCons = tags.LB_NeuroElf_RM_Mconlist;
    ch.MCovs = tags.LB_NeuroElf_RM_Mcovlist;
    ch.CCons = tags.LB_NeuroElf_RM_Cconlist;
    ch.CCovs = tags.LB_NeuroElf_RM_Ccovlist;

    % put list of loaded GLMs into dropdown
    glmfiles = cell(numel(nelfo), 1);
    glmobjs = cell(size(glmfiles));
    for c = 1:numel(glmfiles)
        [gpath, glmfiles{c}] = fileparts(nelfo(c).FilenameOnDisk);
        glmobjs{c} = nelfo(c);
        if isempty(glmfiles{c})
            glmfiles{c} = sprintf('<#%d: %d subs, %d preds>', ...
                glmobjs{c}.Filenumber, glmobjs{c}.NrOfSubjects, ...
                glmobjs{c}.NrOfSubjectPredictors);
        end
    end
    ne_gcfg.fcfg.RM.GLMs = glmobjs;
    ch.GLMs.String = glmfiles(:);

    % which to select
    ch.GLMs.Value = 1;
    if isxff(ne_gcfg.fcfg.StatsVar, 'glm')
        for c = 1:numel(glmobjs)
            if glmobjs{c} == ne_gcfg.fcfg.StatsVar
                ch.GLMs.Value = c;
                break;
            end
        end
    end

    % set callbacks -> GLM selection
    ch.GLMs.Callback = @ne_rm_setglm;
    hFig.UserData.lastglm = glmobjs{ch.GLMs.Value};

    % -> options
    ch.StrBoot.Callback = {@ne_rm_select, 'boot'};
    ch.StrMCMAM.Callback = {@ne_rm_select, 'mcmam'};
    ch.StrSobel.Callback = {@ne_rm_select, 'sobel'};
    ch.Robust.Callback = {@ne_rm_select, 'robust'};
    ch.Ztrans.Callback = {@ne_rm_select, 'ztrans'};
    ch.BootNumSmp.Callback = {@ne_rm_select, 'numsmp'};
    ch.BootReuse.Callback = {@ne_rm_select, 'breuse'};
    ch.BootSEPerc.Callback = {@ne_rm_select, 'perc'};
    ch.BootSEBCa.Callback = {@ne_rm_select, 'bca'};
    ch.BootSEVar.Callback = {@ne_rm_select, 'var'};

    % -> switch X/Y type
    ch.XType.Callback = {@ne_rm_select, 'xtype'};
    ch.YType.Callback = {@ne_rm_select, 'ytype'};

    % -> select single entry
    ch.XList.Callback = {@ne_rm_select, 'xlist'};
    ch.YList.Callback = {@ne_rm_select, 'ylist'};

    % -> select multiple entry
    ch.MCons.Callback = {@ne_rm_select, 'mcons'};
    ch.MCovs.Callback = {@ne_rm_select, 'mcovs'};
    ch.CCons.Callback = {@ne_rm_select, 'ccons'};
    ch.CCovs.Callback = {@ne_rm_select, 'ccovs'};

    % -> do the work
    tags.BT_NeuroElf_RM_cancel.Callback = @ne_rm_closeui;
    tags.BT_NeuroElf_RM_compute.Callback = @ne_rm_compute;

    % set h struct
    ne_gcfg.h.RM.h = ch;

    % load settings from ini
    ci = ne_gcfg.c.ini.Mediation;
    ch.BootNumSmp.String = ...
        sprintf('%d', round(max(500, min(1000000, ci.BootNumSmp))));
    ch.BootReuse.Value = double(ci.BootReuse);
    switch (lower(ci.BootSE))
        case {'perc'}
            hFig.RadioGroupSetOne('bstrap', 1);
        case {'bca'}
            hFig.RadioGroupSetOne('bstrap', 2);
        case {'var'}
            hFig.RadioGroupSetOne('bstrap', 3);
    end
    ch.Robust.Value = double(ci.Robust);
    switch (lower(ci.Strategy))
        case {'boot'}
            hFig.RadioGroupSetOne('abmeth', 1);
            ne_rm_select(0, 0, 'boot');
        case {'mcmam'}
            hFig.RadioGroupSetOne('abmeth', 2);
            ne_rm_select(0, 0, 'mcmam');
        case {'sobel'}
            hFig.RadioGroupSetOne('abmeth', 3);
            ne_rm_select(0, 0, 'sobel');
    end
    ch.Ztrans.Value = double(ci.ZTrans);

    % get position from ini if any good
    try
        lastpos = ne_gcfg.c.ini.Children.RMPosition;
        if any(lastpos ~= -1)
            hFig.Position(1:2) = lastpos;
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end

    % set visible, modal and wait
    hFig.CloseRequestFcn = @ne_rm_closeui;
    hFig.HandleVisibility = 'callback';
    hFig.Visible = 'on';

    % set current GLM, then update (.ListboxTop fields otherwise unheeded!)
    ne_rm_setglm;
    drawnow;

% give warning
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Error loading RFX mediation UI.', 'NeuroElf GUI - error', 'modal'));
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'rm_open')) = [];
end
