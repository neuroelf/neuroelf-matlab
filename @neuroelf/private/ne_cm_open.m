% FUNCTION ne_cm_open: loads the contrast manager (CM)
function varargout = ne_cm_open(varargin)

% Version:  v1.1
% Build:    20031214
% Date:     Mar-12 2020, 2:36 PM EST
% Author:   Jochen Weber, NeuroElf.net
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2020, Jochen Weber
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
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get list of loaded GLMs
nelfo = xff(0, 'objects');
isglm = false(size(nelfo));
for oc = 1:numel(nelfo)
    if strcmpi(nelfo(oc).Filetype, 'glm')
        isglm(oc) = true;
    end
end
nelfo(~isglm) = [];

% already open?
if isfield(ch, 'CM') && isstruct(ch.CM) && isfield(ch.CM, 'CMFig') && ...
    numel(ch.CM.CMFig) == 1 && isxfigure(ch.CM.CMFig, true)

    % update GLMs in list
    glmfiles = cell(numel(nelfo), 1);
    glmobjs = cell(size(glmfiles));
    for c = 1:numel(glmfiles)
        [gpath, glmfiles{c}] = fileparts(nelfo(c).FilenameOnDisk);
        glmobjs{c} = nelfo(c);
        if isempty(glmfiles{c})
            glmfiles{c} = sprintf('<#%d: %d subs, %d preds>', ...
                glmobjs{c}.Filenumber, glmobjs{c}.NrOfSubjects, ...
                numel(glmobjs{c}.SubjectPredictors));
        end
    end
    ne_gcfg.fcfg.CM.GLMs = glmobjs;
    ne_gcfg.h.CM.h.GLMs.String = glmfiles(:);

    % which to select
    ne_gcfg.h.CM.h.GLMs.Value = 1;
    if numel(ne_gcfg.fcfg.CM.GLMs) > 1 && isxff(ne_gcfg.fcfg.StatsVar, 'glm')
        for c = 1:numel(glmobjs)
            if glmobjs{c} == ne_gcfg.fcfg.StatsVar
                ne_gcfg.h.CM.h.GLMs.Value = c;
                ne_cm_setglm();
                break;
            end
        end
    end

    % set visible and bring to the front
    ch.CM.CMFig.Visible = 'on';
    figure(ch.CM.CMFig.MLHandle);
    return;
end

% blocked?
if any(strcmp(ne_gcfg.c.blockcb, 'cm_open'))
    return;
end

% echo
if ne_gcfg.c.echo
    ne_echo('neuroelf_gui(''cm_open'');');
end

% return if no GLMs found
if isempty(nelfo)
    uiwait(msgbox('You must load a GLM file first to use the contrast manager.', ...
        'NeuroElf GUI - info', 'modal'));
    return;
end

% load contrast manager
try
    hFig = xfigure([neuroelf_path('tfg') '/ne_glmcontrasts.tfg']);
    ne_gcfg.h.CM.CMFig = hFig;

    % block further callbacks
    ne_gcfg.c.blockcb{end+1} = 'cm_open';

    % get required controls
    tags = hFig.TagStruct;
    ch = struct;
    ch.GLMs = tags.DD_NeuroElf_CM_GLMs;
    ch.UseGroups = tags.CB_NeuroElf_CM_groups;
    ch.Groups = tags.DD_NeuroElf_CM_groups;
    ch.GroupCmp = tags.RB_NeuroElf_CM_grpcmp;
    ch.GroupSep = tags.RB_NeuroElf_CM_grpsep;
    ch.Subjects = tags.LB_NeuroElf_CM_subjects;
    ch.SubSels = tags.DD_NeuroElf_CM_subsels;
    ch.SubjectsTxt = tags.TX_NeuroElf_CM_subjects;
    ch.Covs = tags.LB_NeuroElf_CM_covariates;
    ch.Contrasts = tags.DD_NeuroElf_CM_contrast;
    ch.Predictors = tags.LB_NeuroElf_CM_predlist;
    ch.PredWeights = tags.LB_NeuroElf_CM_pweights;
    ch.StatType = tags.DD_NeuroElf_CM_stattype;
    ch.RFXstats = tags.CB_NeuroElf_CM_RFXstats;
    ch.StoreFXSize = tags.CB_NeuroElf_CM_sfxsize;
    ch.StoreResVTC = tags.CB_NeuroElf_CM_sresvtc;
    ch.StoreSmoothEst = tags.CB_NeuroElf_CM_ssmooest;
    ch.RobWMaps = tags.CB_NeuroElf_CM_srobinfo;
    ch.SubWMaps = tags.CB_NeuroElf_CM_robwmaps;
    ch.RankTrans = tags.CB_NeuroElf_CM_rnktrans;
    ch.AddGlobalMean = tags.CB_NeuroElf_CM_addGmean;
    ch.AllRegressors = tags.CB_NeuroElf_CM_allregs;
    ch.StoreInVMP = tags.DD_NeuroElf_CM_storeVMP;
    ch.StoreInVMPTxt = tags.TX_NeuroElf_CM_storeVMP;
    ch.SmoothData = tags.CB_NeuroElf_CM_smdata;
    ch.SmoothDataKernel = tags.ED_NeuroElf_CM_smdatak;
    ch.BRange = tags.CB_NeuroElf_CM_brange;
    ch.BRangeMin = tags.ED_NeuroElf_CM_brange1;
    ch.BRangeMax = tags.ED_NeuroElf_CM_brange2;
    ch.InterpMethod = tags.DD_NeuroElf_CM_imeth;

    % put list of loaded GLMs into dropdown
    glmfiles = cell(numel(nelfo), 1);
    glmobjs = cell(size(glmfiles));
    for c = 1:numel(glmfiles)
        [gpath, glmfiles{c}] = fileparts(nelfo(c).FilenameOnDisk);
        glmobjs{c} = nelfo(c);
        if isempty(glmfiles{c})
            glmfiles{c} = sprintf('<#%d: %d subs, %d preds>', ...
                glmobjs{c}.Filenumber, glmobjs{c}.NrOfSubjects, ...
                numel(glmobjs{c}.SubjectPredictors));
        end
    end
    ne_gcfg.fcfg.CM.GLMs = glmobjs;
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

    % put currently selected VMP into struct (if any)
    if strcmpi(ne_gcfg.fcfg.StatsVar.Filetype, 'vmp')
        ne_gcfg.fcfg.CM.VMP = ne_gcfg.fcfg.StatsVar;
    else
        ne_gcfg.fcfg.CM.VMP = [];
    end

    % set callbacks -> GLM selection
    ch.GLMs.Callback = @ne_cm_setglm;

    % -> group/subject selection
    ch.UseGroups.Callback = @ne_cm_usegroups;
    tags.BT_NeuroElf_CM_addgrp.Callback = @ne_cm_addgroup;
    tags.BT_NeuroElf_CM_delgrp.Callback = @ne_cm_delgroup;
    tags.BT_NeuroElf_CM_rengrp.Callback = @ne_cm_rengroup;
    ch.Groups.Callback = @ne_cm_selgroup;
    ch.SubSels.Callback = @ne_cm_selectsubsel;
    tags.BT_NeuroElf_CM_addssel.Callback = @ne_cm_addsubsel;
    tags.BT_NeuroElf_CM_delssel.Callback = @ne_cm_delsubsel;
    ch.Subjects.Callback = @ne_cm_selectsubs;

    % -> covariates
    ch.Covs.Callback = @ne_cm_selectcovs;
    tags.BT_NeuroElf_CM_covload.Callback = @ne_cm_loadcovs;
    tags.BT_NeuroElf_CM_covroi.Callback = @ne_cm_roicovs;
    tags.BT_NeuroElf_CM_coveval.Callback = @ne_cm_coveval;
    tags.BT_NeuroElf_CM_covprops.Callback = @ne_cm_covprops;
    tags.BT_NeuroElf_CM_covadd.Callback = @ne_cm_addcovs;
    tags.BT_NeuroElf_CM_covdel.Callback = @ne_cm_delcovs;

    % -> contrasts
    ch.Contrasts.Callback = @ne_cm_selectcon;
    tags.BT_NeuroElf_CM_addcon.Callback = @ne_cm_addcon;
    tags.BT_NeuroElf_CM_delcon.Callback = @ne_cm_delcon;
    tags.BT_NeuroElf_CM_rencon.Callback = @ne_cm_rencon;
    ch.Predictors.Callback = @ne_cm_selectpreds;
    ch.PredWeights.Callback = @ne_cm_selectpweights;
    tags.BT_NeuroElf_CM_cwp1.Callback = {@ne_cm_setpredweight, '1'};
    tags.BT_NeuroElf_CM_cw0.Callback ={@ne_cm_setpredweight, '0'};
    tags.BT_NeuroElf_CM_cwn1.Callback = {@ne_cm_setpredweight, '-1'};
    tags.BT_NeuroElf_CM_cwval.Callback = {@ne_cm_setpredweight, ''};
    tags.BT_NeuroElf_CM_cwc.Callback = @ne_cm_clearpredweights;
    tags.BT_NeuroElf_CM_cwb.Callback = @ne_cm_balancepredweights;
    tags.BT_NeuroElf_CM_cselect.Callback = {@ne_cm_setpredweight, 'select'};
    tags.BT_NeuroElf_CM_cava.Callback = @ne_ancova_open;

    % -> type of stats
    %tags.RB_NeuroElf_CM_OLSonly.Callback = {@ne_cm_setstype, 'ols'};
    %tags.RB_NeuroElf_CM_Robust.Callback = {@ne_cm_setstype, 'rob'};
    %tags.RB_NeuroElf_CM_OLSRobust.Callback = {@ne_cm_setstype, 'olsrob'};
    tags.DD_NeuroElf_CM_stattype.Callback = @ne_cm_setstype;

    % -> switching RFX off
    tags.CB_NeuroElf_CM_RFXstats.Callback = @ne_cm_rfxswitch;

    % -> switching smoothing for data on/off
    tags.CB_NeuroElf_CM_smdata.Callback = @ne_cm_smoothdata;
    tags.ED_NeuroElf_CM_smdatak.Callback = @ne_cm_smoothdata;

    % -> load/save of contrasts
    tags.BT_NeuroElf_CM_loadcons.Callback = @ne_cm_loadcons;
    tags.BT_NeuroElf_CM_savecons.Callback = @ne_cm_savecons;
    tags.BT_NeuroElf_CM_copycons.Callback = @ne_cm_copycons;

    % -> do the work
    tags.BT_NeuroElf_CM_compute1.Callback = {@ne_cm_compute, 1};
    tags.BT_NeuroElf_CM_compute.Callback = @ne_cm_compute;

    % set h struct
    ne_gcfg.h.CM.h = ch;

    % get position from ini if any good
    try
        lastpos = ne_gcfg.c.ini.Children.CMPosition;
        if any(lastpos ~= -1)
            hFig.Position(1:2) = lastpos;
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end

    % last selected GLM is currently empty!
    hFig.UserData = struct('lastglm', []);

    % set visible, modal and wait
    hFig.CloseRequestFcn = 'set(gcbf, ''Visible'', ''off'');';
    hFig.HandleVisibility = 'callback';
    hFig.Visible = 'on';

    % set current GLM, then update (.ListboxTop fields otherwise unheeded!)
    ne_cm_setglm;
    drawnow;

% give warning
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Error loading contrast manager.', 'NeuroElf GUI - error', 'modal'));
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'cm_open')) = [];
end
