% FUNCTION ne_ancova_open: loads the ANCOVA config UI
function varargout = ne_ancova_open(varargin)

% Version:  v1.1
% Build:    16012618
% Date:     Jan-26 2016, 6:41 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, Jochen Weber
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
if isfield(ne_gcfg.h, 'AC') && ...
    isstruct(ne_gcfg.h.AC) && ...
    isfield(ne_gcfg.h.AC, 'ACFig') && ...
    numel(ne_gcfg.h.AC.ACFig) == 1 && ...
    isxfigure(ne_gcfg.h.AC.ACFig, true)
    ne_gcfg.h.AC.ACFig.Visible = 'on';
    figure(ne_gcfg.h.AC.ACFig.MLHandle);
    return;
end

% blocked?
if any(strcmp(ne_gcfg.c.blockcb, 'ancova_open'))
    return;
end

% echo
if ne_gcfg.c.echo
    ne_echo('neuroelf_gui(''ancova_open'');');
end

% get list of loaded GLMs
nelfo = xff(0, 'objects');
isglm = false(size(nelfo));
for oc = 1:numel(nelfo)
    if strcmpi(nelfo(oc).Filetype, 'glm') && ...
        nelfo(oc).NrOfSubjects > 1 && ...
        nelfo(oc).SeparatePredictors == 2
        isglm(oc) = true;
    end
end
nelfo(~isglm) = [];

% return if no GLMs found
if isempty(nelfo)
    uiwait(msgbox('You must load a multi-subjects GLM file first to use the ANCOVA UI.', ...
        'NeuroElf GUI - info', 'modal'));
    return;
end

% load contrast manager
try
    hFig = 0;
    hFig = xfigure([neuroelf_path('tfg') '/ne_glmancova.tfg']);
    ne_gcfg.h.AC.ACFig = hFig;

    % block further callbacks
    ne_gcfg.c.blockcb{end+1} = 'ancova_open';

    % get required controls
    tags = hFig.TagStruct;
    ch = struct;
    ch.GLMs = tags.DD_NeuroElf_AC_GLMs;
    ch.UseGroups = tags.CB_NeuroElf_AC_groups;
    ch.Groups = tags.LB_NeuroElf_AC_groups;
    ch.Subjects = tags.LB_NeuroElf_AC_subjects;
    ch.Covs = tags.LB_NeuroElf_AC_covariates;
    ch.CovInteractions = tags.CB_NeuroElf_AC_coviact;
    ch.Layout1 = tags.DD_NeuroElf_AC_layoutn;
    ch.LayoutX = tags.DD_NeuroElf_AC_layout;
    ch.CellType = tags.DD_NeuroElf_AC_celltype;
    ch.CellCont = cell(10, 3, 2);
    for c1 = 1:10
        for c2 = 1:3
            for c3 = 1:2
                ch.CellCont{c1, c2, c3} = ...
                    tags.(sprintf('DD_NeuroElf_AC_L%02dL%dL%d', c1, c2, c3));
            end
        end
    end
    ch.FLabelX11 = tags.TX_NeuroElf_AC_FLx11;
    ch.FLabelX21 = tags.TX_NeuroElf_AC_FLx21;
    ch.FLabelX31 = tags.TX_NeuroElf_AC_FLx31;
    ch.FLabelX12 = tags.TX_NeuroElf_AC_FLx12;
    ch.FLabelX22 = tags.TX_NeuroElf_AC_FLx22;
    ch.FLabelX32 = tags.TX_NeuroElf_AC_FLx32;
    ch.AddGlobalMean = tags.CB_NeuroElf_AC_addGmean;
    ch.StoreInVMP = tags.DD_NeuroElf_AC_storeVMP;
    ch.MRCovNone = tags.RB_NeuroElf_AC_mrcnone;
    ch.MRCovGlobal = tags.RB_NeuroElf_AC_mrcglob;
    ch.MRCovGroup = tags.RB_NeuroElf_AC_mrcgroup;
    ch.SmoothData = tags.CB_NeuroElf_AC_smdata;
    ch.SmoothDataKernel = tags.ED_NeuroElf_AC_smdatak;
    ch.BRange = tags.CB_NeuroElf_AC_brange;
    ch.BRangeMin = tags.ED_NeuroElf_AC_brange1;
    ch.BRangeMax = tags.ED_NeuroElf_AC_brange2;
    ch.InterpMethod = tags.DD_NeuroElf_AC_imeth;
    ch.Compute = tags.BT_NeuroElf_AC_compute;
    ch.HDivider = tags.FR_NeuroElf_AC_hdivide;

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
    ne_gcfg.fcfg.AC.GLMs = glmobjs;
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
        ne_gcfg.fcfg.AC.VMP = ne_gcfg.fcfg.StatsVar;
    else
        ne_gcfg.fcfg.AC.VMP = [];
    end

    % set callbacks -> GLM selection
    ch.GLMs.Callback = {@ne_ancova_ui, 'setglm'};

    % -> group/subject selection
    ch.UseGroups.Callback = {@ne_ancova_ui, 'usegroups'};
    ch.Groups.Callback = {@ne_ancova_ui, 'usegroups', true};

    % -> layout and cell type
    ch.Layout1.Callback = {@ne_ancova_ui, 'layout'};
    ch.LayoutX.Callback = {@ne_ancova_ui, 'layout'};
    ch.CellType.Callback = {@ne_ancova_ui, 'celltype'};

    % -> switching smoothing for data on/off
    tags.CB_NeuroElf_AC_smdata.Callback = {@ne_ancova_ui, 'smoothdata'};
    tags.ED_NeuroElf_AC_smdatak.Callback = {@ne_ancova_ui, 'smoothdata'};

    % -> load/save of contrasts
    tags.BT_NeuroElf_AC_loadcfg.Callback = {@ne_ancova_ui, 'loadcfg'};
    tags.BT_NeuroElf_AC_savecfg.Callback = {@ne_ancova_ui, 'savecfg'};

    % -> do the work
    ch.Compute.Callback = @ne_ancova_compute;

    % set h struct
    ne_gcfg.h.AC.h = ch;

    % get position from ini if any good
    try
        lastpos = ne_gcfg.c.ini.Children.ACPosition;
        if any(lastpos ~= -1)
            hFig.Position(1:2) = lastpos;
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end

    % last selected GLM is currently empty!
    hFig.UserData = struct('lastglm', []);

    % set visible, modal and wait
    hFig.CloseRequestFcn = {@ne_ancova_ui, 'closeui'};
    hFig.HandleVisibility = 'callback';
    hFig.Visible = 'on';

    % set current GLM, then update (.ListboxTop fields otherwise unheeded!)
    ne_ancova_ui(0, 0, 'setglm');
    drawnow;

% give warning
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    ne_gcfg.fcfg.AC = [];
    ne_gcfg.h.AC = [];
    if isxfigure(hFig, true)
        hFig.Delete;
    end
    uiwait(warndlg('Error loading ANVOCA config UI.', 'NeuroElf GUI - error', 'modal'));
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'ancova_open')) = [];
end
