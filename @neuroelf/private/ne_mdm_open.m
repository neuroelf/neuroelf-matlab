% PUBLIC FUNCTION ne_mdm_open: open the MDM dialog
function varargout = ne_mdm_open(varargin)

% Version:  v1.1
% Build:    16052509
% Date:     May-25 2016, 9:59 AM EST
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% already open?
if isfield(ne_gcfg.h, 'MDM') && isstruct(ne_gcfg.h.MDM) && isfield(ne_gcfg.h.MDM, 'MDMFig') && ...
    numel(ne_gcfg.h.MDM.MDMFig) == 1 && isxfigure(ne_gcfg.h.MDM.MDMFig, true)
    ne_gcfg.h.MDM.MDMFig.Visible = 'on';
    figure(ne_gcfg.h.MDM.MDMFig.MLHandle);
    return;
end

% blocked?
if any(strcmp(ne_gcfg.c.blockcb, 'mdm_open'))
    return;
end

% load contrast manager
try
    hFig = xfigure([neuroelf_path('tfg') '/ne_mdmconfig.tfg']);
    ne_gcfg.h.MDM.MDMFig = hFig;

    % block further callbacks
    ne_gcfg.c.blockcb{end+1} = 'mdm_open';

    % get required controls
    tags = hFig.TagStruct;
    ch = struct;
    ch.Basefolder = tags.ED_NeuroElf_mdm_basefld;
    ch.FuncPattern = tags.ED_NeuroElf_mdm_xtcpatt;
    ch.DsgnPattern = tags.ED_NeuroElf_mdm_despatt;
    ch.UseMotParms = tags.CB_NeuroElf_mdm_motpatt;
    ch.MParPattern = tags.ED_NeuroElf_mdm_motpatt;
    ch.MParDiff = tags.CB_NeuroElf_mdm_mparmsd;
    ch.MParSquared = tags.CB_NeuroElf_mdm_mparmsq;
    ch.FuncFiles = tags.LB_NeuroElf_mdm_func;
    ch.DsgnFiles = tags.LB_NeuroElf_mdm_design;
    ch.MParFiles = tags.LB_NeuroElf_mdm_mparam;
    ch.Trans0 = tags.RB_NeuroElf_mdm_trnull;
    ch.TransPSC = tags.RB_NeuroElf_mdm_trpsc;
    ch.Transz = tags.RB_NeuroElf_mdm_trz;
    ch.IThresh = tags.ED_NeuroElf_mdm_ithresh;
    ch.MaskFile = tags.DD_NeuroElf_mdm_mask;
    ch.RegrOLS = tags.RB_NeuroElf_mdm_rols;
    ch.RegrRobust = tags.RB_NeuroElf_mdm_rrob;
    ch.Deconv = tags.CB_NeuroElf_mdm_prtdcnv;
    ch.DeconvLags = tags.ED_NeuroElf_mdm_prtlags;
    ch.TFiltNone = tags.RB_NeuroElf_mdm_fltnull;
    ch.TFiltDCT = tags.RB_NeuroElf_mdm_fltdct;
    ch.TFiltFourier = tags.RB_NeuroElf_mdm_fltfour;
    ch.TFiltCutOff = tags.ED_NeuroElf_mdm_fltsecs;
    ch.GlobalSignals = tags.DD_NeuroElf_mdm_gsig;
    ch.ModelFFX = tags.RB_NeuroElf_mdm_mdlffx;
    ch.ModelSPSB = tags.RB_NeuroElf_mdm_mdlspsb;
    ch.ModelSPST = tags.RB_NeuroElf_mdm_mdlspst;
    ch.ModelRFX = tags.RB_NeuroElf_mdm_mdlrfx;
    ch.SingleTrial = tags.CB_NeuroElf_mdm_sngtrl;
    ch.CombineFFX = tags.CB_NeuroElf_mdm_cmbffx;
    ch.VarWeight = tags.CB_NeuroElf_mdm_vweight;
    ch.GLMFile = tags.ED_NeuroElf_mdm_glmfile;
    ch.RedoSubjects = tags.CB_NeuroElf_mdm_redosub;
    ch.ShowDesigns = tags.CB_NeuroElf_mdm_showdsg;
    ch.Derivs0 = tags.RB_NeuroElf_mdm_deriv0;
    ch.Derivs1 = tags.RB_NeuroElf_mdm_deriv1;
    ch.Derivs12 = tags.RB_NeuroElf_mdm_deriv12;
    ch.DerivBoost = tags.DD_NeuroElf_mdm_dboost;
    ch.PRTConds = tags.LB_NeuroElf_mdm_pcond;
    ch.RestConds = tags.LB_NeuroElf_mdm_dropcnd;
    ch.NoDerivConds = tags.LB_NeuroElf_mdm_nodvcnd;
    ch.STSkipConds = tags.LB_NeuroElf_mdm_sstcnd;
    ch.VTCAvgConds = tags.LB_NeuroElf_mdm_acond;
    ch.VTCAvgCondCols = tags.BT_NeuroElf_acond_ccols;
    ch.VTCAvgBaseWin = tags.ED_NeuroElf_avgvtcbwin;
    ch.VTCAvgSampleTR = tags.ED_NeuroElf_avgvtcstr;
    ch.VTCAvgWin = tags.ED_NeuroElf_avgvtcswin;
    ch.VTCAvgNaive = tags.CB_NeuroElf_avgvtcnaive;
    ch.VTCAvgRemGS = tags.CB_NeuroElf_avgvtcremgs;
    ch.VTCAvgSngTrl = tags.CB_NeuroElf_avgvtcrsngt;
    ch.VTCAvgFFX = tags.CB_NeuroElf_mdm_ffxavg;
    ch.VTCAvgRFX = tags.CB_NeuroElf_mdm_rfxavg;
    ch.VTCAvgWRFX = tags.CB_NeuroElf_mdm_wrfxavg;
    ch.VTCAvgRobust = tags.CB_NeuroElf_mdm_robavg;
    ch.VTCAvgRun = tags.BT_NeuroElf_mdm_avgvtcs;

    % pre-set page 2 controls
    ch.PRTConds.String = {'<detecting conditions...>'};
    ch.PRTConds.Value = [];
    ch.PRTConds.UserData = {};
    ch.RestConds.String = {'<no conditions dropped>'};
    ch.RestConds.Value = [];
    ch.NoDerivConds.String = {'<all with derivatives>'};
    ch.NoDerivConds.Value = [];
    ch.STSkipConds.String = {'<all split into trials>'};
    ch.STSkipConds.Value = [];
    ch.VTCAvgConds.String = {'<no conditions selected>'};
    ch.VTCAvgConds.Value = [];
    ch.VTCAvgConds.UserData = cell(0, 2);
    
    % set callbacks as required
    tags.ED_NeuroElf_mdm_basefld.Callback = @ne_mdm_bfedit;
    tags.BT_NeuroElf_mdm_browse.Callback = @ne_mdm_browse;
    tags.CB_NeuroElf_mdm_motpatt.Callback = @ne_mdm_togglerps;
    tags.BT_NeuroElf_mdm_ffiles.Callback = @ne_mdm_ffiles;
    tags.BT_NeuroElf_mdm_condopt.Callback = @ne_mdm_showopt;
    tags.LB_NeuroElf_mdm_func.Callback = {@ne_mdm_select, 'func'};
    tags.LB_NeuroElf_mdm_design.Callback = {@ne_mdm_select, 'design'};
    tags.LB_NeuroElf_mdm_mparam.Callback = {@ne_mdm_select, 'mparam'};
    tags.BT_NeuroElf_mdm_dsc.Callback = @ne_mdm_vtcqa;
    tags.BT_NeuroElf_mdm_pmp.Callback = @ne_mdm_mpplot;
    tags.BT_NeuroElf_mdm_add.Callback = @ne_mdm_addfiles;
    tags.BT_NeuroElf_mdm_del.Callback = @ne_mdm_delfiles;
    tags.RB_NeuroElf_mdm_trnull.Callback = {@ne_mdm_transsetone, 1};
    tags.RB_NeuroElf_mdm_trpsc.Callback = {@ne_mdm_transsetone, 2};
    tags.RB_NeuroElf_mdm_trz.Callback = {@ne_mdm_transsetone, 3};
    tags.DD_NeuroElf_mdm_mask.Callback = @ne_mdm_selectmask;
    tags.RB_NeuroElf_mdm_rols.Callback = {@ne_mdm_regsetone, 1};
    tags.RB_NeuroElf_mdm_rrob.Callback = {@ne_mdm_regsetone, 2};
    tags.RB_NeuroElf_mdm_fltnull.Callback = {@ne_mdm_filtsetone, 1};
    tags.RB_NeuroElf_mdm_fltdct.Callback = {@ne_mdm_filtsetone, 2};
    tags.RB_NeuroElf_mdm_fltfour.Callback = {@ne_mdm_filtsetone, 3};
    tags.RB_NeuroElf_mdm_mdlffx.Callback = {@ne_mdm_modelsetone, 1};
    tags.RB_NeuroElf_mdm_mdlspsb.Callback = {@ne_mdm_modelsetone, 2};
    tags.RB_NeuroElf_mdm_mdlspst.Callback = {@ne_mdm_modelsetone, 3};
    tags.RB_NeuroElf_mdm_mdlrfx.Callback = {@ne_mdm_modelsetone, 4};
    tags.BT_NeuroElf_mdm_glmset.Callback = @ne_mdm_setglmfile;
    tags.BT_NeuroElf_mdm_mainpag.Callback = {@ne_mdm_showopt, false};
    tags.BT_NeuroElf_mdm_load.Callback = @ne_mdm_loadmdm;
    tags.BT_NeuroElf_mdm_save.Callback = @ne_mdm_savemdm;
    tags.BT_NeuroElf_mdm_cancel.Callback = @ne_mdm_closeui;
    tags.BT_NeuroElf_mdm_compvoi.Callback = @ne_mdm_computevoi;
    tags.BT_NeuroElf_mdm_compute.Callback = @ne_mdm_compute;
    tags.BT_NeuroElf_acond_up.Callback = {@ne_mdm_acond, 'moveup'};
    tags.BT_NeuroElf_acond_down.Callback = {@ne_mdm_acond, 'movedown'};
    tags.BT_NeuroElf_acond_addcl.Callback = {@ne_mdm_acond, 'addcoll'};
    tags.BT_NeuroElf_acond_add.Callback = {@ne_mdm_acond, 'addconds'};
    tags.BT_NeuroElf_acond_del.Callback = {@ne_mdm_acond, 'delconds'};
    tags.CB_NeuroElf_avgvtcnaive.Callback = {@ne_mdm_acond, 'naive'};
    tags.BT_NeuroElf_mdm_avgvtcs.Callback = {@ne_mdm_acond, 'average'};

    % create new config
    ne_gcfg.fcfg.MDM = struct('condcols', zeros(0, 3), 'mdm', xff('new:mdm'));
    mdm = ne_gcfg.fcfg.MDM.mdm;
    mdm.TypeOfFunctionalData = 'VTC';
    mdm.RFX_GLM = 1;
    mdm.PSCTransformation = 1;
    mdm.zTransformation = 0;
    mdm.SeparatePredictors = 2;
    mdm.NrOfStudies = 0;
    mdm.XTC_RTC = cell(0, 2);

    % set h struct
    ne_gcfg.h.MDM.h = ch;

    % get position from ini if any good
    try
        lastpos = ne_gcfg.c.ini.Children.MDMPosition;
        if any(lastpos ~= -1)
            hFig.Position(1:2) = lastpos;
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end

    % update UI from object
    ne_mdm_updateui;

    % set visible, modal and wait
    hFig.CloseRequestFcn = @ne_mdm_closeui;
    hFig.ShowPage(1);
    hFig.HandleVisibility = 'callback';
    hFig.Visible = 'on';
    drawnow;

% give warning
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Error loading MDM config dialog.', 'NeuroElf GUI - error', 'modal'));
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mdm_open')) = [];
end
