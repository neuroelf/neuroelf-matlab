% PUBLIC FUNCTION ne_mkda: invoke and handle MKDA dialog
function varargout = ne_mkda(varargin)

% Version:  v0.9d
% Build:    14062015
% Date:     Jun-20 2014, 3:32 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, Jochen Weber
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

% don't do anything while an analysis is being run!
if any(strcmp(ne_gcfg.c.blockcb, 'mkda_run'))
    return;
end

% already open?
if isfield(ne_gcfg.h, 'MKDA') && ...
    isstruct(ne_gcfg.h.MKDA) && ...
    isfield(ne_gcfg.h.MKDA, 'MKDAFig') && ...
    numel(ne_gcfg.h.MKDA.MKDAFig) == 1 && ...
    isxfigure(ne_gcfg.h.MKDA.MKDAFig, true)
    ne_gcfg.h.MKDA.MKDAFig.Visible = 'on';

    % bring up
    figure(ne_gcfg.h.MKDA.MKDAFig.MLHandle);

    % make sure display is up-to-date
    ne_mkda_setplp;

    % return;
    return;
end

% blocked?
if any(strcmp(ne_gcfg.c.blockcb, 'mkda_open'))
    return;
end
ne_gcfg.c.blockcb{end+1} = 'mkda_open';

% load contrast manager
try
    hFig = xfigure([neuroelf_path('tfg') '/ne_mkda.tfg']);
    ne_gcfg.h.MKDA.MKDAFig = hFig;

    % get required controls
    tags = hFig.TagStruct;

    % set callbacks
    tags.UIM_NeuroElf_MKDA_Load.Callback = @ne_mkda_load;
    tags.UIM_NeuroElf_MKDA_Import.Callback = @ne_mkda_import;
    tags.UIM_NeuroElf_MKDA_Save.Callback = @ne_mkda_save;
    tags.UIM_NeuroElf_MKDA_SaveAs.Callback = {@ne_mkda_save, 'saveas'};
    tags.UIM_NeuroElf_MKDA_Run.Callback = @ne_mkda_run;
    tags.UIM_NeuroElf_MKDA_Close.Callback = @ne_mkda_closeui;
    tags.UIM_NeuroElf_MKDA_OPT_mskr.Callback = {@ne_mkda_setoption, 'ApplyMask'};
    tags.UIM_NeuroElf_MKDA_OCNT_diff.Callback = {@ne_mkda_setoption, 'ContrastComp', 'diff'};
    tags.UIM_NeuroElf_MKDA_OCNT_excl.Callback = {@ne_mkda_setoption, 'ContrastComp', 'excl'};
    tags.UIM_NeuroElf_MKDA_OCNT_wex.Callback = {@ne_mkda_setoption, 'ContrastComp', 'wexcl'};
    tags.UIM_NeuroElf_MKDA_OCNT_swex.Callback = {@ne_mkda_setoption, 'ContrastCompExclWeight'};
    tags.UIM_NeuroElf_MKDA_OGRP_sum.Callback = {@ne_mkda_setoption, 'GroupMapComp', 'sum'};
    tags.UIM_NeuroElf_MKDA_OGRP_wsum.Callback = {@ne_mkda_setoption, 'GroupMapComp', 'wsum'};
    tags.UIM_NeuroElf_MKDA_OJBM_max.Callback = {@ne_mkda_setoption, 'JoinBlobComp', 'max'};
    tags.UIM_NeuroElf_MKDA_OJBM_rsum.Callback = {@ne_mkda_setoption, 'JoinBlobComp', 'rsum'};
    tags.UIM_NeuroElf_MKDA_OPT_imap.Callback = {@ne_mkda_setoption, 'KeepIndivMaps'};
    tags.UIM_NeuroElf_MKDA_OWPS_none.Callback = {@ne_mkda_setoption, 'PPSWeighting', 'none'};
    tags.UIM_NeuroElf_MKDA_OWPS_conf.Callback = {@ne_mkda_setoption, 'PPSWeighting', 'confidence'};
    tags.UIM_NeuroElf_MKDA_OWPS_log.Callback = {@ne_mkda_setoption, 'PPSWeighting', 'logpoints'};
    tags.UIM_NeuroElf_MKDA_OWPS_npts.Callback = {@ne_mkda_setoption, 'PPSWeighting', 'points'};
    tags.UIM_NeuroElf_MKDA_OWPS_sqrt.Callback = {@ne_mkda_setoption, 'PPSWeighting', 'sqrtpoints'};
    tags.UIM_NeuroElf_MKDA_OPT_svmp.Callback = {@ne_mkda_setoption, 'SummaryVMP'};
    tags.UIM_NeuroElf_MKDA_OSPN_full.Callback = {@ne_mkda_setoption, 'SpatialNull', 'full'};
    tags.UIM_NeuroElf_MKDA_OSPN_near.Callback = {@ne_mkda_setoption, 'SpatialNull', 'near'};
    tags.UIM_NeuroElf_MKDA_OPT_uniq.Callback = {@ne_mkda_setoption, 'UniqueUnitPoints'};
    tags.UIM_NeuroElf_MKDA_PSlice.Callback = @ne_mkda_plotonslice;
    tags.UIM_NeuroElf_MKDA_PSurf.Callback = @ne_mkda_plotonsurf;
    tags.UIM_NeuroElf_MKDA_PSetColor.Callback = @ne_mkda_plotsetcolor;
    tags.UIM_NeuroElf_MKDA_PSetLabel.Callback = {@ne_mkda_setoption, 'LabelColumn'};
    tags.UIM_NeuroElf_MKDA_PSetLColor.Callback = {@ne_mkda_setoption, 'LabelColor'};
    tags.UIM_NeuroElf_MKDA_PSetBColor.Callback = {@ne_mkda_setoption, 'BorderColor'};
    tags.DD_NeuroElf_MKDA_PLPs.Callback = @ne_mkda_setplp;
    tags.LB_NeuroElf_MKDA_studies.Callback = @ne_mkda_listpoints;
    tags.LB_NeuroElf_MKDA_columns.Callback = @ne_mkda_listpoints;
    tags.DD_NeuroElf_MKDA_anas.Callback = @ne_mkda_setana;
    tags.BT_NeuroElf_MKDA_addana.Callback = @ne_mkda_addana;
    tags.BT_NeuroElf_MKDA_delana.Callback = @ne_mkda_delana;
    tags.BT_NeuroElf_MKDA_renana.Callback = @ne_mkda_renana;
    tags.ED_NeuroElf_MKDA_cont.Callback = @ne_mkda_parsecont;
    tags.DD_NeuroElf_MKDA_contcol.Callback = @ne_mkda_setcontcol;
    tags.LB_NeuroElf_MKDA_conds.Callback = {@ne_mkda_conds, 'sel'};
    tags.BT_NeuroElf_MKDA_addcnd.Callback = {@ne_mkda_conds, 'add'};
    tags.BT_NeuroElf_MKDA_setcnd.Callback = {@ne_mkda_conds, 'set'};
    tags.BT_NeuroElf_MKDA_delcnd.Callback = {@ne_mkda_conds, 'del'};
    tags.BT_NeuroElf_MKDA_addpar.Callback = {@ne_mkda_conds, 'addpar'};
    tags.BT_NeuroElf_MKDA_delpar.Callback = {@ne_mkda_conds, 'delpar'};
    tags.DD_NeuroElf_MKDA_cndcol.Callback = @ne_mkda_setcondcol;
    tags.DD_NeuroElf_MKDA_cndop2s.Callback = @ne_mkda_setcondval;
    tags.ED_NeuroElf_MKDA_weights.Callback = @ne_mkda_listpoints;
    tags.DD_NeuroElf_MKDA_stdcol.Callback = @ne_mkda_listpoints;
    tags.CB_NeuroElf_MKDA_selpts.Callback = @ne_mkda_listpoints;
    tags.BT_NeuroElf_MKDA_addcol.Callback = @ne_mkda_addplpcol;
    tags.BT_NeuroElf_MKDA_delcol.Callback = @ne_mkda_delplpcol;
    tags.BT_NeuroElf_MKDA_editpt.Callback = @ne_mkda_editplppt;
    tags.BT_NeuroElf_MKDA_setcol.Callback = @ne_mkda_setplpcol;
    tags.BT_NeuroElf_MKDA_delpts.Callback = @ne_mkda_delplppts;
    tags.ED_NeuroElf_MKDA_iter.Callback = @ne_mkda_updana;
    tags.DD_NeuroElf_MKDA_mask.Callback = {@ne_mkda_readui, 'mask'};
    tags.DD_NeuroElf_MKDA_res.Callback = @ne_mkda_updana;
    tags.RB_NeuroElf_MKDA_scindic.Callback = {@ne_mkda_readui, 'scindic'};
    tags.RB_NeuroElf_MKDA_scgauss.Callback = {@ne_mkda_readui, 'scgauss'};
    tags.ED_NeuroElf_MKDA_smkern.Callback = @ne_mkda_updana;
    tags.RB_NeuroElf_MKDA_spatnul.Callback = {@ne_mkda_readui, 'nulspat'};
    tags.RB_NeuroElf_MKDA_unitnul.Callback = {@ne_mkda_readui, 'nulunit'};
    tags.ED_NeuroElf_MKDA_smkdist.Callback = @ne_mkda_updana;
    tags.BT_NeuroElf_MKDA_closeui.Callback = @ne_mkda_closeui;
    tags.BT_NeuroElf_MKDA_runana.Callback = @ne_mkda_run;
    tags.BT_NeuroElf_MKDA_runanas.Callback = {@ne_mkda_run, 'runall'};

    % get some shorthands
    ch = struct;
    ch.PLPs = tags.DD_NeuroElf_MKDA_PLPs;
    ch.Studies = tags.LB_NeuroElf_MKDA_studies;
    ch.Columns = tags.LB_NeuroElf_MKDA_columns;
    ch.Analyses = tags.DD_NeuroElf_MKDA_anas;
    ch.Contrast = tags.ED_NeuroElf_MKDA_cont;
    ch.ContColumn = tags.DD_NeuroElf_MKDA_contcol;
    ch.CndParts = tags.LB_NeuroElf_MKDA_conds;
    ch.CndAnd = tags.RB_NeuroElf_MKDA_cndand;
    ch.CndOr = tags.RB_NeuroElf_MKDA_cndor;
    ch.CndColumn = tags.DD_NeuroElf_MKDA_cndcol;
    ch.CndOperator = tags.DD_NeuroElf_MKDA_cndoper;
    ch.CndStaticOp = tags.DD_NeuroElf_MKDA_cndop2s;
    ch.CndFlexiOp = tags.ED_NeuroElf_MKDA_cndop2m;
    ch.Weights = tags.ED_NeuroElf_MKDA_weights;
    ch.StudyColumn = tags.DD_NeuroElf_MKDA_stdcol;
    ch.PointsLabel = tags.TX_NeuroElf_MKDA_points;
    ch.SelectedPts = tags.CB_NeuroElf_MKDA_selpts;
    ch.Points = tags.LB_NeuroElf_MKDA_points;
    ch.Iterations = tags.ED_NeuroElf_MKDA_iter;
    ch.Mask = tags.DD_NeuroElf_MKDA_mask;
    ch.Resolution = tags.DD_NeuroElf_MKDA_res;
    ch.SphereIndicator = tags.RB_NeuroElf_MKDA_scindic;
    ch.SphereGaussian = tags.RB_NeuroElf_MKDA_scgauss;
    ch.SphereSize = tags.ED_NeuroElf_MKDA_smkern;
    ch.SphereTaper = tags.ED_NeuroElf_MKDA_smkdist;
    ch.NullSpatial = tags.RB_NeuroElf_MKDA_spatnul;
    ch.NullStatUnit = tags.RB_NeuroElf_MKDA_unitnul;
    ch.LabelColumnMenu = tags.UIM_NeuroElf_MKDA_PSetLabel;

    % settings
    ch.PLPs.UserData = cell(0, 3);
    ch.CndParts.FontName = 'Courier';
    ch.Points.FontName = 'Courier';

    % set h struct
    ne_gcfg.h.MKDA.h = ch;

    % get position from ini if any good
    try
        lastpos = ne_gcfg.c.ini.Children.MKDAPosition;
        if any(lastpos ~= -1)
            hFig.Position(1:2) = lastpos;
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end

    % initialize configuration
    ne_gcfg.fcfg.MKDA = ne_gcfg.c.ini.MKDA;

    % set visible, modal and wait
    hFig.CloseRequestFcn = @ne_mkda_closeui;
    ne_mkda_setplp;
    hFig.HandleVisibility = 'callback';
    hFig.Visible = 'on';
    drawnow;

% give warning
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Error loading MKDA dialog.', 'NeuroElf GUI - error', 'modal'));
end

% unblock call
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_open')) = [];
