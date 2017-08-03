% FUNCTION ne_mdm_updatemdm: update MDM object from UI
function ne_mdm_updatemdm(varargin)

% Version:  v0.9b
% Build:    12042618
% Date:     Apr-10 2011, 4:52 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% update fields in RunTimeVars
ch = ne_gcfg.h.MDM.h;
mdm = ne_gcfg.fcfg.MDM.mdm;
mdm.TypeOfFunctionalData = 'VTC';
if ch.ModelRFX.Value > 0
    mdm.RFX_GLM = 1;
    mdm.SeparatePredictors = 2;
else
    mdm.RFX_GLM = 0;
    if ch.ModelFFX.Value > 0
        mdm.SeparatePredictors = 0;
    elseif ch.ModelSPSB.Value > 0
        mdm.SeparatePredictors = 2;
    elseif ch.ModelSPST.Value > 0
        mdm.SeparatePredictors = 1;
    else
        mdm.RFX_GLM = 1;
        mdm.SeparatePredictors = 2;
    end
end
mdm.PSCTransformation = double(ch.TransPSC.Value > 0);
mdm.zTransformation = double(ch.Transz.Value > 0);
try
    mdm.XTC_RTC = [ch.FuncFiles.String, ch.DsgnFiles.String];
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    mdm.XTC_RTC = cell(0, 2);
    ch.FuncFiles.Enable = 'off';
    ch.FuncFiles.String = {'<no files selected>'};
    ch.FuncFiles.Value = [];
    ch.DsgnFiles.Enable = 'off';
    ch.DsgnFiles.String = {'<no files selected>'};
    ch.DsgnFiles.Value = [];
    uiwait(warndlg('Internal error with handling filenames!', ...
        'NeuroElf - warning', 'modal'));
end
mdm.NrOfStudies = size(mdm.XTC_RTC, 1);
if isempty(ch.MParFiles.String) || ...
   (numel(ch.MParFiles.String) == 1 && ...
    strcmpi(ch.MParFiles.String, '<no files selected>'))
    mdm.RunTimeVars.MotionParameters = {};
elseif numel(ch.MParFiles.String) == mdm.NrOfStudies
    mdm.RunTimeVars.MotionParameters = ch.MParFiles.String;
else
    mdm.RunTimeVars.MotionParameters = {};
    ch.MParFiles.Enable = 'off';
    ch.MParFiles.String = {'<no files selected>'};
    ch.MParFiles.Value = [];
    uiwait(warndlg('Internal error with handling filenames!', ...
        'NeuroElf - warning', 'modal'));
end
mdm.RunTimeVars.BaseFolder = ch.Basefolder.String;
mdm.RunTimeVars.CombineFFX = (ch.CombineFFX.Value > 0);
mdm.RunTimeVars.Deconvolution = ...
    double(ch.Deconv.Value > 0 && strcmpi(ch.Deconv.Enable, 'on'));
try
    dcls = str2double(ch.DeconvLags.String);
    if numel(dcls) ~= 1 || ...
        isinf(dcls) || ...
        isnan(dcls) || ...
        dcls < 1 || ...
        dcls ~= fix(dcls)
        error('BAD_DECONVLAGS');
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    dcls = 12;
    ch.DeconvLags.String = '12';
end
mdm.RunTimeVars.DeconvolutionLags = dcls;
if ch.Derivs12.Value > 0
    mdm.RunTimeVars.Derivatives = 2;
elseif ch.Derivs1.Value > 0
    mdm.RunTimeVars.Derivatives = 1;
else
    mdm.RunTimeVars.Derivatives = 0;
end
switch (ch.DerivBoost.Value)
    case {1}
        mdm.RunTimeVars.DerivativeBoost = 'auc';
    case {2}
        mdm.RunTimeVars.DerivativeBoost = 'boost';
    case {3}
        mdm.RunTimeVars.DerivativeBoost = 'max';
    case {4}
        mdm.RunTimeVars.DerivativeBoost = 'posauc';
    otherwise
        mdm.RunTimeVars.DerivativeBoost = 'none';
end
mdm.RunTimeVars.DerivativeSkip = ch.NoDerivConds.String;
if isempty(mdm.RunTimeVars.DerivativeSkip)
    mdm.RunTimeVars.DerivativeSkip = {};
elseif ~iscell(mdm.RunTimeVars.DerivativeSkip)
    mdm.RunTimeVars.DerivativeSkip = cellstr(mdm.RunTimeVars.DerivativeSkip);
end
if numel(mdm.RunTimeVars.DerivativeSkip) == 1 && ...
   (isempty(mdm.RunTimeVars.DerivativeSkip{1}) || ...
    (mdm.RunTimeVars.DerivativeSkip{1}(1) == '<' && ...
     mdm.RunTimeVars.DerivativeSkip{1}(end) == '>'))
    mdm.RunTimeVars.DerivativeSkip = {};
end
mdm.RunTimeVars.DesignFilesPattern = ch.DsgnPattern.String;
mdm.RunTimeVars.GlobalSignals = ch.GlobalSignals.Value - 1;
if mdm.RunTimeVars.GlobalSignals > 2
    mdm.RunTimeVars.GlobalSignals = ch.GlobalSignals.String{ch.GlobalSignals.Value};
end
try
    ithr = str2double(ch.IThresh.String);
    if numel(ithr) ~= 1 || ...
        isinf(ithr) || ...
        isnan(ithr) || ...
        ithr < 0
        error('BAD_ITHRESH');
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    ithr = 100;
    ch.IThresh.String = '100';
end
mdm.RunTimeVars.IntensityThreshold = ithr;
mskf = ch.MaskFile.String{ch.MaskFile.Value};
if strcmpi(mskf, 'click to select...')
    mskf = '';
end
mdm.RunTimeVars.MaskFilename = mskf;
mdm.RunTimeVars.MotionParamsDiff = (ch.MParDiff.Value > 0);
mdm.RunTimeVars.MotionParamsSquared = (ch.MParSquared.Value > 0);
mdm.RunTimeVars.MotionParamsUse = (ch.UseMotParms.Value > 0);
mdm.RunTimeVars.MotionParamFilesPattern = ch.MParPattern.String;
mdm.RunTimeVars.GLMOutputFilename = ch.GLMFile.String;
mdm.RunTimeVars.RedoSubjects = (ch.RedoSubjects.Value > 0);
mdm.RunTimeVars.RestConds = ch.RestConds.String;
if isempty(mdm.RunTimeVars.RestConds)
    mdm.RunTimeVars.RestConds = {};
elseif ~iscell(mdm.RunTimeVars.RestConds)
    mdm.RunTimeVars.RestConds = cellstr(mdm.RunTimeVars.RestConds);
end
if numel(mdm.RunTimeVars.RestConds) == 1 && ...
   (isempty(mdm.RunTimeVars.RestConds{1}) || ...
    (mdm.RunTimeVars.RestConds{1}(1) == '<' && ...
     mdm.RunTimeVars.RestConds{1}(end) == '>'))
    mdm.RunTimeVars.RestConds = {};
end
mdm.RunTimeVars.RobustRegression = (ch.RegrRobust.Value > 0);
mdm.RunTimeVars.ShowDesigns = (ch.ShowDesigns.Value > 0);
mdm.RunTimeVars.SingleTrial = (ch.SingleTrial.Value > 0);
mdm.RunTimeVars.SingleTrialSkip = ch.STSkipConds.String;
if isempty(mdm.RunTimeVars.SingleTrialSkip)
    mdm.RunTimeVars.SingleTrialSkip = {};
elseif ~iscell(mdm.RunTimeVars.SingleTrialSkip)
    mdm.RunTimeVars.SingleTrialSkip = cellstr(mdm.RunTimeVars.SingleTrialSkip);
end
if numel(mdm.RunTimeVars.SingleTrialSkip) == 1 && ...
   (isempty(mdm.RunTimeVars.SingleTrialSkip{1}) || ...
    (mdm.RunTimeVars.SingleTrialSkip{1}(1) == '<' && ...
     mdm.RunTimeVars.SingleTrialSkip{1}(end) == '>'))
    mdm.RunTimeVars.SingleTrialSkip = {};
end
try
    tfco = str2double(ch.TFiltCutOff.String);
    if numel(tfco) ~= 1 || ...
        isnan(tfco) || ...
        tfco <= 15
        tfco = Inf;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    tfco = Inf;
end
mdm.RunTimeVars.TempFilterCutoff = tfco;
if ch.TFiltNone.Value > 0
    ftype = 'none';
elseif ch.TFiltDCT.Value > 0
    ftype = 'dct';
else
    ftype = 'fourier';
end
mdm.RunTimeVars.TempFilterType = ftype;
mdm.RunTimeVars.XTCFilesPattern = ch.FuncPattern.String;
