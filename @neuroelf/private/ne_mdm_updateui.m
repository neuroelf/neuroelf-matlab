% FUNCTION ne_mdm_updateui: update UI from MDM object
function ne_mdm_updateui(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:36 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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
cf = ne_gcfg.h.MDM.MDMFig;
ch = ne_gcfg.h.MDM.h;
mdm = ne_gcfg.fcfg.MDM.mdm;
if ~strcmpi(mdm.TypeOfFunctionalData, 'vtc')
    uiwait(warndlg('Currently only VTC-based MDMs are supported.', ...
        'NeuroElf - user information', 'modal'));
    ne_gcfg.fcfg.MDM.mdm.ClearObject;
    ne_gcfg.fcfg.MDM.mdm = xff('new:mdm');
    mdm = ne_gcfg.fcfg.MDM.mdm;
    mdm.TypeOfFunctionalData = 'VTC';
    mdm.RFX_GLM = 1;
    mdm.PSCTransformation = 1;
    mdm.zTransformation = 0;
    mdm.SeparatePredictors = 2;
    mdm.NrOfStudies = 0;
    mdm.XTC_RTC = cell(0, 2);
end
rtv = mdm.RunTimeVars;
if ~isfield(rtv, 'BaseFolder')
    rtv.BaseFolder = '';
end
if ~isfield(rtv, 'CombineFFX')
    rtv.CombineFFX = false;
end
if ~isfield(rtv, 'Deconvolution')
    rtv.Deconvolution = false;
end
if ~isfield(rtv, 'DeconvolutionLags')
    rtv.DeconvolutionLags = 12;
end
if ~isfield(rtv, 'Derivatives')
    rtv.Derivatives = 0;
end
if ~isfield(rtv, 'DerivativeBoost')
    rtv.DerivativeBoost = 'boost';
end
if ~isfield(rtv, 'DerivativeSkip')
    rtv.DerivativeSkip = {};
end
if ~isfield(rtv, 'DesignFilesPattern')
    rtv.DesignFilesPattern = '  subj*/*/fun*/*.prt';
end
if ~isfield(rtv, 'GlobalSignals')
    rtv.GlobalSignals = 0;
end
if ~isfield(rtv, 'IntensityThreshold')
    rtv.IntensityThreshold = 100;
end
if ~isfield(rtv, 'MaskFilename')
    rtv.MaskFilename = '';
end
if ~isfield(rtv, 'MotionParameters')
    rtv.MotionParameters = {};
end
if ~isfield(rtv, 'MotionParamsDiff')
    rtv.MotionParamsDiff = false;
end
if ~isfield(rtv, 'MotionParamsSquared')
    rtv.MotionParamsSquared = false;
end
if ~isfield(rtv, 'MotionParamsUse')
    rtv.MotionParamsUse = false;
end
if ~isfield(rtv, 'MotionParamFilesPattern')
    rtv.MotionParamFilesPattern = '  subj*/*/fun*/*/rp*.txt';
end
if ~isfield(rtv, 'GLMOutputFilename')
    rtv.GLMOutputFilename = '';
end
if ~isfield(rtv, 'RedoSubjects')
    rtv.RedoSubjects = false;
end
if ~isfield(rtv, 'RestConds')
    rtv.RestConds = {};
end
if ~isfield(rtv, 'RobustMaxIter')
    rtv.RobustMaxIter = 30;
end
if ~isfield(rtv, 'RobustRegression')
    rtv.RobustRegression = false;
end
if ~isfield(rtv, 'ShowDesigns')
    rtv.ShowDesigns = false;
end
if ~isfield(rtv, 'SingleTrial')
    rtv.SingleTrial = false;
end
if ~isfield(rtv, 'SingleTrialSkip')
    rtv.SingleTrialSkip = {};
end
if ~isfield(rtv, 'TempFilterCutoff')
    rtv.TempFilterCutoff = Inf;
end
if ~isfield(rtv, 'TempFilterType')
    rtv.TempFilterType = 'none';
end
if ~isfield(rtv, 'VarWeight')
    rtv.VarWeight = true;
end
if ~isfield(rtv, 'XTCFilesPattern')
    rtv.XTCFilesPattern = '  subj*/*/fun*/*.vtc';
end
mdm.RunTimeVars = rtv;

% put information out there
try
    ch.Basefolder.String = rtv.BaseFolder;
    ch.FuncPattern.String = rtv.XTCFilesPattern;
    ch.DsgnPattern.String = rtv.DesignFilesPattern;
    ch.UseMotParms.Value = double(rtv.MotionParamsUse);
    ch.FuncFiles.String = mdm.XTC_RTC(:, 1);
    ch.DsgnFiles.String = mdm.XTC_RTC(:, 2);
    ch.MParPattern.String = rtv.MotionParamFilesPattern;
    ch.MParDiff.Value = double(rtv.MotionParamsDiff);
    ch.MParSquared.Value = double(rtv.MotionParamsSquared);
    if numel(rtv.MotionParameters) == size(mdm.XTC_RTC, 1)
        ch.MParFiles.String = rtv.MotionParameters;
    else
        ch.UseMotParms.Value = 0;
        cf.SetGroupEnabled('MotParm', 'off');
        if ~isempty(mdm.XTC_RTC)
            ch.UseMotParms.Enable = 'off';
        end
    end
    if mdm.PSCTransformation > 0
        ch.TransPSC.Value = 1;
        ch.Transz.Value = 0;
        ch.Trans0.Value = 0;
    elseif mdm.zTransformation > 0
        ch.TransPSC.Value = 0;
        ch.Transz.Value = 1;
        ch.Trans0.Value = 0;
    else
        ch.TransPSC.Value = 0;
        ch.Transz.Value = 0;
        ch.Trans0.Value = 1;
    end
    ch.IThresh.String = sprintf('%g', rtv.IntensityThreshold);
    if ~isempty(rtv.MaskFilename)
        ch.MaskFile.String = {rtv.MaskFilename};
    else
        ch.MaskFile.String = {'Click to select...'};
    end
    if rtv.RobustRegression
        ch.RegrOLS.Value = 0;
        ch.RegrRobust.Value = 1;
    else
        ch.RegrOLS.Value = 1;
        ch.RegrRobust.Value = 0;
    end
    ch.Deconv.Value = double(rtv.Deconvolution);
    ch.DeconvLags.String = sprintf('%d', rtv.DeconvolutionLags);
    if ~isempty(rtv.TempFilterType) && ...
        lower(rtv.TempFilterType(1)) == 'd'
        ch.TFiltDCT.Value = 1;
        ch.TFiltFourier.Value = 0;
        ch.TFiltNone.Value = 0;
    elseif ~isempty(rtv.TempFilterType) && ...
        lower(rtv.TempFilterType(1)) == 'f'
        ch.TFiltDCT.Value = 0;
        ch.TFiltFourier.Value = 1;
        ch.TFiltNone.Value = 0;
    else
        ch.TFiltDCT.Value = 0;
        ch.TFiltFourier.Value = 0;
        ch.TFiltNone.Value = 1;
    end
    ch.TFiltCutOff.String = sprintf('%d', rtv.TempFilterCutoff);
    if mdm.RFX_GLM > 0
        cf.RadioGroupSetOne('Model', 4);
        mdm.RFX_GLM = 1;
        mdm.SeparatePredictors = 2;
    else
        mdm.RFX_GLM = 0;
        switch mdm.SeparatePredictors
            case {0}
                cf.RadioGroupSetOne('Model', 1);
            case {1}
                cf.RadioGroupSetOne('Model', 3);
            case {2}
                cf.RadioGroupSetOne('Model', 2);
            otherwise
                cf.RadioGroupSetOne('Model', 4);
                mdm.RFX_GLM = 1;
                mdm.SeparatePredictors = 2;
        end
    end
    ch.CombineFFX.Value = 0;
    if mdm.RFX_GLM
        ch.CombineFFX.Value = double(rtv.CombineFFX);
    end
    ch.GLMFile.String = rtv.GLMOutputFilename;
    ch.RedoSubjects.Value = double(rtv.RedoSubjects);
    ch.ShowDesigns.Value = double(rtv.ShowDesigns);
    if isempty(mdm.XTC_RTC)
        ch.FuncFiles.Value = [];
        ch.DsgnFiles.Value = [];
    else
        cf.SetGroupEnabled('FFound', 'on');
    end
    if rtv.MotionParamsUse
        cf.SetGroupEnabled('MotParm', 'on');
    else
        cf.SetGroupEnabled('MotParm', 'off');
    end
    ch.SingleTrial.Value = double(rtv.SingleTrial);
    ch.VarWeight.Value = double(rtv.VarWeight);
    bf = rtv.BaseFolder;
    if isempty(bf) || ...
        exist(bf, 'dir') ~= 7 || ...
        numel(dir(bf)) == 0
        cf.SetGroupEnabled('FoldOK', 'off');
    else
        cf.SetGroupEnabled('FoldOK', 'on');
    end
    switch rtv.Derivatives
        case {1}
            cf.SetGroupEnabled('Derivs', 'on');
            cf.RadioGroupSetOne('Deriv', 2);
        case {2}
            cf.SetGroupEnabled('Derivs', 'on');
            cf.RadioGroupSetOne('Deriv', 3);
        otherwise
            cf.SetGroupEnabled('Derivs', 'off');
            cf.RadioGroupSetOne('Deriv', 1);
    end
    cf.ShowPage(1);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Error in MDM object. Please check on the console!', ...
        'NeuroElf - warning', 'modal'));
    return;
end
