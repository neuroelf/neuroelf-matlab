% FUNCTION ne_mdm_compute: compute GLM from MDM
function ne_mdm_compute(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:35 AM EST
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
cc = ne_gcfg.fcfg.MDM;
cf = ne_gcfg.h.MDM.MDMFig;

% get options
ne_mdm_updatemdm;
mdm = cc.mdm;

% create options
rtv = mdm.RunTimeVars;
mgo = struct( ...
    'ithresh',  rtv.IntensityThreshold, ...
    'mask',     rtv.MaskFilename, ...
    'motparsd', rtv.MotionParamsDiff, ...
    'motparsq', rtv.MotionParamsSquared, ...
    'outfile',  rtv.GLMOutputFilename, ...
    'robust',   rtv.RobustRegression);
if mdm.RFX_GLM > 0
    mgo.loadglm = true;
end
if rtv.MotionParamsUse
    mgo.motpars = rtv.MotionParameters;
end
if rtv.Deconvolution
    mgo.ndcreg = rtv.DeconvolutionLags;
end
if any(strcmpi(rtv.TempFilterType, {'dct', 'fourier'})) && ...
   ~isinf(rtv.TempFilterCutoff)
    mgo.tfilter = rtv.TempFilterCutoff;
    mgo.tfilttype = rtv.TempFilterType;
end
if rtv.RedoSubjects
    mgo.redo = true;
end
if rtv.ShowDesigns
    mgo.showsdms = 'plot';
end

% hide window
cf.Visible = 'off';

% set main pointer arrow
mfp = ne_gcfg.h.MainFig.Pointer;
ne_gcfg.h.MainFig.Pointer = 'watch';

% with error handling
try

    % run computation
    glm = cc.mdm.ComputeGLM(mgo);

    % open output in GUI
    neuroelf_gui('openfile', glm);

    % close MDM
    ne_mdm_closeui;

% handle errors
catch ne_eo;

    % make figure visible again
    cf.Visible = 'on';

    % show error
    uiwait(warndlg(['An error occurred:' char([10, 10]) ne_eo.message], ...
        'NeuroElf - error message', 'modal'));
end

% re-set pointer
ne_gcfg.h.MainFig.Pointer = mfp;
