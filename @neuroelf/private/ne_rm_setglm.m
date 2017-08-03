% FUNCTION ne_rm_setglm: change the current GLM
function ne_rm_setglm(varargin)

% Version:  v0.9b
% Build:    11050712
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

% get shortcuts
cc = ne_gcfg.fcfg.RM;
ch = ne_gcfg.h.RM.h;
ci = ne_gcfg.c.ini;
try

    % get currently selected GLM
    glm = cc.GLMs{ch.GLMs.Value};
    ne_gcfg.fcfg.RM.GLM = glm;
    ne_gcfg.h.RM.RMFig.UserData.lastglm = glm;

    % for now, only RFX-GLMs supported
    if glm.ProjectTypeRFX <= 0
        uiwait(msgbox('At the moment, only RFX-GLMs are supported', ...
            'NeuroElf GUI - Info', 'modal'));
        return;
    end

    % set subjects
    subs = glm.Subjects;
    ch.Subjects.String = subs;
    ch.Subjects.Value = 1:numel(subs);
    ch.Subjects.ListboxTop = 1;

    % set X and Y to covariates
    ch.XType.Value = ci.Mediation.XType;
    ch.YType.Value = ci.Mediation.YType;

    % set M and C to first contrast/covariate and empty
    ch.MCons.String = glm.RunTimeVars.Contrasts(:, 1);
    ch.MCons.Value = 1;
    ch.MCons.ListboxTop = 1;
    ch.CCons.String = glm.RunTimeVars.Contrasts(:, 1);
    ch.CCons.Value = 1;
    ch.CCons.ListboxTop = 1;
    ch.MCovs.String = glm.RunTimeVars.CovariatesNames;
    ch.MCovs.Value = 1;
    ch.MCovs.ListboxTop = 1;
    ch.CCovs.String = glm.RunTimeVars.CovariatesNames;
    ch.CCovs.Value = 1;
    ch.CCovs.ListboxTop = 1;
    if ci.Mediation.XType == 1 || ...
        ci.Mediation.YType == 1
        ch.MCons.Value = [];
        ch.MCovs.Value = 1;
    else
        ch.MCons.Value = 1;
        ch.MCovs.Value = [];
    end
    ch.CCons.Value = [];
    ch.CCovs.Value = [];

    % ensure validity
    ne_rm_select(0, 0, 'xtype');
    ne_rm_select(0, 0, 'ytype');

catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    ne_cm_closeui;
    uiwait(warndlg('Error using the selected GLM; closing UI.', 'NeuroElf GUI - error', 'modal'));
    return;
end
