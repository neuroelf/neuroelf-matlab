% FUNCTION ne_cm_setstype: set statistics type
function ne_cm_setstype(varargin)

% Version:  v1.1
% Build:    16052413
% Date:     May-24 2016, 1:13 PM EST
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
ch = ne_gcfg.h.CM.h;
hFig = ne_gcfg.h.CM.CMFig;

% what value
stype = ch.StatType.String{ch.StatType.Value};

% update UI settings
if isempty(regexpi(stype, 'robust'))
    hFig.SetGroupEnabled('Robust', 'off');
else
    hFig.SetGroupEnabled('Robust', 'on');
end
if ~isempty(regexpi(stype, 'spm'))
    ch.RFXstats.Enable = 'off';
    ch.RFXstats.Value = 1;
else
    glm = ne_gcfg.fcfg.CM.glm;
    if isxff(glm, 'glm') && (glm.ProjectTypeRFX > 0 || glm.SeparatePredictors == 2)
        hFig.SetGroupEnabled('RFXGLM', 'on');
        ch.RFXstats.Value = 1;
        if glm.ProjectTypeRFX > 0
            ch.RFXstats.Enable = 'off';
        elseif numel(glm.Subjects) < 3
            ch.RFXstats.Enable = 'off';
            ch.RFXstats.Value = 0;
        end
    elseif isxff(glm, 'glm')
        hFig.SetGroupEnabled('RFXGLM', 'off');
        ch.RFXstats.Value = 0;
    end
end
