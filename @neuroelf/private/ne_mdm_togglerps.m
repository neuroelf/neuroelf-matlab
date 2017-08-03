% FUNCTION ne_mdm_togglerps: toggle motion parameter controls
function ne_mdm_togglerps(varargin)

% Version:  v0.9b
% Build:    11111510
% Date:     Aug-11 2010, 9:00 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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
cf = ne_gcfg.h.MDM.MDMFig;
ch = ne_gcfg.h.MDM.h;

% get current selection
userps = (ch.UseMotParms.Value > 0);

% set controls
cf.SetGroupEnabled('MSelect', 'off');
if userps
    cf.SetGroupEnabled('MotParm', 'on');
    if numel(ch.MParFiles.Value) == 1
        cf.SetGroupEnabled('MSelect', 'on');
    end
else
    cf.SetGroupEnabled('MotParm', 'off');
end

% ensure that FFound is heeded all the same
if isempty(ch.MParFiles.String) || ...
   (numel(ch.MParFiles.String) == 1 && ...
    strcmpi(ch.MParFiles.String, '<no files selected>'))
    ch.MParFiles.Enable = 'off';
    cf.SetGroupEnabled('MSelect', 'off');
end
