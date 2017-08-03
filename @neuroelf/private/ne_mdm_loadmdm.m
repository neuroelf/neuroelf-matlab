% FUNCTION ne_mdm_loadmdm: load MDM from disk
function ne_mdm_loadmdm(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:35 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
ch = ne_gcfg.h.MDM.h;

% try to load MDM
try
    newmdm = xff('*.mdm');
    if isempty(newmdm)
        return;
    end
catch ne_eo;
    uiwait(warndlg(ne_eo.message, 'NeuroElf - error', 'modal'));
    return;
end

% clear old MDM, replace
cc.mdm.ClearObject;
cc.mdm = newmdm;
ne_gcfg.fcfg.MDM.mdm = newmdm;

% update UI
ne_mdm_updateui;
if ~any(cellfun('isempty', regexpi(ch.DsgnFiles.String, '\.prt$')))
    cf.SetGroupEnabled('PRTDsgn', 'on');
else
    cf.SetGroupEnabled('PRTDsgn', 'off');
end
