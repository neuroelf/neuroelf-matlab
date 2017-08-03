% FUNCTION ne_cm_delcon: delete currently selected contrast
function ne_cm_delcon(varargin)

% Version:  v1.1
% Build:    16022911
% Date:     Feb-29 2016, 11:37 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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
cc = ne_gcfg.fcfg.CM;
ch = ne_gcfg.h.CM.h;

% currently selected contrast
scon = ch.Contrasts.Value;

% remove from config
cc.cons(scon, :) = [];

% if no contrasts left
if isempty(cc.cons)

    % set group disabled
    ne_gcfg.h.CM.CMFig.SetGroupEnabled('HasCons', 'off');

    % and set default config again
    ch.Contrasts.Value = 1;
    ch.Contrasts.String = {'<as currently configured>'};

    % then write back to global var and return
    ne_gcfg.fcfg.CM = cc;
    return;
end

% update controls
scon = min(scon, size(cc.cons, 1));
ch.Contrasts.String = cc.cons(:, 1);
ch.Contrasts.Value = scon;
ne_cm_setweights(cc.cons{scon, 2});

% and write to config
ne_gcfg.fcfg.CM = cc;
ne_cm_updatertv;
ne_cm_updateuis(0, 0, cc.glm);
