% FUNCTION ne_cm_rencon: rename currently selected contrast
function ne_cm_rencon(varargin)

% Version:  v1.1
% Build:    16022911
% Date:     Feb-29 2016, 11:42 AM EST
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

% get current name and weights
scon = ch.Contrasts.Value;
oldname = cc.cons{scon, 1};
oldweights = cc.cons{scon, 2};

% if only one weight is set, change name
if sum(oldweights ~= 0) == 1 && ...
    sum(oldweights) > 0
    oldname = cc.preds{find(oldweights ~= 0)};

% if two are set and sum is 0
elseif sum(oldweights ~= 0) == 2 && ...
    sum(oldweights) == 0
    oldname = [cc.preds{find(oldweights > 0)} ' > ' ...
        cc.preds{find(oldweights < 0)}];
end

% ask for new name
newname = inputdlg({'Please enter the contrast''s new name:'}, ...
    'NeuroElf GUI - input', 1, {oldname});
if isequal(newname, 0) || ...
    isempty(newname)
    return;
end
if iscell(newname)
    newname = newname{1};
end

% set back to config and control
ne_gcfg.fcfg.CM.cons{scon, 1} = newname;
ch.Contrasts.String = ne_gcfg.fcfg.CM.cons(:, 1);

% and update current GLM
ne_cm_updatertv;
ne_cm_updateuis(0, 0, cc.glm);
