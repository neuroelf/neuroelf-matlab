% FUNCTION ne_cm_selectsubs: change subject selection
function ne_cm_selectsubs(varargin)

% Version:  v1.1
% Build:    16022911
% Date:     Feb-29 2016, 11:50 AM EST
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
glm = cc.glm;

% get subjects selection
ssel = ch.Subjects.Value;
ssel = ssel(:)';

% if groups are not selected
if ~cc.usegroups

    % apply to SubSel instead
    sstr = ch.Subjects.String;
    if ~iscell(sstr)
        sstr = cellstr(sstr);
    end
    ne_gcfg.fcfg.CM.subsel = sstr(ssel);
    cc = ne_gcfg.fcfg.CM;
    glm.RunTimeVars.SubSel = cc.subsel;
    glm.RunTimeVars.SubSels{findfirst( ...
        strcmpi(glm.RunTimeVars.SubSels(:, 1), ...
        ch.SubSels.String{ch.SubSels.Value})), 2} = cc.subsel;
    return;
end

% apply to current group
cgrp = ch.Groups.Value;
cc.groups{cgrp, 2} = ssel;

% remove from all other groups
ogrp = setdiff(1:size(cc.groups, 1), cgrp);
for gc = ogrp(:)'
    cc.groups{gc, 2} = setdiff(cc.groups{gc, 2}(:)', ssel);
end

% and set back in selection
ne_gcfg.fcfg.CM = cc;

% update in GLM
glm.RunTimeVars.Groups = cc.groups;
