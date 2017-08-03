% FUNCTION ne_cm_delgroup: remove group (and set to other group)
function ne_cm_delgroup(varargin)

% Version:  v0.9c
% Build:    12050819
% Date:     Apr-29 2011, 8:11 PM EST
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
cc = ne_gcfg.fcfg.CM;
ch = ne_gcfg.h.CM.h;

% if this leads to only 1 group
if size(cc.groups, 1) < 3

    % ask if all groups are to be removed
    a = questdlg('Do you wish to disable grouping permanently?', ...
        'NeuroElf - user input', 'Yes', 'Cancel', 'Cancel');
    if ~ischar(a) || ...
       ~strcmpi(a, 'yes')
        return;
    end

    % disable groups
    ch.UseGroups.Value = 0;
    ch.Groups.Value = 1;
    ch.Groups.String = {'<no groups specified>'};
    cc.groups = cell(0, 2);
    ne_gcfg.fcfg.CM = cc;
    ne_cm_usegroups;
    ch.Subjects.Value = cc.subsel;
    ch.Subjects.ListboxTop = 1;
    ne_cm_updatertv;
    return;
end

% remove from list and make sure value is OK
cgrp = ch.Groups.Value;
cc.groups(cgrp, :) = [];
ne_gcfg.fcfg.CM = cc;
ch.Groups.Value = min(ch.Groups.Value, size(cc.groups, 1));
ch.Groups.String = cc.groups(:, 1);

% and refresh selection
ne_cm_selgroup;

% update in GLM
ne_cm_updatertv;
