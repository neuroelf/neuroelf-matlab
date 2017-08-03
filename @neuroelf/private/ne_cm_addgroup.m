% FUNCTION ne_cm_addgroup: add a subject group to the contrast manager
function ne_cm_addgroup(varargin)

% Version:  v0.9b
% Build:    11051315
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

% get handle shortcut
ch = ne_gcfg.h.CM.h;

% ask for new group name
newgroup = inputdlg({'Please enter the group''s name:'}, ...
    'NeuroElf GUI - input', 1, {'new group'});
if isequal(newgroup, 0) || ...
    isempty(newgroup)
    return;
end
if iscell(newgroup)
    newgroup = newgroup{1};
end
if any(strcmpi(newgroup, ne_gcfg.fcfg.CM.groups(:, 1)))
    uiwait(warndlg('Group names must be unique.', 'NeuroElf GUI - info', 'modal'));
    return;
end

% add to list of groups
ne_gcfg.fcfg.CM.groups{end+1, 1} = newgroup;
groups = ne_gcfg.fcfg.CM.groups(:, 1);
ch.Groups.String = groups(:);
ch.Groups.Value = numel(groups);

% then pretend we switched via UI
ne_cm_selgroup;

% update in GLM
ne_cm_updatertv;
