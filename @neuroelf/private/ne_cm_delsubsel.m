% FUNCTION ne_cm_addsubsel: add subject selection to lists
function ne_cm_addsubsel(varargin)

% Version:  v1.0
% Build:    16011019
% Date:     Jan-10 2016, 7:23 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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

% only valid if more than one selection
if size(glm.RunTimeVars.SubSels, 1) < 2
    return;
end

% confirm deletion
vans = questdlg('Are you sure you want to delete this subject selection?', ...
    'NeuroElf - request', 'Yes', 'No', 'Yes');
if ~any(lower(vans) == 'y')
    return;
end

% remove selection from list, then reset
glm.RunTimeVars.SubSels(ch.SubSels.Value, :) = [];
cursels = ch.SubSels.String;
if ~iscell(cursels)
    cursels = cellstr(cursels);
end
cursels(ch.SubSels.Value) = [];
ch.SubSels.Value = 1;
ch.SubSels.String = cursels;
if numel(cursels) < 2
    ch.SubSels.Enable = 'off';
end
ne_cm_selectsubsel;
