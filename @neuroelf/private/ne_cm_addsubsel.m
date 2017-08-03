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

% name
vans = inputdlg({'Please enter a name for a subject selection:'}, ...
    'NeuroElf - input', 1, {'selection'});
if ~iscell(vans) || ...
    numel(vans) ~= 1 || ...
   ~ischar(vans{1}) || ...
    isempty(vans{1})
    return;
end
if any(strcmpi(ddeblank(vans{1}), ch.SubSels.String))
    uiwait(warndlg('Selection already exists.', 'NeuroElf - warning', 'modal'));
    return;
end
newselname = ddeblank(vans{1});

% add selection to list
glm.RunTimeVars.SubSels(end+1, :) = {newselname, glm.Subjects};
ch.SubSels.String = glm.RunTimeVars.SubSels(:, 1);
ch.SubSels.Value = size(glm.RunTimeVars.SubSels, 1);
ne_cm_selectsubsel;
