% FUNCTION ne_glmplotbetaskyr: keep track of released modifiers
function ne_glmplotbetaskyr(varargin)

% Version:  v0.9b
% Build:    11120808
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

% get key
if nargin < 2
    return;
end
ke = varargin{2};

% only valid for valid Tag string
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    numel(varargin{3}) ~= 8 || ...
   ~isfield(ne_gcfg.cc, varargin{3}(:)')
    return;
end
tsat = varargin{3}(:)';

% remove from list
mo = ne_gcfg.cc.(tsat).Config.mods;
if isempty(mo)
    return;
end
mo(strcmpi(mo, ke.Key)) = [];

% and definitely remove command
if ~isempty(mo)
    mo(strcmpi(mo, 'command')) = [];
end
ne_gcfg.cc.(tsat).Config.mods = mo;
