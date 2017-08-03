function s = root_UpdateState(xo, type, N)
% ROOT::UpdateState  - read or make UpdateState setting
%
% FORMAT:       oldvalue = root.UpdateState(type [, value])
%
% Input fields:
%
%       type        3- or 4-letter type for which setting is accessed
%       value       optional new value to be set (boolean)
%
% Output fields:
%
%       oldvalue    old value of settings
%
% TYPES: ROOT

% Version:  v1.1
% Build:    16012316
% Date:     Jan-23 2016, 4:08 PM EST
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

% global factory
global xffsngl;

% requires root object
if numel(xo) ~= 1 || ~strcmpi(xo.S.Extensions{1}, 'root')
    error('neuroelf:xff:badObject', 'Config requires ROOT object');
end

% get extensions
s = xffsngl.CONF.update;
exn = fieldnames(s);

% full set (with single number)
if nargin > 1 && islogical(type) && numel(type) == 1

    % generate full set and possible update values
    for tc = 1:numel(exn)
        xffsngl.CONF.update.(exn{tc}) = type;
    end
    return;

% full set with existing settings
elseif nargin > 1 && isstruct(type) && numel(type) == 1 && numel(fieldnames(type)) == numel(exn)

    % generate full set and possible update values
    ifnames = fieldnames(type);
    for tc = 1:numel(exn)
        bext = exn{tc};
        bexti = find(strcmpi(ifnames, bext));
        if numel(bexti) == 1
            newval = type.(ifnames{bexti});
            if islogical(newval) && numel(newval) == 1
                xffsngl.CONF.update.(exn{tc}) = newval;
            end
        end
    end
    return;
end

% no more valid inputs
if nargin < 2 || ~ischar(type) || ~any(strcmpi(type(:)', exn))
    
    % return full list
    return;
end
type = lower(type(:)');
s = s.(type);

% value given
if nargin > 2 && islogical(N) && numel(N) == 1
    xffsngl.CONF.update.(type) = N;
end
