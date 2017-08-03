function s = root_Config(xo, type, N, V)
% ROOT::Config  - read or make configuration setting
%
% FORMAT:       oldvalue = root.Config(type, setting [, value])
%
% Input fields:
%
%       type        3- or 4-letter type for which setting is accessed
%       setting     which setting
%       value       optional new value to be set
%
% Output fields:
%
%       oldvalue    old value of settings
%
% TYPES: ROOT

% Version:  v1.1
% Build:    16012313
% Date:     Jan-23 2016, 1:22 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2016, Jochen Weber
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

% no more valid inputs
conf = xffsngl.CONF;
ext  = fieldnames(xffsngl.EXT);
if nargin < 2 || ~ischar(type) || (~any(strcmpi(type(:)', ext)) && ~any(strcmp(type(:)', fieldnames(conf))))
    s = conf;
    return;
end
type = lower(type(:)');

% type-setting
if any(strcmp(type, ext))

    % get config
    if isfield(conf.type, type)
        s = conf.type.(type);
    else
        s = struct;
    end
    
    % no more valid inputs
    if nargin < 3 || ~ischar(N) || isempty(N) || (~isfield(s, N(:)') && nargin < 4)
        return;
    end
    if isfield(s, N(:)')
        s = s.(N(:)');
    else
        s = [];
    end

    % valid input
    if nargin > 3 && isvarname(N(:)') && ~isempty(V)

        % set
        xffsngl.CONF.type.(type).(N(:)') = V;
    end

% general config setting
else

    % get config
    s = conf.(type);

    % valid input
    if nargin > 2 && strcmp(class(s), class(N)) && ~isempty(N)
        xffsngl.CONF.(type) = N;
    end
end
