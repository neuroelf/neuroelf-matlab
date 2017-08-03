function s = root_TransIOSize(xo, type, N)
% ROOT::TransIOSize  - read or make TransIOSize setting
%
% FORMAT:       oldvalue = root.TransIOSize(type [, value])
%
% Input fields:
%
%       type        3- or 4-letter type for which setting is accessed
%       value       optional new value to be set
%
% Output fields:
%
%       oldvalue    old value of settings
%
% TYPES: ROOT

% Version:  v1.1
% Build:    16012314
% Date:     Jan-23 2016, 2:52 PM EST
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
bff = xffsngl.BFF;
ext = xffsngl.EXT;
exn = fieldnames(ext);

% full set (with single number)
if nargin > 1 && isa(type, 'double') && numel(type) == 1 && type >= 1

    % generate full set and possible update values
    s = struct;
    type = round(type);
    for tc = 1:numel(bff)
        s.(bff(tc).Extensions{1}) = bff(tc).TransIOSize;
        xffsngl.BFF(tc).TransIOSize = type;
        xffsngl.FF.bff(tc).TransIOSize = type;
    end
    return;

% full set with existing settings
elseif nargin > 1 && isstruct(type) && numel(type) == 1 && numel(fieldnames(type)) == numel(xffsngl.BFF)

    % generate full set and possible update values
    s = struct;
    ifnames = fieldnames(type);
    for tc = 1:numel(bff)
        bext = bff(tc).Extensions{1};
        s.(bext) = bff(tc).TransIOSize;
        bexti = find(strcmpi(ifnames, bext));
        if numel(bexti) == 1
            newval = type.(ifnames{bexti});
            if isa(newval, 'double') && numel(newval) == 1 && ~isnan(newval) && newval >= 1
                xffsngl.BFF(tc).TransIOSize = round(newval);
                xffsngl.FF.bff(tc).TransIOSize = round(newval);
            end
        end
    end
    return;
end

% no more valid inputs
if nargin < 2 || ~ischar(type) || ~any(strcmpi(type(:)', exn))
    
    % compile list
    s = struct;
    for tc = 1:numel(bff)
        s.(bff(tc).Extensions{1}) = bff(tc).TransIOSize;
    end
    return;
end
type = lower(type(:)');

% type-setting
ttest = ext.(type);
if isempty(regexpi(ttest{1}, '\.bff$'))
    warning('neuroelf:xff:badInputType', 'TransIOSize can only be set for BFF types.');
    return;
end
tti = ttest{2};
s = xffsngl.BFF(tti).TransIOSize;

% value given
if nargin > 2 && isa(N, 'double') && numel(N) == 1 && ~isnan(N) && N >= 1
    xffsngl.BFF(tti).TransIOSize = round(N);
    xffsngl.FF.bff(tti).TransIOSize = round(N);
end
