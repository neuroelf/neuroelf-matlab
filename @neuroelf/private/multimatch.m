function m = multimatch(c1, c2, ur, ex)
% multimatch  - match two lists of strings and return matched in first
%
% FORMAT:       matching = multimatch(list1, list2 [, userxi [, exact]])
%
% Input fields:
%
%       list1       Lx1 cell array of strings
%       list2       Cx1 cell array of strings to compare/match with
%       userxi      flag, use regexpi for matching (default: false)
%       exact       flag, use strcmp for matching (default: false)
%
% Output fields:
%
%       m           Lx1 integer, string was matched (by Nth cell, 0 if not)

% Version:  v0.9d
% Build:    14061709
% Date:     Jun-17 2014, 9:51 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

% argument check
if nargin < 2 || ...
   ~iscell(c1) || ...
   ~iscell(c2)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
   ~islogical(ur) || ...
    numel(ur) ~= 1
    ur = false;
end
if nargin < 4 || ...
   ~islogical(ex) || ...
    numel(ex) ~= 1
    ex = false;
end
c1 = c1(:);
for c = 1:numel(c1)
    if ~ischar(c1{c}) || ...
        isempty(c1{c}) || ...
        numel(c1{c}) ~= size(c1{c}, 2)
        c1{c} = '';
    end
end

% prepare output
m = zeros(numel(c1), 1);

% regexpi
if ur

    % use compiled function
    for c2c = 1:numel(c2)
        mi = ~cellfun('isempty', regexpi(c1, c2{c2c}));
        if any(mi)
            m(mi) = c2c;
        end
    end

% strcmp
elseif ex

    % iterate also
    for c2c = 1:numel(c2)
        mi = strcmp(c1, c2{c2c});
        if any(mi)
            m(mi) = c2c;
        end
    end

% strcmpi
else

    % iterate also
    for c2c = 1:numel(c2)
        mi = strcmpi(c1, c2{c2c});
        if any(mi)
            m(mi) = c2c;
        end
    end
end
