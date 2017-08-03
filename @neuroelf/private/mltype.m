function mt = mltype(t)
% mltype  - return Matlab's datatype association
%
% FORMAT:       mt = mltype(t)
%
% Input fields:
%
%       t           either name or numeric type in Matlab
%
% Output fields:
%
%       mt          corresponding type
%
% Note: this function is rather useful for coders, not so much for
%       users...

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
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

% persistent list of datatypes
persistent ml_it;
if isempty(ml_it)
    ml_it = {'cell', 'struct', 'logical', 'char', 'void', ...
        'double', 'single', 'int8', 'uint8', 'int16', 'uint16', ...
        'int32', 'uint32', 'int64', 'uint64'};
end

% argument check
if nargin ~= 1 || ...
   ((~ischar(t) || ...
     ~any(strcmp(t(:)', ml_it))) && ...
    (~isa(t, 'double') || ...
     numel(t) ~= 1 || ...
     isinf(t) || ...
     isnan(t) || ...
     ~any(1:15 == t)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument' ...
    );
end

% char
if ischar(t)
    mt = find(strcmp(t(:)', ml_it));

% double
else
    mt = ml_it{t};
end
