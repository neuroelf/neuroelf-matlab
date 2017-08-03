function ff = findfirst(v, spos)
% findfirst  - find first (or last) element != 0 in input
%
% FORMAT:         ff = findfirst(v [, spos]);
%
% Input fields:
%
%       v             vector
%       spos          start position for search (default: 1),
%                     if negative do search from end to begin
%
% Output fields:
%
%       ff            first occurrance of "true" in v (as of spos)

% Version:  v0.9b
% Build:    11071011
% Date:     Jun-15 2010, 12:37 PM EST
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

% use Matlab code if MEX is unavailable
if nargin < 2 || ...
   ~isa(spos, 'double') || ...
    numel(spos) ~= 1 || ...
    isinf(spos) || ...
    isnan(spos) || ...
    spos <= -numel(v) || ...
    spos > numel(v) || ...
    floor(spos) == 0
    spos = 1;
else
    spos = floor(spos);
end
if spos ~= 1
    if spos > 0
        ff = find(v(spos:end), 1, 'first') + (spos - 1);
    else
        ff = find(v(1:end+1+spos), 1, 'last');
    end
else
    ff = find(v, 1, 'first');
end
