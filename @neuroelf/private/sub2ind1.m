function idx = sub2ind1(sz, s, z)
% sub2ind1  - sub2ind with single argument
%
% FORMAT:       idx = sub2ind1(sz, s [, z])
%
% Input fields:
%
%       sz          1xS size of array
%       s           IxS index expressions
%       z           optional boolean flag (true := zero-based)
%
% Output fields:
%
%       idx         Ix1 index expression to access an array

% Version:  v0.9c
% Build:    13012515
% Date:     Jan-25 2013, 3:31 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2013, Jochen Weber
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
   ~isa(sz, 'double') || ...
    isempty(sz) || ...
    ndims(sz) > 2 || ...
    length(sz) ~= numel(sz) || ...
    any(isinf(sz) | isnan(sz) | sz < 1 | sz ~= fix(sz)) || ...
   ~isnumeric(s) || ...
    ndims(s) > 2 || ...
    size(s, 2) ~= numel(sz)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
sz = sz(:);

% cumprod of size for multiplication
sz = [1; cumprod(sz(1:end-1))];

% 1-based
if nargin < 3 || ...
   ~islogical(z) || ...
    numel(z) ~= 1 || ...
   ~z

    % sum of product
    idx = 1 + (double(s) - 1) * sz;

% 0-based
else
    idx = double(s) * sz;
end
