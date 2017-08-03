function [p, px, py] = pairs(m)
% pairs  - returns the pairs of values in a vector/matrix
%
% FORMAT:       [p, px, py] = pairs(matrix)
%
% Input fields:
%
%       matrix      vector or matrix (pairs taken along column dim)
%
% Output fields
%
%       p           all possible pairs (in third dim)
%       px, py      indices to get unique pairs (length: 0.5 * M * (M-1))

% Version:  v0.9a
% Build:    11043015
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

% argument check
if nargin ~= 1 || ...
   ~isnumeric(m) || ...
    ndims(m) > 2 || ...
    size(m, 1) < 2 || ...
    isempty(m)
    error( ...
        'neuroelf:BadArgument', ...
        'Function pairs requires exactly one numeric, 2-d argument.' ...
    );
end

% create output of same class
s = size(m);
np = 0.5 * s(1) * (s(1) - 1);
p = m(1);
p(1) = 0;
p(np, s(2), 2) = 0;

% get nd-grid for access
[px, py] = meshgrid(1:s(1)-1, 2:s(1));
px(px >= py) = 0;
py = py(px > 0);
px = px(px > 0);

% fill output
p(:, :, 1) = m(px, :);
p(:, :, 2) = m(py, :);
