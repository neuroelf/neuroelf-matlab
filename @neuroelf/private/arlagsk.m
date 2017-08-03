function k = arlagsk(r, n, tr)
% arlagsk  - convert AR lag(s) into smoothing kernel
%
% FORMAT:       k = arlagsk(r [, n [, tr]])
%
% Input fields:
%
%       r           auto-correlation r (or map)
%       n           length of time series
%       tr          TR (to generate k in seconds)
%
% Output fields:
%
%       k           smoothing kernel (or map)

% Version:  v0.9c
% Build:    12121816
% Date:     Dec-18 2012, 4:39 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, Jochen Weber
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
if nargin < 1 || ...
   ~isa(r, 'double') || ...
    any(isinf(r(:)) | isnan(r(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isa(n, 'double') || ...
    numel(n) ~= 1 || ...
    isinf(n) || ...
    isnan(n) || ...
    n <= 1
    n = [];
end
if nargin < 3 || ...
   ~isa(tr, 'double') || ...
    numel(tr) ~= 1 || ...
    isinf(tr) || ...
    isnan(tr) || ...
    tr <= 0
    tr = 1;
end

% fix r
r(r < 0) = 0;
r(r > 1) = 1;

% with or without n (length)
if isempty(n)
    k = (tr * sqrt(8 * log(2))) .* sqrt(-1 ./ (4 .* log(r)));
else
    k = (n * tr * sqrt(8 * log(2)) / (n - 1)) .* sqrt(-1 ./ (4 .* log(r)));
end
