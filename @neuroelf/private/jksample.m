function jks = jksample(n, o)
% jksample  - create jack-knife sampling over N cases
%
% FORMAT:       jks = jksample(n [, o])
%
% Input fields:
%
%       n           number of samples in dataset
%       o           omit, either {1} or 2
%
% Output fields:
%
%       jks         jack-knife samples (n-o)x(n! / (n-o! * o!))

% Version:  v0.9b
% Build:    11111711
% Date:     Aug-31 2010, 1:32 PM EST
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
%       documentation and/or other materials provided with the
%       distribution.
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
   ~isnumeric(n) || ...
    numel(n) ~= 1 || ...
    isinf(n) || ...
    isnan(n) || ...
    n < 2 || ...
    n ~= fix(n)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
n = real(double(n));
if nargin < 2 || ...
   ~isnumeric(o) || ...
    numel(o) ~= 1 || ...
    isinf(1) || ...
    isnan(1) || ...
    o <= 1
    o = 1;
else
    o = 2;
end

% create 1-based sampling
jks = (1:(n-1))' * ones(1, n) + tril(ones(n - 1, n));

% simple case
if o == 1

    % just return
    return;
end

% create array
jkn = zeros(n - 2, 0.5 * n * (n - 1));

% reduce size
n = n - 1;

% create smaller sample
jki = reshape((1:(n-1))' * ones(1, n) + tril(ones(n - 1, n)), [n - 1, 1, n]);

% combine samples
ti = 1;
for nc = 1:n
    if nc < n
        jkn(:, ti:ti+n-nc) = indexarray(jks, jki(:, 1, nc:end), nc);
        ti = ti + (1 + n - nc);
    else
        jkn(:, end) = jks(jki(:, end), end - 1);
    end
end

% reset
jks = jkn;
