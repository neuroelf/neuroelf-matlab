function b = bbubble(b, p, n)
% bbubble  - apply a Brownian-like bubble motion in first non-singleton dim
%
% FORMAT:       data = bbubble(data [, prob [, iter]])
%
% Input fields:
%
%       data        any non-singleton data
%       prob        probability for a motion to occur (default: 0.2)
%       iter        number of iterations (default: 1)
%
% Output fields:
%
%       data        shuffled data

% Version:  v0.9d
% Build:    14061414
% Date:     Jun-14 2014, 2:20 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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
if nargin < 1
    error( ...
        'neuroelf:BadArgument', ...
        'Missing first argument.' ...
    );
end
if numel(b) < 2
    return;
end
s = size(b);
if numel(b) ~= s(1)
    if ndims(b) > 2
        b = squeeze(b);
    else
        b = b(:);
    end
    sb = size(b);
    if numel(b) ~= sb(1)
        b = reshape(b, sb(1), prod(sb(2:end)));
    end
end
s1 = size(b, 1);
s1m = s1 - 1;
if nargin < 2 || ...
   ~isa(p, 'double') || ...
    numel(p) ~= 1 || ...
    isinf(p) || ...
    isnan(p) || ...
    p <= 0 || ...
    p > 1
    p = 0.2;
else
    while p > 0.25
        p = 0.5 * p;
    end
end
p = round(p * (s1 - 1));
if nargin < 3 || ...
   ~isa(n, 'double') || ...
    numel(n) ~= 1 || ...
    isinf(n) || ...
    isnan(n)
    n = 1;
else
    n = max(1, round(n));
end

% iterate
for nc = 1:n

    % shuffle
    sh = lsqueeze(randperm(s1m));
    sh = sort(sh(1:p));
    sh(1 + find(diff(sh) == 1)) = [];
    while numel(sh) < p
        shx = lsqueeze(randperm(s1m));
        shx(any(shx * ones(1, numel(sh)) == ones(numel(shx), 1) * sh(:)', 2)) = [];
        shx = shx(1:(p-numel(sh)));
        sh = sort([sh; shx]);
        sh(1 + find(diff(sh) == 1)) = [];
    end

    % create shuffling argument
    src = [sh'; (sh + 1)'];
    trg = lsqueeze(src([2, 1], :));
    src = src(:);

    % shuffle
    b(trg, :) = b(src, :);
end

% reshape again
b = reshape(b, s);
