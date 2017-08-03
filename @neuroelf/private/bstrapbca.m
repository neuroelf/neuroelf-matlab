function z = bstrapbca(pe, jk, b, n)
% bstrapbca  - compute a z-statistic based on the BCa bootstrap
%
% FORMAT:       z = bstrapbca(pe, jk, b, n)
%
% Input fields:
%
%       pe          point-estimate of the statistic
%       jk          jack-knife of statistic
%       b           boot-strapped values of statistic
%       n           null-value (value against which the pe is tested)
%
% Output fields:
%
%       z           z-score that pe is different from n (signed, 1-tailed)
%
% Notes: this computes the bias-corrected and accelerated bootstrap
%        estimate of  Efron, B., & Tibshirani, R. J. (1993):
%        An Introduction to the Bootstrap, Chapman & Hall/CRC: New York.

% Version:  v0.9b
% Build:    10101414
% Date:     Aug-27 2010, 2:01 PM EST
% Author:   Jared Van Snellenberg, CANlab, Columbia University, NYC, NY, USA
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2009, Jared Van Snellenberg
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
if nargin < 3 || ...
   ~isa(pe, 'double') || ...
    isempty(pe) || ...
    ndims(pe) > 2 || ...
    size(pe, 1) ~= 1 || ...
    any(isinf(pe) | isnan(pe)) || ...
   ~isa(jk, 'double') || ...
    ndims(jk) > 2 || ...
    size(jk, 1) < 3 || ...
    size(jk, 2) ~= size(pe, 2) || ...
    any(isinf(jk(:)) | isnan(jk(:))) || ...
   ~isa(b, 'double') || ...
    ndims(b) > 2 || ...
    size(b, 1) < (2 * size(jk, 1)) || ...
    size(b, 2) ~= size(pe, 2) || ...
    any(isinf(b(:)) | isnan(b(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 4 || ...
   ~isa(n, 'double') || ...
   (numel(n) ~= 1 && ...
    ~isequal(size(n), size(pe))) || ...
    any(isinf(n) | isnan(n))
    n = 0;
end

% number of bootstrap samples
nb = size(b, 1);

% 1 / number of bootstrap samples (required for a few computations)
nbi = 1 / nb;

% compute mean-removed jack-knife sample (in-place)
nj = size(jk, 1);
jk = repmat((1 / nj) .* sum(jk), nj, 1) - jk;

% square sample (used twice in pseudo-skew formula)
jk2 = jk .* jk;

% compute acceleration factor from pseudo-skew of jack-knife statistics
a = (1 / 6) .* (sum(jk2 .* jk) ./ ((sum(jk2) .^ 1.5)));

% get z-score to account for placement of pe within bootstrap distribution
if numel(pe) == 1
    z0 = sdist('norminv', nbi .* sum(b < pe), 0, 1);
else
    z0 = sdist('norminv', nbi .* sum(b < pe(ones(1, nb), :)), 0, 1);
end

% get stochastic p-value for how many bootstrap samples are < null
if numel(n) == 1
    p = nbi .* sum(b < n, 1);
else
    p = nbi .* sum(b < n(ones(1, nb), :));
end

% ensure that p doesn't hit 0 or 1
p(p == 0) = nbi;
p(p == 1) = 1 - nbi;

% compute z-score
z = sdist('norminv', p, 0, 1);

% apply acceleration and correction
za = a .* z0;
z = -(z .* (1 - za) + z0 .* (za - 2)) ./ (1 + a .* z - za);
