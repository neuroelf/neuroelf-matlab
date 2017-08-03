function [b, r, w] = fitrobustlog(X, y, tune, iterlim)
% fitrobustlog  - faster fitting function for robust using log weights
%
% FORMAT:  [b, r, w] = fitrobustlog(X, y [, tune [, iterlim]])
%
% Input fields:
%
%       X           model (will NOT be patched, so must contain confounds)
%       y           data
%       tune        bisquare tuning parameter (default: 0.001)
%       iterlim     iteration limit
%
% Output fields:
%
%       b           beta weights
%       r           full residuals
%       w           weighting vector
%
% Note: other than the robustfit function of the stats toolbox, NaNs
%       must be removed *prior* to regression/estimation, manually

% Version:  v1.1
% Build:    16050513
% Date:     May-05 2016, 1:01 PM EST
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

% check arguments
if nargin < 2 || ~isa(X, 'double') || ~isa(y, 'double') || ...
    isempty(X) || ndims(X) > 2 || size(X, 2) > size(X, 1) || ...
    any(isinf(X(:)) | isnan(X(:))) || numel(y) ~= size(X, 1) || ...
    numel(y) ~= max(size(y)) || any(isinf(y) | isnan(y))
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
if nargin < 3 || ~isa(tune, 'double') || numel(tune) ~= 1 || ...
    isinf(tune) || isnan(tune) || tune <= 0 || tune >= 0.5
    tune = 0.001;
end
tune = -sdist('norminv', tune, 0, 1);
if nargin < 4 || ~isa(iterlim, 'double') || numel(iterlim) ~= 1 || ...
    isinf(iterlim) || isnan(iterlim) || iterlim < 1 || iterlim > 1000
    iterlim = 50;
end

% follow initially same logic
[n, p] = size(X);
op = ones(1, p);

% find the least squares solution.
[Q, R, perm] = qr(X, 0);
tol = abs(R(1)) * n * eps(class(R));
xrank = sum(abs(diag(R)) > tol);
if xrank ~= p
    error('neuroelf:robustfit:rankDeficient', ...
        'X is rank deficient, rank = %d', xrank);
end

% fill beta estimates
b(perm, :) = R \ (Q' * y);
b0 = zeros(size(b));
mb0 = max(abs(b0));

% handle case of "almost perfect fit"
D = sqrt(eps(class(X)));

% perform iteratively reweighted least squares to get coefficient estimates
iter = 0;
strect = struct('RECT', true);
sw = 1;
ssw = n;
while any(abs(b - b0) > D * max(abs(b), mb0)) || iter == 0
    iter = iter + 1;
    if (iter > iterlim)
        break;
    end

    % compute residuals from previous fit, then compute scale estimate
    r = y - X * b;
    rw = sw .* r;
    rs = sum(rw);
    rv = 1 / sqrt((sum(rw .* rw) - rs * rs / ssw) / (ssw - xrank));

    % compute new weights from these residuals, then re-fit
    w = max(0, 1 - log(1 + abs(rv .* r) ./ tune) .^ 2);
    b0 = b;
    sw = sqrt(w);
    ssw = sum(w);
    yw = y .* sw;
    xw = X(:, perm) .* sw(:, op);
    b(perm) = linsolve(xw, yw, strect);
end

% compute final (raw) residual
r = y - X * b;
