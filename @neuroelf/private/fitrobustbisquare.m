function [b, r, w, stats] = fitrobustbisquare(X, y, tune, iterlim)
% fitrobustbisquare  - faster fitting function for robust bisquare
%
% FORMAT:  [b, r, w, stats] = fitrobustbisquare(X, y [, tune [, iterlim]])
%
% Input fields:
%
%       X           model (will NOT be patched, so must contain confounds)
%       y           data
%       tune        bisquare tuning parameter (default: 4.685)
%       iterlim     iteration limit
%
% Output fields:
%
%       b           beta weights
%       r           full residuals
%       w           weighting vector
%       stats       struct with fields as robustfit
%
% Note: other than the robustfit function of the stats toolbox, NaNs
%       must be removed *prior* to regression/estimation, manually

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

% check arguments
if nargin < 2 || ...
   ~isa(X, 'double') || ...
   ~isa(y, 'double') || ...
    isempty(X) || ...
    ndims(X) > 2 || ...
    size(X, 2) > size(X, 1) || ...
    any(isinf(X(:)) | isnan(X(:))) || ...
    numel(y) ~= size(X, 1) || ...
    numel(y) ~= max(size(y)) || ...
    any(isinf(y) | isnan(y))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
   ~isa(tune, 'double') || ...
    numel(tune) ~= 1 || ...
    isinf(tune) || ...
    isnan(tune) || ...
    tune <= 0
    tune = 4.685;
end
if nargin < 4 || ...
   ~isa(iterlim, 'double') || ...
    numel(iterlim) ~= 1 || ...
    isinf(iterlim) || ...
    isnan(iterlim) || ...
    iterlim < 1 || ...
    iterlim > 400
    iterlim = 50;
end
[n, p] = size(X);
op = ones(1, p);

% find the least squares solution.
[Q, R, perm] = qr(X, 0);
tol = abs(R(1)) * n * eps(class(R));
xrank = sum(abs(diag(R)) > tol);
if xrank ~= p
    error( ...
        'neuroelf:RankDeficient', ...
        'X is rank deficient, rank = %d', ...
        xrank ...
    );
end

% fill beta estimates
b(perm, :) = R \ (Q' * y);
b0 = zeros(size(b));

% adjust residuals using leverage, as advised by DuMouchel & O'Brien
E = X(:, perm) / R;
h = min(.9999, sum(E .* E, 2));
adjfactor = 1 ./ sqrt(1 - h);
dfe = n - xrank;
ols_s = norm(y - X * b) / sqrt(dfe);

% handle case of "almost perfect fit"
D = sqrt(eps(class(X)));
tiny_s = 1e-6 * sqrt(varc(y));
if tiny_s == 0
    tiny_s = D;
end

% perform iteratively reweighted least squares to get coefficient estimates
iter = 0;
wxrank = xrank;    % rank of weighted version of x
strect = struct('RECT', true);
while any(abs(b - b0) > D * max(abs(b), abs(b0))) || iter == 0
    iter = iter + 1;
    if (iter > iterlim)
        break;
    end

    % compute residuals from previous fit, then compute scale estimate
    r = y - X * b;
    radj = r .* adjfactor;
    s = madsigma(radj, wxrank);

    % compute new weights from these residuals, then re-fit
    wa = radj / (max(s, tiny_s) * tune);
    w = (abs(wa) < 1) .* (1 - wa .^ 2) .^ 2;
    b0 = b;
    sw = sqrt(w);
    yw = y .* sw;
    xw = X(:, perm) .* sw(:, op);
    [b(perm), wxrank] = linsolve(xw, yw, strect);
end

% compute final (raw) residual
r = y - X * b;

if (nargout > 3)
    radj = r .* adjfactor;
    mad_s = madsigma(radj, xrank);

    % compute a robust estimate of s
    if all(w < D | w > 1 - D)

        % all weights 0 or 1, this amounts to ols using a subset of the data
        included = (w > 1 - D);
        robust_s = norm(r(included)) / sqrt(sum(included) - xrank);
    else
        % compute robust mse according to DuMouchel & O'Brien (1989)
        robust_s = myrobustsigma(radj, xrank, max(mad_s, tiny_s), tune, h);
    end

    sigma = max(robust_s, ...
        sqrt((ols_s ^ 2 * xrank ^ 2 + robust_s ^ 2 * n) / (xrank ^ 2 + n)));

    % get coefficient standard errors and related quantities
    RI = R \ eye(xrank);
    tempC = (RI * RI') * sigma ^ 2;
    tempse = sqrt(max(eps(class(tempC)), diag(tempC)));
    C = nan(p);
    se = zeros(p, 1);
    covb(perm, perm) = tempC;
    C(perm, perm) = tempC ./ (tempse * tempse');
    se(perm) = tempse;

    % save everything
    stats.ols_s = ols_s;
    stats.robust_s = robust_s;
    stats.mad_s = mad_s;
    stats.s = sigma;
    stats.resid = r;
    stats.rstud = r .* adjfactor / sigma;
    stats.se = se;
    stats.covb = covb;
    stats.coeffcorr = C;
    stats.t = nan(size(b));
    stats.t(se > 0) = b(se > 0) ./ se(se > 0);
    stats.p = 2 * sdist('tcdf', -abs(stats.t), dfe);
    stats.w = w;
    z = zeros(p);
    z(perm,perm) = R(1:xrank, 1:xrank);
    stats.R = z;
    stats.dfe = dfe;
    stats.h = h;
    stats.iter = iter;
end

% -----------------------------

function s = madsigma(r, p)
%MADSIGMA    Compute sigma estimate using MAD of residuals from 0
rs = sort(abs(r));
ns = numel(rs) + p;
if mod(ns, 2) == 0
    s = rs(ns / 2) / 0.6745;
else
    s = (rs(ns / 2 - 0.5) + rs(ns / 2 + 0.5)) / 1.349;
end

function s = myrobustsigma(r, p, s, t, h)

% include tuning constant in sigma value
st = s * t;

% get standardized residuals
n = length(r);
u = r ./ st;

% compute derivative of phi function
phi = u .* (abs(u) < 1) .* (1 - u .^ 2) .^ 2;
delta = 0.0001;
u1 = u - delta;
phi0 = u1 .* (abs(u1) < 1) .* (1 - u1 .^ 2) .^ 2;
u1 = u + delta;
phi1 = u1 .* (abs(u1) < 1) .* (1 - u1 .^ 2) .^ 2;
dphi = (phi1 - phi0) ./ (2 * delta);

% compute means of dphi and phi^2
m1 = mean(dphi);
m2 = sum((1 - h) .* phi .^ 2) / (n - p);

% compute factor that is called K by Huber and O'Brien
K = 1 + (p / n) * (1 - m1) / m1;

% compute final sigma estimate
s = K * sqrt(m2) * st / m1;
