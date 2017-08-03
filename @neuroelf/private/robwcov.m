function [wcov, w, md, wmean] = robwcov(x, maxiter)
%ROBWCOV  Robust estimate of covariance matrix
%   WCOV = ROBWCOV(X) computes the robustly (iterative least-squares)
%   estimated covariance matrix of X (with up to 50 iterations).
%
%   [WCOV, W] = ROBWCOV(X) also returns the weights vector W applied to X.
%
%   [WCOV, W, WMD] = ROBWCOV(X) also returns the Mahalanobis distance,
%   estimated using the (inverse of the) weighted covariance matrix.
%
%   [WCOV, W, WMD, WMEAN] = ROBWCOV(X) also returns the weighted mean.
%
%   WCOV = ROBWCOV(X, NITER) sets the maximum number of iterations.

% Version:  v1.1
% Build:    16032917
% Date:     Mar-29 2016, 5:08 PM EST
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
% SOFTWARE,` EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% argument check
if nargin < 1 || (~isa(x, 'double') && ~isa(x, 'single')) || ndims(x) > 2 || isempty(x) || size(x, 1) <= size(x, 2)
    error('neuroelf:general:badArgument', 'Invalid data argument.');
end
if nargin < 2 || ~isa(maxiter, 'double') || numel(maxiter) ~= 1 || isinf(maxiter) || isnan(maxiter) || maxiter < 1
    maxiter = 50;
end

% make sure data is double (cost close to 0 if x is already double)
x = double(x);

% get number of rows
sd = size(x, 1);
nd = size(x, 2);
spn = (1:sd)';

% normalization factor, indices for MAD, safety margin, and tuning
sdp = 0.5 * (sd + 1);
sdf = 1 / 0.6745;
if sdp ~= fix(sdp)
    sdp = round(sdp-0.5);
    sdp2 = sdp + 1;
    sdf = 0.5 * sdf;
    sdo = false;
else
    sdo = true;
end
tiny_s = 1e-6 * median(std(x));
tune = 4.685;

% set weights to one
o = ones(sd, 1);
on = ones(1, nd);
w = o;

% discrepancy
omd = zeros(sd, 1);
disc = Inf .* w;

% compute mean
wmean = (1 / sum(w)) .* sum(w(:, on) .* x, 1);

% compute unbiased sample covariance
xbar = x - o * wmean;
wcov = (1 / (sum(w) - 1)) .* (xbar' * sparse(spn, spn, w, sd, sd, sd) * xbar);

% repeat while iterations not run out and any discrepancy > sqrt(eps)
seps = sqrt(eps);
while maxiter > 0 && any(abs(disc) > seps)

    % compute distance
    md = sum((xbar / wcov) .* xbar, 2);

    % discrepancy
    disc = w .* abs(omd - md);
    omd = md;

    % recompute weights
    resmd = abs(md - mean(md));
    sresmd = sort(resmd);
    if sdo
        ms = sdf * sresmd(sdp);
    else
        ms = sdf * (sresmd(sdp) + sresmd(sdp2));
    end
    wa = resmd ./ (max(ms, tiny_s) * tune);
    w = (abs(wa) < 1) .* (1 - wa .* wa) .^ 2;
    w(isnan(w)) = 0;

    % compute weighted mean
    wmean = (1 / sum(w)) .* sum(w(:, on) .* x, 1);

    % compute reweighted, unbiased sample covariance
    xbar = x - o * wmean;
    wcov = (1 / (sum(w) - 1)) .* (xbar' * sparse(spn, spn, w, sd, sd, sd) * xbar);

    % counter
    maxiter = maxiter - 1;
end
