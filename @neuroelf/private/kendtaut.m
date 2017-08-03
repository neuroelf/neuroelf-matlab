function t = kendtaut(v, tau, opts)
% kendtaut  - give approximate t-score for kendall tau (via monte-carlo)
%
% FORMAT:       t = kendtaut(v, tau [, opts])
%
% Input fields:
%
%       v           vector for which tau was computed
%       tau         N-dim array with tau values
%       opts        optional settings
%        .dist      type of distribution to draw random samples from
%                   either of 'binary', 'lognormal', {'normal'}, 'uniform'
%        .pacc      p accuracy factor (default: 20)
%        .pmin      min. p for which t should be accurate (default: 0.001)
%
% Output fields:
%
%       t           approximate t-score with d.f. = numel(v) - 2
%
% Notes: the d.f. corresponds to a correlation-like linear regression
%        to simplify direct comparison of t-scores and/or conjunctions;
%        the number of iterations is computed as 2 * pacc / pmin;
%        additionally, pmin is also the "stepsize" in which t's will be
%        read from the histogram!

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
if nargin < 2 || ...
   ~isnumeric(v) || ...
   ~isa(tau, 'double') || ...
    isempty(v) || ...
    max(size(v)) ~= numel(v)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument.' ...
    );
end
n = numel(v);
if n < 8
    error( ...
        'neuroelf:BadArgument', ...
        'The minimal number of observations for a t estimation is 8.' ...
    );
end
tau(isinf(tau) | isnan(tau)) = 0;
tau(tau < -1) = -1;
tau(tau > 1) = 1;
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dist') || ...
   ~ischar(opts.dist) || ...
    isempty(opts.dist) || ...
   ~any(lower(opts.dist(1)) == 'blnu')
    opts.dist = 'n';
else
    opts.dist = lower(opts.dist(1));
end
if ~isfield(opts, 'pacc') || ...
   ~isa(opts.pacc, 'double') || ...
    numel(opts.pacc) ~= 1 || ...
    isinf(opts.pacc) || ...
    isnan(opts.pacc) || ...
    opts.pacc < 1
    opts.pacc = 20;
else
    opts.pacc = min(1000, opts.pacc);
end
if ~isfield(opts, 'pmin') || ...
   ~isa(opts.pmin, 'double') || ...
    numel(opts.pmin) ~= 1 || ...
    isinf(opts.pmin) || ...
    isnan(opts.pmin) || ...
    opts.pmin <= 0 || ...
    opts.pmin > 16 || ...
   (opts.pmin > 0.1 && ...
    opts.pmin < 1.63)
    opts.pmin = 0.001;
elseif opts.pmin > 1
    opts.pmin = sdist('tcdf', -opts.pmin, n - 2);
end
np = ceil(2 * opts.pacc / opts.pmin);

% perform monte-carlo simulation
v = double(v(:));
switch (opts.dist)
    case {'b'}
        mt = kendtau(v, double(rand(n, np) >= 0.5));
    case {'l'}
        mt = kendtau(v, exp(randn(n, np)));
    case {'n'}
        mt = kendtau(v, randn(n, np));
    case {'u'}
        mt = kendtau(v, rand(n, np));
end

% build histogram in bins of 1 / np
pm = 1 / np;
mth = hist(mt, (-1+pm/2):pm:(1-pm/2));
ch = round(0.5 * numel(mth));
mtp = cumsum(mth) ./ np;
mtp(mtp == 0) = 0.5 * opts.pmin;
mtp(mtp >= (1 - eps)) = 1 - 0.5 * opts.pmin;
mtp(1:ch) = mtp(1:ch) + pm;
mtp(ch+1:end) = mtp(ch+1:end) - pm;

% then resolve
tb = tau < 0;
tau(tb) = ch + ceil((1 / pm) .* tau(tb));
tau(~tb) = (ch + 1) + floor((1 / pm) .* tau(~tb));
t = sdist('tinv', mtp(tau), n - 2);
