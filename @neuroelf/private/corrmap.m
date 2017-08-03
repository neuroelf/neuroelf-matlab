function cm = corrmap(r, m, opts)
% corrmap  - correlation map
%
% FORMAT:       cm = corrmap(r, m [, opts]);
%
% Input fields:
%
%       r           regressor (covariate)
%       m           maps (to be regressed)
%       opts        optional settings struct
%        .dim       regression dim (otherwise autodetect, last-to-first)
%        .nomean    adds mean of each map to model
%        .output    either {'r'} or 't'
%        .type      either of 'kendall', {'ols'}, 'robust'
%
% Output fields:
%
%       cm          either r or t correlation map with N-2 (or N-3) d.f.

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

% argument check
if nargin < 2 || ...
   ~isnumeric(r) || ...
    numel(r) ~= max(size(r)) || ...
    any(isinf(r) | isnan(r)) || ...
   ~isnumeric(m) || ...
   ~any(size(m) == numel(r))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
r = double(r(:));
nr = numel(r);
if ~isa(m, 'double')
    m = double(m);
end
ms = size(m);
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dim') || ...
    numel(opts.dim) ~= 1 || ...
   ~isa(opts.dim, 'double') || ...
    isinf(opts.dim) || ...
    isnan(opts.dim) || ...
    opts.dim < 1 || ...
    opts.dim > numel(ms) || ...
    ms(fix(opts.dim)) ~= nr
    opts.dim = findfirst(ms == nr, -1);
else
    opts.dim = fix(opts.dim);
end
rs = ms;
rs(opts.dim) = [];
if numel(rs) < 2
    rs(2) = 1;
end
if ~isfield(opts, 'nomean') || ...
   ~islogical(opts.nomean') || ...
    numel(opts.nomean) ~= 1
    opts.nomean = false;
end
if ~isfield(opts, 'output') || ...
   ~ischar(opts.output) || ...
    numel(opts.output) ~= 1 || ...
   ~any('rt' == lower(opts.output))
    opts.output = 'r';
else
    opts.output = lower(opts.output);
end
if ~isfield(opts, 'type') || ...
   ~ischar(opts.type) || ...
    isempty(opts.type) || ...
   ~any('kor' == lower(opts.type(1)))
    opts.type = 'o';
else
    opts.type = lower(opts.type(1));
end

% permute and reshape m as necessary
m = reshape(permute(m, [opts.dim, setdiff(1:numel(ms), opts.dim)]), [nr, prod(rs)]);

% build design
if opts.nomean
    mv = zeros(nr, 1);
    for mc = 1:nr
        mp = m(mc, :);
        mv(mc) = mean(mp(~isinf(mp) & ~isnan(mp) & (mp ~= 0)));
    end
    mv(isinf(mv) | isnan(mv)) = 0;
else
    mv = zeros(nr, 0);
end
X = [ztrans(r), mv, ones(nr, 1)];
Xc = [1, zeros(1, size(mv, 2)), 0];
df = nr - size(X, 2);

% don't allow inf/nans, and compute stats only where at least 2/3 is ~= 0
in = ~any(isinf(m) | isnan(m), 1) & (sum(m == 0, 1) < (nr / 3));

% prepare output
cm = zeros(size(m, 2), 1);

% switch depending on type
switch (opts.type)

    % Kendall's tau t-test
    case {'k'}

        % if mean removed, yet remove mean first
        if opts.nomean
            z = m(:, in) - mv(:, ones(1, sum(in)));
            z(m(:, in) == 0) = 0;
            m(:, in) = z;
        end

        % then compute tauz
        z = kendtauz(r, m(:, in));

        % and convert to t
        cm(in) = sdist('tinv', sdist('tcdf', z, 1e6), df);

    % OLS correlation
    case {'o'}

        % fit to data
        iXX = pinv(X' * X);
        b = iXX * X' * m(:, in);

        % compute RSS
        rss = sum((m(:, in) - X * b) .^ 2, 1);

        % don't allow perfect fit
        rss(rss < eps) = eps;

        % and t-statistic
        cm(in) = b(1, :) ./ sqrt((iXX(1) / df) .* rss);

    % robust correlation
    case {'r'}

        % robust fit
        [b, w] = fitrobustbisquare_img(X, m(:, in)');

        % and contrast computation
        cm(in) = robustt(X, m(:, in)', b, w, Xc);

end

% r map requested?
if opts.output ~= 't'
    cm = correlinvtstat(cm, df + 2);
end

% reshape
cm = reshape(cm, rs);
