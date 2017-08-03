function [t, f, w] = regress_coords(c, tc, opts)
% regress_coords  - find best matching quaternion t to match coordinates
%
% FORMAT:       t = regress_coords(c, tc [, opts])
%
% Input fields:
%
%       c           input coordinates (set 1)
%       tc          transformed coordinates (set 2)
%       opts        optional settings
%        .robust    perform robust fit, default: false
%
% Output fields:
%
%       t       transformation matrix so that tc = t * c;
%
% Note: naturally, at least 4 coordinates must be in each set; for the
%       robust regression, at least 6 coordinates must be given

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
   ~isnumeric(c) || ...
   ~isnumeric(tc) || ...
    ndims(c) > 2 || ...
    ndims(tc) > 2 || ...
   ~any(size(c) == 3) || ...
   ~any(size(tc) == 3) || ...
   ~all(sort(size(c)) == sort(size(tc))) || ...
    max(size(c)) < 4
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isa(c, 'double')
    c = double(c(:, :));
end
if ~isa(tc, 'double')
    tc = double(tc(:, :));
end
if size(c, 1) == 3 && ...
    size(c, 2) ~= 3
    c = c';
end
if size(tc, 1) == 2 && ...
    size(tc, 2) ~= 3
    tc = tc';
end
tc = tc(:);
if ~isfield(opts, 'noresize') || ...
   ~islogical(opts.noresize) || ...
    numel(opts.noresize) ~= 1
    opts.noresize = false;
end
if ~isfield(opts, 'noshear') || ...
   ~islogical(opts.noshear) || ...
    numel(opts.noshear) ~= 1
    opts.noshear = false;
end
if ~isfield(opts, 'robust') || ...
   ~islogical(opts.robust) || ...
    numel(opts.robust) ~= 1 || ...
    size(c, 1) < 6
    opts.robust = false;
end

% prepare regression matrix
o = ones(size(c, 1), 1);
z = zeros(size(c, 1), 1);
X = [[c(:, 1), z, z, c(:, 2), z, z, c(:, 3), z, z, o, z, z]; ...
     [z, c(:, 1), z, z, c(:, 2), z, z, c(:, 3), z, z, o, z]; ...
     [z, z, c(:, 1), z, z, c(:, 2), z, z, c(:, 3), z, z, o]];

% add no-resize contraints ? that means that
% b1 * b1 + b2 * b2 + b3 * b3 := 1
% b4 * b4 + b5 * b5 + b6 * b6 := 1
% b7 * b7 + b8 * b8 + b9 * b9 := 1
% but not implemented yet...

% add no-shearing contraints ? that means that
% b1 * b4 + b2 * b5 + b3 * b6 := 0
% b1 * b7 + b2 * b8 + b3 * b9 := 0
% b4 * b7 + b5 * b8 + b6 * b9 := 0
% but not implemented yet...

% normal or robust regression
if opts.robust

    % find the least squares solution.
    [n, p] = size(X);
    op = ones(1, p);
    [Q, h, perm] = qr(X, 0);
    tol = abs(h(1)) * n * eps(class(h));
    xrank = sum(abs(diag(h)) > tol);
    if xrank ~= p
        error( ...
            'neuroelf:RankDeficient', ...
            'X is rank deficient, rank = %d', ...
            xrank ...
        );
    end

    % fill beta estimates
    X = X(:, perm);
    [pp, perm] = sort(perm);
    Q = sqrt(eps(class(X)));
    h = X / h;
    h = min(.9999, sum(h .* h, 2));
    h = 1 ./ sqrt(1 - h);
    strect = struct('RECT', true);

    % adjust residuals using leverage, as advised by DuMouchel & O'Brien
    b = linsolve(X, tc, strect);
    b0 = zsz(b);
    wxrank = xrank;

    % handle case of "almost perfect fit"
    tiny_s = 1e-6 * std(tc);
    if tiny_s == 0
        tiny_s = Q;
    end

    % perform iteratively reweighted least squares to get coefficient estimates
    iter = 1;
    while any(abs(b - b0) > Q * max(abs(b), abs(b0))) || iter == 1

        % break if too many iterations
        if (iter > 250)
            break;
        end

        % compute residuals from previous fit, then compute scale estimate
        r = tc - X * b;
        radj = r .* h;
        s = madsigma(radj, wxrank);

        % compute new weights from these residuals, then re-fit
        w = radj / (max(s, tiny_s) * 4.685);
        w = (abs(w) < 1) .* (1 - w .^ 2) .^ 2;
        w = repmat(min(reshape(w, size(c)), [], 2), 3, 1);
        b0 = b;
        sw = sqrt(w);
        [b, wxrank] = linsolve(X .* sw(:, op), tc .* sw, strect);
        iter = iter + 1;
    end

    % put betas in correct positions in t
    t = b(perm);
    w = w(1:size(c, 1));

    % get tc (hat) to compute estimate of variance and explained variance
    tc = sw .* tc;
    tch = tc - sw .* (X * t);
    t = [reshape(t, 3, 4); [0, 0, 0, 1]];
    tv = sum(var(reshape(tc, size(c))));
    rv = sum(var(reshape(tch, size(c))));

% regular least squares fit
else
    t = pinv(X' * X) * X' * tc(:);

    % get tc (hat) to compute estimate of variance and explained variance
    tch = tc(:) - X * t;
    t = [reshape(t, 3, 4); [0, 0, 0, 1]];
    tv = sum(var(reshape(tc, size(c))));
    rv = sum(var(reshape(tch, size(c))));

    % fill w if necessary
    if nargout > 2
        w = ones(size(c, 1), 1);
    end
end

% compute fit if requested
f = 1 - rv / tv;



function s = madsigma(r, p)
%MADSIGMA    Compute sigma estimate using MAD of residuals from 0
rs = sort(abs(r));
ns = numel(rs) + p;
if mod(ns, 2) == 0
    s = rs(ns / 2) / 0.6745;
else
    s = (rs(ns / 2 - 0.5) + rs(ns / 2 + 0.5)) / 1.349;
end
