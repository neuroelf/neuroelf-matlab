function t = depcorrt(X, Y, Z, nZ, robflag)
% depcorrt  - dependent correlation comparison t statistic
%
% FORMAT:       t = depcorr(X, Y, Z [, nZ])
%
% Input fields:
%
%       X           vector (or matrix/array) X
%       Y           vector (or matrix/array) Y
%       Z           vector (or matrix/array) Z
%       nZ          optional vector (or matrix/array) nZ
%
% Output fields:
%
%       t           t-statistic for difference of correlations between
%                   rXZ and rYZ, d.f. = N - 3
% or
%                   z-statistic for difference of correlations between
%                   rXY and rZnZ
%
% Note: the correlation is computed along the LAST non-singleton dimension
%
% See also: http://ssc.utexas.edu/consulting/answers/general/gen28.html
%     and:  http://luna.cas.usf.edu/~mbrannic/files/regression/corr1.html

% Version:  v0.9a
% Build:    10101312
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
if nargin < 3 || ...
   ~isnumeric(X) || ...
   ~isnumeric(Y) || ...
   ~isnumeric(Z) || ...
    isempty(X) || ...
   ~isequal(size(X), size(Y)) || ...
   ~isequal(size(X), size(Z))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument (combination) given.' ...
    );
end
robust = false;
v4 = false;
if nargin > 3 && ...
    strcmp(class(nZ), class(X)) && ...
    isequal(size(nZ), size(X))
    v4 = true;
    if nargin > 4 && ...
        ischar(robflag) && ...
       ~isempty(robflag) && ...
        lower(robflag(1)) == 'r'
        robust = true;
    end
elseif nargin > 3 && ...
    ischar(nZ) && ...
   ~isempty(nZ) && ...
    lower(nZ(1)) == 'r'
    robust = true;
end

% translate?
sz = size(X);
if numel(sz) == 2 && ...
    sz(2) == 1
    ts = [1, 1];
    X = X';
    Y = Y';
    Z = Z';
    if v4
        nZ = nZ';
    end
    dm = 1;
    n = sz(1);
else
    ts = sz(1:end-1);
    n = sz(end);
    dm = prod(ts);
    X = reshape(X, [dm, n]);
    Y = reshape(Y, [dm, n]);
    Z = reshape(Z, [dm, n]);
    if v4
        nZ = reshape(nZ, [dm, n]);
    end
end
n2 = floor(0.5 * n);

% three variables
if ~v4

    % compute the three correlations
    if robust
        X = X';
        Y = Y';
        Z = Z';
        msk = ~any(isinf(X) | isnan(X) | isinf(Y) | isnan(Z) | isinf(Z) | isnan(Z));
        msk = msk & ((sum(X ~= 0) > n2) & (sum(Y ~= 0) > n2) & (sum(Z ~= 0) > n2));
        X = ztrans(X);
        Y = ztrans(Y);
        Z = ztrans(Z);
        msk = msk & ~any(isinf(X) | isnan(X) | isinf(Y) | isnan(Y) | isinf(Z) | isnan(Z));
        smsk = sum(msk);
        [b, w, ix, se] = fitrobustbisquare_multi(cat(2, ...
            reshape(X(:, msk), [n, 1, smsk]), ones([n, 1, smsk])), ...
            reshape(Y(:, msk), [n, 1, smsk]));
        rxy = correlinvtstat(lsqueeze(b(1, 1, :) ./ se(1, 1, :)), n);
        [b, w, ix, se] = fitrobustbisquare_multi(cat(2, ...
            reshape(X(:, msk), [n, 1, smsk]), ones([n, 1, smsk])), ...
            reshape(Z(:, msk), [n, 1, smsk]));
        rxz = correlinvtstat(lsqueeze(b(1, 1, :) ./ se(1, 1, :)), n);
        [b, w, ix, se] = fitrobustbisquare_multi(cat(2, ...
            reshape(Y(:, msk), [n, 1, smsk]), ones([n, 1, smsk])), ...
            reshape(Z(:, msk), [n, 1, smsk]));
        ryz = correlinvtstat(lsqueeze(b(1, 1, :) ./ se(1, 1, :)), n);

        % compute determinant
        % |R| =  (1 - rxy ^2 - rxz^2 - ryz^2 + (2*rxy*rxv*rvy)), the determinant of
        % the correlation matrix for X, Y, and V.
        R = 1 + (2 .* rxy .* rxz .* ryz) - (rxy .* rxy + rxz .* rxz + ryz .* ryz);

        % compute dependent t-score
        t = zeros(dm, 1);
        t(msk) = (rxz - ryz) .* sqrt((n-1) .* (1 + rxy)) ./ ...
            (sqrt(2 .* ((n-1) / (n-3)) .* R + ((rxz + ryz) ./2 ) .^2 .* (1 - rxy) .^ 3));
    else
        [rxy{1:2}] = cov_nd(X, Y);
        rxy = rxy{2};
        [rxz{1:2}] = cov_nd(X, Z);
        rxz = rxz{2};
        [ryz{1:2}] = cov_nd(Y, Z);
        ryz = ryz{2};
        R = 1 + (2 .* rxy .* rxz .* ryz) - (rxy .* rxy + rxz .* rxz + ryz .* ryz);
        t = (rxz - ryz) .* sqrt((n-1) .* (1 + rxy)) ./ ...
            (sqrt(2 .* ((n-1) / (n-3)) .* R + ((rxz + ryz) ./2 ) .^2 .* (1 - rxy) .^ 3));
    end

% 4 variables
else

    % compute the required
    if robust
        X = ztrans(X');
        Y = ztrans(Y');
        Z = ztrans(Z');
        nZ = ztrans(nZ');
        msk = ~any(isinf(X) | isnan(X) | isinf(Y) | isnan(Y) | ...
            isinf(Z) | isnan(Z) | isinf(nZ) | isnan(nZ));
        smsk = sum(msk);
        [b, w, ix, se] = fitrobustbisquare_multi(cat(2, ...
            reshape(X(:, msk), [n, 1, smsk]), ones([n, 1, smsk])), ...
            reshape(Y(:, msk), [n, 1, smsk]));
        rxy = correlinvtstat(lsqueeze(b(1, 1, :) ./ se(1, 1, :)), n);
        [b, w, ix, se] = fitrobustbisquare_multi(cat(2, ...
            reshape(X(:, msk), [n, 1, smsk]), ones([n, 1, smsk])), ...
            reshape(Z(:, msk), [n, 1, smsk]));
        rxz = correlinvtstat(lsqueeze(b(1, 1, :) ./ se(1, 1, :)), n);
        [b, w, ix, se] = fitrobustbisquare_multi(cat(2, ...
            reshape(Y(:, msk), [n, 1, smsk]), ones([n, 1, smsk])), ...
            reshape(Z(:, msk), [n, 1, smsk]));
        ryz = correlinvtstat(lsqueeze(b(1, 1, :) ./ se(1, 1, :)), n);
        [b, w, ix, se] = fitrobustbisquare_multi(cat(2, ...
            reshape(X(:, msk), [n, 1, smsk]), ones([n, 1, smsk])), ...
            reshape(nZ(:, msk), [n, 1, smsk]));
        rxn = correlinvtstat(lsqueeze(b(1, 1, :) ./ se(1, 1, :)), n);
        [b, w, ix, se] = fitrobustbisquare_multi(cat(2, ...
            reshape(Y(:, msk), [n, 1, smsk]), ones([n, 1, smsk])), ...
            reshape(nZ(:, msk), [n, 1, smsk]));
        ryn = correlinvtstat(lsqueeze(b(1, 1, :) ./ se(1, 1, :)), n);
        [b, w, ix, se] = fitrobustbisquare_multi(cat(2, ...
            reshape(Z(:, msk), [n, 1, smsk]), ones([n, 1, smsk])), ...
            reshape(nZ(:, msk), [n, 1, smsk]));
        rzn = correlinvtstat(lsqueeze(b(1, 1, :) ./ se(1, 1, :)), n);
        rxyzn = (rxy + rzn) ./ 2;
        hxyzn = .5 .* ( ...
            ((rxz - rxy .* ryz) .* (ryn - ryz .* rzn)) + ...
            ((rxn - rxz .* rzn) .* (ryz - rxy .* rxz)) + ...
            ((rxz - rxn .* rzn) .* (ryn - rxy .* rxn)) + ...
            ((rxn - rxy .* ryn) .* (ryz - ryn .* rzn)));
        sxyzn = hxyzn ./ ((1 - rxyzn .* rxyzn) .^ 2);
        rxy(isnan(rxy)) = 0;
        rzn(isnan(rzn)) = 0;
        zxyzn = fisherr2z(rxy) - fisherr2z(rzn);
        t = zeros(dm, 1);
        t(msk) = sqrt(n - 3) .* zxyzn ./ sqrt(2 - 2 .* sxyzn);
    else
        [rxy{1:2}] = cov_nd(X, Y);
        rxy = rxy{2};
        [rxz{1:2}] = cov_nd(X, Z);
        rxz = rxz{2};
        [ryz{1:2}] = cov_nd(Y, Z);
        ryz = ryz{2};
        [rxn{1:2}] = cov_nd(X, nZ);
        rxn = rxn{2};
        [ryn{1:2}] = cov_nd(Y, nZ);
        ryn = ryn{2};
        [rzn{1:2}] = cov_nd(Z, nZ);
        rzn = rzn{2};
        rxyzn = (rxy + rzn) ./ 2;
        hxyzn = .5 .* ( ...
            ((rxz - rxy .* ryz) .* (ryn - ryz .* rzn)) + ...
            ((rxn - rxz .* rzn) .* (ryz - rxy .* rxz)) + ...
            ((rxz - rxn .* rzn) .* (ryn - rxy .* rxn)) + ...
            ((rxn - rxy .* ryn) .* (ryz - ryn .* rzn)));
        sxyzn = hxyzn ./ ((1 - rxyzn .* rxyzn) .^ 2);
        rxy(isnan(rxy)) = 0;
        rzn(isnan(rzn)) = 0;
        zxyzn = fisherr2z(rxy) - fisherr2z(rzn);
        t = sqrt(n - 3) .* zxyzn ./ sqrt(2 - 2 .* sxyzn);
    end
end

% set invalid values to 0
t(isinf(t) | isnan(t)) = 0;

% reshape to result
if numel(ts) < 2
    ts(2) = 1;
end
t = reshape(t, ts);
