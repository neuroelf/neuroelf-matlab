function F = robustF(X, y, b, w, c)
% robustF  - compute F scores given robust stats
%
% FORMAT:       F = robustF(X, y, b, w, c)
%
% Input fields:
%
%       X           TxP model
%       y           ...xT data (n-dim supported)
%       b           ...xP beta estimates
%       w           ...xT weights (e.g. from fitrobustbisquare_img)
%       c           contrast vectors CxP for C contrasts with P weights
%
% Output fields:
%
%       F           ...xC F-contrasts
%
% Note: all values in c that are non-zero will be examined as set!

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-08 2011, 9:16 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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
if nargin < 5 || ...
   ~isa(X, 'double') || ...
    ndims(X) ~= 2 || ...
    size(X, 1) <= size(X, 2) || ...
    any(isinf(X(:))) || ...
    any(isnan(X(:))) || ...
   ~isa(y, 'double') || ...
   (numel(y) ~= size(X, 1) && ...
    size(y, ndims(y)) ~= size(X, 1)) || ...
    any(isinf(y(:)) | isnan(y(:))) || ...
   ~isa(b, 'double') || ...
    ndims(b) ~= ndims(y) || ...
   (numel(b) ~= size(X, 2) && ...
    size(b, ndims(b)) ~= size(X, 2)) || ...
    any(isinf(b(:)) | isnan(b(:))) || ...
   ~isa(w, 'double') || ...
    ndims(w) ~= ndims(y) || ...
    any(size(w) ~= size(y)) || ...
    any(isinf(w(:)) | isnan(w(:)) | w(:) < 0 | w(:) > 1) || ...
   ~isa(c, 'double') || ...
    ndims(c) ~= 2 || ...
    size(c, 2) ~= size(X, 2) || ...
    any(isnan(c(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% treat single-stats case
if ndims(y) == 2 && ...
    ndims(w) == 2 && ...
    size(y, 2) == 1 && ...
    size(w, 2) == 1
    y = y';
    w = w';
end

% reshape arguments
[n, p] = size(X);
nc = size(c, 1);
sy = size(y);
sy(end) = [];
ny = prod(sy);
y = reshape(y, ny, n);
b = reshape(b, ny, p);
w = reshape(w, ny, n);
op = ones(1, p);

% simply get non-zero elements in c
ci = cell(1, nc);
for cc = 1:nc
    ci{nc} = find(c(cc, :) ~= 0);
end

% compute (raw) residual
r = y - b * X';

% compute sqrt of weights
sw = sqrt(w);

% and sum (to get idea of real d.f.)
ws = sum(w .* sw, 2);

% compute weighted residual
wr = sw .* r;
rs = sum(wr .* wr, 2);

% compute penalizing factor
pF = (1 / (n - p)) * (ws - p);

% create output data
F = zeros(ny, nc);

% test xprogress
if ny > 5000
    try
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 200, 640, 36]);
        xprogress(pbar, 'settitle', sprintf('Computing %d contrasts...', nc));
        xprogress(pbar, 0, sprintf('Running %d samples...', ny), 'visible', 0, ny);
        vcs = ceil(max(2000, ny / 100));
        vcn = vcs;
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        pbar = [];
        vcn = Inf;
    end
else
    pbar = [];
    vcn = Inf;
end

% iterate over data elements (voxels)
for vc = 1:ny

    % create weighted matrix
    wX = sw(vc(op), :)' .* X;

    % iterate over contrasts
    for cc = 1:nc

        % compute variance explained by set elements in c
        xSS = sum((wX(:, ci{cc}) * b(vc, ci{cc})') .^ 2);

        % compute F test (without d.f. terms)
        F(vc, cc) = xSS / rs(vc);
    end

    % update progress bar
    if vc >= vcn && ...
       ~isempty(pbar)
        xprogress(pbar, vc);
        vcn = vcn + vcs;
    end
end

% convert to nominal F
F(isinf(F) | isnan(F)) = 0;
for cc = 1:nc
    F(:, cc) = sdist('finv', sdist('fcdf', ...
        ((1 / numel(ci{cc})) * (ws - p)) .* (pF .* F(:, cc)), numel(ci{cc}), ws - p, true), ...
        numel(ci{cc}), n - p, true);
end
F = reshape(F, [sy, nc]);
F(isinf(F) | isnan(F)) = 0;

% clear progress bar
if ~isempty(pbar)
    closebar(pbar);
end
