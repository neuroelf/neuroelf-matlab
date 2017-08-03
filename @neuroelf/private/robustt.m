function t = robustt(X, y, b, w, c)
% robustt  - compute t scores given robust stats
%
% FORMAT:       t = robustt(X, y, b, w, c)
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
%       t           ...xC t-contrasts

% Version:  v1.1
% Build:    16061514
% Date:     Jun-15 2016, 2:01 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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
   (ndims(b) ~= ndims(y) && ...
    ndims(b) ~= (ndims(y) - 1)) || ...
   (numel(b) ~= size(X, 2) && ...
    size(b, ndims(y)) ~= size(X, 2)) || ...
    any(isinf(b(:)) | isnan(b(:))) || ...
   ~isa(w, 'double') || ...
    ndims(w) ~= ndims(y) || ...
    any(size(w) ~= size(y)) || ...
    any(isinf(w(:)) | isnan(w(:)) | w(:) < 0 | w(:) > 1) || ...
   ~isa(c, 'double') || ...
    ndims(c) ~= 2 || ...
    size(c, 2) ~= size(X, 2) || ...
    any(isinf(c(:)) | isnan(c(:)))
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

% compute (raw) residual
r = y - b * X';

% compute sqrt of weights
sw = sqrt(w);

% and sum (to get idea of real d.f.)
ws = sum(w, 2);

% compute weighted residual, penalizing the loss of d.f.
wr = sw .* r .* (n ./ ws(:, ones(1, n)));
rs = sum(wr .* wr, 2);

% compute t-stat se-factor
sf = sqrt(rs ./ (ws - p));

% create output data
t = zeros(ny, nc);

% test xprogress
if ny > 2000
    try
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 200, 640, 36]);
        xprogress(pbar, 'settitle', sprintf('Computing %d contrasts...', nc));
        xprogress(pbar, 0, sprintf('Running %d samples...', ny), 'visible', 0, ny);
        vcs = ceil(max(500, ny / 100));
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

% only use elements for which statistics are ~= 0
nb = all(b == 0, 2);

% single contrast
if nc == 1

    % iterate over data elements (voxels)
    for vc = 1:ny

        % move on if nothing to do
        if nb(vc)
            continue;
        end

        % create weighted matrix
        wX = sw(vc(op), :)' .* X;

        % compute
        t(vc) = b(vc, :) * c' / sqrt(c * pinv(wX' * wX) * c');

        % update progress bar
        if vc >= vcn && ...
           ~isempty(pbar)
            xprogress(pbar, vc);
            vcn = vcn + vcs;
        end
    end

% multi contrast
else

    % iterate over data elements (voxels)
    for vc = 1:ny

        % move on if nothing to do
        if nb(vc)
            continue;
        end

        % create weighted matrix
        wX = sw(vc(op), :)' .* X;

        % invert weighted matrix
        iX = pinv(wX' * wX);

        % iterate over contrasts
        for cc = 1:nc
            t(vc, cc) = b(vc, :) * c(cc, :)' / sqrt(c(cc, :) * iX * c(cc,:)');
        end

        % update progress bar
        if vc >= vcn && ...
           ~isempty(pbar)
            xprogress(pbar, vc);
            vcn = vcn + vcs;
        end
    end
end

% convert to nominal t
t = t ./ sf(:, ones(1, nc));
t(isinf(t) | isnan(t)) = 0;
if (n - p) < 2000
    ts = -sign(t);
    t = reshape(ts .* sdist('tinv', sdist('tcdf', -abs(t), ws(:, ones(1, nc)) - p), n - p), [sy, nc]);
    t(isinf(t) | isnan(t)) = 0;
end

% clear progress bar
if ~isempty(pbar)
    closebar(pbar);
end
