function [bi, wi, ixxi, sei] = fitrobustbisquare_img(X, yi, tune, olt, opts, hp, iXX, wi)
% fitrobustbisquare_img  - faster fitting function for robust bisquare
%
% FORMAT:  [bi, wi] = fitrobustbisquare_img(X, yi [, tune [, olt]])
%
% Input fields:
%
%       X           model (will NOT be patched, so must contain confounds)
%       yi          data
%       tune        bisquare tuning parameter (default: 4.685)
%       olt         outlier threshold (value from 0...1, default sqrt(N))
%
% Output fields:
%
%       bi          beta weights
%       wi          weighting vector
%
% Note: other than the robustfit function of the stats toolbox, NaNs
%       must be removed *prior* to regression/estimation, manually
%       also, if the data contain more than olt zeros (per case), those
%       will be removed first and the weights then set to zero after
%       regression

% Version:  v1.1
% Build:    16061510
% Date:     Jun-15 2016, 10:26 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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
if nargin < 2 || ~isa(X, 'double') || ~isnumeric(yi) || isempty(X) || ...
    ndims(X) > 2 || size(X, 2) > size(X, 1) || any(isinf(X(:)) | isnan(X(:))) || ...
   (size(yi, ndims(yi)) ~= size(X, 1) && size(yi, 1) ~= size(X, 1)) || ...
    any(isinf(yi(:)) | isnan(yi(:)))
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
sqeps = sqrt(eps);
if ~isa(yi, 'double')
    yi = double(yi);
end
if nargin < 3 || ~isa(tune, 'double') || numel(tune) ~= 1 || ...
    isinf(tune) || isnan(tune) || tune <= 0
    tune = 4.685;
end
[n, p] = size(X);
op = ones(1, p);
if nargin < 4 || ~isa(olt, 'double') || numel(olt) ~= 1 || ...
    isinf(olt) || isnan(olt) || olt <= 0
    olt = ceil(n - (n .^ 0.75));
elseif olt > 1
    olt = n;
else
    olt = n * olt;
end
if nargin < 5 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dim') || ~isa(opts.dim, 'double') || numel(opts.dim) ~= 1 || ...
    isinf(opts.dim) || isnan(opts.dim) || opts.dim < 1 || opts.dim ~= fix(opts.dim) || ...
    opts.dim > ndims(yi) || size(yi, opts.dim) ~= size(X, 1)
    if size(yi, 1) == size(X, 1)
        opts.dim = 1;
    else
        opts.dim = ndims(yi);
    end
end
yd1 = (opts.dim == 1);
if ~isfield(opts, 'maxiter') || ~isa(opts.maxiter, 'double') || numel(opts.maxiter) ~= 1 || ...
    isinf(opts.maxiter) || isnan(opts.maxiter) || opts.maxiter < 2
    maxiter = 50;
else
    maxiter = min(100, round(opts.maxiter));
end
if ~isfield(opts, 'pbar') || numel(opts.pbar) ~= 1 || ...
   ~any(strcmpi(class(opts.pbar), {'xfigure', 'xprogress'}))
    opts.pbar = [];
end
if ~isfield(opts, 'prange') || ~isa(opts.prange, 'double') || numel(opts.prange) ~= 2 || ...
    any(isinf(opts.prange) | isnan(opts.prange) | opts.prange < 0 | opts.prange > 1)
    opts.prange = [0, 1];
else
    opts.prange = sort(opts.prange(:))';
    if opts.prange(2) == opts.prange(1)
        opts.prange(2) = opts.prange(1) + sqeps;
    end
end

% find the least squares solution.
if (sum(X(:) ~= 0) / numel(X)) <= 0.25 && ~issparse(X)
    X = sparse(X);
end
if nargin < 6 || ~iscell(hp) || numel(hp) ~= 2 || ~isa(hp{1}, 'double') || ...
   ~isequal(size(hp{1}), [p, p]) || any(isinf(hp{1}(:)) | isnan(hp{1}(:))) || ...
   ~isa(hp{2}, 'double') || ~isequal(size(hp{2}), [1, p]) || ~isequal(sort(hp{2}(:)), (1:p)')
    [Q, h, perm] = qr(X, 0);
else
    h = hp{1};
    perm = hp{2};
end
if nargin < 7 || ~isa(iXX, 'double') || ~isequal(size(iXX), [p, p]) || ...
    any(isinf(iXX(:)) | isnan(iXX(:)))
    iXX = [];
else
    iXX = iXX(perm, perm);
end
tol = abs(h(1)) * n * eps(class(h));
xrank = sum(abs(diag(h)) > tol);
if xrank ~= p
    error('neuroelf:general:rankDeficient', 'X is rank deficient, rank = %d', xrank);
end

% preparation values
X = X(:, perm);
[pp, perm] = sort(perm);
Q = sqrt(eps(class(X)));
h = X / h;
h = min(.9999, sum(h .* h, 2));
h = 1 ./ sqrt(1 - h);

% input/output size
sy = size(yi);
if yd1
    sy(1) = [];
else
    sy(end) = [];
end
ny = prod(sy);

% bring data into better format
if yd1
    yi = reshape(yi, n, ny);
else
    yi = reshape(yi, ny, n);
end

% create output "images"
bi = zeros(ny, p);
if nargin < 8 || ~any(numel(wi) == [n, n * ny]) || (~isa(wi, 'double') && ~isa(wi, 'single'))
    if yd1
        wi = ones(n, ny);
    else
        wi = ones(ny, n);
    end
    iw = false;
else
    wi = limitrangec(double(wi), 0, 1, 1);
    if yd1
        if numel(wi) == n
            wi = repmat(wi(:), 1, ny);
        else
            wi = reshape(wi, [n, ny]);
        end
    else
        if numel(wi) == n
            wi = repmat(wi(:)', ny, 1);
        else
            wi = reshape(wi, [ny, n]);
        end
    end
    iw = true;
end
dix = false;
dis = false;
if nargout > 2
    dix = true;
    ixxi = zeros([p, p, ny]);
    if nargout > 3
        dis = true;
        sei = ones(size(bi));
        wsi = ones([1, 1, ny]);
    end
end
b = zeros(p, 1);
b0 = b;

% test xprogress
pbar = [];
vcn = Inf;
pbmin = opts.prange(1);
pbdif = opts.prange(2) - pbmin;
if ny > 500
    vcs = round(max(100, ny / 100));
    vcn = vcs;
    if isempty(opts.pbar)
        try
            pbar = xprogress;
            xprogress(pbar, 'setposition', [80, 200, 640, 36]);
            xprogress(pbar, 'settitle', sprintf('Robustly fitting %d-x-%d matrix...', n, p));
            xprogress(pbar, pbmin, sprintf('Estimating %d samples...',  ...
                ny), 'visible', pbmin, opts.prange(2));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            pbar = [];
            vcn = Inf;
        end
    else
        pbar = opts.pbar;
        pbar.Progress(pbmin, pbmin);
    end
end

% iterate over voxels
pbdif = pbdif / ny;
for vc = 1:ny

    % fill beta estimates
    if yd1
        y = yi(:, vc);
    else
        y = yi(vc, :)';
    end

    % how many values
    y0 = y ~= 0;
    y0s = sum(y0);

    % too few values
    if y0s <= p

        % continue immediately
        continue;

    % enough values to use full data (hopefully rejecting outliers)
    elseif y0s >= olt

        % first figure out if X is sparse and get indices
        if issparse(X)
            xsp = true;
            [xsi, xsj, xsv] = find(X);
        else
            xsp = false;
        end

        % compute iXX if necessary
        if isempty(iXX)
            iXX = inv(X' * X);
            ixxt = (iXX ~= 0 & abs(iXX) < sqeps);
            iXX(ixxt) = 0;
            if sum(iXX(:) ~= 0) < (0.25 * numel(iXX))
                iXX = sparse(iXX);
            end
        end

        % compute intial beta estimate with initial weights
        if iw

            % get weights
            if yd1
                w = wi(:, vc);
            else
                w = wi(vc, :)';
            end

            % follow path (see comments below)
            sw = sqrt(w);
            if xsp
                Xsw = sparse(xsi, xsj, xsv .* sw(xsi), n, p);
            else
                Xsw = X .* sw(:, op);
            end
            iXXw = inv(Xsw' * Xsw);
            ixxt = (iXXw ~= 0 & abs(iXXw) < sqeps);
            iXXw(ixxt) = 0;
            if issparse(iXX)
                iXXw = sparse(iXXw);
            end
            b = iXXw * (Xsw' * (y .* sw));

        % un-biased (unweighted) estimate
        else
            b = iXX * (X' * y);
        end

        % and set test beta to 0s
        b0(:) = 0;

        % handle cases of "almost perfect fit"
        tiny_s = 1e-6 * sqrt(varc(y));
        if tiny_s == 0
            tiny_s = Q;
        end

        % initialize interation counter
        iter = 1;

        % loop until solution is found
        while any(abs(b - b0) > Q * max(abs(b), abs(b0))) || iter == 1

            % break if iterations are over
            if (iter > maxiter)
                break;
            end

            % compute residual
            r = y - X * b;

            % and adjust by factor
            radj = h .* r;

            % then compute sigma
            s = madsigma(radj, xrank);

            % and weights
            w = radj / (max(s, tiny_s) * tune);
            w = (abs(w) < 1) .* (1 - w .* w) .^ 2;

            % set test beta
            b0 = b;

            % and compute reweighted design matrix
            sw = sqrt(w);
            if xsp
                Xsw = sparse(xsi, xsj, xsv .* sw(xsi), n, p);
            else
                Xsw = X .* sw(:, op);
            end

            % follow same logic
            iXXw = inv(Xsw' * Xsw);
            ixxt = (iXXw ~= 0 & abs(iXXw) < sqeps);
            iXXw(ixxt) = 0;
            if issparse(iXX)
                iXXw = sparse(iXXw);
            end

            % then re-compute betas
            b = iXXw * (Xsw' * (y .* sw));

            % and increase iterations counter
            iter = iter + 1;
        end

        % store values
        if yd1
            wi(:, vc) = w;
        else
            wi(vc, :) = w';
        end
        bi(vc, :) = b;
        if dix
            ixxi(:, :, vc) = iXXw;

            % compute std-errors?
            if dis

                % compute some factors
                ws = sum(w);
                wn = n / ws;
                wr = wn .* (sw .* r);
                rs = sum(wr .* wr);
                sf = 1 / sqrt(rs ./ (ws - p));

                % compute t-score and store correction factor
                sei(vc, :) = sf .* (b ./ sqrt(diag(iXXw)))';
                wsi(vc) = ws;
            end
        end

    % not enough data
    else

        % get reduced model
        Xr = X(y0, :);
        Xrc = any(Xr ~= 0, 1);

        % use function on partial data
        try
            [bex, wex] = fitrobustbisquare_img(Xr(:, Xrc), y(y0)', tune, olt, opts);

        % yet error, bad design matrix -> unestimable
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            continue;
        end

        % store values into good portion
        if yd1
            wi(y0, vc) = wex(:);
            wi(~y0, vc) = 0;
        else
            wi(vc, y0) = wex(:)';
            wi(vc, ~y0) = 0;
        end
        bi(vc, Xrc) = bex(:)';
        if dix
            sw = sqrt(wex(:));
            Xsw = Xr(:, Xrc) .* repmat(sw, 1, sum(Xrc));
            ixxi(Xrc, Xrc, vc) = inv(Xsw' * Xsw);
            if dis
                r = y(y0) - Xr(:, Xrc) * bex(:);
                ws = sum(wex);
                wn = n / ws;
                wr = wn .* (sw .* r);
                rs = sum(wr .* wr);
                sf = 1 / sqrt(rs ./ (ws - p));
                sei(vc, Xrc) = sf .* (bex(:) ./ sqrt(diag(ixxi(Xrc, Xrc, vc))))';
                wsi(vc) = ws;
            end
        end
    end

    % update progress bar
    if vc >= vcn && ...
       ~isempty(pbar)
        pbar.Progress(pbmin + pbdif * vc);
        vcn = vcn + vcs;
    end

end

% reshape outputs
bi = reshape(bi(:, perm), [sy, p]);
if yd1
    wi = reshape(wi, [n, sy]);
else
    wi = reshape(wi, [sy, n]);
end
if dis
    sei = reshape(sei(:, perm), [sy, p]);
end

% clear progress bar
if ~isempty(pbar) && ...
    isempty(opts.pbar)
    closebar(pbar);
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
