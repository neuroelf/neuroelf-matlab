function [bi, wi, ixxi, sei] = fitrobustbisquare_multi(Xi, yi, tune, olt, opts)
% fitrobustbisquare_img  - faster fitting function for robust bisquare
%
% FORMAT:  [bi, wi, ixxi, sei] = fitrobustbisquare_img(Xi, yi [, tune [, olt]])
%
% Input fields:
%
%       Xi          models (will NOT be patched, so must contain confounds)
%       yi          data
%       tune        bisquare tuning parameter (default: 4.685)
%       olt         outlier threshold (value from 0...1, default sqrt(N))
%
% Output fields:
%
%       bi          beta weights
%       wi          weighting vector
%       ixxi        inverse matrices
%       sei         standard error
%
% Note: other than the robustfit function of the stats toolbox, NaNs
%       must be removed *prior* to regression/estimation, manually
%       also, if the data contain more than olt zeros (per case), those
%       will be removed first and the weights then set to zero after
%       regression

% Version:  v0.9c
% Build:    13022210
% Date:     Feb-22 2013, 10:46 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2013, Jochen Weber
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
   ~isa(Xi, 'double') || ...
   ~isnumeric(yi) || ...
    isempty(Xi) || ...
    size(Xi, 2) > size(Xi, 1) || ...
    any(isinf(Xi(:)) | isnan(Xi(:))) || ...
    size(yi, 2) ~= 1 || ...
    any(isinf(yi(:)) | isnan(yi(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
xs = size(Xi);
ys = size(yi);
nx = prod(xs(3:end));
ny = prod(ys(3:end));
if ~any(ny == [1, nx])
    error( ...
        'neuroelf:BadArgument', ...
        'Xi and yi must match in upper dims.' ...
    );
end
n = xs(1);
p = xs(2);
op = ones(1, p);
if ~isa(yi, 'double')
    yi = double(yi);
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
   ~isa(olt, 'double') || ...
    numel(olt) ~= 1 || ...
    isinf(olt) || ...
    isnan(olt) || ...
    olt <= 0
    olt = ceil(n - sqrt(n));
elseif olt > 1
    olt = n;
else
    olt = n * olt;
end
if nargin < 5 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'pbar') || ...
    numel(opts.pbar) ~= 1 || ...
   ~any(strcmpi(class(opts.pbar), {'xfigure', 'xprogress'}))
    opts.pbar = [];
end
if ~isfield(opts, 'prange') || ...
   ~isa(opts.prange, 'double') || ...
    numel(opts.prange) ~= 2 || ...
    any(isinf(opts.prange) | isnan(opts.prange) | opts.prange < 0 | opts.prange > 1)
    opts.prange = [0, 1];
else
    opts.prange = sort(opts.prange(:))';
    if opts.prange(2) == opts.prange(1)
        opts.prange(2) = opts.prange(1) + sqrt(eps);
    end
end

% find the least squares solution
[ixxi, bi] = mmregress(Xi, yi);
b0 = zeros(size(ixxi, 1), 1);
wxrank = size(ixxi, 1);

% pre-allocation wi and sei and get significance threshold
if ny > 1
    wi = ones(size(yi));
else
    wi = ones([n, 1, xs(3:end)]);
end
if nargout > 3
    sei = ones(size(bi));
    wsi = ones([1, 1, xs(3:end)]);
end
Q = sqrt(eps(class(Xi)));
strect = struct('RECT', true);

% test xprogress
pbar = [];
vcn = Inf;
pbmin = opts.prange(1);
pbdif = opts.prange(2) - pbmin;
if nx > 500
    vcs = round(max(100, nx / 100));
    vcn = vcs;
    if isempty(opts.pbar)
        try
            pbar = xprogress;
            xprogress(pbar, 'setposition', [80, 200, 640, 36]);
            xprogress(pbar, 'settitle', sprintf('Robustly fitting %d-x-%d matrix...', n, p));
            xprogress(pbar, pbmin, sprintf('Estimating %d samples...',  ...
                nx), 'visible', pbmin, opts.prange(2));
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

% fix values
if ny == 1
    y = yi;
    y0 = y ~= 0;
    y0s = sum(y0);
    if y0s < n
        error( ...
            'neuroelf:BadArgument', ...
            'Data must not contain zeros.' ...
        );
    end
end

% iterate over voxels
pbdif = pbdif / nx;
for vc = 1:nx

    % get data
    if ny > 1
        y = yi(:, 1, vc);
        y0 = y ~= 0;
        y0s = sum(y0);
    end

    % enough values to use full data (hopefully rejecting outliers)
    if y0s >= olt

        % get initial model and estimate
        X = Xi(:, :, vc);
        w = wi(:, 1, vc);
        b = bi(:, 1, vc);
        b0(:) = 0;

        % compute residual adjustment factor
        [qrr, h, perm] = qr(X, 0);
        h = X(:, perm) / h;
        h = min(.9999, sum(h .* h, 2));
        h = 1 ./ sqrt(1 - h);

        % handle case of "almost perfect fit"
        tiny_s = 1e-6 * sqrt(varc(y));
        if tiny_s == 0
            tiny_s = Q;
        end

        % perform iteratively reweighted least squares to get coefficient estimates
        iter = 1;
        while any(abs(b - b0) > Q * max(abs(b), abs(b0))) || iter == 1

            % break if too many iterations
            if (iter > 50)
                break;
            end

            % compute residuals from previous fit, then compute scale estimate
            r = y - X * b;
            radj = r .* h;
            s = madsigma(radj, wxrank);

            % compute new weights from these residuals, then re-fit
            w = radj / (max(s, tiny_s) * tune);
            w = (abs(w) < 1) .* (1 - w .* w) .^ 2;
            b0 = b;
            sw = sqrt(w);
            [b, wxrank] = linsolve(X .* sw(:, op), y .* sw, strect);
            iter = iter + 1;
        end

    % too few values
    elseif y0s <= p

        % continue immediately
        continue;

    % not enough data
    else

        % get reduced model
        Xr = Xi(y0, :);
        Xrc = any(Xr ~= 0, 1);
        b = 0 .* b0;
        r = 0 .* y;
        w = double(y0);

        % use external function on partial data
        try
            [bex, rex, wex] = fitrobustbisquare(Xr(:, Xrc), y(y0), tune);

        % yet error, bad design matrix -> unestimable
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            continue;
        end

        % store in actually required values
        b(Xrc) = bex;
        r(y0) = rex;
        w(y0) = wex;
    end

    % store values in output
    wi(:, 1, vc) = w;
    bi(:, 1, vc) = b;
    ixxi(:, :, vc) = inv(transmul(X .* sw(:, op)));

    % compute std-errors?
    if nargout > 3

        % compute some factors
        sw = sqrt(w);
        ws = sum(w);
        wn = n / ws;
        wr = wn .* (sw .* r);
        rs = sum(wr .* wr);
        sf = 1 / sqrt(rs ./ (ws - p));

        % compute t-score and store correction factor
        sei(:, 1, vc) = sf .* (b ./ sqrt(diag(ixxi(:, :, vc))));
        wsi(vc) = ws;
    end

    % update progress bar
    if vc >= vcn && ...
       ~isempty(pbar)
        pbar.Progress(pbmin + pbdif * vc);
        vcn = vcn + vcs;
    end
end

% compute actual se
sei(isinf(sei) | isnan(sei)) = 0;
sei = bi ./ (-sign(sei) .* sdist('tinv', ...
    sdist('tcdf', -abs(sei), repmat(wsi - p, size(ixxi, 1), 1)), n - p));
sei(isinf(sei) | isnan(sei)) = 0;

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
