function [ax, st] = robscatter(x, y, opts)
% robscatter  - scatter data with option of robustfit
%
% FORMAT:       [ax, st] = robscatter(x, y [, opts])
%
% Input fields:
%
%       x           Nx1 regressor (model)
%       y           NxD dependent variable (data)
%       opts        optional 1x1 struct with settings
%        .addaxes   boolean flag, add X and Y axes to the graph (false)
%        .ax        axes object, if not given create new figure/axes
%        .cbands    boolean flag, also plot confidence bands (false)
%        .cbandz    z-level for which confidence bands are given (1.96)
%        .color     Dx3 color (scaled 0 .. 1)
%        .cweight   color weighting (0 .. 1, default 0.5)
%        .errors    plot errors between points and regression line
%        .hold      enable hold for axes object, default false
%        .iterlim   iteration limit (default: 30)
%        .lwidth    line width (default, don't set)
%        .marker    Marker type [see set(SCATTER, 'Marker')]
%        .msize     Marker size (default: Matlab default)
%        .ocolor    1x3 outlier color
%        .regrline  plot regression line, default true
%        .xrange    X-axis range
%        .yrange    Y-axis range
%
% Output fields:
%
%       ax          axes object handle
%       st          robust stats output

% Version:  v1.1
% Build:    16041910
% Date:     Apr-19 2016, 10:02 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2015, 2016, Jochen Weber
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
   (~isa(x, 'double') && ...
    ~isa(x, 'single')) || ...
   (~isa(y, 'double') && ...
    ~isa(y, 'single')) || ...
    size(x, 1) ~= size(y, 1) || ...
    isempty(x) || ...
    isempty(y) || ...
    size(x, 1) ~= numel(x) || ...
    ndims(y) ~= 2 || ...
    any(isinf(x) | isnan(x)) || ...
    any(isinf(y(:)) | isnan(y(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument provided.' ...
    );
end
nd = size(y, 2);
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'addaxes') || ...
   ~islogical(opts.addaxes) || ...
    numel(opts.addaxes) ~= 1
    opts.addaxes = false;
end
if ~isfield(opts, 'ax') || ...
   (~isa(opts.ax, 'double') && ...
    ~isa(opts.ax, 'matlab.graphics.axis.Axes')) || ...
    numel(opts.ax) ~= 1 || ...
   ~ishandle(opts.ax) || ...
   ~strcmpi(get(opts.ax, 'Type'), 'axes')
    opts.ax = [];
end
if ~isfield(opts, 'cbands') || ...
   ~islogical(opts.cbands) || ...
    numel(opts.cbands) ~= 1
    opts.cbands = false;
end
if ~isfield(opts, 'cbandz') || ...
   ~isa(opts.cbandz, 'double') || ...
    numel(opts.cbandz) ~= 1 || ...
    isinf(opts.cbandz) || ...
    isnan(opts.cbandz) || ...
    opts.cbandz < 1
    opts.cbandz = sdist('norminv', 0.975, 0, 1);
end
if ~isfield(opts, 'color') || ...
   ~isa(opts.color, 'double') || ...
    ndims(opts.color) ~= 2 || ...
    size(opts.color, 1) ~= nd || ...
    any(isinf(opts.color(:)) | isnan(opts.color(:)) | opts.color(:) < 0 | opts.color(:) > 1)
    opts.color = repmat([0, 0, 1], [nd, 1]);
end
if ~isfield(opts, 'cweight') || ...
   ~isa(opts.cweight, 'double') || ...
    numel(opts.cweight) ~= 1 || ...
    isinf(opts.cweight) || ...
    isnan(opts.cweight) || ...
    opts.cweight < 0 || ...
    opts.cweight > 1
    opts.cweight = 0.5;
end
if ~isfield(opts, 'errors') || ...
   ~islogical(opts.errors) || ...
    numel(opts.errors) ~= 1
    opts.errors = false;
end
if ~isfield(opts, 'hold') || ...
   ~islogical(opts.hold) || ...
    numel(opts.hold) ~= 1
    opts.hold = false;
end
if ~isfield(opts, 'iterlim') || ...
   ~isa(opts.iterlim, 'double') || ...
    numel(opts.iterlim) ~= 1 || ...
    isinf(opts.iterlim) || ...
    isnan(opts.iterlim) || ...
    opts.iterlim < 0 || ...
    opts.iterlim > 400
    opts.iterlim = 30;
else
    opts.iterlim = round(opts.iterlim);
end
if ~isfield(opts, 'lwidth') || ...
   ~isa(opts.lwidth, 'double') || ...
    numel(opts.lwidth) ~= 1 || ...
    isinf(opts.lwidth) || ...
    isnan(opts.lwidth) || ...
    opts.lwidth < 0.2 || ...
    opts.lwidth > 5
    opts.lwidth = [];
end
if ~isfield(opts, 'marker') || ...
   ~ischar(opts.marker) || ...
   ~any(strcmpi(opts.marker(:)', ...
        {'+', 'o', '*', '.', 'x', 'square', 'diamond', ...
         'v', '^', '>', '<', 'pentagram', 'hexagram', 'none'}))
    opts.marker = 'o';
else
    opts.marker = lower(opts.marker(:)');
end
if ~isfield(opts, 'msize') || ...
   ~isa(opts.msize, 'double') || ...
    numel(opts.msize) ~= 1 || ...
    isinf(opts.msize) || ...
    isnan(opts.msize) || ...
    opts.msize < 1 || ...
    opts.msize > 100
    opts.msize = [];
else
    opts.msize = round(opts.msize);
end
if ~isfield(opts, 'ocolor') || ...
   ~isa(opts.ocolor, 'double') || ...
    numel(opts.ocolor) ~= 3 || ...
    any(isinf(opts.ocolor) | isnan(opts.ocolor) | opts.ocolor < 0 | opts.ocolor > 1)
    opts.ocolor = [1, 1, 1];
else
    opts.ocolor = opts.ocolor(:)';
end
if ~isfield(opts, 'regrline') || ...
   ~islogical(opts.regrline) || ...
    numel(opts.regrline) ~= 1
    opts.regrline = true;
end
if ~isfield(opts, 'xrange') || ...
   ~isa(opts.xrange, 'double') || ...
    numel(opts.xrange) ~= 2
    opts.xrange = [-Inf, Inf];
else
    opts.xrange = opts.xrange(:)';
    if isnan(opts.xrange(1))
        opts.xrange(1) = -Inf;
    end
    if isnan(opts.xrange(2))
        opts.xrange(2) = Inf;
    end
end
if ~isfield(opts, 'yrange') || ...
   ~isa(opts.yrange, 'double') || ...
    numel(opts.yrange) ~= 2
    opts.yrange = [-Inf, Inf];
else
    opts.yrange = opts.yrange(:)';
    if isnan(opts.yrange(1))
        opts.yrange(1) = -Inf;
    end
    if isnan(opts.yrange(2))
        opts.yrange(2) = Inf;
    end
end

% make sure to plot and bring to front
ax = opts.ax;
if isempty(ax)
    f = figure;
    figure(f);
    ax = axes;
else
    f = get(ax, 'Parent');
    while ~strcmpi(get(f, 'Type'), 'figure')
        f = get(f, 'Parent');
        if strcmpi(get(f, 'Type'), 'root')
            error( ...
                'neuroelf:UIError', ...
                'Invalid UI object nesting chain detected.' ...
            );
        end
    end
    set(0, 'CurrentFigure', f);
    set(f, 'CurrentAxes', ax);
end
if opts.hold
    hold(ax, 'on');
else
    hold(ax, 'off');
end

% prepare stats
mx = mean(x);
nx = numel(x);
mnx = min(x);
mxx = max(x);
X = [x - mx, ones(nx, 1)];
st = cell(1, nd);

% compute some additional numbers for confidence bands
if opts.cbands
    df = nx - 2;
    dfx = eps + mxx - mnx;
    mnx = mnx - 0.1 * dfx;
    mxx = mxx + 0.1 * dfx;
    bandx = mnx:(0.001 * dfx):mxx;
    bandxm = bandx - mx;
    bandxm2 = bandxm .* bandxm;
    yh = cell(4, nd);
    sxx = sum(X(:, 1) .* X(:, 1));
    cbandt = sdist('tinv', normcdfc(opts.cbandz), df);
end

% perform robust fit
if opts.iterlim > 0
    b = zeros(nd, 2);
    for nc = 1:nd
        [b(nc, :), br, bw, st{nc}] = ...
            fitrobustbisquare(X, y(:, nc), [], opts.iterlim);
        if opts.cbands
            dfe = nx / sum(bw .* bw);
            resi = bw .* (y(:, nc) - X * b(nc, :)');
            yh{4, nc} = (1 / df) .* dfe .* sum(resi .* resi);
        end
    end

% or OLS fit (if iterlim == 0)
else
    b = (pinv(X'*X) * X' * y)';
    for nc = 1:nd
        st{nc} = struct('w', ones(size(y, 1), 1));
        if opts.cbands
            resi = y(:, nc) - X * b(nc, :)';
            yh{4, nc} = (1 / df) .* sum(resi .* resi);
        end
    end
end

% confidence bands
if opts.cbands
    yh{1, nc} = b(nc, 2) + b(nc, 1) .* bandxm;
    yh{2, nc} = cbandt .* sqrt(yh{4, nc} .* (1 / nx + (1 / sxx) .* bandxm2));
    yh{3, nc} = [yh{1, nc} - yh{2, nc}; yh{1, nc} + yh{2, nc}]';
end

% scatter plot(s)
for nc = 1:nd
    s = scatter(x, y(:, nc));
    set(s, 'Marker', opts.marker);
    sc = get(s, 'Children');
    if ~isempty(opts.msize)
        set(sc, 'MarkerSize', opts.msize);
    end
    sc = sc(end:-1:1);
    for scc = 1:numel(sc)
        set(sc(scc), 'MarkerEdgeColor', (1 - opts.cweight) * opts.color(nc, :) + ...
            opts.cweight * st{nc}.w(scc) * opts.color(nc, :) + ...
            opts.cweight * (1 - st{nc}.w(scc)) * opts.ocolor);
    end
    if opts.regrline
        l = line([mnx; mxx], ...
            [b(nc, 2) + (mnx - mx) * b(nc, 1); b(nc, 2) + (mxx - mx) * b(nc, 1)], ...
            'Parent', ax);
        set(l, 'Color', opts.color(nc, :));
        if ~isempty(opts.lwidth)
            set(l, 'LineWidth', opts.lwidth);
        end
    end
    hold(ax, 'on');
    if opts.errors
        el = zeros(size(y, 1), 1);
        for elc = 1:numel(el)
            el(elc) = line([x(elc), x(elc)], ...
                [y(elc, nc), b(nc, 2) + (x(elc) - mx) * b(nc, 1)], 'Parent', ax);
            set(el(elc), 'Color', (1 - opts.cweight) * opts.color(nc, :) + ...
                opts.cweight * st{nc}.w(elc) * opts.color(nc, :) + ...
                opts.cweight * (1 - st{nc}.w(elc)) * opts.ocolor);
        end
        set(el, 'LineStyle', '-.');
        if ~isempty(opts.lwidth)
            set(el, 'LineWidth', opts.lwidth);
        end
    end
    if opts.cbands
        el = plot(bandx, yh{3, nc}(:, 1), 'r--', bandx, yh{3, nc}(:, 2), 'r--');
        set(el, 'Color', opts.color(nc, :));
    end
end

% set range
xl = get(ax, 'XLim');
yl = get(ax, 'YLim');
if ~isinf(opts.xrange(1))
    xl(1) = opts.xrange(1);
end
if ~isinf(opts.xrange(2))
    xl(2) = opts.xrange(2);
end
if ~isinf(opts.yrange(1))
    yl(1) = opts.yrange(1);
end
if ~isinf(opts.yrange(2))
    yl(2) = opts.yrange(2);
end
set(ax, 'XLim', xl, 'YLim', yl);

% add axes ?
if opts.addaxes
    xd = xl(2) - xl(1);
    yd = yl(2) - yl(1);
    xl = line([xl(1) + 0.05 * xd; xl(2) - 0.05 * xd], [0; 0], 'Parent', ax);
    yl = line([0; 0], [yl(1) + 0.05 * yd; yl(2) - 0.05 * yd], 'Parent', ax);
    set(xl, 'Color', [0, 0, 0]);
    set(yl, 'Color', [0, 0, 0]);
    c = get(ax, 'Children');
    set(ax, 'Children', c([end-1:end, 1:end-2]));
end

% continue hold
if ~opts.hold
    hold(ax, 'off');
end
