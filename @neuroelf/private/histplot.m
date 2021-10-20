function p = histplot(x, opts)
if nargin < 1 || ~isa(x, 'double')
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
x = x(:);
x(isinf(x) | isnan(x)) = [];
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'ax') || numel(opts.ax) ~= 1 || (~ishandle(opts.ax) && ~isa(opts.ax, 'double'))
    opts.ax = gca;
end
p = opts.ax;
if isempty(x)
    return;
end
if ~isfield(opts, 'bins') || numel(opts.bins) ~= 1 || ~isa(opts.bins, 'double') || ...
    isinf(opts.bins) || isnan(opts.bins) || opts.bins < 1
    opts.bins = 10;
else
    opts.bins = round(opts.bins);
end
if ~isfield(opts, 'color') || isempty(opts.color)
    opts.color = [0.25, 0.25, 0.25];
end
if ~isfield(opts, 'density') || ~islogical(opts.density) || numel(opts.density) ~= 1
    opts.density = false;
end
if ~isfield(opts, 'max') || numel(opts.max) ~= 1 || ~isa(opts.max, 'double') || ...
    isinf(opts.bins) || isnan(opts.bins) || opts.bins < 1
    opts.max = max(x);
end
if ~isfield(opts, 'min') || numel(opts.min) ~= 1 || ~isa(opts.min, 'double')
    opts.min = min(x);
end
if opts.max <= opts.min
    opts.max = opts.min + 1;
end
if ~isfield(opts, 'medline') || numel(opts.medline) ~= 1 || ~isa(opts.medline, 'double') || ...
    isinf(opts.medline) || isnan(opts.medline) || opts.medline < 1
    opts.medline = [];
end
if ~isfield(opts, 'medlinec') || isempty(opts.medlinec)
    opts.medlinec = [0.5, 0.5, 0.5];
end
if ~isfield(opts, 'sline') || numel(opts.sline) ~= 1 || ~isa(opts.sline, 'double') || ...
    isinf(opts.sline) || isnan(opts.sline) || opts.sline < 1
    opts.sline = [];
end
if ~isfield(opts, 'slinew') || numel(opts.slinew) ~= 1 || ~isa(opts.slinew, 'double') || ...
    isinf(opts.slinew) || isnan(opts.slinew) || opts.slinew <= 0
    opts.slinew = 2.5;
end
if ~isfield(opts, 'xbins') || ~islogical(opts.xbins) || numel(opts.xbins) ~= 1 || opts.xbins
    x(x < opts.min | x > opts.max) = [];
end
if ~isfield(opts, 'ymax') || numel(opts.ymax) ~= 1 || ~isa(opts.ymax, 'double')
    opts.ymax = [];
end

% histogram
d = opts.max - opts.min;
s = d / opts.bins;
e = opts.min:s:(opts.max + s/2);
h = histcount(x, opts.min, opts.max - d / (2 * opts.bins), s);
if opts.density
    h = h ./ sum(h);
end
if isempty(opts.ymax)
    opts.ymax = max(h);
end
hh = histogram(p, 'BinEdges', e, 'BinCounts', h);
set(hh, 'FaceColor', opts.color);
hold(p, 'on');
if ~isempty(opts.sline)
    hx = repmat(h(:)', 10, 1);
    hx = [hx(:); 0];
    sk = smoothkern(10 * opts.sline);
    sk(sk < 1e-4) = [];
    lsk = floor(0.5 * numel(sk)) - 4;
    hc = conv(hx, sk, 'full');
    hc = hc(lsk:end+1-lsk);
    xe = (opts.min - s/2):s/10:(opts.max + s/2);
    sl = plot(p, xe(:), hc);
    set(sl, 'LineWidth', opts.slinew);
end
if ~isempty(opts.medline)
    m = median(x);
    ml = plot(p, [m; m], [0; opts.ymax]);
    set(ml, 'LineWidth', opts.medline, 'LineStyle', '-.');
    if ~isempty(opts.medlinec)
        try
            set(ml, 'Color', opts.medlinec);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    end
end
set(p, 'XLim', [opts.min, opts.max], 'YLim', [0, opts.ymax]);
