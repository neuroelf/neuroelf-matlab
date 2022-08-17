function tile_svs(svsfile, opts)
% TILE_SVS  Split SVS image into tiles
%
% FORMAT: TILE_SVS(svsfile) splits svsfile with default options.
%
%         TILE_SVS(svsfile, options) splits svsfile with specified
%         options (1x1 struct), valid settings are
%         .infoint  1x1 info printout interval, default: 5 (seconds)
%         .jpgqual  1x1 JPG writing quality, default: 95
%         .minvar   1x1 minimum variance in RGB data, default: 72
%         .outdir   char (string) with target folder, default: 'tiles'
%         .prefix   char (string) with tile filename prefix
%         .tilesize 1x2 double of [rows, columns]; default: auto-detect
%

% check inputs
if nargin < 1 || ~ischar(svsfile) || isempty(svsfile)
    error('neuroelf:badArgument', 'Bad or missing SVSFILE input.');
end
svsfile = svsfile(:)';
if exist(svsfile, 'file') ~= 2
    error('neuroelf:fileNotFound', 'SVSFILE not found.');
end
try
    svsinfo = imfinfo(svsfile);
    ih = svsinfo(1).Height;
    iw = svsinfo(1).Width;
    th = svsinfo(1).TileLength;
    tw = svsinfo(1).TileWidth;
    xt = ceil(iw / tw);
    yt = ceil(ih / th);
    if numel(svsinfo(1).TileOffsets) ~= (xt * yt)
        error('neureolf:badSVSFile', 'Invalid number of tiles in SVSFILE.');
    end
catch ne_eo;
    rethrow(ne_eo);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'infoint') || ~isa(opts.infoint, 'double') || ...
    numel(opts.infoint) ~= 1 || isinf(opts.infoint) || isnan(opts.infoint) || ...
    opts.infoint < 1
    infoint = 5;
else
    infoint = opts.infoint;
end
if ~isfield(opts, 'jpgqual') || ~isa(opts.jpgqual, 'double') || ...
    numel(opts.jpgqual) ~= 1 || isinf(opts.jpgqual) || isnan(opts.jpgqual) || ...
    opts.jpgqual < 1 || opts.jpgqual > 100
    jpgqual = 95;
else
    jpgqual = round(opts.jpgqual);
end
if ~isfield(opts, 'minvar') || ~isa(opts.minvar, 'double') || ...
    numel(opts.minvar) ~= 1 || isinf(opts.minvar) || isnan(opts.minvar) || ...
    opts.minvar < 0 || opts.minvar > 255
    minvar = 72;
else
    minvar = opts.minvar;
end
if ~isfield(opts, 'outdir') || ~ischar(opts.outdir)
    opts.outdir = 'tiles';
elseif isempty(opts.outdir)
    opts.outdir = '.'
elseif exist(opts.outdir(:)', 'dir') ~= 7
    error('neuroelf:badArgument', 'Invalid OUTDIR folder in OPTS.');
end
outdir = opts.outdir(:)';
if exist(outdir, 'dir') ~= 7
    try
        mkdir(outdir);
    catch ne_eo;
        rethrow(ne_eo);
    end
end
if outdir(end) ~= filesep
    outdir(end+1) = filesep;
end
if ~isfield(opts, 'prefix') || ~ischar(opts.prefix) || isempty(opts.prefix)
    [~, prefix] = fileparts(svsfile);
else
    prefix = opts.prefix(:)';
end
if ~isfield(opts, 'tilesize') || ~isa(opts.tilesize, 'double') || numel(opts.tilesize) ~= 2 || ...
    any(isinf(opts.tilesize) | isnan(opts.tilesize) | opts.tilesize < 16)
    tilesize = [th, tw];
else
    tilesize = ceil(opts.tilesize(:)');
end

% iterate over width and height
y1 = 1;
rt = 1;
tc = 1;
tt = ceil(ih / tilesize(1)) * ceil(iw / tilesize(2));
fprintf('Tiling SVS %s into %d tiles...\n', svsfile, tt);
tp = 0.005 * tt;
tn = tp - 0.05;
tm1 = now;
tmn = tm1 + (infoint / 86400);
wt = 0;
while y1 <= ih
    y2 = min(y1 + tilesize(1) - 1, ih);
    x1 = 1;
    ct = 1;
    while x1 <= iw
        x2 = min(x1 + tilesize(2) - 1, iw);
        tdata = imread(svsfile, 'PixelRegion', {[y1, y2], [x1, x2]});
        tvar = var(double(tdata(:)));
        if tvar > minvar
            imwrite(tdata, sprintf('%s%s-tile-%04d-%04d.jpg', ...
                outdir, prefix, rt, ct), 'Quality', jpgqual);
            wt = wt + 1;
        end
        x1 = x1 + tilesize(2);
        ct = ct + 1;
        tc = tc + 1;
        tm2 = now;
        if tm2 > tmn
            tmn = tm2 + (infoint / 86400);
            tm2 = 86400 * (tm2 - tm1);
            fprintf('Written: %5d/%5d/%d tiles, Elapsed: %.1fsec, Remaining: %.1fsec\n', ...
                wt, tc, tt, tm2, tm2 * (tt - tc) / tc);
            tn = tn + tp;
        end
    end
    y1 = y1 + tilesize(1);
    rt = rt + 1;
end
fprintf('\n');
