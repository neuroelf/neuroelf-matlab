function [im, ima] = image_font(text, font, fsize, opts)
%IMAGE_FONT  Set text in font and return as image.
%   FIMAGE = IMAGE_FONT(TEXT) sets the text in TEXT in the default font.
%
%   FIMAGE = IMAGE_FONT(TEXT, FONT) uses font named FONT instead (if the
%   named font isn't found, the default font is used).
%
%   FIMAGE = IMAGE_FONT(TEXT, FONT, SIZE) uses the specified font size
%   instead.
%
%   FIMAGE = IMAGE_FONT(TEXT, FONT, SIZE, OPTS) allows to specify some
%   additional options in a 1x1 struct with supported fields:
%      .bcolor      1x3 background color (double, default: [255, 255, 255])
%      .color       1x3 font color (double, default: [0, 0, 0])
%      .color2      1x3 font color (double, default: []), which can be
%                   used to generate a color gradient
%      .colorgrad   either 1x4 double (relative or absolute xy1 / xy2) or
%                   1x3 (point + radius); relative is [0, 1] in both dims
%                   and radius is relative to the diagonal of the image
%      .emboss      1x2 double (smoothing strength and lighting change)
%      .gcolor      1x3 double (default: half-way towards [255, 255, 255])     
%      .gcolor2     1x3 double (default: half-way towards [255, 255, 255])     
%      .glow        1x3 double (smoothing strength, sharpness threshold,
%                   default: [0, 0, 0] = no glow)
%      .halign      horizontal alignment one of {'center'}, 'left', 'right'
%      .invert      invert the output of font setting
%      .outsize     1x2 size of output image (default: determined by TEXT)
%      .padding     1x1 double (only used if outsize not given, default: 0)
%      .scolor      1x3 double for shadow color (default: [0, 0, 0])
%      .shadow      1x4 double with x/y distance, smoothing kernel and
%                   strength (e.g. [8, 4, 12, 0.5], default: [0, 0, 0, 0])
%      .tablines    table-based lines either of {'all'}, 'inner', or 'none'
%
%   Example:
%
%   im = IMAGE_FONT(['This is a test...' char(10) 'And another line.']);
%   image(im);

% Version:  v1.1
% Build:    20061213
% Date:     Jun-12 2020, 1:15 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, 2020, Jochen Weber
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

% persistent fonts
persistent imfonts;
if numel(imfonts) ~= 1 || ~isstruct(imfonts)

    % load fonts
    imfonts = struct;
    fontnames = findfiles(neuroelf_path('fonts'), '*.mat', '-d1');
    for fc = 1:numel(fontnames)
        [fontfolder, fontname] = fileparts(fontnames{fc});
        imfonts.(fontname) = load(fontnames{fc});
        imfonts.(fontname) = imfonts.(fontname).font;
        timname = [tempname '.png'];
        binwrite(timname, imfonts.(fontname).image);
        imfonts.(fontname).image = imread(timname);
        delete(timname);
        imfonts.(fontname).lmap = zeros(1, 1024);
        imfonts.(fontname).lmap(imfonts.(fontname).letters) = 1:numel(imfonts.(fontname).letters);
        imfonts.(fontname).xstart = 1 + [0; imfonts.(fontname).xstop(1:end-1, 1)];
        imfonts.(fontname).yret = imfonts.(fontname).size - size(imfonts.(fontname).image, 1);
        imfonts.(fontname) = fontaddkerning(imfonts.(fontname));
    end

    % settings
    imfonts.defaults = struct('bcolor', [255, 255, 255], 'color', [0, 0, 0], ...
        'font', 'calibri', 'fsize', 24);
end

% check arguments
if nargin > 0 && ischar(text)

    % parse into lines
    if ~isempty(strfind(text, char([13, 10])))
        text = splittocellc(text, char([13, 10]), false, false);
    else
        text = splittocellc(text, char([10, 13]), false, true);
    end
    text = text(:);
end
if nargin < 1 || ~iscell(text) || ~all(cellfun(@ischar, text(:)))
    error('neuroelf:general:badArgument', 'Invalid text argument.');
end
tabsize = size(text);
text = text(:);
ntext = numel(text);
if nargin < 2 || isempty(font)
    font = {imfonts.(imfonts.defaults.font)};
elseif ischar(font) && ~isempty(font) && isfield(imfonts, lower(font(:)'))
    font = {imfonts.(lower(font(:)'))};
elseif ~iscell(font) || ~any(numel(font) == [1, ntext]) || ~all(cellfun(@ischar, font(:)))
    font = {imfonts.(imfonts.defaults.font)};
else
    font = font(:);
    for fc = 1:numel(font)
        if ~isfield(imfonts, lower(font{fc}(:)'))
            font{fc} = imfonts.(imfonts.defaults.font);
        else
            font{fc} = imfonts.(lower(font{fc}(:)'));
        end
    end
end
if ntext > 1 && numel(font) == 1
    font = repmat(font, ntext, 1);
end
if nargin < 3 || ~isa(fsize, 'double') || ~any(numel(fsize) == [1, ntext])
    fsize = imfonts.defaults.fsize;
else
    fsize(isinf(fsize(:)) | isnan(fsize(:)) | fsize(:) <= 0) = imfonts.defaults.fsize;
end
fsize = ceil(fsize);
if ntext > 1 && numel(fsize) == 1
    fsize = fsize .* ones(ntext, 1);
end
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'bcolor') || ~isnumeric(opts.bcolor) || numel(opts.bcolor) ~= 3 || ...
    any(isinf(opts.bcolor) | isnan(opts.bcolor) | opts.bcolor < 0 | opts.bcolor > 255)
    opts.bcolor = imfonts.defaults.bcolor;
elseif all(opts.bcolor <= 1) && ~all(opts.bcolor == round(opts.bcolor))
    opts.bcolor = round(255 .* opts.bcolor(:)');
else
    opts.bcolor = round(opts.bcolor(:)');
end
if ~isfield(opts, 'color') || ~isnumeric(opts.color) || numel(opts.color) ~= 3 || ...
    any(isinf(opts.color) | isnan(opts.color) | opts.color < 0 | opts.color > 255)
    opts.color = imfonts.defaults.color;
elseif all(opts.color <= 1) && ~all(opts.color == round(opts.color))
    opts.color = round(255 .* opts.color(:)');
else
    opts.color = round(opts.color(:)');
end
if ~isfield(opts, 'color2') || ~isnumeric(opts.color2) || numel(opts.color2) ~= 3 || ...
    any(isinf(opts.color2) | isnan(opts.color2) | opts.color2 < 0 | opts.color2 > 255)
    opts.color2 = 255 - imfonts.defaults.color;
elseif all(opts.color2 <= 1) && ~all(opts.color2 == round(opts.color))
    opts.color2 = round(255 .* opts.color2(:)');
else
    opts.color2 = round(opts.color2(:)');
end
if ~isfield(opts, 'colorgrad') || ~isa(opts.colorgrad, 'double') || ~any(numel(opts.colorgrad) == [3, 4]) || ...
    any(isinf(opts.colorgrad(:)) | isnan(opts.colorgrad(:)) | opts.colorgrad(:) < 0)
    opts.colorgrad = [];
else
    opts.colorgrad = opts.colorgrad(:)';
end
if ~isfield(opts, 'emboss') || ~isa(opts.emboss, 'double') || numel(opts.emboss) ~= 2 || ...
    any(isinf(opts.emboss) | isnan(opts.emboss) | opts.emboss < -1) || opts.emboss(1) <= 0
    opts.emboss = [0, 0];
else
    opts.emboss = opts.emboss(:)';
    opts.emboss(2) = min(1, opts.emboss(2));
end
if ~isfield(opts, 'gcolor') || ~isnumeric(opts.gcolor) || numel(opts.gcolor) ~= 3 || ...
    any(isinf(opts.gcolor) | isnan(opts.gcolor) | opts.gcolor < 0 | opts.gcolor > 255)
    opts.gcolor = round(0.5 .* (opts.color + [255, 255, 255]));
elseif all(opts.gcolor <= 1) && ~all(opts.gcolor == round(opts.gcolor))
    opts.gcolor = round(255 .* opts.gcolor(:)');
else
    opts.gcolor = round(opts.gcolor(:)');
end
if ~isfield(opts, 'gcolor2') || ~isnumeric(opts.gcolor2) || numel(opts.gcolor2) ~= 3 || ...
    any(isinf(opts.gcolor2) | isnan(opts.gcolor2) | opts.gcolor2 < 0 | opts.gcolor2 > 255)
    opts.gcolor2 = round(0.5 .* (opts.color2 + [255, 255, 255]));
elseif all(opts.gcolor2 <= 1) && ~all(opts.gcolor2 == round(opts.gcolor2))
    opts.gcolor2 = round(255 .* opts.gcolor2(:)');
else
    opts.gcolor2 = round(opts.gcolor2(:)');
end
if ~isfield(opts, 'glow') || ~isa(opts.glow, 'double') || numel(opts.glow) ~= 3 || ...
    any(isinf(opts.glow) | isnan(opts.glow) | opts.glow < 0) || opts.glow(3) > opts.glow(2)
    opts.glow = [0, 0, 0];
else
    opts.glow = opts.glow(:)';
end
if ~isfield(opts, 'halign') || ~ischar(opts.halign) || isempty(opts.halign)
    opts.halign = 'center';
else
    opts.halign = lower(opts.halign(:)');
    if ~any(strcmp(opts.halign, {'center', 'left', 'right'}))
        opts.halign = 'center';
    end
end
if ~isfield(opts, 'invert') || ~islogical(opts.invert) || numel(opts.invert) ~= 1
    opts.invert = false;
end
if ~isfield(opts, 'outsize') || ~isa(opts.outsize, 'double') || numel(opts.outsize) ~= 2 || ...
    any(isinf(opts.outsize) | opts.outsize < 0)
    opts.ooutsize = [];
else
    opts.ooutsize = round(opts.outsize(:)');
end
opts.outsize = [NaN, NaN];
if ~isfield(opts, 'padding') || ~isa(opts.padding, 'double') || numel(opts.padding) ~= 1 || ...
    isinf(opts.padding) || isnan(opts.padding) || opts.padding < 0
    opts.padding = 0;
end
if ~isfield(opts, 'shadow') || ~isa(opts.shadow, 'double') || numel(opts.shadow) ~= 4 || ...
    any(isinf(opts.shadow(:)) | isnan(opts.shadow(:))) || opts.shadow(3) < 0 || ...
    opts.shadow(4) <= 0 || opts.shadow(4) > 1
    opts.shadow = [0, 0, 0, 0];
    doshadow = false;
else
    opts.shadow = opts.shadow(:)';
    doshadow = true;
end
if ~isfield(opts, 'scolor') || ~isnumeric(opts.scolor) || numel(opts.scolor) ~= 3 || ...
    any(isinf(opts.scolor) | isnan(opts.scolor) | opts.scolor < 0 | opts.scolor > 255)
    opts.scolor = [0, 0, 0];
elseif all(opts.scolor <= 1) && ~all(opts.scolor == round(opts.scolor))
    opts.scolor = round(255 .* opts.scolor(:)');
else
    opts.scolor = round(opts.scolor(:)');
end
if ~isfield(opts, 'spkern') || ~isa(opts.spkern, 'double') || numel(opts.spkern) ~= 1 || ...
    isinf(opts.spkern) || isnan(opts.spkern)
    opts.spkern = 0;
else
    opts.spkern = round(opts.spkern);
end
if ~isfield(opts, 'tablines') || ~ischar(opts.tablines) || isempty(opts.tablines)
    opts.tablines = 'all';
else
    opts.tablines = lower(opts.tablines(:)');
    if ~any(strcmp(opts.tablines, {'all', 'cols', 'icols', 'inner', 'irows', 'none', 'rows'}))
        opts.tablines = 'all';
    end
end
if ~strcmp(opts.tablines, 'none') && any(tabsize > 1)
    opts.padding = max(1, opts.padding);
end
if ~isfield(opts, 'xkern') || ~isa(opts.xkern, 'double') || numel(opts.xkern) ~= 1 || ...
    isinf(opts.xkern) || isnan(opts.xkern)
    opts.xkern = 0;
else
    opts.xkern = round(opts.xkern);
end
for fc = 1:numel(font)
    font{fc}.spkern = opts.spkern;
    font{fc}.xkern = opts.xkern;
end

% set each line with current settings
lines = fontsetlines(text, font, fsize);

% get sizes
fsize = cellfun('size', lines, 1);
fsize(cellfun('isempty', lines)) = 0;
fwid = cellfun('size', lines, 2);
fwid(cellfun('isempty', lines)) = 0;

% tabular option
if tabsize(2) > 1
    lines = reshape(lines, tabsize);
    fsize = reshape(fsize, tabsize);
    fwid = reshape(fwid, tabsize);
end
mfsize = max(fsize, [], 2);
mfwid = max(fwid, [], 1);

% outsize needs to be determined
if isnan(opts.outsize(1))
    opts.outsize(1) = sum(mfsize) + tabsize(1) * round(2 * opts.padding) + abs(opts.shadow(2));
    xpad = opts.padding;
else
    xpad = 0;
end
if isnan(opts.outsize(2))
    opts.outsize(2) = sum(mfwid) + tabsize(2) * (round(2 * opts.padding) + abs(opts.shadow(1)));
    ypad = opts.padding;
else
    ypad = 0;
end

% create image
im = uint8(0);
im(opts.outsize(1), opts.outsize(2), 3) = 0;
for pc = 1:3
    im(:, :, pc) = opts.bcolor(pc);
end
ima = zeros(opts.outsize);

% color gradient
if numel(opts.colorgrad) == 4
    if all(opts.colorgrad <= 1)
        psize = 2 .* [xpad, ypad, xpad, ypad];
        opts.colorgrad = 1 + round((opts.outsize([2, 1, 2, 1]) - (psize + 1)) .* opts.colorgrad);
    else
        opts.colorgrad = max(ones(1, 4), min(opts.outsize([2, 1, 2, 1]), opts.colorgrad));
    end
    xf = opts.colorgrad(1);
    xd = opts.colorgrad(3) - xf;
    xf = xf + xpad;
    yf = opts.colorgrad(2);
    yd = opts.colorgrad(4) - yf;
    yf = yf + ypad;
    xyd = 1 / sqrt(xd * xd + yd * yd);
    xya = atan2(yd, xd);
    if isinf(xyd)
        cg = zeros(opts.outsize);
    else
        [cgy, cgx] = ndgrid(1:opts.outsize(1), 1:opts.outsize(2));
        cg = max(0, min(1, (xyd * cos(xya)) .* (cgx - xf) + (xyd * sin(xya)) .* (cgy - yf)));
    end
elseif numel(opts.colorgrad) == 3
    cg = zeros(opts.outsize);
else
    cg = [];
end

% store font pieces into
for cc = 1:tabsize(2)
    xfrom = (2 * (cc - 1) + 1) * xpad + sum(mfwid(1:cc-1));
    yfrom = ypad + 1;
    for lc = 1:tabsize(1)
        lim = (1 / 255) .* double(lines{lc, cc});
        if size(lim, 2) > opts.outsize(2)
            lim = lim(:, 1:opts.outsize(2));
        end
        xafrom = xfrom;
        tsize = size(lim, 2);
        if tsize < mfwid(cc)
            if opts.halign(1) == 'c'
                xafrom = xafrom + floor(0.5 * (mfwid(cc) - tsize));
            elseif opts.halign(1) == 'r'
                xafrom = xafrom + mfwid(cc) - tsize;
            end
        end
        if opts.invert
            lim = 1 - lim;
        end
        ima(yfrom:yfrom+mfsize(lc)-1, xafrom+1:xafrom+tsize) = lim;
        lim = double(lim > 0);
        if isempty(cg)
            cgim = [];
        else
            cgim = cg(yfrom:yfrom+mfsize(lc)-1, xafrom+1:xafrom+tsize);
            cgimn = 1 - cgim;
        end
        for pc = 1:3
            cim = double(im(yfrom:yfrom+mfsize(lc)-1, xafrom+1:xafrom+tsize, pc));
            if isempty(cgim) && ~isempty(cim)
                im(yfrom:yfrom+mfsize(lc)-1, xafrom+1:xafrom+tsize, pc) = ...
                    round((1 - lim) .* cim + lim .* opts.color(pc));
            elseif ~isempty(cgim) && ~isempty(cim)
                im(yfrom:yfrom+mfsize(lc)-1, xafrom+1:xafrom+tsize, pc) = round((1 - lim) .* cim + ...
                    lim .* (cgimn .* opts.color(pc) + cgim .* opts.color2(pc)));
            end
        end
        yfrom = yfrom + mfsize(lc) + 2 * ypad;
    end
end

% table lines
if tabsize(2) > 1 && any(strcmp(opts.tablines, {'all', 'cols', 'icols', 'inner'}))
    for cc = 1:tabsize(2)
        if cc == 1 && opts.tablines(1) == 'i'
            continue;
        end
        xfrom = 1 + (2 * (cc - 1)) * xpad + sum(mfwid(1:cc-1));
        ima(:, xfrom) = 1;
        if isempty(cg)
            cgim = [];
        else
            cgim = cg(:, xfrom);
            cgimn = 1 - cgim;
        end
        for pc = 1:3
            if isempty(cgim)
                im(:, xfrom, pc) = round(opts.color(pc));
            else
                im(:, xfrom, pc) = round(...
                    (cgimn .* opts.color(pc) + cgim .* opts.color2(pc)));
            end
        end
    end
    if opts.tablines(1) ~= 'i'
        ima(:, end) = 1;
        if isempty(cg)
            cgim = [];
        else
            cgim = cg(:, end);
            cgimn = 1 - cgim;
        end
        for pc = 1:3
            if isempty(cgim)
                im(:, end, pc) = round(opts.color(pc));
            else
                im(:, end, pc) = round(...
                    (cgimn .* opts.color(pc) + cgim .* opts.color2(pc)));
            end
        end
    end
end
if tabsize(1) > 1 && any(strcmp(opts.tablines, {'all', 'irows', 'inner', 'rows'}))
    for lc = 1:tabsize(1)
        if lc == 1 && opts.tablines(1) == 'i'
            continue;
        end
        yfrom = 1 + (2 * (lc - 1)) * ypad + sum(mfsize(1:lc-1));
        ima(yfrom, :) = 1;
        if isempty(cg)
            cgim = [];
        else
            cgim = cg(yfrom, :);
            cgimn = 1 - cgim;
        end
        for pc = 1:3
            if isempty(cgim)
                im(yfrom, :, pc) = round(opts.color(pc));
            else
                im(yfrom, :, pc) = round(...
                    (cgimn .* opts.color(pc) + cgim .* opts.color2(pc)));
            end
        end
    end
    if opts.tablines(1) ~= 'i'
        ima(end, :) = 1;
        if isempty(cg)
            cgim = [];
        else
            cgim = cg(end, :);
            cgimn = 1 - cgim;
        end
        for pc = 1:3
            if isempty(cgim)
                im(end, :, pc) = round(opts.color(pc));
            else
                im(end, :, pc) = round(...
                    (cgimn .* opts.color(pc) + cgim .* opts.color2(pc)));
            end
        end
    end
end

% build output image
t = transimg(size(ima, 2), size(ima, 1), opts.bcolor);

% compute shadow
if doshadow

    % get cubic kernel
    [nullipd, ckern] = flexinterpn_method(zeros(5, 1), [Inf; 2; 1; 4], 'cubic');
    ckern{1} = conv(ckern{1}(1:256:end, :), smoothkern(opts.shadow(3) * ckern{2} / 256));
    ckern{2} = ckern{2} / 256;

    % transform alpha channel
    simas = [1, 1; size(ima)] - ones(2, 1) * opts.shadow(1, [2, 1]);
    sima = min(1, max(0, opts.shadow(4) .* ...
        flexinterpn(ima, [inf, inf; simas(1, :); 1, 1; simas(2, :)], ckern{1}, ckern{2})));

    % build output image
    rgbp = uint8(ones(size(ima)));
    addlayer(t, cat(3, uint8(opts.scolor(1)) .* rgbp, uint8(opts.scolor(2)) .* rgbp, ...
        uint8(opts.scolor(3)) .* rgbp), sima);
end

% glow
if opts.glow(1) > 0
    gima = limitrangec(2.5 .* ...
        flexinterpn(ima, [inf, inf; ones(2, 2); size(ima)], smoothkern(2.5), 1), 0, 1, 0);
    gima = limitrangec(-opts.glow(3) + opts.glow(2) .* ...
        flexinterpn(gima, [inf, inf; ones(2, 2); size(ima)], smoothkern(opts.glow(1)), 1), 0, 1, 0);
    if any(gima(:) > 0)
        grgb = repmat(reshape(uint8(opts.gcolor), [1, 1, 3]), size(ima));
        if ~isempty(cg)
            cgim = repmat(cg, [1, 1, 3]);
            cgimn = 1 - cgim;
            grgb2 = repmat(reshape(uint8(opts.gcolor2), [1, 1, 3]), size(ima));
            grgb = uint8(round(cgimn .* double(grgb) + cgim .* double(grgb2)));
        end
        addlayer(t, grgb, gima);
    end
end

% emboss
if opts.emboss(1) ~= 0
    eima = ima - flexinterpn(ima, [inf, inf; ones(2, 2); size(ima)], smoothkern(opts.emboss(1)), 1);
    hsv = hsvconv(im, 2);
    hsv(:, :, 3) = limitrangec(hsv(:, :, 3) + opts.emboss(2) .* eima, 0, 1, 0);
    im = hsvconv(hsv, 1);
end

% add regular layer
addlayer(t, im, ima);

% render, and get rendered and alpha image
nlay = numel(t.Layer);
if nlay > 1
    render(t);
    im = t.Rendered;
    joinlayers(t, 1:nlay);
    tl = t.Layer(1);
    ima = limitrangec(double(tl.Alpha), 0, 1, 0);
else
    im = t.Layer(1).Pixel;
end
delete(t);

% join alpha information
if nargout < 2
    im = double(im);
    for pc = 1:3
        im(:, :, pc) = opts.bcolor(pc) .* (1 - ima) + ima .* im(:, :, pc);
    end
    im = uint8(round(im));
end



% internal function, set single lines into images
function lim = fontsetlines(lines, font, fsize)
lim = cell(numel(lines), 1);
for lc = 1:numel(lim)

    % font size
    lfont = font{lc};
    ifsize = size(lfont.image);

    % get all letter images first
    line = double(lines{lc});
    nline = numel(line);
    if nline == 0
        lim{lc} = uint8(zeros(ifsize(1), 0));
        continue;
    end
    letims = cell(nline, 1);
    letspc = zeros(size(letims));
    for letc = 1:nline
        leti = lfont.lmap(line(letc));
        if leti == 0
            leti = lfont.lmap(63);
        end
        letims{letc} = lfont.image(:, lfont.xstart(leti):lfont.xstop(leti));

        % kerning
        if letc < nline
            nleti = lfont.lmap(line(letc+1));
            if nleti == 1
                letspc(letc) = lfont.kerning(leti, nleti) + font{lc}.spkern;
            elseif nleti > 1
                letspc(letc) = lfont.kerning(leti, nleti) + font{lc}.xkern;
            end
        end
    end

    % total size
    xsims = cellfun('size', letims, 2);

    % combine total image
    xstot = sum(xsims) + sum(letspc);
    lineimage = uint8(0);
    lineimage(ifsize(1), xstot) = 0;
    lii = 1;
    for letc = 1:numel(letims)
        lineimage(:, lii:lii+xsims(letc)-1) = ...
            max(lineimage(:, lii:lii+xsims(letc)-1), letims{letc});
        lii = lii + xsims(letc) + letspc(letc);
    end

    % resize
    if fsize ~= font{lc}.size
        ffact = fsize(lc) / font{lc}.size;
        if ffact < 1
            lineimage = image_resize(lineimage, ffact);
        else
            lineimage = image_resize(lineimage, ...
                round(ffact * size(lineimage, 1)), round(ffact * size(lineimage, 2)));
        end
        mli = double(max(lineimage(:)));
        if mli < 255
            lineimage = uint8((255 / mli) .* double(lineimage));
        end 
    end

    % store
    lim{lc} = lineimage;
end


% add kerning
function font = fontaddkerning(font)

% for each letter
nl = numel(font.letters);
fh = size(font.image, 1);
fhi = ([[1, 1:fh-1]; (1:fh); [2:fh, fh]])';
lsf = NaN .* zeros(fh, nl);
rsf = lsf;
for lc = 1:nl

    % get the letter image (masked)
    lmi = (font.image(:, font.xstart(lc):font.xstop(lc)) >= 128);
    lms = size(lmi, 2);
    lmi = squeeze(any(reshape(lmi(fhi, :), [fh, 3, lms]), 2));

    % find the first pixel that is not background
    cpix = find(sum(lmi, 2) > 1);
    for cc = 1:numel(cpix)
        cpc = cpix(cc);
        lsf(cpc, lc) = findfirst(lmi(cpc, :)) - 1;
        rsf(cpc, lc) = lms - findfirst(lmi(cpc, :), -1);
    end
end

% next compute the median for each pair
ktab = NaN .* zeros(nl, nl);
for rc = 1:nl
    rsfc = rsf(:, rc);
    if all(isnan(rsfc))
        continue;
    end
    nrf = sum(~isnan(rsfc));
    for lc = 1:nl
        rsflsf = rsfc + lsf(:, lc);
        if all(isnan(rsflsf))
            continue;
        end
        nlf = sum(~isnan(lsf(:, lc)));
        rsflsf = rsflsf(~isnan(rsflsf));
        nrl = numel(rsflsf);
        rlmin = min(rsflsf);
        rlmed = 1 + median(rsflsf);
        minw = (rlmed - rlmin) / rlmed;
        ktab(rc, lc) = (minw * rlmin + (1 - minw) * rlmed) * (nrl * nrl / (nrf * nlf));
    end
end

% add any existing "kerning additions"
if isfield(font, 'xktab') && isequal(size(ktab), size(font.xktab))
    ktab = ktab - full(font.xktab);
end

% overall median
ktmed = median(ktab(~isnan(ktab(:)))) + 0.1 * size(font.image, 2) / nl;
ktab = ktmed - ktab;
ktab(isnan(ktab)) = 0;

% adjust "space" kerning
ktabsp = ceil(mean(mean(ktab(2:end, 2:end))));
ktab(1, :) = -ktabsp;
ktab(:, 1) = -ktabsp;

% store table
font.kerning = fix(ktab);
