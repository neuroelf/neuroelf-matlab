function [font, glyphs, glyphw] = createfontmat(fontname, fontsize, letters)
%CREATEFONTMAT  Store MAT file with FONT glyph image and information
%   FONT = CREATEFONTMAT(FONTNAME) creates a MAT file (content) for the
%   font FONTNAME.
%
%   FONT = CREATEFONTMAT(FONTNAME, FONTSIZE) overrides the font size.
%
%   FONT = CREATEFONTMAT(FONTNAME, FONTSIZE, LETTERS) also overrides the
%   default letters selection ([32:127, 161:255]), allowing UniCode as well.

% Version:  v1.1
% Build:    16060110
% Date:     Jun-01 2016, 10:21 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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

% check inputs
if nargin < 1 || ~ischar(fontname) || isempty(fontname)
    error('neuroelf:general:badArgument', 'Invalid fontname argument.');
end
if nargin < 2 || ~isa(fontsize, 'double') || numel(fontsize) ~= 1
    fontsize = 72;
end
figsize = ceil(fontsize .* [5, 4]);
fs8 = 2 * figsize(2);

% generate figure and axes
f = figure;
set(f, 'Units', 'pixels');
set(f, 'Position', [10, 10, figsize]);
set(f, 'Units', 'normalized');
set(f, 'Visible', 'off');
limage = [tempname '.png'];
ne_screenshot(0, 0, f, limage);
lim = imread(limage);
delete(limage);
if size(lim, 1) == fs8
    set(f, 'Units', 'pixels');
    set(f, 'Position', [10, 10, ceil(0.5 .* figsize)]);
    set(f, 'Units', 'normalized');
    ufontsize = 0.5 * fontsize;
else
    ufontsize = fontsize;
end
a = axes('Parent', f);
set(a, 'Color', [0, 0, 0]);
set(a, 'Position', [0, 0, 1, 1]);
set(f, 'Units', 'pixels');

% add text element
t = text(a, 0.5, 0.5, 'A');
set(t, 'Interpreter', 'none');
set(t, 'Color', [1, 1, 1]);
set(t, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
set(t, 'FontSize', ufontsize);
set(t, 'FontName', fontname(:)');
set(t, 'Position', [0.5, 0.4, 0]);
set(a, 'XLim', [0, 1], 'YLim', [0, 1]);
set(a, 'XTick', [], 'YTick', []);
set(f, 'Visible', 'off');

% progress bar
if nargin < 3 || ~isa(letters, 'double') || isempty(letters) || ...
    any(isinf(letters(:)) | isnan(letters(:)) | letters (:) < 1 | letters(:) > 65535)
    letters = [32:127, 161:255];
else
    letters = unique(fix(letters(:)))';
end
pbar = xprogress;
xprogress(pbar, 'settitle', 'Creating font (glyph) images...');
xprogress(pbar, 0, 'Letter images...', 'visible', 0, numel(letters));

% now get to work
glyphs = cell(max(letters), 1);
glyphw = zeros(size(glyphs));
for lc = 1:numel(letters)

    % set text
    letter = letters(lc);
    if letter == 32
        set(t, 'String', '|| |');
    else
        set(t, 'String', char(letter));
    end

    % take and load image
    limage = [tempname '.png'];
    ne_screenshot(0, 0, f, limage);
    lim = imread(limage);
    delete(limage);

    % get middle portion
    lim([1, end], :, :) = 0;
    lim(:, [1, end], :) = 0;
    lim = max(lim, [], 3);
    limset = (lim > 0);
    limx = find(any(limset, 1));
    if letter == 32
        limx = find(limset(ceil(size(limset, 1) / 2), :));
        limxd = find(diff(limx(:)) > 1);
        limx = (limx(limxd(end)+1) - limx(limxd(end))) - ...
            ceil(0.67 * (limx(limxd(1)+1) - limx(limxd(1))));
        lim = uint8(zeros(size(lim, 1), limx));
        limx = 1:limx;
    elseif isempty(limx)
        continue;
    end
    glyphs{letter} = lim(:, limx);
    glyphw(letter) = numel(limx);

    % progress
    xprogress(pbar, lc);
end

% clean up
closebar(pbar);
delete(f);

% compile image
letters = find(glyphw > 0);
letters = letters(:);
xstop = cumsum(glyphw(letters, 1));
lim = cat(2, glyphs{letters});
maxlim = max(lim, [], 2);
usey = find(maxlim > 0);
ybase = 1 + round(0.6 .* size(lim, 1) - usey(1));
lim = lim(usey, :);
timname = [tempname '.png'];
imwrite(lim, timname);
lim = binread(timname);
delete(timname);
font = struct('image', lim, 'letters', letters, 'size', fontsize, ...
    'xktab', [], 'xstop', xstop, 'ybase', ybase);
save([neuroelf_path('fonts') '/' fontname(:)' '.mat'], 'font');
