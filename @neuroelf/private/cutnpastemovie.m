function cutnpastemovie(bmps, pos, config)
% cutnpastemovie  - create an AVI movie from bitmap files
%
% FORMAT:       cutnpastemovie(bmps, pos, config)
%
% Input fields:
%
%       bmps        1xL cell array of cell arrays with bitmap filenames
%                   each cell must have either one or a common N filenames
%       pos         Lx6 position information (xf, yf, xt, yt, w, h), where
%                       xf, yf   - X and Y source coordinate in bitmaps
%                       xt, yt   - target coordinate within final frame
%                       w, h     - width and height of clipping (Inf: all)
%       config      1x1 struct with mandatory fields
%        .color     background color (1x3 RGB, e.g. [0, 0, 0])
%        .filename  target AVI filename
%        .size      output size (e.g. [720, 576])
%                   and optional fields
%        .compress  Compression setting (default: None)
%        .fps       FPS setting (default: 30)
%        .kps       Keyframe (per second) setting (default: 2)
%        .quality   Quality setting (default: 75, but unused with None)
%        .vidname   Videoname setting (default: 'Bitmap clipping')
%
% No output fields.

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 1:55 PM EST
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
if nargin ~= 3 || ...
   ~iscell(bmps) || ...
    isempty(bmps) || ...
   ~isa(pos, 'double') || ...
    isempty(pos) || ...
    size(pos, 2) > 6 || ...
    any(isnan(pos(:)) | pos(:) < 1) || ...
   ~isstruct(config) || ...
    numel(config) ~= 1 || ...
   ~isfield(config, 'color') || ...
   ~isa(config.color, 'double') || ...
    numel(config.color) ~= 3 || ...
    any(isinf(config.color) | isnan(config.color) | config.color < 0 | config.color > 255) || ...
   ~isfield(config, 'filename') || ...
   ~ischar(config.filename) || ...
    isempty(config.filename) || ...
   ~isfield(config, 'size') || ...
   ~isa(config.size, 'double') || ...
    numel(config.size) ~= 2 || ...
    any(isinf(config.size) | isnan(config.size) | config.size < 16 | config.size > 1680)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
bgc = round(config.color(:)');
vsz = 16 * ceil(config.size(:) ./ 16)';

% create AVI file options
aoo = {config.filename(:)'};
if isfield(config, 'compress') && ...
    ischar(config.compress)
    aoo(end+1:end+2) = {'Compression', config.compress(:)'};
else
    aoo(end+1:end+2) = {'Compression', 'None'};
end
if isfield(config, 'fps') && ...
    isa(config.fps, 'double') && ...
    numel(config.fps) == 1
    aoo(end+1:end+2) = {'FPS', config.fps};
else
    aoo(end+1:end+2) = {'FPS', 30};
end
if isfield(config, 'kps') && ...
    isa(config.kps, 'double') && ...
    numel(config.kps) == 1
    aoo(end+1:end+2) = {'Keyframe', config.kps};
else
    aoo(end+1:end+2) = {'Keyframe', 2};
end
if isfield(config, 'quality') && ...
    isa(config.quality, 'double') && ...
    numel(config.quality) == 1
    aoo(end+1:end+2) = {'Quality', config.quality};
else
    aoo(end+1:end+2) = {'Quality', 75};
end
if isfield(config, 'vidname') && ...
    ischar(config.vidname)
    aoo(end+1:end+2) = {'Videoname', config.vidname(:)'};
else
    aoo(end+1:end+2) = {'Videoname', 'Bitmap clipping'};
end

% check cell arrays
na = numel(bmps);
ni = 1;
for cc = 1:na
    if isempty(bmps{cc})
        error( ...
            'neuroelf:BadArgument', ...
            'Empty bitmap array #%d.', ...
            cc ...
        );
    end
    for bc = 1:numel(bmps{cc})
        if ~ischar(bmps{cc}{bc}) || ...
            exist(bmps{cc}{bc}(:)', 'file') ~= 2
            error( ...
                'neuroelf:BadArgument', ...
                'Invalid file or file not found in array %d.', ...
                cc ...
            );
        end
    end
    if ni > 1
        if numel(bmps{cc}) > 1 && ...
            numel(bmps{cc}) ~= ni
            error( ...
                'neuroelf:BadArgument', ...
                'Invalid number of bitmaps in array %d.', ...
                cc ...
            );
        end
    else
        ni = max(ni, numel(bmps{cc}));
    end
end
if ni > 1
    for cc = 1:na
        if numel(bmps{cc}) == 1
            bmps{cc} = bmps{cc}(ones(1, ni));
        end
    end
end

% try to create AVI file object
try
    ao = avifile(aoo{:});
catch ne_eo;
    rethrow(ne_eo);
end

% prepare target frame
emv = uint8(zeros([vsz(2), vsz(1), 3]));
emv(:, :, 1) = bgc(1);
emv(:, :, 2) = bgc(2);
emv(:, :, 3) = bgc(3);

% try progress bar
try
    pbar = xprogress;
    xprogress(pbar, 'setposition', [80, 200, 640, 36]);
    xprogress(pbar, 'settitle', 'Creating clipped movie...');
    xprogress(pbar, 0, 'Processing frames...', 'visible', 0, ni);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    pbar = [];
end

% load bitmaps along time axis
for ic = 1:ni
    try
        ipos = pos;
        for bc = 1:na
            limg = imread(bmps{bc}{ic});
            if isinf(ipos(bc, 5))
                ipos(bc, 5) = size(limg, 2);
            end
            if isinf(ipos(bc, 6))
                ipos(bc, 6) = size(limg, 1);
            end
            srcd = limg(pos(bc, 2):pos(bc, 2) + ipos(bc, 6) - 1, pos(bc, 1):pos(bc, 1) + ipos(bc, 5) - 1, :);
            [srcx, srcy] = ndgrid(1:ipos(bc, 6), 1:ipos(bc, 5));
            [trgx, trgy] = ndgrid(pos(bc, 4):pos(bc, 4) + ipos(bc, 6) - 1, pos(bc, 3):pos(bc, 3) + ipos(bc, 5) - 1);
            srci = sub2ind(size(srcd), srcx(:), srcy(:));
            trgi = sub2ind(size(emv), trgx(:), trgy(:));
            srco = any(srcd > 0, 3);
            srci = srci(srco);
            trgi = trgi(srco);
            srcs = size(srcd, 1) * size(srcd, 2);
            trgs = size(emv, 1) * size(emv, 2);
            emv(trgi) = srcd(srci);
            emv(trgi + trgs) = srcd(srci + srcs);
            emv(trgi + 2 * trgs) = srcd(srci + 2 * srcs);
        end
    catch ne_eo;
        if ~isempty(pbar)
            closebar(pbar);
        end
        close(ao);
        delete(aoo{1});
        rethrow(ne_eo);
    end

    % cut frame if needed
    if size(emv, 1) > vsz(2)
        emv(vsz(2)+1:end, :, :) = [];
    end
    if size(emv, 2) > vsz(1)
        emv(:, vsz(1)+1:end, :) = [];
    end

    % add frame to video
    ao = addframe(ao, emv);
    if ~isempty(pbar)
        xprogress(pbar, ic);
    end

    % clear vid frame again for next use
    emv(:, :, 1) = bgc(1);
    emv(:, :, 2) = bgc(2);
    emv(:, :, 3) = bgc(3);
end

% close progress bar
if ~isempty(pbar)
    closebar(pbar);
end

% close video
try
    close(ao);
catch ne_eo;
    rethrow(ne_eo);
end
