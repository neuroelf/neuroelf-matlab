function rim = image_colrotate(im, r, opts)
% image_colrotate  - rotate hue portion in image
%
% FORMAT:       rim = image_colrotate(im, r [, opts])
%
% Input fields:
%
%       im          HxWx3 image data
%       r           rotation spec [ifrom, ito, ofrom, oto]
%       opts        optional settings
%        .grayout   gray out the non-used portion (default: false)
%        .hsv       boolean flag, input/output HSV (default: false)
%
% Output fields:
%
%       rim         color-rotated image data

% Version:  v1.1
% Build:    21072614
% Date:     Jul-26 2021, 2:34 PM EST
% Author:   Jochen Weber, Memorial Sloan Kettering Cancer Center, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2021, Jochen Weber
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
if nargin < 2 || isempty(im) || ~isnumeric(im) || ndims(im) ~= 3 || size(im, 3) ~= 3 || ...
    numel(r) ~= 4 || ~isa(r, 'double') || any(isinf(r(:)) | isnan(r(:)) | r(:) < 0)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
r = r(:);
if any(r > 2)
    r = (1/360) .* r;
end
r = mod(r, 1);
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'grayout') || ~islogical(opts.grayout) || numel(opts.grayout) ~= 1
    opts.grayout = false;
end
if ~isfield(opts, 'hsv') || ~islogical(opts.hsv) || numel(opts.hsv) ~= 1
    opts.hsv = false;
end

% loop required?
if numel(im) > 4e8
    
    % create placeholder
    if opts.hsv
        rim = zeros(size(im));
    else
        rim = uint8(0);
        rim(size(im, 1), size(im, 2), 3) = 0;
    end
    
    % setup columns ranges
    k = max(1, floor(5e7 / size(im, 1)));
    cols = 1:k:size(im, 2);
    if cols(end) == size(im, 2)
        cols(end) = [];
    end
    colx = [cols(2:end)-1, size(im, 2)];
    
    % iterate
    for c = 1:numel(cols)
        rim(:, cols(c):colx(c), :) = image_colrotate(im(:, cols(c):colx(c), :), r, opts);
    end
    
    % early return
    return;
end

% convert to HSV
if ~opts.hsv
    im = hsvconv(im, 2, 1);
end

% prepare rotation arguments
if r(2) <= r(1)
    r(2) = r(2) + 1;
end
if r(4) <= r(3)
    r(4) = r(4) + 1;
end
if r(3) < r(1)
    r(3:4) = r(3:4) + 1;
end

% read and prepare hue
hue = im(:, :, 1);
hue(hue < 0) = 0;
hue360 = false;
if any(hue(:) > 2)
    hue = (1/360) .* hue;
    hue360 = true;
else
    hue = mod(hue, 1);
end

% make sure all hue is above cutoff
msk = hue < r(1);
hue(msk) = hue(msk) + 1;
hue = hue - r(1);
r = r - r(1);

% gray out saturation if hue > threshold
if opts.grayout
    msk = hue > r(2);
    sat = im(:, :, 2);
    sat(msk) = 0;
    im(:, :, 2) = sat;
end

% apply rotation
hd1 = mod(r(2) - r(1), 1);
if hd1 == 0
    hd1 = 1;
end
hd2 = mod(r(4) - r(3), 1);
if hd2 == 0
    hd2 = 1;
end
hue = r(3) + (hd2 / hd1) .* hue;

% if not grayed out, rotate inverse portion
if ~opts.grayout
    hd1r = 1 - hd1;
    hd2r = 1 - hd2;
    if hd2r == 0
        hd2r = 1e-9;
    end
    msk = hue > r(4);
    hue(msk) = r(4) + (hd2r / hd1r) .* hue(msk);
end

% fix hue back into [0..1] range
hue = mod(hue, 1);

% replace hue
im(:, :, 1) = hue;

% and convert back if needed
if ~opts.hsv
    rim = hsvconv(im, 1, 1);
else
    if hue360
        im(:, :, 1) = 360 .* im(:, :, 1);
    end
    rim = im;
end
