function image_matchprops(ims, opts)
% image_matchprops  - match list of images in properties
%
% FORMAT:       image_matchprops(ims [, opts])
%
% Input fields:
%
%       ims         list of image filenames
%       opts        optional settings
%        .gray      grayscaling ('intensity', 'luma', {'none'}, 'value')
%        .hsmooth   histogram smoothing kernel (default: 2.5)
%        .hue       match on H (default: false)
%        .huexsat   interact hue with saturation (improve estimate, true)
%        .outpat    output pattern (default: '%s_matched')
%        .outtype   override output type (default: '', equals keep type)
%        .refimage  single reference image to match (default: generate)
%        .sat       match on S (default: false)
%        .sfreq     spatial frequency/hull (default: false)
%        .size      resize images (default: true, to smallest image size)
%        .type      matching type, one of 'hist', {'meanstd'}
%        .value     match on V/gray-scale lightness (default: true)
%        .writeref  write reference image (first + _ref.png, default: true)
%
% No output fields.
%
% Note: if one of the images is gray-scale, all images will be gray-scaled
%       using the luma-based calculation

% Version:  v0.9c
% Build:    12121219
% Date:     Nov-21 2012, 2:48 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, Jochen Weber
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
if nargin < 1 || ...
   ~iscell(ims) || ...
    isempty(ims)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument.' ...
    );
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'gray') || ...
   ~ischar(opts.gray) || ...
   ~any(strcmpi(opts.gray(:)', {'i', 'intensity', 'l', 'luma', 'n', 'none', 'v', 'value'}))
    opts.gray = 'n';
else
    opts.gray = lower(opts.gray(1));
end
if ~isfield(opts, 'hsmooth') || ...
   ~isa(opts.hsmooth, 'double') || ...
    numel(opts.hsmooth) ~= 1 || ...
    isinf(opts.hsmooth) || ...
    isnan(opts.hsmooth) || ...
    opts.hsmooth < 0
    opts.hsmooth = 2.5;
else
    opts.hsmooth = min(16, opts.hsmooth);
end
if ~isfield(opts, 'hue') || ...
   ~islogical(opts.hue) || ...
    numel(opts.hue) ~= 1
    opts.hue = false;
end
if ~isfield(opts, 'huexsat') || ...
   ~islogical(opts.huexsat) || ...
    numel(opts.huexsat) ~= 1
    opts.huexsat = true;
end
if ~isfield(opts, 'outpat') || ...
   ~ischar(opts.outpat) || ...
    isempty(opts.outpat) || ...
    strcmp(opts.outpat(:)', '%s') || ...
    isempty(strfind(opts.outpat(:)', '%s'))
    opts.outpat = '%s_matched';
else
    opts.outpat = opts.outpat(:)';
end
if ~isfield(opts, 'outtype') || ...
   ~ischar(opts.outtype) || ...
   ~any(strcmpi(opts.outtype(:)', {'', 'bmp', 'jpg', 'png'}))
    opts.outtype = '';
else
    opts.outtype = opts.outtype(:)';
end
if ~isfield(opts, 'refimage') || ...
    isempty(opts.refimage) || ...
   ((~ischar(opts.refimage) || ...
     exist(opts.refimage(:)', 'file') ~= 2) && ...
    ~isa(opts.refimage, 'uint8'))
    opts.refimage = [];
elseif ischar(opts.refimage)
    try
        opts.refimage = imread(opts.refimage(:)');
    catch ne_eo;
        rethrow(ne_eo);
    end
end
refimage = opts.refimage;
if ~isfield(opts, 'sat') || ...
   ~islogical(opts.sat) || ...
    numel(opts.sat) ~= 1
    opts.sat = false;
end
if ~isfield(opts, 'sfreq') || ...
   ~islogical(opts.sfreq) || ...
    numel(opts.sfreq) ~= 1
    opts.sfreq = false;
end
resize = [];
if ~isfield(opts, 'size') || ...
   ((~isa(opts.size, 'double') || ...
     numel(opts.size) ~= 2 || ...
     any(isinf(opts.size) | isnan(opts.size) | opts.size < 1)) && ...
    (~islogical(opts.size) || ...
     numel(opts.size) ~= 1))
    opts.size = true;
elseif isa(opts.size, 'double')
    resize = opts.size;
    opts.size = true;
end
if ~isfield(opts, 'type') || ...
   ~ischar(opts.type) || ...
    isempty(opts.type) || ...
   ~any(strcmpi(opts.type(:)', {'h', 'hist', 'm', 'meanstd'}))
    opts.type = 'm';
else
    opts.type = lower(opts.type(1));
end
if ~isfield(opts, 'value') || ...
   ~islogical(opts.value) || ...
    numel(opts.value) ~= 1
    opts.value = true;
end
if ~isfield(opts, 'writeref') || ...
   ~islogical(opts.writeref) || ...
    numel(opts.writeref) ~= 1
    opts.writeref = true;
end

% input check
ims = ims(:);
nim = numel(ims);
imsz = zeros(nim, 3);
for ic = 1:nim
    if ~ischar(ims{ic}) || ...
        isempty(ims{ic}) || ...
        exist(ims{ic}(:)', 'file') ~= 2
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid image filename for image %d.', ...
            ic ...
        );
    end
    try
        ims{ic} = ims{ic}(:)';
        sz = size(imread(ims{ic}));
        imsz(ic, 1:numel(sz)) = sz;
    catch ne_eo;
        rethrow(ne_eo);
    end
end
if all(imsz(:, 1) == imsz(1, 1)) && ...
    all(imsz(:, 2) == imsz(1, 2))
    opts.size = true;
    resize = imsz(1, 1:2);
end
if any(imsz(:, 3) < 3) && ...
    opts.gray == 'n'
    opts.gray = 'l';
end
if isempty(resize) && ...
    opts.size
    resize = min(imsz(:, 1:2), [], 1);
    presize = sqrt(prod(resize));

    % but keep ratio as constant as possible
    rx = harmmean(imsz(:, 1) ./ imsz(:, 2));
    resize = round([presize * sqrt(rx), presize / sqrt(rx)]);

    % allow format specifics
    if all(all(rem(imsz(:, 1:2) ./ 16, 1) == 0))
        resize = 16 .* round(resize ./ 16);
    elseif all(all(rem(imsz(:, 1:2) ./ 8, 1) == 0))
        resize = 8 .* round(resize ./ 8);
    end
end

% gray-scaling
if opts.gray ~= 'n'

    % for histogram matching
    if opts.type == 'h'

        % histogram smoothing kernel
        hsmooth = smoothkern(3 * opts.hsmooth, 1/4096);

        % generate data for images
        imh = zeros(nim, 765);

    % mean/std matching
    elseif opts.type == 'm'

        % generate data for images
        imh = zeros(nim, 2);
    end

    % for each image
    for ic = 1:nim

        % read in image
        im = imread(ims{ic});

        % does this image need gray-scaling
        if size(im, 3) > 1

            % type of gray-scaling
            switch (opts.gray)

                % intensity
                case {'i'}

                    % average
                    im = mean(double(im(:, :, 1:3)), 3);

                % luma
                case {'l'}

                    % sRGB formula (see wikipedia: HSL_and_HSV)
                    im = 0.21 .* double(im(:, :, 1)) + ...
                         0.72 .* double(im(:, :, 2)) + ...
                         0.07 .* double(im(:, :, 3));

                % value (HSV)
                case {'v'}

                    % HSV conversion
                    im = hsvconv(im, 2);
                    im = 255 .* im(:, :, 3);
            end
        end

        % histogram
        if opts.type == 'h'
            imh(ic, :) = histcount(im(:), 1 / 6, 255, 1 / 3) ./ numel(im);
        elseif opts.type == 'm'
            imh(ic, :) = [mean(im(:)), std(im(:))];
        end
    end

    % average histograms
    if opts.type == 'h'

        % regular average
        imha = (1 / nim) .* sum(imh, 1);

        % and smooth histogram
        imha = flexinterpn(imha(:), [inf; 1; 1; numel(imha)], hsmooth, 1);

    % average means and compute std to match
    elseif opts.type == 'm'

        imha = [mean(imh(:, 1)), (1 / sqrt(nim)) .* sqrt(sum(imh(:, 2) .* imh(:, 2)))];
    end

    % match each image
    for ic = 1:nim

        % read image
        im = imread(ims{ic});

        % generate patterns
        [impath, imfile, imext] = fileparts(ims{ic});

        % resize?
        if opts.size && ...
           (resize(1) ~= size(im, 1) || ...
            resize(2) ~= size(im, 2))
            im = image_resize(im, resize(1), resize(2));
        end

        % does this image need gray-scaling, follow same path
        if size(im, 3) > 1
            switch (opts.gray)
                case {'i'}
                    im = mean(double(im(:, :, 1:3)), 3);
                case {'l'}
                    im = 0.21 .* double(im(:, :, 1)) + ...
                         0.72 .* double(im(:, :, 2)) + ...
                         0.07 .* double(im(:, :, 3));
                case {'v'}
                    im = hsvconv(im, 2);
                    im = 255 .* im(:, :, 3);
            end
        end

        % generate matched image
        if opts.type == 'h'

            % using histscale for histogram-based matching
            im = histscale(imha, im, struct('bins', 0:1/3:255.001));
        else

            % recalculating image pixel values otherwise
            im = limitrangec( ...
                imha(2) .* ((im - imh(ic, 1)) ./ imh(ic, 2)) + imha(1), ...
                0, 255, 0);
        end

        % write reference image
        if opts.writeref
            if ic == 1
                newref = double(im);
            else
                newref = newref + double(im);
            end
        end

        % write image
        imf = sprintf('%s/%s%s', impath, sprintf(opts.outpat, imfile), imext);
        disp(sprintf('Writing out ''%s''...', imf));
        pause(0.01);
        imwrite(uint8(im), imf);
    end

% in color
else

    % for histogram matching
    if opts.type == 'h'

        % histogram smoothing kernel
        hsmooth = smoothkern(3 * opts.hsmooth, 1/4096);

        % generate data for images
        imh = zeros(nim, 765, 3);

    % mean/std matching
    elseif opts.type == 'm'

        % generate data for images
        imh = zeros(nim, 2, 3);
    end

    % for each image
    for ic = 1:nim

        % read in image
        im = imread(ims{ic});

        % reject gray-scale images for hue-based matching
        if size(im, 3) < 3 && ...
            opts.hue
            error( ...
                'neuroelf:BadArgument', ...
                'Color images required for hue-based matching.' ...
            );
        elseif size(im, 3) < 3
            im = im(:, :, [1, 1, 1]);
        end

        % compute HSV from RGB
        im = hsvconv(im, 2, 1);
        nim = size(im, 1) * size(im, 2);

        % histogram
        if opts.type == 'h'

            % hue
            if opts.hue
                if opts.sat && ...
                    opts.huexsat
                    imh(ic, :, 1) = histcount(lsqueeze(im(:, :, 1)), ...
                        1 / 1530, 1, 1 / 765, lsqueeze(im(:, :, 2)), []) ./ ...
                        sum(lsqueeze(im(:, :, 2)));
                else
                    imh(ic, :, 1) = histcount(lsqueeze(im(:, :, 1)), ...
                        1 / 1530, 1, 1 / 765) ./ nim;
                end
            end
            if opts.sat
                imh(ic, :, 2) = histcount(lsqueeze(im(:, :, 2)), ...
                    1 / 1530, 1, 1 / 765) ./ nim;
            end
            if opts.value
                imh(ic, :, 3) = histcount(lsqueeze(im(:, :, 3)), ...
                    1 / 1530, 1, 1 / 765) ./ nim;
            end

        % means
        elseif opts.type == 'm'
            imh(ic, :) = [mean(im(:)), std(im(:))];
        end
    end

    % error out for now
    error( ...
        'neuroelf:NotYetImplemented', ...
        'HSV-based matching not yet implemented.' ...
    );
end
