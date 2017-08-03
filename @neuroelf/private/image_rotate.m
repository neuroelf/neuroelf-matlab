function [ir, ira] = image_rotate(im, r, opts)
% image_resize  - resize a HxWxD numeric image to new measures
%
% FORMAT:       [ir, ira] = image_resize(im, r [, opts])
%
% Input fields:
%
%       im          HxWxD image or image filename
%       r           rotation angle
%       opts        optional settings
%        .around    1x2 [x,y] around which rotation is done (default center)
%        .bgcolor   1xD background information or either 'mean' or {'mode'}
%        .maskpost  post-rotational mask (HxW [0 .. 1] scaled data)
%        .maskpre   pre-rotational mask (HxW [0 .. 1] scaled data)
%        .method    interpolation method (see flexinterpn_method)
%        .outfile   output file (either fixed name or pattern with one %)
%        .quality   if output format is JPG, quality setting
%        .resize    flag, if true, make larger to fit (default: false)
%
% Output fields:
%
%       ir          rotated image data
%       ira         rotated alpha channel (if in image)

% Version:  v0.9b
% Build:    11110216
% Date:     Apr-09 2011, 2:11 PM EST
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
if nargin < 2 || ...
    isempty(im) || ...
   (~isnumeric(im) && ...
    ~ischar(im)) || ...
    numel(r) ~= 1 || ...
   ~isa(r, 'double') || ...
    isinf(r) || ...
    isnan(r) || ...
    r < -360 || ...
    r > 360
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
fromfile = '';
if ischar(im)
    fromfile = im(:)';
    try
        [im, imap, ima] = imread(fromfile);
    catch ne_eo;
        rethrow(ne_eo);
    end
    im = single(im);
    if ~isempty(ima)
        im(:, :, end + 1) = ima;
    end
else
    im = single(im);
end
ims = size(im);
if numel(ims) < 3
    ims(3) = 1;
end
imp = ims(3);
ims(3) = [];
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'around') || ...
   ~isa(opts.around, 'double') || ...
    numel(opts.around) ~= 2 || ...
    any(isinf(opts.around) | isnan(opts.around) | ...
        opts.around < -ims | opts.around > 2 * ims)
    opts.around = 0.5 .* (ims + 1);
elseif all(opts.around > 0 & opts.around < 1)
    opts.around = opts.around .* ims;
end
if ~isfield(opts, 'bgcolor') || ...
    isempty(opts.bgcolor) || ...
   (~isa(opts.bgcolor, 'double') && ...
    ~ischar(opts.bgcolor))
    opts.bgcolor = 'mode';
end
if ischar(opts.bgcolor)
    bgrgb = [im(2:end-1, 1, :);  im(2:end-1, end, :); ...
           permute(im(1, 2:end-1, :), [2, 1, 3]); ...
           permute(im(end, 2:end-1, :), [2, 1, 3])];
    if strcmpi(opts.bgcolor(:)', 'mode')
        if imp > 2
            bgrgb = sum((ones(size(bgrgb, 1), 1) * [65536, 256, 1]) .* ...
                reshape(round(bgrgb(:, :, 1:3)), size(bgrgb, 1), 3), 2);
            rgbmode = mode(bgrgb);
            rgb = [floor(rgbmode / 65536), ...
                floor(mod(rgbmode, 65536) / 256), mod(rgbmode, 256)];
        else
            bgrgb = round(bgrgb(:, 1, 1));
            rgb = mode(bgrgb);
        end
    else
        if imp > 2
            rgb = single(round(mean(bgrgb(:, :, 1:3))));
        else
            rgb = single(round(mean(bgrgb(:, :, 1))));
        end
    end
else
    rgb = opts.bgcolor(:)';
end
if (imp > 2 && ...
    numel(rgb) ~= 3) || ...
   (imp < 3 && ...
    numel(rgb) ~= 1)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid bgcolor option.' ...
    );
end
if isfield(opts, 'maskpost') && ...
    ischar(opts.maskpost) && ...
    strcmpi(opts.maskpost(:)', 'ellipse')
    ctr = 0.5 .* (ims + 1);
    [er1, er2] = ndgrid(1:ims(1), 1:ims(2));
    er1 = [(1 / ims(1)) .* (er1(:) - ctr(1)), (1 / ims(2)) .* (er2(:) - ctr(2))];
    er2 = (min(ims) / max(ims)) * 0.5 - 1 / min(ims);
    opts.maskpost = flexinterpn(single(reshape(sqrt(sum(er1 .* er1, 2)) < er2, ims)), ...
        [Inf, Inf; 1, 1; 1, 1; ims], smoothkern(2, 0.001), 1, 0);
end
if ~isfield(opts, 'maskpost') || ...
   ~isnumeric(opts.maskpost) || ...
   ~isequal(size(opts.maskpost), ims)
    opts.maskpost = [];
end
if isfield(opts, 'maskpre') && ...
    ischar(opts.maskpre) && ...
    strcmpi(opts.maskpre(:)', 'ellipse')
    ctr = 0.5 .* (ims + 1);
    [er1, er2] = ndgrid(1:ims(1), 1:ims(2));
    er1 = [(1 / ims(1)) .* (er1(:) - ctr(1)), (1 / ims(2)) .* (er2(:) - ctr(2))];
    er2 = (min(ims) / max(ims)) * 0.5 - 1 / min(ims);
    opts.maskpre = flexinterpn(single(reshape(sqrt(sum(er1 .* er1, 2)) < er2, ims)), ...
        [Inf, Inf; 1, 1; 1, 1; ims], smoothkern(2, 0.001), 1, 0);
end
if ~isfield(opts, 'maskpre') || ...
   ~isnumeric(opts.maskpre) || ...
   ~isequal(size(opts.maskpre), ims)
    opts.maskpre = [];
end
if ~isfield(opts, 'method') || ...
   ~ischar(opts.method)
    opts.method = 'cubic';
else
    opts.method = lower(opts.method(:)');
end
if ~isfield(opts, 'outfile') || ...
   ~ischar(opts.outfile) || ...
    numel(opts.outfile) < 5 || ...
   ~any(opts.outfile(:) == '.')
    opts.outfile = '';
else
    opts.outfile = opts.outfile(:)';
    ofp = find(opts.outfile == '%');
    if numel(ofp) > 1 || ...
       (~isempty(strfind(opts.outfile, '%s')) && ...
         isempty(fromfile))
        error( ...
            'neuroelf:BadArgument', ...
            'Bad or missing outfile option.' ...
        );
    end
    if ~isempty(ofp) && ...
        isempty(regexpi(opts.outfile(ofp+1:end), '^([+\-]?\d+[\.0-9]*)?[dfgs]'))
        error( ...
            'neuroelf:BadArgument', ...
            'Bad or missing outfile option.' ...
        );
    end
end
if ~isfield(opts, 'quality') || ...
   ~isa(opts.quality, 'double') || ...
    numel(opts.quality) ~= 1 || ...
    isinf(opts.quality) || ...
    isnan(opts.quality) || ...
    opts.quality < 25
    opts.quality = 80;
else
    opts.quality = min(100, opts.quality);
end
if ~isfield(opts, 'resize') || ...
   ~islogical(opts.resize) || ...
    numel(opts.resize) ~= 1
    opts.resize = false;
end

% pre-mask
if ~isempty(opts.maskpre)
    mp = opts.maskpre;
    if ~islogical(mp)
        mp = single(limitrangec(mp(:, :), 0, 1, 0));
    end
    for imc = 1:imp
        if ~islogical(mp)
            im(:, :, imc) = im(:, :, imc) .* mp + rgb(imc) .* (1 - mp);
        else
            im(:, :, imc) = im(:, :, imc) .* mp + rgb(imc) .* (~mp);
        end
    end
end

% permute to have last dim first
im = permute(im, [3, 1, 2]);

% resize?
if opts.resize

    % implementation later
    ir = [];
    return;
end

% build rotation matrix
trf = spmtrf([0, opts.around]) * spmtrf([0, 0, 0], [pi * r / 180, 0, 0]) * ...
    spmtrf([0, - opts.around]);

% rotate each plane
ir = cell(1, imp);
for imc = 1:imp
    ir{imc} = flexinterpn_method(im(imc, :, :), ...
        [Inf, Inf, Inf; 1, 1, 1; 1, 1, 1; 1, ims], ...
        double(rgb(imp)), trf, opts.method);
end

% concatenate back in first dim
ir = cat(1, ir{:});

% and re-permute
ir = permute(ir, [2, 3, 1]);

% post-mask
if ~isempty(opts.maskpost)
    mp = opts.maskpost;
    if ~islogical(mp)
        mp = single(limitrangec(mp(:, :), 0, 1, 0));
    end
    for imc = 1:imp
        if ~islogical(mp)
            ir(:, :, imc) = ir(:, :, imc) .* mp + rgb(imc) .* (1 - mp);
        else
            ir(:, :, imc) = ir(:, :, imc) .* mp + rgb(imc) .* (~mp);
        end
    end
end

% and ensure the datatype makes sense
if any(imp == [2, 4])
    ira = double(ir(:, :, end));
    ir = uint8(round(ir(:, :, 1:end-1)));
else
    ira = [];
    ir = uint8(round(ir));
end

% save data
if ~isempty(opts.outfile)

    % extension
    imx = regexp(opts.outfile, '\.([^\.]+)$', 'tokens');
    if isempty(imx) || ...
        isempty(imx{1})
        imx = 'bmp';
    else
        imx = lower(imx{1}{1});
    end

    % pattern
    if ~isempty(ofp)

        % generate filename with old filename
        if ~isempty(regexpi(opts.outfile, '\%([+\-]?\d+[\.0-9]*)?s'))

            % without extension though!
            [fp, fn] = fileparts(fromfile);
            opts.outfile = sprintf(opts.outfile, fn);

        % or degrees
        else
            r = mod(r + 720, 360);
            opts.outfile = sprintf(opts.outfile, r);
        end

    end

    % JPG/JPEG
    if any(strcmp(imx, {'jpeg', 'jpg'}))
        imwrite(ir, opts.outfile, 'Quality', opts.quality);

    % other format
    else

        % transparency?
        if ~isempty(ira)
            imwrite(ir, opts.outfile, 'Alpha', ira);
        else
            imwrite(ir, opts.outfile);
        end
    end
end
