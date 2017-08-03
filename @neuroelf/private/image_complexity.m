function cmi = image_complexity(im, isz)
% image_complexity  - give image complexity estimate
%
% FORMAT:       cmi = image_complexity(im [, isz])
%
% Input fields:
%
%       im          HxWxC image (where C is either 1 for gray or 3 for RGB!)
%       isz         resampling size for the estimation of complexity
%                   default is [64, 64]
%
% Output fields:
%
%       cmi         complexity estimate

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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
   ~isnumeric(im) || ...
    isempty(im)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument.' ...
    );
end
im = double(im);
if nargin < 2 || ...
   ~isa(isz, 'double') || ...
    numel(isz) ~= 2 || ...
    any(isinf(isz) | isnan(isz) | isz < 4 | isz > 1024)
    isz = [64, 64];
else
    isz = round(isz(:)');
end

% grayscale image
if ndims(im) == 2

    % compute gradient of image
    [gx, gy] = gradient(im);
    gr = sqrt(gx .* gx + gy .* gy);

    % resample image to isz and back
    nim = image_resize(image_resize(im, isz(1), isz(2)), size(im, 1), size(im, 2));

    % recompute gradient
    [gx, gy] = gradient(nim);
    ngr = sqrt(gx .* gx + gy .* gy);

    % compute estimate
    cmi = mean(max(0, noinfnan((gr(:) - ngr(:)) ./ gr(:))));

% RGB image
else

    % make sure to comply with RGB
    mmm = minmaxmean(im);
    uim = uint8(floor(255.999 .* ((1 / (mmm(2) - mmm(1)) .* (im - mmm(1))))));

    % compute gray scale version first
    gim = double(rgb2gray(uim));

    % compute gradient of image
    [gx, gy] = gradient(gim);
    gr = sqrt(gx .* gx + gy .* gy);

    % resample image to isz and back
    nim = image_resize(image_resize(gim, isz(1), isz(2)), size(im, 1), size(im, 2));

    % recompute gradient
    [gx, gy] = gradient(nim);
    ngr = sqrt(gx .* gx + gy .* gy);
    dgr = max(0, noinfnan((gr(:) - ngr(:)) ./ gr(:)));

    % perform the same for each color component
    [gx, gy] = gradient(im(:, :, 1));
    cgr = sqrt(gx .* gx + gy .* gy);

    % resample image to isz and back
    nim = image_resize(image_resize(im(:, :, 1), isz(1), isz(2)), size(im, 1), size(im, 2));

    % recompute gradient
    [gx, gy] = gradient(nim);
    ngr = sqrt(gx .* gx + gy .* gy);
    dcgr = max(0, noinfnan((cgr(:) - ngr(:)) ./ cgr(:)));

    % 2nd component
    [gx, gy] = gradient(im(:, :, 2));
    cgr = sqrt(gx .* gx + gy .* gy);

    % resample image to isz and back
    nim = image_resize(image_resize(im(:, :, 2), isz(1), isz(2)), size(im, 1), size(im, 2));

    % recompute gradient
    [gx, gy] = gradient(nim);
    ngr = sqrt(gx .* gx + gy .* gy);
    dcgr = max(dcgr, max(0, noinfnan((cgr(:) - ngr(:)) ./ cgr(:))));

    % 3rd component
    [gx, gy] = gradient(im(:, :, 3));
    cgr = sqrt(gx .* gx + gy .* gy);

    % resample image to isz and back
    nim = image_resize(image_resize(im(:, :, 3), isz(1), isz(2)), size(im, 1), size(im, 2));

    % recompute gradient
    [gx, gy] = gradient(nim);
    ngr = sqrt(gx .* gx + gy .* gy);
    dcgr = max(dcgr, max(0, noinfnan((cgr(:) - ngr(:)) ./ cgr(:))));

    % mixed estimate
    cmi = (1 / (size(im, 1) * size(im, 2))) * ...
        (0.5 * sum(dgr(:)) + 0.5 * sum(dcgr(:)));
end
