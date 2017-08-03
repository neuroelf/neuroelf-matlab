function is = image_smooth(im, sk)
% image_smooth  - smooth an image with a gaussian kernel
%
% FORMAT:       is = image_smooth(im, sk)
%
% Input fields:
%
%       im          HxWxD image (height x width x color depth)
%       sk          1x1 or 1x2 smoothing kernel(s) in pixel
%
% Output fields:
%
%       is          smoothed image

% Version:  v1.0
% Build:    14121621
% Date:     Dec-16 2014, 9:45 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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
   ~isnumeric(im) || ...
    isempty(im) || ...
   ~isa(sk, 'double') || ...
    isempty(sk) || ...
    numel(sk) > 2 || ...
    any(isinf(sk) | isnan(sk) | sk < 0)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
skx = smoothkern(sk(1), 0, true, 'linear');
skx(skx < 0.004 .* max(skx)) = [];
if numel(sk) == 1
    sky = skx;
else
    sky = smoothkern(sk(2), 0, true, 'linear');
end
skn = [0; 1; 0];

% get image dims
imc = class(im);
ims = size(im);
if numel(ims) < 3
    ims(3) = 1;
end
imi = [Inf; 1; 1; 1] * ones(1, numel(ims));
imi(end, :) = ims;
if numel(ims) > 2
    sk3 = {skn};
else
    sk3 = {};
end

% interpolate
is = flexinterpn(flexinterpn(im, imi, {skx, skn, sk3{:}}, {1, 1, 1}, 0), ...
    imi, {skn, sky, sk3{:}}, {1, 1, 1}, 0);

% make class correct
switch (lower(imc))
    case {'int16'}
        is = int16(fix(max(min(is, 32767), -32768)));
    case {'int32'}
        is = int32(fix(max(min(is, 2^31 - 1), -(2^31))));
    case {'int8'}
        is = int8(fix(max(min(is, 127), -128)));
    case {'single'}
        is = single(is);
    case {'uint16'}
        is = uint16(fix(max(min(is, 65535), 0)));
    case {'uint32'}
        is = uint32(fix(max(min(is, 2^32 - 1), 0)));
    case {'uint8'}
        is = uint8(fix(max(min(is, 255), 0)));
end
