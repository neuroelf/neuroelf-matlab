function rgb = hsvconv(col, dir, scale)
% hsvconv  - convert from and to HSV (RGB) color scheme
%
% FORMAT:       col = hsvconv(col, dir, scale)
%
% Input fields:
%
%       col         Xx3 or X*Y*3 or X*Y*Z*3 data (uint8, single, double)
%       dir         either {1} (from HSV to RGB) or 2 (from RGB to HSV)
%       scale       scaling for hue value, {1} ([0 .. 1]) or 2 ([0 .. 360))
%
% Output fields:
%
%       col         converted data (uint8->double for hsv, uint8 for rgb)

% Version:  v0.9c
% Build:    11042919
% Date:     Apr-29 2011, 8:11 PM EST
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

% allow basic use (HSV->RGB)
if nargin < 1 || ...
   (~isa(col, 'single') && ...
    ~isa(col, 'double'))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isa(dir, 'double') || ...
    numel(dir) ~= 1 || ...
    isinf(dir) || ...
    isnan(dir) || ...
   ~any(dir == [1, 2])
    dir = 1;
end
if nargin < 3 || ...
   ~isa(scale, 'double') || ...
    numel(scale) ~= 1 || ...
    isinf(scale) || ...
    isnan(scale) || ...
   ~any(scale == [1, 2])
    scale = 1;
end
if dir ~= 1 || ...
    scale ~= 1
    error( ...
        'neuroelf:Unsupported', ...
        'This is a compiled function, M-file with limited capabilities.' ...
    );
end

% perform compuation, first, get dim and size
hd = ndims(col);
hs = size(col);
if hs(hd) ~= 3
    error( ...
        'neuroelf:BadArgument', ...
        'Last dimension of col input must be of size 3.' ...
    );
end
hs(hd) = [];
col = reshape(col, prod(hs), 3);

% get hue
h = mod(6 .* col(:, 1), 6);

% produce output buffer
rgb = zeros(prod(hs), 3);

% fill R, G, B, assuming S = 1, V = 1
rgb(:, 1) = max(0, min(1, abs(h - 3) - 1));
rgb(:, 2) = max(0, min(1, 2 - abs(2 - h)));
rgb(:, 3) = max(0, min(1, 2 - abs(4 - h)));

% compute S/V scaling and offset of RGB values
rgs = col(:, 2) .* col(:, 3);
rgo = col(:, 3) - rgs;

% compute full RGB as uint8
rgb = uint8(255 .* reshape( ...
    rgo(:, [1, 1, 1]) + rgs(:, [1, 1, 1]) .* rgb, [hs, 3]));
