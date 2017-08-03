function ir = image_resize(im, h, w)
% image_resize  - resize a HxWxD numeric image to new measures
%
% FORMAT:       ir = image_resize(im, h [, w])
%
% Input fields:
%
%       im          HxWxD image (height x width x color depth)
%       h           new height (or maximum dimension)
%       w           new width
%
% Output fields:
%
%       ir          resized image
%
% Note: if both h and w are given, one of the two can be set
%       to 0 which will lead to auto-detection of the second
%       parameter; alternatively only h can be given which will
%       set the larger of the two dimensions to this size

% Version:  v1.1
% Build:    16051917
% Date:     May-19 2016, 5:25 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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
if nargin < 2 || ~isnumeric(im) || isempty(im) || ...
    numel(h) ~= 1 || ~isa(h, 'double') || isinf(h) || isnan(h) || ...
    h <= 0 || (h > 1 && h ~= fix(h)) || ...
   (nargin < 3 && h > 1 && h < 4)
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
imc = class(im);
if nargin < 3
    maxsize = true;
    w = 0;
else
    maxsize = false;
end
if h < 1
    if nargin < 3
        h = round(h * max(size(im, 1), size(im, 2)));
    else
        h = round(h * size(im, 1));
    end
end
if w < 1 && w > 0
    w = w * size(im, 2);
end
if numel(w) ~= 1 || ~isa(w, 'double') || isinf(w) || isnan(w) || ...
    w < 0 || w ~= fix(w) || all([h, w] < 4)
    error('neuroelf:general:badArgument', 'Bad argument.');
end

% get image dims
ims = size(im);
if numel(ims) < 3
    ims(3) = 1;
end

% get requested size
if maxsize
    [ms, msp] = max(ims(1:2));
    if msp == 1
        newheight = h;
        newwidth = round(h * ims(2) / ms);
    else
        newheight = round(h * ims(1) / ms);
        newwidth = h;
    end
else
    newheight = h;
    newwidth = w;
    if h == 0
        newheight = round(w * ims(1) / ims(2));
    elseif w == 0
        newwidth = round(h * ims(2) / ims(1));
    end
end

% get interpolation points
imht = ims(1) / newheight;
imhs = 0.5 * (imht + 1);
imhe = (ims(1) + 1) - imhs;
imwt = ims(2) / newwidth;
imws = 0.5 * (imwt + 1);
imwe = (ims(2) + 1) - imws;

% get interpolation kernels
[nulld, hk] = resampleaa([0; 0], imht);
if imwt ~= imht
    [nulld, wk] = resampleaa([0; 0], imwt);
else
    wk = hk;
end

% interpolate
ir = flexinterpn(im, ...
    [Inf, Inf, Inf; imhs, imws, 1; imht, imwt, 1; imhe, imwe, ims(3)], ...
    {hk{1}, wk{1}, [0; 1; 0]}, {hk{2}, wk{2}, 1}, 0);

% make class correct
switch (lower(imc))
    case 'uint8'
        ir = uint8(fix(max(min(ir, 255), 0)));
    case 'int16'
        ir = int16(fix(max(min(ir, 32767), -32768)));
    case 'int32'
        ir = int32(fix(max(min(ir, 2^31 - 1), -(2^31))));
    case 'int8'
        ir = int8(fix(max(min(ir, 127), -128)));
    case 'single'
        ir = single(ir);
    case 'uint16'
        ir = uint16(fix(max(min(ir, 65535), 0)));
    case 'uint32'
        ir = uint32(fix(max(min(ir, 2^32 - 1), 0)));
end
