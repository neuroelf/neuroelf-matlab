function ti = setlayerpixel(ti, l, p)
% transimg::setlayerpixel  - set the pixels of one layer of an image
%
% Using: limitrangec.

% Version:  v0.9d
% Build:    14082216
% Date:     Aug-22 2014, 4:24 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2012, 2014, Jochen Weber
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

% global variables for storage
global tiobj ...
       tiobjlup ...
       ne_methods;

% check arguments
lup = find(tiobjlup == ti.L);
if numel(lup) ~= 1
    error( ...
        'transimg:ObjectRemoved', ...
        'Object removed from global storage.' ...
    );
end
tio = tiobj(lup);
if nargin < 3 || ...
   ~isa(l, 'double') || ...
    numel(l) ~= 1 || ...
    isinf(l) || ...
    isnan(l) || ...
    l < 1 || ...
    l > (numel(tio.Layer) + 1) || ...
    l ~= fix(l) || ...
   (~isa(p, 'uint8') && ...
    ~isa(p, 'double') && ...
    ~isa(p, 'single') && ...
    ~isa(p, 'transio')) || ...
   (l <= numel(tio.Layer) && ...
    (size(p, 1) ~= size(tio.Layer(l).Pixel, 1) || ...
     size(p, 2) ~= size(tio.Layer(l).Pixel, 2))) || ...
   (l > numel(tio.Layer) && ...
    (size(p, 1) ~= tio.Height || ...
     size(p, 2) ~= tio.Width)) || ...
   ~any([1, 3] == size(p, 3))
    error( ...
        'transimg:InvalidCall', ...
        'Invalid call to transimg::addlayer.' ...
    );
end
l = floor(real(l));

% convert pixel data if necessary
if ~isa(p, 'uint8')
    p = ne_methods.limitrangec(single(p(:, :, :)), 0, 255, 0);
end

% reset IsRendered flag
tio.IsRendered = false;

% add layer ?
if l > numel(tio.Layer)
    tio.Layer(l).Type = 'f';
    tio.Layer(l).Alpha = single(1);
end

% set pixel data
tio.Layer(l).Pixel = p;

% re-render layer as well
if tio.Layer(l).Type ~= 'f'
    tio.Layer(l).IsRendered = false;
end

% set back in global storage
tiobj(lup) = tio;
