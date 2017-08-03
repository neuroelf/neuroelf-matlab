function ti = setlayeralpha(ti, l, a)
% transimg::setlayeralpha  - set alpha channel of a layer of an image
%
% Using: limitrangec.

% Version:  v0.9d
% Build:    14082216
% Date:     Aug-22 2014, 4:24 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2012, 2014, Jochen Weber
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
if nargin < 3 || ...
   ~isa(l, 'double') || ...
    numel(l) ~= 1 || ...
    isinf(l) || ...
    isnan(l) || ...
    l < 1 || ...
    l > numel(tiobj(lup).Layer) || ...
    l ~= fix(l) || ...
   (~isa(a, 'double') && ...
    ~isa(a, 'single') && ...
    ~isa(a, 'uint8')) || ...
   (numel(a) ~= 1 && ...
    (ndims(a) ~= 2 || ...
     (l <= numel(tiobj(lup).Layer) && ...
      (size(a, 1) ~= size(tiobj(lup).Layer(l).Pixel, 1) || ...
       size(a, 2) ~= size(tiobj(lup).Layer(l).Pixel, 2))) || ...
     (l > numel(tiobj(lup).Layer) && ...
      (size(a, 1) ~= tiobj(lup).Height || ...
       size(a, 2) ~= tiobj(lup).Width))))
    error( ...
        'transimg:InvalidCall', ...
        'Invalid call to transimg::addlayer.' ...
    );
end
if isa(a, 'uint8')
    a = single(1 / 255) .* single(a);
else
    a = ne_methods.limitrangec(single(a), 0, 1, 0);
end
l = floor(real(l));

% reset IsRendered flag
tiobj(lup).IsRendered = false;

% re-render layer as well
if tiobj(lup).Layer(l).Type ~= 'f'
    tiobj(lup).Layer(l).IsRendered = false;
end

% set layer alpha
tiobj(lup).Layer(l).Alpha = a;
