function r = renderlayers(tio)
% renderlayers  - rendering of transimg layers into a HxWx3 uint8 image
%
% FORMAT:       r = renderlayers(tio)
%
% Input fields:
%
%       tio         1x1 struct of transimg content
%
% Output fields:
%
%       r           HxWx3 uint8 image
%
% Note: this has an alternative MEX (compiled) function implementation!

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

% Matlab implementation in case MEX file is missing
if nargin ~= 1 || ...
    numel(tio) ~= 1 || ...
   ~isstruct(tio) || ...
   ~isfield(tio, 'Layer') || ...
   ~isfield(tio, 'Background') || ...
   ~isfield(tio, 'Height') || ...
   ~isfield(tio, 'Width') || ...
   ~isstruct(tio.Layer) || ...
   ~isfield(tio.Layer, 'Pixel') || ...
   ~isfield(tio.Layer, 'Alpha')
    error( ...
        'Bad or missing input argument.' ...
    );
end

% get settings
h = tio.Height;
w = tio.Width;

% create output buffer
r = single(0);
r(h, w, 3) = 0;
r(:, :, 1) = tio.Background(1);
r(:, :, 2) = tio.Background(2);
r(:, :, 3) = tio.Background(3);

% iterate over layers
for lc = 1:numel(tio.Layer)

    % get alpha
    a = single(tio.Layer(lc).Alpha);

    % nothing to do
    if numel(a) == 1 && ...
        a == 0
        continue;
    end

    % get pixel data
    p = single(tio.Layer(lc).Pixel);

    % extend pixel
    if size(p, 3) == 1
        p = p(:, :, [1, 1, 1]);
    end

    % extend alpha
    if numel(a) > 1
        a = a(:, :, [1, 1, 1]);
    end

    % add to buffer
    r = (1 - a) .* r + a .* p;
end

% set type
r = uint8(min(single(255), max(single(0), r)));
