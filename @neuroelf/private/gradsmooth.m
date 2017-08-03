function d = gradsmooth(d)
% gradsmooth  - apply gradient based smoothing
%
% FORMAT:       v = gradsmooth(v)
%
% Input fields:
%
%       v           3D volumetric data
%
% Output fields:
%
%       v           data with gradient based smoothing applied

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
   ~isnumeric(d) || ...
    ndims(d) ~= 3
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing input argument.' ...
    );
end

% go to double if necessary
if ~isa(d, 'double')
    d = double(d);
end

% create three convolution operators
cox = zeros(3, 3, 3);
coy = zeros(3, 3, 3);
coz = zeros(3, 3, 3);
cox([1, 3], 2, 2) = [-0.5, 0.5];
coy(2, [1, 3], 2) = [-0.5, 0.5];
coz(2, 2, [1, 3]) = [-0.5, 0.5];

% compute X-gradient
g = conv3d(d, cox, 2);
gr = g .* g;

% compute Y-gradient
g = conv3d(d, coy, 2);
g = g .* g;
gr = gr + g;

% compute Z-gradient
g = conv3d(d, coz, 2);
g = g .* g;
gr = gr + g;
gr = sqrt(gr);

% compute direct X-gradient
cox(:, 2, 2) = [-0.25, 0.25, 0];
g = abs(conv3d(d, cox, 2));
cox(:, 2, 2) = [0, -0.25, 0.25];
g = g + abs(conv3d(d, cox, 2));
gd = g .* g;

% compute direct Y-gradient
coy(2, :, 2) = [-0.25, 0.25, 0];
g = abs(conv3d(d, coy, 2));
coy(2, :, 2) = [0, -0.25, 0.25];
g = g + abs(conv3d(d, coy, 2));
g = g .* g;
gd = gd + g;

% compute direct Z-gradient
coz(2, 2, :) = [-0.25, 0.25, 0];
g = abs(conv3d(d, coz, 2));
coz(2, 2, :) = [0, -0.25, 0.25];
g = g + abs(conv3d(d, coz, 2));
g = g .* g;
gd = gd + g;

% compute ratio
mgr = minmaxmean(gr);
mgd = minmaxmean(gd);
gr = (1 / mgr(2)) .* gr;
gd = (1 / mgd(2)) .* gd;
gr = gr + 1;
gd = gd + 1;
gd = gd ./ gr;
mgd = minmaxmean(gd);
gd = gd - mgd(1);
gd = (1 / (mgd(2) - mgd(1))) .* gd;
gr = -gd;
gr = gr + 1;

% get smoothed data
sd = smoothdata3(d, [2, 2, 2]);

% compute sum
d = gd .* sd + gr .* d;
