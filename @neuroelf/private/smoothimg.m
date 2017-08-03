function s = smoothimg(d, f, t)
% smoothimg  - apply gaussian smoothing to 2D data (layers, over time)
%
% FORMAT:       s = smoothimg(d, f [, t])
%
% Input fields:
%
%       d           2D/3D data (for >2 dim, every plane is smoothed)
%       f           FWHM kernel sizes (1x1 or 1x2 double)
%       t           smoothing kernel threshold (default: 0.001)
%
% Output fields:
%
%       s           smoothed data

% Version:  v0.9c
% Build:    11092413
% Date:     Sep-24 2011, 1:46 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
    ~isnumeric(d) || ...
   ~isa(f, 'double') || ...
   ~any(numel(f) == [1, 2]) || ...
    any(isinf(f) | isnan(f) | f < 0 | f > 32)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
sd = size(d);
sd2 = sd(1:2);
sd3 = floor(0.5 + numel(d) / prod(sd2));
fd2 = [Inf, Inf; 1, 1; 1, 1; sd2];
if nargin < 3 || ...
   ~isa(t, 'double') || ...
    numel(t) ~= 1 || ...
    isinf(t) || ...
    isnan(t) || ...
    t < 0
    t = 0.001;
elseif t > 0.25
    t = 0.25;
end
if numel(f) == 1
    f = f([1, 1]);
end

% build 2 separate smoothing kernels for each dim
k1 = smoothkern(f(1), t);
k2 = smoothkern(f(2), t);

% create double-based copy of data (also resolves transio)
s = double(d);

% iterate over third+ dim
for d3 = 1:sd3
    s(:, :, d3) = flexinterpn(flexinterpn(d(:, :, d3), ...
        fd2, {[0; 1; 0]; k2}, {1, 1}, 0), fd2, {k1, [0; 1; 0]}, {1, 1}, 0);
end
