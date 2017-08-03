function g = gradmag(x, m)
% gradmag  - gradient magnitude of 3D data
%
% FORMAT:       g = gradmag(x [, imeth])
%
% Input fields:
%
%       x           3D data to compute gradient over
%       imeth       interpolation method (for sub-index scanning)
%
% Output fields:
%
%       g           gradient

% Version:  v0.9c
% Build:    13020213
% Date:     Feb-02 2013, 1:10 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, 2013, Jochen Weber
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
   ~isnumeric(x) || ...
    ndims(x) ~= 3
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~ischar(m) || ...
    isempty(m) || ...
    isempty(regexpi(m(:)', '^(cubic|lanczos[2-9]|linear|nearest|spline[23])$'))
    m = 'lanczos3';
else
    m = lower(m(:)');
end

% resolve data if necessary
if istransio(x)
    x = x(:, :, :);
end
xs = size(x);

% interpolation arguments
fi = [Inf, Inf, Inf; 1, 1, 1; 1, 1, 1; xs];
xi = zeros(4, 3);
xi([2, 4], 1) = 0.25;
yi = zeros(4, 3);
yi([2, 4], 2) = 0.25;
zi = zeros(4, 3);
zi([2, 4], 3) = 0.25;

% compute gradient magnitude
g = 4 .* sqrt( ...
    (flexinterpn_method(x, fi - xi, m) - flexinterpn_method(x, fi + xi, m)) .^ 2 + ...
    (flexinterpn_method(x, fi - yi, m) - flexinterpn_method(x, fi + yi, m)) .^ 2 + ...
    (flexinterpn_method(x, fi - zi, m) - flexinterpn_method(x, fi + zi, m)) .^ 2);
