function [ex, ey, a, b, t] = cellipse(x, y, opts)
% cellipse  - compute confidence-ellipse of bivariate data
%
% FORMAT:       [ex, ey, rx, ry, t] = cellipse(x, y [, opts])
%
% Input fields:
%
%       x           X-variate
%       y           Y-variate
%       opts        optional settings
%        .dim       dim along to compute (default: first non-singleton)
%        .nonull    reject 0-values in x and y (default: false)
%        .npts      number of points for the ellipse (default: 721)
%        .pa        p-value (bivariate confidence, default: 0.05)
%
% Output fields:
%
%       ex          ellipse X-coordinates (variate 1)
%       ey          ellipse Y-coordinates (variate 2)
%       rx          ellipse X-radius (1-SD length)
%       ry          ellipse Y-radius (1-SD length)
%       t           theta angle of ellipse

% Version:  v0.9b
% Build:    11040113
% Date:     Apr-01 2011, 1:21 PM EST
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
    numel(x) < 2 || ...
   (~isa(x, 'double') && ...
    ~isa(x, 'single')) || ...
   (~isa(y, 'double') && ...
    ~isa(y, 'single')) || ...
   ~isequal(size(x), size(y))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dim') || ...
   ~isa(opts.dim, 'double') || ...
    numel(opts.dim) ~= 1 || ...
    isinf(opts.dim) || ...
    isnan(opts.dim) || ...
    opts.dim ~= fix(opts.dim) || ...
    opts.dim < 1 || ...
    opts.dim > ndims(x) || ...
    size(x, opts.dim) == 1
    opts.dim = findfirst(size(x) > 1);
end
dim = opts.dim;
if ~isfield(opts, 'nonull') || ...
   ~islogical(opts.nonull) || ...
    numel(opts.nonull) ~= 1
    opts.nonull = false;
end
if ~isfield(opts, 'npts') || ...
   ~isa(opts.npts, 'double') || ...
    numel(opts.npts) ~= 1 || ...
    isinf(opts.npts) || ...
    isnan(opts.npts) || ...
    opts.npts ~= fix(opts.npts) || ...
    opts.npts < 1
    opts.npts = 721;
else
    opts.npts = max(19, min(7201, opts.npts));
end
if ~isfield(opts, 'pa') || ...
   ~isa(opts.pa, 'double') || ...
    numel(opts.pa) ~= 1 || ...
    isinf(opts.pa) || ...
    isnan(opts.pa) || ...
    opts.pa <= 0 || ...
    opts.pa >= 0.5
    opts.pa = 0.05;
end

% get number of elements (for subsref argument)
n = size(x, dim);
raa = repmat({':'}, 1, ndims(x));
rac = raa;
raa{dim} = ones(1, n);
rab = size(x);
rab(dim) = 1;

% remove mean from data
[mx, ex] = meannoinfnan(x, dim, opts.nonull);
[my, ey] = meannoinfnan(y, dim, opts.nonull);
ex = (ex & ey);
ey = sum(ex, dim);
x = double(x) - mx(raa{:});
y = double(y) - my(raa{:});

% and set invalid values to 0
x(~ex) = 0;
y(~ex) = 0;

% compute variance and covariance -> angle t(heta)
t = 0.5 * atan(2 .* sum(x .* y, dim) ./ (sum(y .* y, dim) - sum(x .* x, dim)));

% rotate data (into b, y)
tc = cos(t);
ts = sin(t);
b = tc(raa{:}) .* x - ts(raa{:}) .* y;
y = ts(raa{:}) .* x + tc(raa{:}) .* y;

% compute radii
a = sqrt(sum(b .* b, dim) ./ ey);
b = sqrt(sum(y .* y, dim) ./ ey);

% create list of points and update subsref argument
y = (0:(2 * pi / opts.npts):(2 * pi))';
rac{dim} = ones(1, numel(y));

% compute factor for ellipse
pz = sdist('norminv', 1 - 0.5 * opts.pa, 0, 1);

% compute ellipse data
x = pz .* a(rac{:}) .* repmat(cos(y), rab);
y = pz .* b(rac{:}) .* repmat(sin(y), rab);
ex = mx(rac{:}) + tc(rac{:}) .* x + ts(rac{:}) .* y;
ey = my(rac{:}) - ts(rac{:}) .* x + tc(rac{:}) .* y;
