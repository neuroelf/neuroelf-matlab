function [hc, bins, v] = histcountn(v, from, to, step, w)
% histcountn  - ND histogram across equidistant bins
%
% FORMAT:       hc = histcount(v, from, to, [step [, w]]);
%
% Input fields:
%
%       v           PxD values (points by dimensions)
%       from        1x1 or 1xD from value(s)
%       to          1x1 or 1xD range definition
%       step        1x1 or 1xD step size(s) (take as number of steps of <0)
%       w           optional Px1 weighting information (otherwise uniform)
%
% Output fields:
%
%       hc          ND histogram count
%
% Example:
%
%    % create random numbers (along 3 dims)
%    r = randn(1000, 3);
%
%    % create random weights
%    w = rand(size(r, 1), 1);
%
%    % compute weighted histogram in 3D (61 bins in each dim)
%    wh = histcountn(r, -3, 3, 0.1, w);
%
%    % smooth histogram
%    swh = smoothdata3(wh, [2, 2, 2]);

% Version:  v1.1
% Build:    16032312
% Date:     Mar-23 2016, 12:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2016, Jochen Weber
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

% check arguments
if nargin < 1 || ...
   ~isnumeric(v) || ...
    ndims(v) > 2
    eror( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
np = size(v, 1);
nc = size(v, 2);
ns = ceil(np ^ (1 / nc));
if (ns ^ nc) > (10 * np)
    ns = ns - 1;
end
ss = false;
if nargin < 2 || ...
   (nargin == 2 && ...
    isa(from, 'double') && ...
    numel(from) == 1 && ...
    ~isinf(from) && ...
    ~isnan(from) && ...
    from > 3 && ...
    from <= 512 && ...
    from == fix(from))
    ns = from(1, ones(1, nc));
    mmm = zeros(nc, 6);
    for cc = 1:nc
        mmm(cc, :) = minmaxmean(v(:, cc));
    end
    from = mmm(:, 1)';
    to = mmm(:, 2)';
    ss = true;
    step = (to - from) ./ ns;
elseif nargin < 3 || ...
   ~isa(from, 'double') || ...
   (numel(from) ~= 1 && ...
    numel(from) ~= nc) || ...
    any(isinf(from(:)) | isnan(from(:))) || ...
   ~isa(to, 'double') || ...
   (numel(to) ~= 1 && ...
    numel(to) ~= nc) || ...
    any(isinf(to(:)) | isnan(to(:))) || ...
    any(to(:) - from(:) <= 0)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
else
    from = from(:)';
    if numel(from) == 1
        from = from(1, ones(1, nc));
    end
    to = to(:)';
    if numel(to) == 1
        to = to(1, ones(1, nc));
    end
end
if nargin < 4 && ...
   ~ss
    step = (to - from) ./ ns;
elseif ~isa(step, 'double') || ...
   (numel(step) ~= 1 && ...
    numel(step) ~= nc) || ...
    any(isinf(step(:)) | isnan(step(:)) | step(:) <= 0)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
else
    step = step(:)';
    if numel(step) == 1
        step = step(1, ones(1, nc));
    end
end
if nargin < 5 || ...
   (~isa(w, 'double') && ...
    ~isa(w, 'single')) || ...
    numel(w) ~= np
    w = [];
else
    w = w(:);
end

% compute number of bins in each dimension
nsteps = ceil((0.25 .* step + to - from) ./ step) - 1;

% compute indices
v = double(v);
v = v - repmat(from, np, 1);
hc = 1 + repmat(nsteps ./ (to - from), np, 1) .* v;

% compute integer indices
hc = floor(hc);

% remove any index that is invalid
hc(any(hc < 1, 2) | any(hc > repmat(nsteps, np, 1), 2), :) = [];

% multiply
hc = 1 + sum((hc - 1) .* repmat(cumprod([1, nsteps(1, 1:end-1)], 2), np, 1), 2);

% count
if isempty(w)
    hc = reshape(histcount(hc, 0.5, prod(nsteps) - 0.5, 1), [nsteps, 1]);
else
    hc = reshape(histcount(hc, 0.5, prod(nsteps) - 0.5, 1, w, []), [nsteps, 1]);
end

% compute bins (edges)
if nargout > 1
    bins = cell(1, nc);
    for cc = 1:nc
        bins{cc} = (from(cc)+0.5*step(cc)):step(cc):to(cc);
    end
end
