function t = signpermt(x, n, d, pbo)
% signpermt  - compute t scores from sign permutation test
%
% FORMAT:       t = signpermt(x [, n [, d]])
%
% Input fields:
%
%       x           N-dim data
%       n           number of permutations (default: 10000)
%       d           dimension to run across (default: last)
%
% Output fields:
%
%       t           t-stat with d.f. size(x, d) - 1

% Version:  v0.9c
% Build:    12041015
% Date:     Apr-10 2012, 2:16 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, Jochen Weber
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
    numel(x) < 3 || ...
   ((isa(x, 'single') || ...
     isa(x, 'double')) && ...
    any(isnan(x(:))))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isa(n, 'double') || ...
    numel(n) ~= 1 || ...
    isinf(n) || ...
    isnan(n) || ...
    n < 1
    n = 10000;
else
    n = max(80, min(1e6, n));
end
if nargin < 3 || ...
   ~isa(d, 'double') || ...
    numel(d) ~= 1 || ...
    isinf(d) || ...
    isnan(d) || ...
    d < 1 || ...
    d > ndims(x) || ...
    d ~= fix(d)
    d = findfirst(size(x) > 1, -1);
end
if nargin < 4 || ...
    numel(pbo) ~= 1 || ...
   (~isa(pbo, 'xprogress') && ...
    ~isxfigure(pbo, true))
    pbo = [];
end

% convert to double
if ~isa(x, 'double')
    x = double(x);
end

% get number of values along dimension
nx = size(x, d);

% full and half number of possible sign permutations
fnx = 2 ^ nx;
hnx = 2 ^ (nx - 1);

% compute sum (point-estimate for testing)
sx = sum(x, d);

% if number of iterations is larger or equal to available permutations
if fnx <= n

    % create permutation array as all possible combinations
    pa = ones(nx, fnx);

    % fill the first row
    pa(1, :) = repmat([1, -1], 1, hnx);

    % fill in other rows
    for rc = 2:nx
        pa(rc, 1:2:end) = pa(rc - 1, 1:hnx);
        pa(rc, 2:2:end) = pa(rc - 1, hnx+1:fnx);
    end

% otherwise
else

    % create random sign values
    pa = -1 + 2 .* (rand(nx, n) >= 0.5);
end

% create counter variable
cx = zeros(size(sx));

% create replication and reshaping arguments
rm = size(x);
rm(d) = 1;
rs = ones(1, ndims(x));
rs(d) = nx;

% test xprogress
np = size(pa, 2);
vcn = Inf;
vcd = 2 / 86400;
if (np * numel(x)) > 2e8
    try
        if isempty(pbo)
            pbar = xprogress;
            xprogress(pbar, 'setposition', [80, 200, 640, 36]);
            xprogress(pbar, 'settitle', sprintf('Computing %d permutations...', np));
            xprogress(pbar, 0, 'Permutation 1', 'visible', 0, 1);
        else
            pbar = pbo;
            pbar.Progress(0, 'Permutation 1');
        end
        vcn = now + vcd;
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        pbar = [];
    end
else
    pbar = [];
end

% iterate over permutations
for pc = 1:np;

    % compute sum over product of x with signs
    tsx = sum(x .* repmat(reshape(pa(:, pc), rs), rm), d);

    % add to counter
    cx = cx + (sx > tsx) + 0.5 .* (sx == tsx);

    % update progress bar
    if now >= vcn
        pbar.Progress(pc / np, sprintf('Permutation %d / %d', pc, np));
        vcn = now + vcd;
    end
end

% convert from count to p
p = squeeze(cx ./ np);
ph = (p > 0.5);
p(ph) = p(ph) - 0.5 / np;
ph = (p < 0.5);
p(ph) = p(ph) + 0.5 / np;

% convert to nominal t
t = sdist('tinv', p, nx - 1);
t(isinf(t) | isnan(t)) = 0;

% clear progress bar
if ~isempty(pbar) && ...
    isempty(pbo)
    closebar(pbar);
end
