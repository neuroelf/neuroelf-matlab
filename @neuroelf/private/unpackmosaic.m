function d = unpackmosaic(md, sd, td, n)
% unpackmosaic  - unpack mosaic from two dimensions into third
%
% FORMAT:       md = packmosaic(d [, sd [, td [, n]]]);
%
% Input fields:
%
%       md          data to be unpacked (must be numerical)
%       sd          source (mosaiced) dimensions, default: [2, 1]
%       td          target dimension(s) (default: after last)
%       n           packing scheme (default: prod of any factor ~= 2)
%
% Output fields:
%
%       d           un-mosaic'ed data

% Version:  v0.9a
% Build:    10062206
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
   ~isnumeric(md) || ...
    isempty(md)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing data argument.' ...
    );
end
ds = size(md);
if nargin < 2 || ...
   ~isa(sd, 'double') || ...
    numel(sd) ~= 2 || ...
    numel(unique(round(sd))) ~= 2 || ...
    any(isinf(sd) | isnan(sd) | sd < 1 | sd > numel(ds))
    sd = [2, 1];
else
    sd = round(sd(:))';
end
if nargin < 3 || ...
   ~isa(td, 'double') || ...
    numel(td) ~= 1 || ...
    isinf(td) || ...
    isnan(td) || ...
    td < 1 || ...
    td > (numel(ds) + 1) || ...
    any(sd == round(td))
    td = numel(ds) + 1;
else
    td = round(td);
end
if nargin < 4 || ...
   ~isa(n, 'double') || ...
    numel(n) ~= 1 || ...
    isinf(n) || ...
    isnan(n) || ...
    n < 1 || ...
    any(n > ds(sd))
    nf1 = factor(ds(sd(1)));
    nf2 = factor(ds(sd(2)));
    nf = intersect(nf1, nf2);
    n = prod(nf(nf ~= 2));
end
if any(mod(ds(sd), n) ~= 0)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid mosaic factor.' ...
    );
end
nn = ds(sd) ./ n;

% generate output
us = ds;
us(sd) = us(sd) ./ n;
us(td) = n * n;
d = repmat(md(1) - md(1), us);

% prepare subsref and subsasgn arguments
sr = struct('type', '()', 'subs', {repmat({':'}, 1, numel(us))});
sa = sr;
sr.subs{sd(1)} = 1;
sr.subs{sd(2)} = 1;
sa.subs{td} = 1;

% iterate over slices
sc = 1;
for c2 = 1:n

    % set outer dim in subsasgn
    sr.subs{sd(2)} = ((c2 - 1) * nn(2) + 1):(c2 * nn(2));
    for c1 = 1:n

        % set value
        sr.subs{sd(1)} = ((c1 - 1) * nn(1) + 1):(c1 * nn(1));
        sa.subs{td} = sc;
        ud = subsref(md, sr);
        if any(ud(:) ~= 0)
            lsc = sc;
        end
        d = subsasgn(d, sa, ud);
        sc = sc + 1;
    end
end

% remove empty slices at end
if (lsc + 1) < sc
    sa.subs{td} = (lsc + 1):(sc - 1);
    d = subsasgn(d, sa, []);
end
