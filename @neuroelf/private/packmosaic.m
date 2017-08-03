function md = packmosaic(d, sd, td, n)
% packmosaic  - pack one dimension into one/two other(s)
%
% FORMAT:       md = packmosaic(d [, sd [, td [, n]]]);
%
% Input fields:
%
%       d           data to be packed (must be numerical)
%       sd          source (to-mosaic) dimension, default last
%       td          target dimension(s) (default: [2, 1])
%       n           packing scheme (default: ceil(sqrt(size(sd))) )
%
% Output fields:
%
%       md          mosaic'ed data
%
% Note: the result will be squeezed to prevent singleton dims

% Version:  v1.0
% Build:    16010821
% Date:     Jan-08 2016, 9:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
    isempty(d)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing data argument.' ...
    );
end
ds = size(d);
if sum(ds ~= 1) < 3
    md = squeeze(d);
    return;
end
if nargin < 2 || ...
   ~isa(sd, 'double') || ...
    numel(sd) ~= 1 || ...
    isinf(sd) || ...
    isnan(sd) || ...
    sd < 1 || ...
    sd > numel(ds)
    sd = numel(ds);
end
if ds(sd) < 2
    md = squeeze(d);
    return;
end
if nargin < 3 || ...
   ~isa(td, 'double') || ...
    isempty(td) || ...
    numel(td) > 2 || ...
    any(isinf(td) | isnan(td) | td < 1 | td > numel(ds) | td == sd)
    td = [2, 1];
end
if nargin < 4 || ...
   ~isa(n, 'double') || ...
    numel(n) ~= numel(td) || ...
    any(isinf(n) | isnan(n) | n < 1 | n > ds(sd))
    if numel(td) == 1
        n = ds(sd);
    else
        n = ceil(sqrt(ds(sd))) * ones(1, 2);
    end
end
if numel(td) == 1
    n = ds(sd);
end

% generate output
ms = ds;
ms(sd) = 1;
ms(td) = ms(td) .* n;
md = repmat(d(1) - d(1), ms);

% prepare subsref and subsasgn arguments
sr = struct('type', '()', 'subs', {repmat({':'}, 1, numel(ds))});
sa = sr;
sa.subs{sd} = 1;

% pack 1-d
if numel(td) == 1

    % iterate over slices
    ss = ds(td);
    sw = ss - 1;
    fe = 1;
    for sc = 1:n
        sr.subs{sd} = sc;
        sa.subs{td} = fe:(fe+sw);
        md = subsasgn(md, sa, subsref(d, sr));
        fe = fe + ss;
    end

% pack 2-d
else

    % iterate over slices
    sc = 1;
    ss = ds(td);
    sw = ss - 1;
    f2 = 1;
    for c2 = 1:n(2)

        % set outer dim in subsasgn
        sa.subs{td(2)} = f2:(f2+sw(2));
        f1 = 1;
        for c1 = 1:n(1)

            % set value
            sr.subs{sd} = sc;
            sa.subs{td(1)} = f1:(f1+sw(1));
            md = subsasgn(md, sa, subsref(d, sr));
            f1 = f1 + ss(1);
            sc = sc + 1;
            if sc > ds(sd)
                break;
            end
        end

        % check sc
        if sc > ds(sd)
            break;
        end

        % otherwise continue with f2
        f2 = f2 + ss(2);
    end
end

% squeeze output
md = squeeze(md);
