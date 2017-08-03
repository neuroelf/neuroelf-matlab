function h = histnd(v, range, hin, dim, mask)
% histnd  - ND histogram function
%
% FORMAT:       h = histnd(v, range, hin [, dim [, mask]])
%
% Input fields:
%
%       v           ND value array to compute histogram over
%       range       1x2 histogram range
%       hin         1x1 histogram bin edge steps or existing histogram
%       dim         histogram dimension, either {1} or 2
%       mask        optional mask for data
%
% Output fields:
%
%       h           SxV or VxS histogram (created or updated)
%
% Note: this function can either be called to 1) create an initial
%       histogram or to 2) update an existing histogram

% Version:  v0.9c
% Build:    11051913
% Date:     May-19 2011, 12:25 PM EST
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
if nargin < 3 || ...
   ~isnumeric(v) || ...
    isempty(v) || ...
   ~isa(range, 'double') || ...
    numel(range) ~= 2 || ...
    any(isinf(range) | isnan(range)) || ...
    range(2) <= range(1) || ...
   ~isnumeric(hin) || ...
    isempty(hin)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
rd = range(2) - range(1);
if nargin < 4 || ...
   ~isa(dim, 'double') || ...
    numel(dim) ~= 1 || ...
   (~isequal(dim, 1) && ...
    ~isequal(dim, 2))
    dim = [];
end
if nargin < 5 || ...
   ~islogical(mask) || ...
    numel(mask) ~= numel(v)
    mask = [];
end
if ~isempty(mask)
    v = double(lsqueeze(v(mask)));
else
    v = double(v(:));
end
nv = numel(v);

% second round of checks
if numel(hin) == 1
    if ~isa(hin, 'double') || ...
        isinf(hin) || ...
        isnan(hin) || ...
        hin ~= fix(hin) || ...
        hin < 1 || ...
        hin > 1e6
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid number of histogram bin edge steps.' ...
        );
    end
    hin = rd / hin;
    if isempty(dim)
        dim = 1;
    end
    h = histcount(reshape(v, 1, nv), range(1), range(2), hin, 1);
    if dim > 1
        h = h';
    end

% not a valid histogram
elseif ~any(size(hin) == nv) || ...
   (~isempty(dim) && ...
    size(hin, 3 - dim) ~= nv)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid input histogram size.' ...
    );

% otherwise
else

    % copy over
    h = hin;

    % fill empty dim
    if isempty(dim)
        dim = 3 - findfirst(size(h) == nv);
    end

    % get number of bins
    hin = size(h, dim);

    % compute correct bin
    cb = 1 + floor((hin - 1) .* limitrangec((v - range(1)) ./ rd, 0, 1, 0));

    % update histogram
    if dim == 1
        cb = cb + reshape(hin .* (0:(nv-1)), nv, 1);
    else
        cb = nv .* cb + reshape(0:(nv-1), nv, 1);
    end
    h(cb) = h(cb) + 1;
end
