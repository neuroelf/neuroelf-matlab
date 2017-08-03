function z = stoufferz(z, dim, w)
% stoufferz  - combine z scores into Stouffer's Z
%
% FORMAT:       z = stoufferz(z [, dim [, w]])
%
% Input fields:
%
%       z           z-scores to combine
%       dim         dimension along to combine (default: last)
%       w           weighting vector, number must match size(z, dim)
%
% Output fields:
%
%       z           combined z-score

% Version:  v0.9b
% Build:    10102710
% Date:     Oct-27 2010, 10:30 AM EST
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

% check arguments
if nargin < 1 || ...
    numel(z) < 2 || ...
   (~isa(z, 'single') && ...
    ~isa(z, 'double'))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isa(dim, 'double') || ...
    numel(dim) ~= 1 || ...
    isinf(dim) || ...
    isnan(dim) || ...
    dim < 1 || ...
    dim > ndims(z) || ...
    size(z, fix(dim)) < 2
    dim = findfirst(size(z) > 1, -1);
else
    dim = fix(dim);
end
nw = size(z, dim);
if nargin < 3 || ...
   ~isa(w, 'double') || ...
    numel(w) ~= nw || ...
    numel(w) ~= max(size(w)) || ...
    any(isinf(w) | isnan(w) | w < 0)
    w = [];
end

% replace bad entries in z
z(isinf(z) | isnan(z)) = 0;

% no weighting
if isempty(w)

    % simple formula
    z = sqrt(1 / nw) .* sum(z, dim, 'double');

% weighting
else

    % create access arguments
    a1 = repmat({':'}, 1, ndims(z));
    an = a1;
    a1{dim} = 1;

    % multiply first score
    zs = w(1) .* double(z(a1{:}));

    % iterate along w
    for wc = 2:nw

        % add to sum
        an{dim} = wc;
        zs = zs + w(wc) .* double(z(an{:}));
    end

    % return sum
    z = sqrt(1 / sum(w .* w)) .* zs;
end
