function ties = numberofties(d, dim)
% numberofties  - returns the number of ties in a vector or matrix
%
% FORMAT:       ties = numberofties(d [, dim])
%
% Input fields:
%
%       d           data vector or matrix
%       dim         dim along which to count ties, default: last > 1
%
% Output fields
%
%       ties        number of ties
%
% Note: see http://www.statsdirect.com/help/nonparametric_methods/kend.htm

% Version:  v0.9c
% Build:    13012810
% Date:     Jan-28 2013, 10:18 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2013, Jochen Weber
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
   (~isa(d, 'double') && ...
    ~isa(d, 'single')) || ...
    numel(d) < 2
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
    dim ~= round(dim) || ...
    dim < 1 || ...
    size(d, dim) < 2
    dim = findfirst(size(d) > 1, -1);
end

% get size of input
isz = size(d);

% set output size
osz = isz;
osz(dim) = 1;

% sort data, then diff data
dsd = double(diff(sort(d, dim), 1, dim) == 0);

% permute to make first dim
if numel(d) ~= size(d, 1)
    if dim > 1
        dsd = permute(dsd, [dim, setdiff(1:ndims(d), dim)]);
    end
    dsd = reshape(dsd, isz(dim) - 1, prod(osz));
else
    dsd = dsd(:);
end

% prepend and append one 0, then re-diff
dsd = diff(cat(1, zeros(1, size(dsd, 2)), dsd, zeros(1, size(dsd, 2))));

% generate output
ties = zeros(osz);

% for along columns
for c = 1:size(dsd, 2)

    % find first and last in ties
    fties = find(dsd > 0);
    lties = find(dsd < 0);

    % combine
    nties = lties - fties;

    % and compute
    ties(c) = sum(0.5 .* nties .* (nties + 1));
end
