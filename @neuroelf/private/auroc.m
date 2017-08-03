function [a, hr, fa] = auroc(v, l)
% auroc -  area under ROC computation
%
% FORMAT:       [a, hr, fa] = auroc(v, l)
%
% Input fields:
%
%       v           values by which labels are discriminated
%       l           labels (vector with two unique values)
%
% Output fields:
%
%       a           area-under-ROC estimate (using trapezoid integral)
%       hr          hit-rate (along dimension)
%       fa          false-alarms (along dimension)

% Version:  v0.9b
% Build:    10111819
% Date:     Nov-18 2010, 7:17 PM EST
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
if nargin < 2 || ...
   ~isnumeric(v) || ...
    isempty(v) || ...
   (~isnumeric(l) && ...
    ~islogical(l)) || ...
    numel(l) ~= max(size(l))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
l = double(l(:));
ul = unique(l);
if any(isinf(l) | isnan(l)) || ...
    numel(ul) ~= 2
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid labels argument.' ...
    );
end
l = (l == ul(2));

% assume dim as the one in which matches number of labels
dim = findfirst(size(v) == numel(l));
if isempty(dim)
    error( ...
        'neuroelf:BadArgument', ...
        'Number of labels must match one of size(v).' ...
    );
end

% get group sizes
n = numel(l);
n2 = sum(l);
n1 = n - n2;

% sort values to get label sorting
[v, vi] = sort(v, dim, 'descend');
l = l(vi);

% generate hr and fa arrays (needed for computation of area anyway!)
sv = size(v);
sv(dim) = 1;
hr = cat(dim, zeros(sv), (1 / n2) .* cumsum(l, dim));
fa = cat(dim, zeros(sv), (1 / n1) .* cumsum(~l, dim));

% compute area (see trapz.m)
sr1 = repmat({':'}, 1, numel(sv));
sr1{dim} = 1:n;
sr2 = sr1;
sr2{dim} = 2:n+1;
a = squeeze(sum(diff(fa, 1, dim) .* (0.5 .* (hr(sr1{:}) + hr(sr2{:}))), dim));
