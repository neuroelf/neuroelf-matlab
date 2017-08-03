function a = autocorr(d, dim)
% autocorr  - estimate auto-correlation in data
%
% FORMAT:       a = autocorr(d [, dim])
%
% Input fields:
%
%       d           data matrix (or vector)
%       dim         dimension to work along (default: first non-singleton)
%
% Output fields:
%
%       a           auto-correlation value ([-1 .. 1])

% Version:  v0.9a
% Build:    10051716
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
   ~isnumeric(d) || ...
    numel(d) < 3 || ...
    any(isinf(d(:)) | isnan(d(:)))
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
   ~any(1:ndims(d) == dim)
    dim = findfirst(size(d) > 1);
end
if ~isa(d, 'double')
    d = double(d);
end
nd = size(d, dim);

% build subsref arguments
sr = {':'};
sr = sr(ones(1, ndims(d)));
sr1 = sr;
sr2 = sr;
sr{dim} = ones(1, nd);
sr1{dim} = 1:(nd - 1);
sr2{dim} = 2:nd;

% remove mean from signal
md = (1 / nd) .* sum(d, dim);
d = d - md(sr{:});

% compute correlation
a = (nd / (nd - 1)) .* (sum(d(sr1{:}) .* d(sr2{:}), dim) ./ sum(d .* d, dim));
