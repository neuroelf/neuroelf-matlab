function k = kurt(data, dim, bc)
% kurt  - compute kurtosis of data
%
% FORMAT:       k = kurt(data [, dim [, bc]])
%
% Input fields:
%
%       data        N-dimensional data (single/double)
%       dim         dimension to work along (default: last)
%       bc          bias-correct kurtosis (default: false)
%
% Output fields:
%
%       k           kurtosis
%
% Note: instead of computing the 4th power/moment of the data, the
%       squared data is squared again (faster implementation), results
%       can thus differ from Matlab's output!
%
% See also  KURTOSIS.

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
   ~isnumeric(data) || ...
    numel(data) < 4
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing data argument to compute kurtosis.' ...
    );
end
if nargin < 2 || ...
   ~isa(dim, 'double') || ...
    numel(dim) ~= 1 || ...
   ~any(1:ndims(data) == dim)
    dim = findfirst(size(data) > 1, -1);
end
if size(data, dim) < 4
    error( ...
        'neuroelf:BadArgument', ...
        'Kurtosis can only be estimated for at least 4 values.' ...
    );
end
if nargin < 3 || ...
   ~islogical(bc) || ...
    numel(bc) ~= 1
    bc = false;
end

% force class (also resolves transio)
if ~isa(data, 'double')
    data = double(data);
end

% compute mean (without inf/nan values)
[mdata, ge, ges] = meannoinfnan(data, dim);

% remove mean
sa = {':'};
sa = sa(ones(1, ndims(data)));
sa{dim} = ones(1, size(data, dim));
data = data - mdata(sa{:});

% set illegal values to 0
data(~ge) = 0;

% compute squared data
dsq = data .* data;

% compute square of square (faster than .^ 4!)
ds4 = dsq .* dsq;

% compute means (without inf/nan values) of ds4
mds4 = sum(ds4, dim) ./ ges;

% then compute that mean also (as 1 / mean)
mdsqi = ges ./ sum(dsq, dim);

% finally compute kurtosis
k = mds4 .* mdsqi .* mdsqi;

% correct for bias?
if bc

    % correct k (see kurtosis.m)
    k = ((ges + 1) .* k - 3 .* (ges - 1)) .* ...
         (ges - 1) ./ ((ges - 2) .* (ges - 3)) + 3;
end
