function s = skew(data, dim, bc)
% skew  - compute skewness of data
%
% FORMAT:       s = skew(data [, dim [, bc]])
%
% Input fields:
%
%       data        N-dimensional data (single/double)
%       dim         dimension to work along (default: last)
%       bc          bias-correct skewness (default: false)
%
% Output fields:
%
%       s           skewness
%
% Note: instead of computing the 3rd power/moment of the data, the
%       squared data is multiplied again (faster implementation), results
%       can thus differ from Matlab's output!
%
% See also  SKEWNESS.

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
    numel(data) < 3
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing data argument to compute skewness.' ...
    );
end
if nargin < 2 || ...
   ~isa(dim, 'double') || ...
    numel(dim) ~= 1 || ...
   ~any(1:ndims(data) == dim)
    dim = findfirst(size(data) > 1, -1);
end
if size(data, dim) < 3
    error( ...
        'neuroelf:BadArgument', ...
        'Skewness can only be estimated for at least 4 values.' ...
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

% compute 3rd power (faster than .^ 3!)
ds3 = dsq .* data;

% compute means (without inf/nan values) of ds3
mds3 = sum(ds3, dim) ./ ges;

% then compute that mean also (as 1 / mean)
mdsqi = ges ./ sum(dsq, dim);

% finally compute skewness (and don't use .^ 1.5 but multiply by sqrt!)
s = mds3 .* mdsqi .* sqrt(mdsqi);

% correct for bias?
if bc

    % correct s (see skewness.m)
    s = s .* sqrt((ges - 1) ./ ges) .* ges ./ (ges - 2);
end
