function [m, ge, ges, sd] = meannoinfnan(x, dim, nonull)
% meannoinfnan  - build mean of elements other than Inf/NaN
%
% FORMAT:       [m, ge, ges, sd] = meannoinfnan(x [, dim [, nonull]])
%
% Input fields:
%
%       x           data
%       dim         optional dimension, default: 1st non-singleton
%       nonull      if given and true, also take out 0 values from average
%
% Output fields:
%
%       m           mean over dim (if all elements Inf/NaN := 0)
%       ge          good entries in x
%       ges         sum of ge over dim (number of good entries)
%       sd          standard deviation over same entries

% Version:  v0.9b
% Build:    11041216
% Date:     Apr-12 2011, 4:01 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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
   ~isnumeric(x)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing first argument.' ...
    );
end
if nargin < 2 || ...
   ~isa(dim, 'double') || ...
    numel(dim) ~= 1 || ...
    isinf(dim) || ...
    isnan(dim) || ...
    dim < 1 || ...
    dim > ndims(x)
    dim = findfirst(size(x) > 1);
    if isempty(dim)
        dim = 1;
    end
end
if nargin < 3 || ...
   ~islogical(nonull) || ...
    numel(nonull) ~= 1
    nonull = false;
end
if ~isa(x, 'double') && ...
   ~isa(x, 'single') && ...
   ~nonull
    m = sum(x, dim, 'double') ./ size(x, dim);
    ge = true(size(x));
    ges = size(x, dim) .* ones(size(x));
    if nargout > 3
        sd = sqrt(varc(x, dim));
    end
    return;
end

% build sum of good elements
if nonull
    ge = (~isinf(x) & ~isnan(x) & x ~= 0);
else
    ge = (~isinf(x) & ~isnan(x));
end
x(~ge) = 0;
ges = sum(ge, dim, 'double');
m = sum(x, dim) ./ ges;
m(isnan(m)) = 0;

% also compute standard deviation
if nargout > 3
    x(~ge) = NaN;
    sd = sqrt(varc(x, dim, true));
end
