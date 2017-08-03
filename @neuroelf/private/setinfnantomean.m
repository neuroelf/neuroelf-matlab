function x = setinfnantomean(x, dim, t)
% setinfnantomean  - set Inf/NaN in input to mean (over dim)
%
% FORMAT:       x = setinfnantomean(x [, dim [, t]])
%
% Input fields:
%
%       x           data
%       dim         optional dimension, default: 1st non-singleton
%       t           threshold (if more than t Inf/NaN values, set all to 0)
%                   whereas t can be 0...1 (relative) 1...N absolute or
%                   special value 's' for >= sqrt(size(x,dim))
%
% Output fields:
%
%       x           data with Inf/NaN values replaced by mean over dim

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

% get mean over dim
[m, ge] = meannoinfnan(x, dim);

% replicate along dim
xnd = ndims(x);
rsz = ones(1, xnd);
rsz(dim) = size(x, dim);
m = repmat(m, rsz);
x(~ge) = m(~ge);

% threshold number of Inf/NaN values?
if nargin > 2 && ...
    numel(t) == 1 && ...
   ((isa(t, 'double') && ...
    ~isinf(t) && ...
    ~isnan(t) && ...
     t >= 0 && ...
     t <= size(x, dim)) || ...
    (ischar(t) && ...
     lower(t) == 's'))
    if ischar(t)
        t = sqrt(size(x, dim));
    elseif t < 1
        t = t * size(x, dim);
    end
    ge = repmat(sum(~ge, dim, 'double') >= t, rsz);
    x(ge) = 0;
end
