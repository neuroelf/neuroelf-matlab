function [t, df] = onesamplet(x, dim, nonull, nomt)
% onesamplet  - compute onesample t-test over dimension
%
% FORMAT:       [t, df] = onesamplet(x [, dim [, nonull [, nomt]]])
%
% Input fields:
%
%       x           data
%       dim         optional dimension, default: 1st non-singleton
%       nonull      if given and true, also take out 0 values from test
%       nomt        if given and true, nominalize to d.f. size(x,dim)-1
%
% Output fields:
%
%       t           t-test over dimension
%       df          degrees of freedom (skipped values)

% Version:  v0.9b
% Build:    11041216
% Date:     Apr-12 2011, 4:06 PM EST
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
if nargin < 4 || ...
   ~islogical(nomt) || ...
    numel(nomt) ~= 1
    nomt = false;
end

% get mean and SD
[mx, ge, ges, sdx] = meannoinfnan(x, dim, nonull);

% compute t-test
t = sqrt(ges) .* (mx ./ sdx);

% replace bad values
t(isinf(t) | isnan(t) | ges < 2) = 0;

% nominalize
if nomt
    df = size(x, dim) - 1;
    t = sdist('tinv', sdist('tcdf', t, max(0, ges - 1)), df);
    t(isinf(t)) = 0;

% otherwise return correct df if requested
elseif nargout > 1
    df = ges - 1;
end
