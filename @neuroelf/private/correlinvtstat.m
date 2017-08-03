function r = correlinvtstat(t, n)
% correlinvtstat  - convert t-statistic into correlation value
%
% FORMAT:       r = correlinvtstat(t, n)
%
% Input fields:
%
%       t           t statistic
%       n           number of points in correlated series
%
% Output fields:
%
%       r           matching correlation r
%
% See also correltstat

% Version:  v0.9a
% Build:    10062205
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
if nargin < 2 || ...
   ~isnumeric(t) || ...
   ~isnumeric(n) || ...
    isempty(t) || ...
    isempty(n) || ...
    any(isnan(t(:)) | isinf(t(:))) || ...
    any(isinf(n(:)) | isnan(n(:)) | n(:) < 2)
    error( ...
        'neuroelf:BadArgument', ...
        'Missing or invalid argument given.' ...
    );
end
osize = size(t);
if numel(t) == 1 && ...
    numel(n) > 1
    osize = size(t);
end
if ~isa(t, 'double')
    t = double(t(:));
else
    t = t(:);
end
if ~isa(n, 'double')
    n = double(n(:));
else
    n = n(:);
end

% formula
r = reshape(t ./ sqrt(t .* t + n - 2), osize);
