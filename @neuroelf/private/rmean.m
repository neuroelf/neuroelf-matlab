function [rm, w] = rmean(d, dim)
% rmean  - re-weighted mean
%
% FORMAT:       [rm, w] = rmean(d [, dim])
%
% Input fields:
%
%       d           data array
%       dim         optional dim (otherwise first non-singleton)
%
% Output fields
%
%       rm          re-weighted mean
%       w           weight

% Version:  v0.9a
% Build:    11043015
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
   ~isnumeric(d)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument supplied.' ...
    );
end
if nargin < 2 || ...
   ~isa(dim, 'double') || ...
    numel(dim) ~= 1 || ...
    isinf(dim) || ...
    isnan(dim) || ...
   ~any(dim == 1:ndims(d))
    dim = findfirst(size(d) ~= 1);
    if isempty(dim)
        dim = 1;
    end
end

% compute general mean first
sz = size(d);
ne = sz(dim);
sz(dim) = 1;
rs = osz(sz);
rs(dim) = ne;
ed = 1 / ne;
rm = ed .* sum(d, dim);
ow = ones(rs);

% now iterate until difference to last iteration gets small
mi = 50;
mb = sqrt(1 / ne);
while true
    w = (d - repmat(rm, rs)) .^ 2;
    for dc = 1:ndims(d)
        if dc ~= dim
            w = sum(w, dc);
        end
    end
    w = 1 ./ sqrt(w);
    w = min((1 / mean(w)) .* w, 2);
    nm = (1 / sum(w)) .* sum(repmat(w, sz) .* d, dim);
    wd = abs(w - ow);
    if mi < 1 || ...
        sum(wd) < mb
        rm = nm;
        break;
    end
    rm = nm;
    mi = mi - 1;
end
