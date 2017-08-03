function x = setinfnantoval(x, v, dim, z)
% setinfnantoval  - set Inf/NaN in input to specific values(s)
%
% FORMAT:       x = setinfnantoval(x, v [, dim [, z]])
%
% Input fields:
%
%       x           data
%       v           target value(s, if multiple, alternate along dim)
%       dim         optional dimension, default: 1st non-singleton
%       z           also replace zero values
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
if nargin < 2 || ...
   ~isnumeric(x) || ...
   ~isnumeric(v) || ...
    isempty(v)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing first argument.' ...
    );
end
if nargin < 3 || ...
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
if nargin < 4 || ...
   ~islogical(z) || ...
    numel(z) ~= 1
    z = false;
end

% make v a 1d array
v = v(:);

% find Inf/NaN values
xf = (isinf(x) | isnan(x));

% if requested add 0 values
if z
    xf = (xf | (x == 0));
end

% special case
if numel(x) < 2
    x(xf) = v(xf);
    return;
end

% only one value in v
if numel(v) == 1
    x(xf) = v;

% no permutation needed
elseif dim == findfirst(size(x) > 1)

    % set values in x alternativel from v
    x(xf) = v(1 + mod(1:sum(xf(:)), numel(v)));

% permutation required
else

    % get permutation right
    p = 1:ndims(x);
    p(dim) = [];
    p = [dim, p];
    r = 2:ndims(x);
    r = [r(1:dim-1), 1, r(dim:end)];
    xf = permute(xf, p);

    % depending on numel v
    if numel(v) <= 255
        xi = uint8(xf);
    elseif numel(v) <= 65535
        xi = uint16(xf);
    else
        xi = uint32(xf);
    end

    % put in numbers
    xi(xf) = 1 + mod(1:sum(xf(:)), numel(v));

    % re-permute
    xi = permute(xi, r);
    xf = (xi > 0);

    % set in x
    x(xf) = v(xi(xf));
end
