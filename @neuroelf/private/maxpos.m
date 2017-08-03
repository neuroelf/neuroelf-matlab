function [mi, ms] = maxpos(varargin)
% maxpos  - return not the max but only the pos
%
% FORMAT:       [mi, ms] = maxpos(x [, y, dim])
%
% Input fields:
%
%       x       vector or array
%       y       if given and not empty must match x in type and size
%       dim     dimension to work along if y is empty
%
% Output fields:
%
%       mi      maximum index
%       ms      maximum subscript (cell array)
%
% Note: uses the internal max() function of Matlab.

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 11:08 PM EST
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

% no argument check, simply pass on
try
    [mi{1:2}] = max(varargin{:});
    mi = mi{2};
catch ne_eo;
    rethrow(ne_eo);
end

if nargout > 1
    sz = size(varargin{1});
    if nargin > 2
        dim = varargin{3};
    elseif sz(1) > 1
        dim = 1;
    else
        dim = findfirst(sz > 1);
    end
    if dim > 1
        mi = 1 + (mi - 1) .* prod(sz(1:dim-1));
    end
    [ms{1:ndims(varargin{1})}] = ind2sub(sz, mi);
    if numel(mi) > 1
        sz(dim) = 1;
        for c = setdiff(1:ndims(varargin{1}), dim)
            rm = sz;
            rm([c, dim]) = 1;
            ms{c} = repmat(reshape(1:sz(c), ...
                [ones(1, c - 1), sz(c), ones(1, numel(sz) - c)]), rm);
        end
    end
end
