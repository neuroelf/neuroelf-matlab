function x = indexarraynb(a, varargin)
% indexarraynb  - additional ways to index into an array (incl. NN interp)
%
% FORMAT:       x = indexarraynb(a, i1, ...)
%
% Input fields:
%
%       a           N-D array (up to 4D)
%       i1, ...     indexing expressions
%
% Output fields:
%
%       x           returned value
%
% Note: supported indexing expressions (either double or uint32) are:
%       > single argument
%         - NxD number-of-values -by- number-of-input-dimensions
%       > multiple arguments (as many as dimensions)
%         - Nx1 (or 1xN) specific indexing for this dimension (ordered list)
%         - N-D argument where as input has lower dimensionality
%           e.g. indexing into a 10x10 array a with a(s, :) where s is a
%           10x1x10 array with numbers
%           note: for now, only first dimension can be expanded this way!!
%       > other than indexarray, this function does not perform
%         boundary checks, return 0 values for index out of bounds

% Version:  v0.9c
% Build:    14013111
% Date:     Nov-21 2011, 3:38 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% pass on to indexarray for now
try
    if size(varargin{1}, 2) == ndims(a)
        rn = ~(any(varargin{1} < 0.5, 2) | ...
            any(varargin{1} > (ones(size(varargin{1}, 1), 1) * size(a)), 2));
        x = zeros(size(varargin{1}, 1), 1);
        if any(rn)
            x(rn) = indexarray(a, varargin{1}(rn, :), varargin{2:end});
        end
    else
        x = indexarray(a, varargin{:});
    end
catch ne_eo;
    rethrow(ne_eo);
end
