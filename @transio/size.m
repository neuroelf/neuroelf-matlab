function [varargout] = size(A,N)
% transio::size  - overloaded method

% Version:  v0.9c
% Build:    11042919
% Date:     Apr-29 2011, 8:11 PM EST
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

% get size
sA = struct(A);
csz = sA.DataDims;
csl = length(csz);

% without dim
if nargin < 2
    sz = csz;
elseif ~isa(N, 'double') || ...
    numel(N) ~= 1 || ...
    isinf(N) || ...
    isnan(N) || ...
    N < 1 || ...
    N > 63
    error( ...
        'transio:BadArgument', ...
        'Invalid call to size.' ...
    );
else
    N = fix(N);
    if N > csl
        sz = 1;
    else
        sz = csz(N);
    end
end

% outputs
if nargout <= 1
    varargout{1} = sz;
else
    for ac = 1:nargout
        varargout{ac} = sz(ac);
    end
    varargout{nargout} = prod(sz(nargout:end));
end
