function es = end(A,K,N)
% transio::end  - overloaded method
%
% FORMAT:       ES = end(obj, K, N)
%
% Input fields:
%
%       obj         1x1 transio object
%       K           end-dimension
%       N           number of total dimensions in expression
%
% Output fields:
%
%       ES          end (size) of (remaining) dimension(s)

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

% check arguments
if nargin < 3 || ...
   ~isa(A, 'transio') || ...
   ~isa(N, 'double') || ...
   ~isa(K, 'double') || ...
    numel(N) ~= 1 || ...
    numel(K) ~= 1 || ...
    isinf(K) || ...
    isnan(K) || ...
    K < 1 || ...
    N < K
    error( ...
        'transio:BadArgument', ...
        'Bad or too few arguments for transio::end.' ...
    );
end
K = floor(K);

% get size
csz = A.DataDims;
csl = numel(csz);

% if end requested beyond defined size
if K > csl

    % return 1
    es = 1;

% otherwise
else

    % make sure N is valid
    N = max(1, ceil(N));

    % collapse if N < number of desired dims
    if N < csl && ...
        K == N
        csz(N) = prod(csz(N:end));
    end

    es = csz(K);
end
