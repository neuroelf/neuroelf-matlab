function m = mtimesnd(m1, m2, iflag)
%MTIMESND  Matrix multiplication with support along higher dimensions.
%   M = MTIMESND(M1, M2) multiplies the matrices in M1 with size AxMx...
%   with matrices in M2 with size MxBx... into matrices with size AxBx...
%   whereas the additional dimensions must match.
%
%   M = MTIMESND(M1, M2, PERMUTEFLAG) will permute either M1 and/or M2
%   such that the first two dimensions will be transposed for the following
%   values of PERMUTEFLAG:
%
%   0  - do not permute (sizes must be AxMx... and MxBx...)
%   1  - permute dims 1 and 2 in M1 only, for sizes MxAx... and MxBx...
%   2  - permute dims 1 and 2 in M2 only, for sizes AxMx... and BxMx...
%   3  - permute dims 1 and 2 in M1 and M2, for sizes MxAx... and BxMx...
%
%   Example: To compute the covariance matrices of a large set of 2D
%   matrices in a single array variable X (with size AxBx...) with output
%   size AxAx..., this would be the syntax
%
%   X = randn(500, 10, 100, 100);
%   XTX = MTIMESND(X, X, 1);
%   size(XTX)
%
%   ans =
%
%       10    10   100   100
%
%   See also MTIMES.

% Version:  v1.1
% Build:    16031110
% Date:     Mar-11 2016, 10:11 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014 - 2016, Jochen Weber
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
if nargin < 2 || nargin > 3 || ...
   (~isa(m1, 'double') && ~isa(m1, 'single')) || ...
   (~isa(m2, 'double') && ~isa(m2, 'single'))
    error('neuroelf:mtimesnd:badArgument', 'Bad M1 or M2 argument.');
end
s1 = size(m1);
s2 = size(m2);

% check dims
if (~any(s1 == s2(1)) && ~any(s1 == s2(2))) || ~isequal(s1(3:end), s2(3:end))
    error('neuroelf:mtimesnd:dimsMismatch', 'Arguments must match in size.');
end
osx = prod(s1(3:end));

% permute ?
flag = 0;
if nargin > 2
    if ~isnumeric(iflag) || numel(iflag) ~= 1 || ...
        isinf(iflag) || isnan(iflag) || iflag < 0 || iflag > 3
        error('neuroelf:mtimesnd:badArgument', 'Bad PERMUTEFLAG argument.');
    end
    flag = floor(real(double(iflag)));
end

% determine output class
if isa(m1, 'single') && isa(m2, 'single')
    m = single(0);
else
    m = 0;
end

% flag value
switch (flag)
    case 0

        % check sizes
        if size(m1, 2) ~= size(m2, 1)
            error('neuroelf:mtimesnd:dimsMismatch', 'Arguments must match in size.');
        end
        os = [s1(1), s2(2), s1(3:end)];

        % grow output
        m(prod(os)) = 0;
        m = reshape(m, os);

        % classical MTIMES
        for c1 = 1:osx
            m(:, :, c1) = m1(:, :, c1) * m2(:, :, c1);
        end

    % flip dims 1/2 for M1 (follow same path with correct transpose)
    case 1
        if size(m1, 1) ~= size(m2, 1)
            error('neuroelf:mtimesnd:dimsMismatch', 'Arguments must match in size.');
        end
        os = [s1(2), s2(2), s1(3:end)];
        m(prod(os)) = 0;
        m = reshape(m, os);
        for c1 = 1:osx
            m(:, :, c1) = m1(:, :, c1)' * m2(:, :, c1);
        end

    % flip dims 1/2 for M2
    case 2
        if size(m1, 2) ~= size(m2, 2)
            error('neuroelf:mtimesnd:dimsMismatch', 'Arguments must match in size.');
        end
        os = [s1(1), s2(1), s1(3:end)];
        m(prod(os)) = 0;
        m = reshape(m, os);
        for c1 = 1:osx
            m(:, :, c1) = m1(:, :, c1) * m2(:, :, c1)';
        end

    % flip dims 1/2 for both M1/M2
    case 3
        if size(m1, 1) ~= size(m2, 2)
            error('neuroelf:mtimesnd:dimsMismatch', 'Arguments must match in size.');
        end
        os = [s1(2), s2(1), s1(3:end)];
        m(prod(os)) = 0;
        m = reshape(m, os);
        for c1 = 1:osx
            m(:, :, c1) = m1(:, :, c1)' * m2(:, :, c1)';
        end
end
