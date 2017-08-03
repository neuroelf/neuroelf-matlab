function t = slttest(v, varargin)
%SLTTEST  Searchlight-based (single vector) t-test.
%   T = SLTTEST(V) computes a one-sample t-test along the first dimension
%   in vector V.
%
%   T = SLTTEST(V, G) computes a two-sample t-test along the first
%   dimension if (and only if) G is a vector with only 1's and 2's. The
%   test doesn't assume equal variances.
%
%   T = SLTTEST(V, W) computes a weighted one-sample t-test if (and only
%   if) W is a vector of weights between 0 and 1 (inclusive)
%
%   T = SLTTEST(PV) computes a paired t-test if PV is a Vx2 array.
%
%   T = SLTTEST(V, X) computes a GLM regression V ~ X, computing the
%   t-statistic for the first column (e.g. to account for nuisance on V)
%
%   T = SLTTEST(V, X, C) computes a GLM regression, and computes the
%   t-statistic for contrast C.
%
%   T = SLTTEST(V, X, W) computes a GLM regression applying the weights W.
%
%   T = SLTTEST(V, X, W, C) computes contrast C for a weighted regression.

% Version:  v1.1
% Build:    16032400
% Date:     Mar-24 2016, 12:07 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2016, Jochen Weber
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

% arguments
if nargin < 1 || (~isa(v, 'double') && ~isa(v, 'single')) || isempty(v) || size(v, 2) > 2
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end

% simple t-test
if nargin == 1
    if size(v, 2) > 1
        error('neuroelf:general:badArgument', 'Bad or missing argument.');
    end
    gv = (~isnan(v) & ~isinf(v));
    t = sqrt(sum(gv)) * mean(v(gv)) / std(v(gv));
else
    error('neuroelf:general:notYetImplemented', 'Not yet implemented.');
end
