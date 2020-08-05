function [k, p, vs] = cohenskappa(v1, v2)
%COHENSKAPPA  Compute Cohen's kappa coefficient of reliability.
%   KAPPA = COHENSKAPPA(V1, V2) computes the kappa between (binary) vectors
%   V1 and V2. Both must have the same number of elements and only 0 or 1
%   values.
%
%   KAPPA = COHENSKAPPA(V) computes the average kappa between all pairs of
%   vectors, V(:,1) with V(:,2) through V(:,end-1) with V(:,end).
%
%   KAPPA = COHENSKAPPA(V1, V) computes the average kappa between all pairs
%   of vectors V1 with V(:,1) through V1 with V(:,end).
%
%   Example:
%
%   % data matrix
%   D = double(rand(100, 5) >= 0.5)
%   KAPPA = COHENSKAPPA(D)
%
%   See: https://en.wikipedia.org/wiki/Cohen%27s_kappa

% Version:  v1.1
% Build:    20061018
% Date:     Jun-10 2020, 6:08 PM EST
% Author:   Jochen Weber, Memorial Sloan Kettering Cancer Center, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2020, Jochen Weber
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
if nargin < 2 && numel(v1) == length(v1)
    error('neuroelf:cohenskappa:badArgument', 'Missing second variable.');
end

% multiple variables in a matrix
if nargin < 2 || (numel(v2) == 1 && isa(v2, double) && ~isinf(v2) && ~isnan(v2) && v2 >= 1 && v2 <= size(v1, 2))
    
    % with single column reference
    nv = size(v1, 2);
    if nargin > 1
        ks = zeros(nv, 1);
        for c1 = 1:nv
            if c1 ~= v2
                ks(c1) = cohenskappa(v1(:, v2), v1(:, c1));
            end
        end
        ks(v2) = [];
        k = mean(ks);
        return;
    end
    
    % iterate over columns
    ks = zeros(round(0.5 * nv * (nv - 1)), 1);
    ki = 1;
    for c1 = 1:nv-1
        for c2 = c1+1:nv
            ks(ki) = cohenskappa(v1(:,c1), v1(:,c2));
            ki = ki + 1;
        end
    end
    k = mean(ks);
    return
end

% multiple columns in second matrix
nv = numel(v1);
v1 = double(v1(:));
if nv == size(v2, 1) && size(v2, 2) > 1
    ks = zeros(size(v2, 2), 1);
    for c2 = 1:numel(ks)
        ks(c2) = cohenskappa(v1, v2(:, c2));
    end
    k = mean(ks);
    return;
end

% must have same size
if nv ~= numel(v2)
    error('neuroelf:cohenskappa:badArgument', 'Vector length mismatch.');
end
v2 = double(v2(:));

% values
if any(v1 ~= 0 & v1 ~= 1) || any(v2 ~= 0 & v2 ~= 1)
    error('neuroelf:cohenskappa:badArgument', 'Bad values in vectors.');
end

% compute elements
v11 = (v1 == 0 & v2 == 0);
v12 = (v1 == 0 & v2 == 1);
v21 = (v1 == 1 & v2 == 0);
v22 = (v1 == 1 & v2 == 1);
p11 = sum(v11) / nv;
p12 = sum(v12) / nv;
p21 = sum(v21) / nv;
p22 = sum(v22) / nv;

% combinations
p1_ = p11 + p12;
p2_ = p21 + p22;
p_1 = p11 + p21;
p_2 = p12 + p22;

% helper variable
pe = p_1 * p1_ + p_2 * p2_;

% and kappa
k = (p11 + p22 - pe) / (1 - pe);

% second output
if nargout > 1
    p = [p11, p12; p21, p22];
    if nargout > 2
        vs = {v11, v12; v21, v22};
    end
end
