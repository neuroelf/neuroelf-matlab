function vif = vifactor(m)
% vifactor  - variance inflation factor computation
%
% FORMAT:       vif = vifactor(m)
%
% Input fields:
%
%       m           design matrix, may contain constant
%
% Output fields:
%
%       vif         variance inflation factors for regressors

% Version:  v0.9c
% Build:    13052018
% Date:     May-01 2013, 11:33 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2013, Jochen Weber
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
   ~isa(m, 'double') || ...
    isempty(m) || ...
    ndims(m) > 2 || ...
    size(m, 2) >= size(m, 1) || ...
    any(isinf(m(:)) | isnan(m(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% don't consider constants
cv = (all(diff(m) == 0, 1));
cm = size(m, 2);

% requires extension
t = ones(size(m, 1), 1);
b = inv(m' * m) * m' * t;
t = t - m * b;
if sum(abs(t)) > (sqrt(eps) * size(m, 1))
    m(:, end+1) = 1;
    cv(end+1) = false;
end
cv(~all(m == 0)) = false;

% initialize vif
vif = ones(1, cm);

% compute VIF for each regressor
vm = var(m);
for c = 1:cm

    % only if regressor is valid
    if cv(c)
        continue;
    end

    % which regressors to use as matrix
    uc = ~cv;
    uc(c) = false;

    % regress
    b = inv(m(:, uc)' * m(:, uc)) * m(:, uc)' * m(:, c);

    % residual
    r = m(:, c) - m(:, uc) * b;

    % variance inflation factor
    vif(c) = 1 / (1 - (vm(c) - var(r)) / vm(c));
end
