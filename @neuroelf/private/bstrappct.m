function z = bstrappct(pe, b, n)
% bstrappct  - compute a z-statistic from the percentile bootstrap
%
% FORMAT:       z = bstrapbca(pe, b, n)
%
% Input fields:
%
%       pe          point-estimate of the statistic (used for sign)
%       b           boot-strapped values of statistic
%       n           null-value (value against which the pe is tested)
%
% Output fields:
%
%       z           z-score that pe is different from n (signed, 1-tailed)
%
% Notes: this computes the bias-corrected and accelerated bootstrap
%        estimate of  Efron, B., & Tibshirani, R. J. (1993):
%        An Introduction to the Bootstrap, Chapman & Hall/CRC: New York.

% Version:  v0.9b
% Build:    10082911
% Date:     Aug-27 2010, 2:01 PM EST
% Author:   Jared van Snellenberg, CANlab, Columbia University, NYC, NY, USA
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2009, Jared van Snellenberg
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
   ~isa(pe, 'double') || ...
    isempty(pe) || ...
    ndims(pe) > 2 || ...
    size(pe, 1) ~= 1 || ...
    any(isinf(pe) | isnan(pe)) || ...
   ~isa(b, 'double') || ...
    ndims(b) > 2 || ...
    size(b, 2) ~= size(pe, 2) || ...
    any(isinf(b(:)) | isnan(b(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
   ~isa(n, 'double') || ...
   (numel(n) ~= 1 && ...
    ~isequal(size(n), size(pe))) || ...
    any(isinf(n) | isnan(n))
    n = 0;
end

% number of bootstrap samples
nb = size(b, 1);

% correction factor (even if all b < n, still remain 1!)
pel = 2 .* double(pe < n) - 1;

% get z-score as norminv of fraction of bootstraps < n
if numel(n) == 1
    z = sdist('norminv', (1 / nb) .* (pel + sum(b > n)), 0, 1);
else
    z = sdist('norminv', (1 / nb) .* (1 + sum(b > n(ones(1, nb), :))), 0, 1);
end
