function z = bstrapz(pe, b)
% histztrans  - compute a z-statistic from a bootstrapped NULL
%
% FORMAT:       z = bstrapbca(pe, b)
%
% Input fields:
%
%       pe          point-estimate(s) of the statistic
%       b           boot-strapped NULL values of statistic
%
% Output fields:
%
%       z           z-score(s)

% Version:  v1.1
% Build:    19072313
% Date:     Jul-23 2019, 1:42 PM EST
% Author:   Jochen Weber, Memorial Sloan Kettering Cancer Center, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2019, Jochen Weber
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
if nargin < 2 || ~isa(pe, 'double') || isempty(pe) || any(isinf(pe(:)) | isnan(pe(:))) || ...
   ~isa(b, 'double') || isempty(b) || any(isinf(b(:)) | isnan(b(:)))
    error('neuroelf:BadArgument', 'Bad or missing argument.');
end

% compute according Z stats
b = b(:);
nb = size(b, 1);
zt = (0:(1/nb):(1 - 0.5/nb))';
zt = sdist('norminv', zt + 0.5/nb, 0, 1);
zt(end+1) = Inf;
nh = 0.5 * nb;

% sort p
[ps, pidx] = sort(pe(:));

% iterate over values
zi = 1;
z = zeros(size(pe));
for pc = 1:numel(ps)
    
    % find next transfer value
    pv = ps(pc);
    while pv > zt(zi)
        zi = zi + 1;
    end
    
    % store
    if zi <= nh
        z(pidx(pc)) = zt(zi);
    else
        z(pidx(pc)) = zt(zi-1);
    end
end
