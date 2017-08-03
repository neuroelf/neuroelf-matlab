function z = indepcorrz(r1, n1, r2, n2)
% indepcorrz  - independent correlation comparison z statistic
%
% FORMAT:       z = indepcorrz(r1, n1, r2, n2)
%
% Input fields:
%
%       r1, r2      N-D correlation values to compare
%       n1, n2      number of cases for r1, r2
%
% Output fields:
%
%       z           z-statistic for difference of correlations between
%                   r1 and r2

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

% argument check
if nargin ~= 4 || ...
   ~isnumeric(r1) || ...
   ~isnumeric(r2) || ...
   ~isnumeric(n1) || ...
   ~isnumeric(n2) || ...
    isempty(r1) || ...
   ~isequal(size(r1), size(r2)) || ...
    numel(n1) ~= 1 || ...
    numel(n2) ~= 1 || ...
    isinf(n1) || ...
    isnan(n1) || ...
    n1 < 3 || ...
    n1 ~= round(n1) || ...
    isinf(n2) || ...
    isnan(n2) || ...
    n2 < 3 || ...
    n2 ~= round(n2)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument (combination) given.' ...
    );
end

% set values to 0 where either is invalid
invval = (isinf(r1) | isnan(r1) | isinf(r2) | isnan(r2));
r1(invval) = 0;
r2(invval) = 0;

% transform to fisher z
r1 = fisherr2z(r1);
r2 = fisherr2z(r2);

% compute statistic
z = (1 / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))) .* (r1 - r2);
