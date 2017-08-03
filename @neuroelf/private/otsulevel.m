function l = otsulevel(d)
% otsulevel  - compute Otsu-based threshold level for data
%
% FORMAT:       l = otsulevel(d)
%
% Input fields:
%
%       d           numeric data
%
% Output fields:
%
%       l           level which minimizes the intra-class variance

% Version:  v0.9d
% Build:    14071118
% Date:     Jul-11 2014, 6:54 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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
   ~isnumeric(d) || ...
    numel(d) < 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% sort data
d = single(d);
d = sort(d(:));

% build running sum of squares and mean for variance computation
dnum = 1 ./ ((1:numel(d))');
dnumm = 1 ./ ((0:numel(d)-1)');
rmean = cumsum(d);
rmeansq = cumsum(d .* d);
rvarc1 = (rmeansq - rmean .* rmean .* dnum) .* dnumm;
d = d(end:-1:1, 1);
rmean = cumsum(d);
rmeansq = cumsum(d .* d);
rvarc2 = (rmeansq - rmean .* rmean .* dnum) .* dnumm;

% find level at the minimum of the sum of variances
rvarc1 = rvarc1 + rvarc2(end:-1:1, 1);
rvarc1(isnan(rvarc1)) = Inf;
l = 0.5 * sum(d([0, 1] + (numel(d) - minpos(rvarc1))));
