function al = skarlags(sk, n, o)
% skarlags  - convert temporal smoothing kernel into AR lags
%
% FORMAT:       al = skarlags(sk [, n [, o]])
%
% Input fields:
%
%       sk          smoothing kernel (1x1 double FWHM in TRs)
%       n           number of time points
%       o           if given and true, create V auto-correlation matrix
%
% Output fields:
%
%       al          auto-correlation lags (lags 1 through L) or V matrix

% Version:  v0.9c
% Build:    12121816
% Date:     Dec-18 2012, 4:37 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, Jochen Weber
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
   ~isa(sk, 'double') || ...
    numel(sk) ~= 1 || ...
    isinf(sk) || ...
    isnan(sk) || ...
    sk < 0
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isa(n, 'double') || ...
    numel(n) ~= 1 || ...
    isinf(n) || ...
    isnan(n) || ...
    n < 2
    n = [];
end

% generate lag array and set index to 1
al = zeros(ceil(8 * sk), 1);
ai = 1;

% iterate until lag r is below threshold
while true

    % compute lag value
    if ~isempty(n)
        al(ai) = exp(-1 / (((n - 1) * (sk / ai) / (n * sqrt(2 * log(2)))) .^ 2));
    else
        al(ai) = exp(-1 / (((sk / ai) / (sqrt(2 * log(2)))) .^ 2));
    end

    % threshold reached
    if al(ai) < eps

        % remove value
        ai = ai - 1;

        % break
        break;
    end

    % increase index
    ai = ai + 1;
end

% only keep required lags
al = al(1:ai);

% don't create matrix of size n
if nargin < 3 || ...
   ~islogical(o) || ...
    numel(o) ~= 1 || ...
   ~o
    return;
end

% create matrix
na = numel(al);
a = [al(end:-1:1); 1; al];
al = zeros(n, n);

% fill in lags
for c = 1:n
    mf = max(1, c - na);
    mt = min(n, c + na);
    af = max(1, 1 + na - (c - mf));
    at = af + (mt - mf);
    al(mf:mt, c) = a(af:at);
end
