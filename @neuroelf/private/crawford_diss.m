function [tval, varargout] = crawford_diss(x,gx,y,gy,opts)
% crawford_diss -  calculate Crawford's single case dissociation t-score
%
% FORMAT:       tval = crawford_diss(x, gx, y, gy [,opts])
%
% Input fields:
%
%       x           score observed in single case for task X
%       gx          scores observed in group for task X
%       y           score observed in single case for task Y
%       gy          scores observed in group for task Y
%       opts        optional 1x1 struct with additional settings
%        .permute   make a permutation test
%                   if requested, the next output argument will be
%                   a matrix with the size of gx, where for each
%                   group's sample data this dissociation test
%                   will be performed on its own
%        .spermute  make a permutation test
%                   if requested, the next output argument will be
%                   a matrix with the size of gx, where for each
%                   group's sample data the test will be performed
%                   if this sample is omited in the calculus
%
% Output fields:
%
%       tval        t-values for (X > Y) ^ (S > G)
%
% References:
%
% [1] Crawford, J.R., Garthwaite, P.H. (2005a):
%     Testing for suspected impairments and dissociations in
%     single-case studies in neuropsychology: Evaluation of
%     alternatives using Monte Carlo simulations and revised
%     tests for dissociations. Neuropsychology, 19, 318-331.
%
% This function works over the last dimension, so x and y may each
% have a size of X-by-Y-by-Z, where gx and gy must then have a
% size of X-by-Y-by-Z-by-N each, with N being > 2!
%
% See also crawford_diff, crawford_abnorm

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
if nargin < 4 || ...
   ~isa(x,  'double') || ...
   ~isa(gx, 'double') || ...
   ~isa(y,  'double') || ...
   ~isa(gy, 'double') || ...
    isempty(x)  || ...
    isempty(gx) || ...
    isempty(y)  || ...
    isempty(gy)
    error( ...
        'neuroelf:BadArgument', ...
        'Too few or invalid argument.' ...
    );
end

% get sizes
sx  = size(x);
psx = prod(sx);
sgx = size(gx);
sy  = size(y);
psy = prod(sy);
sgy = size(gy);
if psx == 1 && ...
    psy == 1 && ...
    sgx(1) == prod(sgx) && ...
    sgy(1) == prod(sgy)
    gx  = gx';
    gy  = gy';
    sgx = size(gx);
    sgy = size(gy);
end
N = sgx(end);
if all([numel(sy), numel(sgx), numel(sgy)] == numel(sx)) && ...
    sx(end) == 1 && ...
    sy(end) == 1
    sx(end) = [];
    sy(end) = [];
end

% check sizes
if numel(sx) ~= numel(sy) || ...
    any(sx ~= sy) || ...
    numel(sgx) ~= numel(sgy) || ...
    any(sgx ~= sgy) || ...
    numel(sx) ~= (numel(sgx) - 1) || ...
    any(sx ~= sgx(1:end-1)) || ...
    N < 3
    error( ...
        'neuroelf:BadArgumentSize', ...
        'Invalid sized argument (N must be > 2)' ...
    );
end

% prepare options
if nargin < 5 || ...
   ~isstruct(opts) || ...
    isempty(opts)
    opts = struct;
else
    opts = opts(1);
end
nout = 0;

% calculation
%
% -> group mean, std, and correlation for formula a, b, and c terms
%
[r{1:2}] = cov_nd(gx, gy);
r = r{2};
r(r < -0.99999) = -0.99999;
r(r > 0.99999) = 0.99999;
mx = mean(gx, numel(sgx));
my = mean(gy, numel(sgy));
dx = std(gx, 0, numel(sgx));
dy = std(gy, 0, numel(sgy));

% revised criteria formulated by Crawford in [1]:
%
% To obtain a p value, we solve phi = y, which is a quadratic equation
% in y^2. Choosing the positive root gives...
%
%     ( -b + sqrt(b^2 - 4ac) ) 1/2
% y = ( -------------------- )
%     (          2a          )
%
% where
%
% a = (1 + r)(1 - r^2)
%
% b = (1 - r){4 (n - 1)^2 + 4(1 + r)(n - 1) + (1 + r)(5 + r)}
%
%          [ x - E(X)         y - E(Y) ] 2   ( n*(n-1)^2 )
% c = -2 * [ --------    -    -------- ]   * ( --------- )
%          [   S(X)             S(Y)   ]     (   n + 1   )
%
% The p value equals Prob(t > y), where t has a standard t distribution
% on n-1 degrees of freedom.
%
% -> term in (c)'s square brackets computed first, to get sign of term
%
mc = (x - mx) ./ dx - (y - my) ./ dy;
mS = 1 - 2 * (mc < 0);

% calculation ob formula terms a, b, and c
%
fa = (1 + r) .* (1 - r.^2);
fb = (1-r) .* (4 * (N-1)^2 + 4 * (1+r) .* (N-1) + (1+r) .* (5+r));
fc = -2 * mc.^2 * (N * (N-1)^2 / (N+1));

% calculation of t-value(s)
%
tsqr = (-fb + sqrt(fb .* fb - 4 * fa .* fc)) ./ (2 * fa);
tval = mS .* real(sqrt(tsqr));

% permute over group
if isfield(opts, 'permute') && ...
    nargout > 1
    tgs = zeros(sgx);
    sro = struct;
    sro.type = '()';
    sro.subs = {};
    for cs = sgx(:)'
        sro.subs{end+1} = 1:cs;
    end
    sri = sro;
    for cN = 1:N
        sri.subs{end} = cN;
        sro.subs{end} = setdiff(1:N, cN);
        tgs = subsasgn(tgs, sri, ...
           crawford_diss( ...
               subsref(gx, sri), subsref(gx, sro), ...
               subsref(gy, sri), subsref(gy, sro)));
    end
    nout = nout + 1;
    varargout{nout} = tgs;
    clear tgs;
end

% single versus permuted group
if isfield(opts, 'spermute') && ...
    nargout > 1
    tgs = zeros(sgx);
    sro = struct;
    sro.type = '()';
    sro.subs = {};
    for cs = sgx(:)'
        sro.subs{end+1} = 1:cs;
    end
    sri = sro;
    for cN = 1:N
        sri.subs{end} = cN;
        sro.subs{end} = setdiff(1:N, cN);
        tgs = subsasgn(tgs, sri, ...
            crawford_diss(x, subsref(gx, sro), y, subsref(gy, sro)));
    end
    nout = nout + 1;
    varargout{nout} = tgs;
    clear tgs;
end
