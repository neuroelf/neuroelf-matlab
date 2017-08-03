function [tval, varargout] = crawford_abnorm(x,gx,opts)
% crawford_abnorm -  calculate Crawford's score abnormality measure
%
% FORMAT:       tval = crawford_abnorm(x, gx, [,opts])
%
% Input fields:
%
%       x           score observed in single case for task X
%       gx          scores observed in group for task X
%       opts        optional 1x1 struct with additional settings
%        .permute   make a permutation test
%                   if requested, the next output argument will be
%                   a matrix with the size of gx, where for each
%                   group's sample data this abnormality test
%                   will be performed on its own
%        .spermute  make a permutation test
%                   if requested, the next output argument will be
%                   a matrix with the size of gx, where for each
%                   group's sample data the test will be performed
%                   if this sample is omited in the calculus
%
% Output fields:
%
%       tval        t-values for (S > G)
%
% References:
%
% [1] Crawford, J.R., Howell, D.C. (1998):
%     Comparing an Individual's Test Score Against Norms Derived
%     from Small Samples, The Clinical Neuropsychologist 12, p. 482-486.
%
% [2] Crawford, J.R., Howell, D.C., Garthwaite, P.H. (1998):
%     Payne and Jones Revisited: Estimating the Abnormality of
%     Test Score Differences Using a Midified Paired Samples
%     t Test, Journal of Clinical and Experimental Neuropsychology,
%     Vol. 20, p. 898-905.
%
% [3] Sokal, R.R., Rohlf, J.F. (1995):
%     Biometry (3rd ed.), San Francisco, CA: W.H. Freeman
%
% This function works over the last dimension, so x and y may each
% have a size of X-by-Y-by-Z, where gx and gy must then have a
% size of X-by-Y-by-Z-by-N each, with N being > 2!
%
% See also crawford_diff, crawford_diss

% Version:  v0.9a
% Build:    11122811
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
if nargin < 2 || ...
   ~isnumeric(x) || ...
   ~isnumeric(gx) || ...
    isempty(x) || ...
    isempty(gx)
    error( ...
        'neuroelf:BadArgument', ...
        'Too few or invalid argument.' ...
    );
end

% get sizes
sx  = size(x);
psx = prod(sx);
sgx = size(gx);
if psx == 1 && ...
    sgx(1) == prod(sgx)
    gx  = double(gx)';
    sgx = size(gx);
end
N = sgx(end);
if numel(sgx) == numel(sx) && ...
    sx(end) == 1
    sx(end) = [];
end

% check sizes
if numel(sx) ~= (numel(sgx) - 1) || ...
    any(sx ~= sgx(1:end-1)) || ...
    N < 3
    error( ...
        'neuroelf:BadArgumentSize', ...
        'Invalid sized argument (N must be > 2)' ...
    );
end

% prepare options
if nargin < 3 || ...
   ~isstruct(opts) || ...
    isempty(opts)
    opts = struct;
else
    opts = opts(1);
end
nout = 0;

% ensure double type
x = double(x);
gx = double(gx);

% calculation
%
% -> group mean and std for formula
%
mx = mean(gx, numel(sgx));
dx = std(gx, 0, numel(sgx));


% Crawford's formula (following Sokal and Rohlf, 1995): [1], p. 483
%
%     (        x - E(X)       )
% t = ( --------------------- )
%     ( S(X) * sqrt((N+1)/N)) )
%
tval = (x - mx) ./ (dx * sqrt((N + 1) /N));

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
            crawford_abnorm(subsref(gx, sri), subsref(gx, sro)));
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
            crawford_abnorm(x, subsref(gx, sro)));
    end
    nout = nout + 1;
    varargout{nout} = tgs;
    clear tgs;
end
