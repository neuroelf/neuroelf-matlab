function cfunc = convones(f, cnum, w)
%CONVONES  Convolves a vector with a block (a number of consecutive) ones.
%   C = CONVONES(A, NUMONES) convolves vector A with the vector
%   ONES(NUMONES, 1).
%
%   C = CONVONES(A, NUMONES, true) convolves the vector and divides by the
%   number of ones. This forms a (block-based) moving average.
%
%   Example:
%
%   r = randn(100000, 1);
%   rc = convones(r, 2000);
%   plot(rc);
%
%   See also CONV.

% Version:  v1.1
% Build:    16031110
% Date:     Mar-11 2016, 10:18 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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
if nargin < 2 || nargin > 3 || ~isa(f, 'double') || ~isa(cnum, 'double') || ...
    numel(cnum) ~= 1 || isinf(cnum) || isnan(cnum) || cnum < 1 || cnum ~= fix(cnum)
    error('neuroelf:convones:badArgument', 'Bad or missing argument.');
end

% special case: cnum == 1
if cnum == 1
    cfunc = f;
    return;
end

% get size and dim of f
fs = size(f);
ofs = fs;
if numel(f) == max(fs)
    f = f(:);
    ofs(ofs == numel(f)) = numel(f) + cnum - 1;
    fs = size(f);
else
    ofs(1) = fs(1) + cnum - 1;
end
fn = fs(1);
fs(1) = [];
if numel(fs) < 2
    fs(2) = 1;
end
fsp = prod(fs);

% produce output and intermediate vector
lf = fn + cnum - 1;
cfunc = zeros(lf, fsp);
cimed = cat(1, f, zeros(cnum - 1, fsp));

% weight
if nargin > 2 && islogical(w) && numel(w) == 1 && w
    wfunc = convones(ones(fn, 1), cnum);
else
    wfunc = [];
end

% continue as long as cnum is not 0 !
cpos = 0;
mfac = fix(cnum / 2);
cfac = 1;
while cnum > 0
    if mod(cnum, 2)
        cfunc = cfunc + ...
            [zeros(cpos, fsp); cimed(1:lf-cpos, :)];
        cpos = cpos + cfac;
    end
    if cfac > mfac
        break;
    end
    cimed = [cimed(1:lf-cfac, :) ; zeros(cfac, fsp)] + ...
            [zeros(cfac, fsp) ; cimed(1:lf-cfac, :)];
    cfac = cfac * 2;
    cnum = fix(cnum / 2);
end

% weight now
if ~isempty(wfunc)
    cfunc = cfunc ./ repmat(wfunc, [1, fs]);
end

% shift to f dim
cfunc = reshape(cfunc, ofs);
