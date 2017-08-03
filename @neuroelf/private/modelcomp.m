function [f, df1, df2, bf, ptcf, iXXf, ser, sef, smest] = modelcomp(ftc, rtc, tcd, dim, tol, w, wr)
% modelcomp  - compare two models (F-test)
%
% FORMAT:       [f, df1, df2, bf, pf, iXXf] = modelcomp(x1, x2, y [, tdim [, tol [, w]]])
%
% Input fields:
%
%       x1, x2      design matrices (TxR), x1 = full, x2 = reduced
%       tcd         time course data (N-dim)
%       tdim        temporal dimension (default: 1)
%       tol         tolerance value to set iXX to 0 (default: 4 * eps)
%       w, wr       regression weights (Tx1 or N-dim-x-T)
%
% Output fields:
%
%       f           F-statistic with difference in model fit
%       df1, df2    appropriate degrees of freedom
%       bf          betas of the full model
%       pf          predicted responses of the full model (for resestsmooth)
%       iXXf        inv(X'X) of the full model
%
% Note: this function uses calcbetas and requires Matlab's dropped
%       output functionality.

% Version:  v1.0
% Build:    14120216
% Date:     Dec-02 2014, 4:42 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012 - 2014, Jochen Weber
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
if nargin < 3 || ...
   ~isa(ftc, 'double') || ...
    ndims(ftc) ~= 2 || ...
    size(ftc, 2) > size(ftc, 1) || ...
    isempty(ftc) || ...
    any(isinf(ftc(:)) | isnan(ftc(:))) || ...
    size(ftc, 1) ~= size(rtc, 1) || ...
   ~isa(rtc, 'double') || ...
    ndims(rtc) ~= 2 || ...
    isempty(rtc) || ...
    size(rtc, 2) >= size(ftc, 2) || ...
    any(isinf(rtc(:)) | isnan(rtc(:))) || ...
   ~isnumeric(tcd) || ...
   ~any(size(tcd) == size(rtc, 1)) || ...
    isempty(tcd)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad design matrix or data supplied.' ...
    );
end
numrows = size(ftc, 1);
if nargin < 4 || ...
   ~isa(dim, 'double') || ...
    numel(dim) ~= 1 || ...
    isinf(dim) || ...
    isnan(dim) || ...
   ~any(dim == 1:ndims(tcd)) || ...
    size(tcd, dim) ~= numrows
    dim = findfirst(size(tcd) == numrows);
    if isempty(dim)
        error( ...
            'neuroelf:BadArgument', ...
            'Data and models must match in number of data points.' ...
        );
    end
end
if nargin < 5 || ...
   ~isa(tol, 'double') || ...
    isempty(tol)
    tol = 4 * eps;
elseif numel(tol) ~= 1 || ...
    isnan(tol) || ...
    abs(tol) > sqrt(eps)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad tol argument supplied.' ...
    );
else
    tol = abs(tol);
end
if nargin < 6 || ...
   ~isa(w, 'double') || ...
   (numel(w) ~= numrows && ...
    size(w, dim) ~= numrows) || ...
   (size(w, dim) == numrows && ...
   ~all(size(w) == size(tcd)))
    w = [];
end
if nargin < 7 || ...
   ~isa(wr, 'double') || ...
   ~isequal(size(wr), size(w)) || ...
    any(isinf(wr(:)) | isnan(wr(:)))
    wr = w;
end

% compute stats
[bf, iXXf, ptcf, ser] = calcbetas(rtc, tcd, dim, tol, wr);
[bf, iXXf, ptcf, sef] = calcbetas(ftc, tcd, dim, tol, w);
nsef = (isinf(sef) | isnan(sef));
if any(nsef(:))
    sef(nsef) = max(ser(nsef), 0);
end

% compute sum of w for d.f. (or use numrows)
if isempty(w)
    sw = numrows;
    swr = numrows;
elseif numel(w) == numrows
    sw = sum(w(:));
    swr = sum(wr(:));
else
    sw = sum(w, dim);
    swr = sum(wr, dim);
end

% compute F-stats
f = modelcompse(sef, ser, sw, size(ftc, 2), size(rtc, 2), swr);

% degrees of freedom
df1 = size(ftc, 2) - size(rtc, 2);
df2 = sw - size(ftc, 2);

% smoothness estimate
if nargout > 8
    smest = {[], []};
    if ndims(f) ~= 3 || ...
       ~any([1, 4] == dim)
        return;
    end
    [smest{1}, smest{2}] = resestsmooth(res, 1);
end
