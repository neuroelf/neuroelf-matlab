function [psctc, pscf] = psctrans(tc, dim, tp)
%PSCTRANS  Perform percent signal change (PSC) transformation.
%   P = PSCTRANS(S) transforms the signal S to percentage change (with a
%   mean value of 100) along the first non-singleton dimension. This is a
%   typical application in functional magnetic resonance imaging (fMRI) to
%   account for scanner hardward differences that lead to highly different
%   scales of raw signals in arbitrary units (a.u.), compared to % change.
%
%   The computation is P = 100 .* (P ./ MEAN(P)) -- with correct expansion.
%   The result will be of type double-precision float regardless of input
%   type.
%
%   P = PSCTRANS(S, DIM) specifies the dimension along with to compute the
%   transform. SIZE(S, DIM) must be greater than 1.
%
%   P = PSCTRANS(S, DIM, TP) only uses the data points indexed by TP to
%   compute the mean necessary for the transform.
%
%   [P, PFACTOR] = PSCTRANS(S, ...) also returns the scaling factor that
%   was applied.
%
%   Note: if S is more than 2-D *and* DIM is neither the first or last
%   dimension, this function requires more memory.
%
%   Class support for input S: all numeric types
%   Class support for input TP: 1xN double or 1xP logical
%
%   See also ZTRANS.

% Version:  v1.1
% Build:    16031110
% Date:     Mar-11 2016, 10:09 AM EST
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

% check arguments
if nargin < 1 || nargin > 3 || ~isnumeric(tc) || numel(tc) < 2
    error('neuroelf:psctrans:badArgument', 'Invalid time course data given.');
end
ts = size(tc);

% convert tc to double if necessary
if ~isa(tc, 'double')
    tc = double(tc);
end

% dimension
if nargin < 2 || ~isa(dim, 'double') || numel(dim) ~= 1 || ...
    isinf(dim) || isnan(dim) || fix(dim) ~= dim || ...
    dim < 1 || dim > numel(ts) || ts(round(dim)) < 2
    dim = find(ts ~= 1);
    dim = dim(1);
end
td = ts(dim);

% reshape for a smaller number of dimensions
if dim == 1

    % merge dims 2-N
    tc = reshape(tc, ts(1), prod(ts(2:end)));

% merge dims flexibly
elseif dim > 2 || numel(ts) > 2

    % keep last dimension
    if dim == numel(ts)
        tc = reshape(tc, prod(ts(1:dim-1)), td);

    % keep single second (middle) dim
    else
        tc = reshape(tc, prod(ts(1:dim-1)), td, prod(ts(dim+1:end)));
    end

    % either way, from now on it's dim 2
    dim = 2;
end

% prepare ()-subsref argument
if dim == 1
    sref = {ones(1, td), ':'};
elseif numel(ts) == 2 % ISMATRIX not available in R2007b
    sref = {':', ones(1, td)};
else
    sref = {':', ones(1, td), ':'};
end

% for all points along dim
if nargin < 3 || (~islogical(tp) || numel(tp) == td) && ...
   (~isa(tp, 'double') || numel(tp) ~= max(size(tp)) || ...
    any(isinf(tp) | isnan(tp) | tp < 1 | tp > td))

    % build PSC factor
    pscf = (100 * td) ./ sum(tc, dim);

% only select tp from dim
else

    % round for not logicals
    if ~islogical(tp)
        tp = round(tp);
        np = numel(tp);
    else
        np = sum(tp);
    end

    % build input subsref argument
    if dim == 1
        dref = {tp(:), ':'};
    elseif numel(ts) == 2 % ISMATRIX not available in R2007b
        dref = {':', tp(:)};
    else
        dref = {':', tp(:), ':'};
    end

    % compute mean and PSC factor over given dim(tp) only
    pscf = (100 * np) ./ sum(tc(dref{:}), dim);
end

% clear illegal factors
pscf(isinf(pscf) | isnan(pscf)) = 0;

% multiplication
if numel(ts) == 2 % ISMATRIX not available in R2007b

    % create diagonal sparse factor matrix
    pscff = sparse(1:numel(pscf), 1:numel(pscf), pscf(:)', ...
        numel(pscf), numel(pscf), numel(pscf));

    % along first dim
    if dim == 1

        psctc = tc * pscff;

    % second dim
    else
        psctc = pscff * tc;
    end

    % work around MATLAB "bug" for sparse 1x1 factors leading to sparse
    if issparse(psctc)
        psctc = full(psctc);
    end

% more complicated array
else

    % compute PSC (psctc = 100 * tc / mean(tc)) (requires more memory!)
    psctc = tc .* pscf(sref{:});
end

% reshape
psctc = reshape(psctc, ts);
