function [ztc, zf, zsh] = ztrans(tc, dim, tp)
%ZTRANS  Perform z-transformation (setting STD to 1 and MEAN to 0).
%   Z = ZTRANS(S) transforms the signal S to a scaled version with MEAN of
%   0 and STD of 1. This is a typical application in functional magnetic
%   resonance imaging (fMRI) to account for scanner hardward differences
%   that lead to highly different scales of raw signals in arbitrary units
%   (a.u.), compared to z scores.
%
%   The function performs the function ZSCORE in the Statistics toolbox,
%   with the added feature of allowing a subset of data points to be used
%   to estimate the MEAN and STD over. The result will be of type
%   double-precision float regardless of input type.
%
%   The computation is Z = (P - MEAN(P)) ./ STD(P) with correct expansion.
%
%   Z = ZTRANS(S, DIM) specifies the dimension along with to compute the
%   transform. SIZE(S, DIM) must be greater than 2.
%
%   Z = ZTRANS(S, DIM, TP) only uses the data points indexed by TP to
%   compute the MEAN and STD necessary for the transform.
%
%   TP can be set to a 1x1 double, which then, in turn, represents the
%   Z-threshold for a two-pass transform, such that only points within
%   +/- Z will be taken to transform the data.
%
%   Example:
%
%   % create a random (time) series with mean 30 and STD 10
%   r = 30 + 10 .* randn(1000, 1);
%
%   % set 10 data points to another (further) distribution
%   r(ceil(1000 .* rand(1, 10))) = 120 + 200 .* randn(10, 1);
%
%   % transform (without and with +/i Z=3 limits)
%   z1 = ZTRANS(r);
%   z2 = ZTRANS(r, 1, 3);
%
%   % plot (in series z2, the properties of the median are correct)
%   plot([z1, z2]);
%
%   [P, ZFACTOR, ZSHIFT] = ZTRANS(S, ...) also returns the scaling factor
%   and removed mean that was applied.
%
%   Class support for input S: all numeric types
%   Class support for input TP: 1xN double or 1xP logical
%
%   See also PSCTRANS, ZSCORE.

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
if nargin < 1 || nargin > 3 || ~isnumeric(tc) || numel(tc) < 3
    error('neuroelf:ztrans:badArgument', 'Invalid time course data given.');
end
ts = size(tc);

% convert tc to double if necessary
if ~isa(tc, 'double')
    tc = double(tc);
end

% dimension
if nargin < 2 || ~isa(dim, 'double') || numel(dim) ~= 1 || ...
    isinf(dim) || isnan(dim) || fix(dim) ~= dim || ...
    dim < 1 || dim > numel(ts) || ts(round(dim)) < 3
    dim = find(ts > 2);
    if isempty(dim)
        error('neuroelf:ztrans:badArgument', 'Invalid time course data given.');
    end
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
if nargin < 3 || (~islogical(tp) || numel(tp) ~= td) && ...
   (~isa(tp, 'double') || numel(tp) ~= max(size(tp)) || ...
     any(isinf(tp) | isnan(tp) | tp < 1 | tp > td))

    % compute mean and 1/std factor over given dim first
    zsh = sum(tc, dim) ./ td;
    ztc = tc - zsh(sref{:});
    zf = 1 ./ sqrt(sum((1 / (td - 1)) .* (ztc .* ztc), dim));

% only select tp from dim
else

    % double time-points
    if ~islogical(tp)
        
        % special case, single number: two-pass z-trans within +/- bounds
        if numel(tp) == 1
            
            % first pass (full data)
            ztc = ztrans(tc, dim);
            
            % then create masked data
            msk = double(abs(ztc) <= tp);
            
            % and compute again (compute weighted sum and sum-of-squares)
            mtc = msk .* ztc;
            stc = sum(mtc .* mtc, dim);
            mtc = sum(mtc, dim);
            msk = sum(msk, dim);
            
            % compute STD
            zf = 1 ./ sqrt((stc - ((mtc .* mtc) ./ msk)) ./ (msk - 1));
            zf(isinf(zf) | isnan(zf) | stc <= 0) = Inf;
            
            % and mean
            zsh = mtc ./ msk;
            zsh(isinf(zsh) | isnan(zsh)) = 0;
            
            % then compute z-transform
            if dim == 1
                rma = [td, 1];
            else
                rma = [1, td];
            end
            ztc = reshape((ztc - repmat(zsh, rma)) .* repmat(zf, rma), ts);
            
            % and return early
            return;
        end
        
        % round for not logicals
        tp = unique(round(tp(:)));
        np = numel(tp);
    else
        np = sum(tp);
    end

    % build input subsref argument
    if dim == 1
        dref = {tp(:), ':'};
    elseif  numel(ts) == 2 % ISMATRIX not available in R2007b
        dref = {':', tp(:)};
    else
        dref = {':', tp(:), ':'};
    end

    % compute mean and 1/std factor over given dim(tp) first
    zsh = (1 / np) .* sum(tc(dref{:}), dim);
    ztc = tc - zsh(sref{:});
    zf = 1 ./ sqrt(sum((1 / (np - 1)) .* (ztc(dref{:}) .^ 2), dim));
end

% no illegal factors
zf(isinf(zf) | isnan(zf)) = 0;

% compute z transform (ztc = (tc - mean(tc)) / std(tc)) with subsref
ztc = reshape(ztc .* zf(sref{:}), ts);
