function [x, c] = lalign(ref, l, opts)
% lalign  - align two lines (1D data)
%
% FORMAT:       [x, c] = lalign(ref, line [, opts])
%
% Input fields:
%
%       ref         reference line (1D data at 1:N)
%       line        to-be-aligned data (at 1:M)
%       opts        optional settings
%        .bfset     basis function set, one of {'dct'}, 'poly'
%        .bfsetn    number of basis functions (default: min(sqrt(M), sqrt(N)))
%        .interpe   estim interpolation, {'cubic'}, 'lanczos3', 'linear'
%
% Output fields:
%
%       x           basis function set
%       c           coefficients

% Version:  v1.0
% Build:    15032216
% Date:     Mar-22 2015, 4:07 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015, Jochen Weber
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
   ~isnumeric(ref) || ...
    size(ref, 1) ~= numel(ref) || ...
    numel(ref) < 4 || ...
   ~isnumeric(l) || ...
    size(l, 1) ~= numel(l) || ...
    numel(l) < 4
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
n = numel(ref);
m = numel(l);
if nargin < 3 || ...
    numel(opts) ~= 1 || ...
   ~isstruct(opts)
    opts = struct;
end
if ~isfield(opts, 'bfset') || ...
   ~ischar(opts.bfset) || ...
    isempty(opts.bfset) || ...
   ~any(strcmpi(opts.bfset, {'dct', 'poly'}))
    opts.bfset = 'd';
else
    opts.bfset = lower(opts.bfset(1));
end
if ~isfield(opts, 'bfsetn') || ...
   ~isa(opts.bfsetn, 'double') || ...
    numel(opts.bfsetn) ~= 1 || ...
    isinf(opts.bfsetn) || ...
    isnan(opts.bfsetn) || ...
    opts.bfsetn < 1
    opts.bfsetn = round(min(sqrt([m, n])));
else
    opts.bfsetn = round(opts.bfsetn);
end
if ~isfield(opts, 'interpe') || ...
   ~ischar(opts.interpe) || ...
   ~any(strcmpi(opts.interpe(:)', {'linear', 'cubic', 'lanczos3', 'lanczos8'}))
    opts.interpe = 'cubic';
else
    opts.interpe = lower(opts.interpe(:)');
end
if ~isfield(opts, 'robust') || ...
   ~islogical(opts.robust) || ...
    numel(opts.robust) ~= 1
    opts.robust = false;
end

% smpling range
ns = n + 0.25;
nh = ns / 2;

% compute gradient components
dgx = 4 * (flexinterpn_method(ref, [Inf; 1.125; 1; ns], 0, opts.interpe) - ...
    flexinterpn_method(ref, [Inf; 0.875; 1; ns], 0, opts.interpe));

% get sample coordinates
cxyz = (1:1:ns)';

% get data and gradients at coordinates
d = double(ref(:));

% put into matrix for covariance
x = zeros(n, opts.bfsetn);
v1c = zeros(n, opts.bfsetn);
x(:, 1) = 1;
if opts.bfset == 'd'
    for pc = 2:opts.bfsetn
        x(:, pc) = cos(((pc - 1) * pi / (n - 1)) .* (0:(n - 1)))';
    end
else
    for pc = 2:opts.bfsetn
        x(:, pc) = ztrans(((0:(n - 1))') .^ (pc - 1));
        x(:, pc) = ztrans(x(:, pc) - x(:, 1:pc-1) * ...
            (((x(:, 1:pc-1)' * x(:, 1:pc-1)) \ x(:, 1:pc-1)') * x(:, pc)));
    end
end
for pc = 1:opts.bfsetn
    v1c(:, pc) = x(:, pc) .* dgx;
end

% get smoothed version of data
v2s = double(l(:));
maxiter = max(101, 8 * opts.bfsetn);
ss = Inf * ones(1, maxiter);
pss = Inf;
stablec = 0;
bc = 1;
c = zeros(opts.bfsetn, maxiter);
while maxiter > 1
    dxyz = cxyz + x(:, 1:bc) * c(1:bc, maxiter);
    msk = (dxyz >= 1 & dxyz <= m);
    if sum(msk) < nh
        error( ...
            'neuroelf:InternalError', ...
            'Overlap too small.' ...
        );
    end
    f = flexinterpn_method(v2s, dxyz(msk, :), 0, opts.interpe);
    cm1 = v1c(msk, 1:bc);
    dm1 = d(msk);
    sc = sum(dm1) / sum(f);
    dm1 = dm1 - f * sc;

    % which kind of regression
    if opts.robust
        c(1:bc, maxiter - 1) = fitrobustbisquare(cm1, dm1) + c(1:bc, maxiter);
    else
        c(1:bc, maxiter - 1) = transmul(cm1) \ (cm1' * dm1) + c(1:bc, maxiter);
    end
    maxiter = maxiter - 1;
    ss(maxiter) = sum(dm1 .^ 2) / numel(dm1);
    if (pss - ss(maxiter)) / pss < 1e-6
        stablec = stablec + 1;
        if stablec > 2
            bc = bc + 1;
            if bc > opts.bfsetn
                break;
            end
            stablec = 0;
        end
    else
        stablec = 0;
    end
    pss = min(ss);
end

% pick best solution
c = c(:, minpos(ss));
