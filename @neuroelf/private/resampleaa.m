function [rsd, k] = resampleaa(d, f, rd, p, k)
% resampleaa  - anti-aliased resampling of data
%
% FORMAT:       [d, k] = resampleaa(d, f [, rd [, p [, k]]])
%
% Input fields:
%
%       d           data
%       f           resampling factor
%       rd          resampling dimension (default: first non-singleton)
%       p           phase (default: 0)
%       k           resampling kernel (default: conv of cubic + gauss)
%
% Output fields:
%
%       d           resampled data
%       k           resampling kernel (useful for multiple calls!)
%
% Note: this function requires flexinterpn to be available.

% Version:  v1.1
% Build:    16030622
% Date:     Mar-06 2016, 10:59 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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

% persistent memory
persistent ne_rsk;
if ~isstruct(ne_rsk)
    ne_rsk = struct;
end

% argument check
if nargin < 2 || ...
   ~isnumeric(d) || ...
    numel(d) < 2 || ...
   ~isa(f, 'double') || ...
    numel(f) ~= 1 || ...
    isinf(f) || ...
    isnan(f) || ...
    f <= 0
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
   ~isa(rd, 'double') || ...
    numel(rd) ~= 1 || ...
    isinf(rd) || ...
    isnan(rd) || ...
    rd < 1 || ...
    rd > ndims(d)
    rd = findfirst(size(d) > 1);
else
    rd = round(rd);
end
if nargin < 4 || ...
   ~isa(p, 'double') || ...
    numel(p) ~= 1 || ...
    isinf(p) || ...
    isnan(p)
    p = 0;
end

% nothing to do
if f == 1

    % copy data (forced to double)
    rsd = double(d);

    % kernel?
    if nargout > 1
        if isfield(ne_rsk, 'k1_0')
            k = {ne_rsk.k1_0, 4096};
        else
            k = smoothkern(0, 1e-6);
            [knull, ck] = flexinterpn_method(0, 1, 'cubic');
            k = flexinterpn(k, [Inf; 1; 1 / 4096; numel(k)], ck{:});
            k = conv(k, ck{1});
            k = 0.5 .* (k + k(end:-1:1));
            k([1, numel(k)]) = 0;
            ne_rsk.k1_0 = k;
            k = {k, 4096};
        end
    end
    return;
end

% reshape accordingly
ods = size(d);
if rd == 1
    d = reshape(d, ods(1), prod(ods(2:end)));
else
    d = reshape(d, ...
        [prod(ods(1:rd-1)), ods(rd), prod(ods(rd+1:end))]);
    rd = 2;
end
ds = size(d);
if ds(1) == 1
    ds(1) = [];
    rd = 1;
end
ds(ds == 1) = [];
if numel(ds) == 1
    d = reshape(d, ds, 1);
elseif ~isequal(ds, size(d))
    d = reshape(d, ds);
end

% check kernel
if nargin < 5 || ...
   ~iscell(k) || ...
    numel(k) ~= 2 || ...
   ~isa(k{1}, 'double') || ...
   ~isa(k{2}, 'double') || ...
    numel(k{1}) ~= size(k{1}, 1) || ...
    numel(k{2}) ~= 1 || ...
    isinf(k{2}) || ...
    isnan(k{2}) || ...
    k{2} < 1 || ...
    k{2} ~= fix(k{2}) || ...
    mod(numel(k{1}) - 1, k{2}) ~= 0

    % create permanent kernel for ratios between 1 and 3
    if f > 1 && ...
        f <= 3
        skname = sprintf('k%.1f', f);
        skname(3) = '_';
        if isfield(ne_rsk, skname)
            k = {ne_rsk.(skname), 4096};
        else
            k = smoothkern(sqrt(2) * (f - 1), 1e-6);
            [knull, ck] = flexinterpn_method(0, 1, 'cubic');
            k = flexinterpn(k, [Inf; 1; 1 / 4096; numel(k)], ck{:});
            k = conv(k, ck{1});
            k = 0.5 .* (k + k(end:-1:1));
            k([1, numel(k)]) = 0;
            ne_rsk.(skname) = k;
            k = {k, 4096};
        end

    % sample gaussian kernel for ratios > 3
    elseif f > 3
        k = smoothkern(f, 1.42e-5 / f);
        kfac = max(1, 2 ^ (12 - ceil(log2(ceil((numel(k) - 1) / 12)))));
        k = {flexinterpn_method(k, [Inf; 1; 1 / kfac; numel(k)], 'cubic'), kfac};
        k{1}([1, end]) = 0;
    % for smaller ratios, take cubic kernel
    else
        [knull, k] = flexinterpn_method(0, 1, 'cubic');
    end
end

% bounds
lb = 1 + p;
if f < 1
    ub = size(d, rd) + p + 1 - f;
elseif f > 1
    ub = size(d, rd) + p + (f - 1) / f;
else
    rsd = d;
end

% one dim
if numel(ds) == 1
    rsd = flexinterpn(d, [Inf; lb; f; ub], k{:}, 0);

% two or more dims
else

    % determine output size in rdim first!
    bb = 1 + floor((ub - lb) / f);

    % resampling in first dim
    if rd == 1
        rsd = zeros(bb, ds(2));
        for vc = 1:ds(2)
            rsd(:, vc) = flexinterpn(d(:, vc), [Inf; lb; f; ub], k{:}, 0);
        end

    % resampling in second dim with two dims
    elseif numel(ds) == 2
        rsd = zeros(ds(1), bb);
        for vc = 1:ds(1)
            rsd(vc, :) = flexinterpn(d(vc, :), ...
                [Inf, Inf; 1, lb; 1, f; 1, ub], k{:}, 0);
        end

    % in second dim with third dim
    else
        rsd = zeros(ds(1), bb, ds(3));
        for tc = 1:ds(3)
            for vc = 1:ds(1)
                rsd(vc, :, tc) = flexinterpn(d(vc, :, tc), ...
                    [Inf, Inf; 1, lb; 1, f; 1, ub], k{:}, 0);
            end
        end
    end
end

% reshape if necessary
if ndims(rsd) ~= numel(ods)
    ods(rd) = size(rsd, rd);
    rsd = reshape(rsd, ods);
end
