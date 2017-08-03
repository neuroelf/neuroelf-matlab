function [d, w, fcount] = prefilter(d, opts)
% prefilter  - spike/shift/resonance filter data
%
% FORMAT:       [fd, w, fcount] = prefilter(d, opts)
%
% Input fields:
%
%       d           data to filter
%       opts        mandatory struct but with optional fields
%        .dim       filter dimension (default: 1)
%        .freq      signal frequency (default: 100)
%        .fwin      filtering window seconds (default: 2)
%        .fwincut   filtering frequency cutoff in seconds (default: 0.5)
%        .fwslide   filtering window slide in seconds (default: fwin/16)
%        .kern      smoothing kernel (default: 2)
%        .post      post-diff smoothing kernel (default: 0)
%        .stdt      std(abs(diff)) change threshold (default: 0.02)
%
% Output fields:
%
%       fd          filtered data (forced to double datatype)
%       w           weights of smoothing requirement
%       fcount      number of windows in which filtering occurred

% Version:  v0.9b
% Build:    10062015
% Date:     Jun-20 2010, 9:57 AM EST
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
%       derived from this software without specific prior written
%       permission.
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
   (~isnumeric(d) && ...
    ~istransio(d)) || ...
    numel(d) < 4
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument in call.' ...
    );
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dim') || ...
   ~isa(opts.dim, 'double') || ...
    numel(opts.dim) ~= 1 || ...
   ~any(opts.dim == (1:ndims(d)))
    opts.dim = findfirst(size(d) > 1);
end
if ~isfield(opts, 'freq') || ...
   ~isa(opts.freq, 'double') || ...
    numel(opts.freq) ~= 1 || ...
    isinf(opts.freq) || ...
    isnan(opts.freq) || ...
    opts.freq < 0.1 || ...
    opts.freq > 1e6
    opts.freq = 100;
end
if ~isfield(opts, 'fwin') || ...
   ~isa(opts.fwin, 'double') || ...
    numel(opts.fwin) ~= 1 || ...
    isinf(opts.fwin) || ...
    isnan(opts.fwin) || ...
    opts.fwin < 0.05
    opts.fwin = 2;
else
    opts.fwin = min(30, opts.fwin);
end
if ~isfield(opts, 'fwincut') || ...
   ~isa(opts.fwincut, 'double') || ...
    numel(opts.fwincut) ~= 1 || ...
    isinf(opts.fwincut) || ...
    isnan(opts.fwincut) || ...
    opts.fwincut < 0.05
    opts.fwincut = min(0.5 * opts.fwin, 0.5);
else
    opts.fwincut = min(0.5 * opts.fwin, opts.fwincut);
end
if ~isfield(opts, 'fwslide') || ...
   ~isa(opts.fwslide, 'double') || ...
    numel(opts.fwslide) ~= 1 || ...
    isinf(opts.fwslide) || ...
    isnan(opts.fwslide)
    opts.fwslide = max(1 / opts.freq, opts.fwin / 16);
else
    opts.fwslide = max(1 / opts.freq, min(0.25 * opts.fwin, opts.fwslide));
end
if ~isfield(opts, 'kern') || ...
   ~isa(opts.kern, 'double') || ...
    numel(opts.kern) ~= 1 || ...
    isinf(opts.kern) || ...
    isnan(opts.kern) || ...
    opts.kern < 0
    opts.kern = 2;
end
if ~isfield(opts, 'post') || ...
   ~isa(opts.post, 'double') || ...
    numel(opts.post) ~= 1 || ...
    isinf(opts.post) || ...
    isnan(opts.post) || ...
    opts.post < 0
    opts.post = 0;
end
if ~isfield(opts, 'stdt') || ...
   ~isa(opts.stdt, 'double') || ...
    numel(opts.stdt) ~= 1 || ...
    isinf(opts.stdt) || ...
    isnan(opts.stdt) || ...
    opts.stdt <= 0
    stdt = 1.01;
else
    stdt = 1 + max(1e-6, opts.stdt - floor(opts.stdt));
end

% get data
dt = class(d);
if istransio(d)
    if ~strcmp(dt, 'double')
        d = single(resolve(d));
    else
        d = resolve(d);
    end
elseif ~strcmp(dt, 'double')
    d = single(d);
end
ods = size(d);
if opts.dim == 1
    d = reshape(d, ods(1), prod(ods(2:end)));
    nd = 2;
else
    d = reshape(d, ...
        [prod(ods(1:opts.dim-1)), ods(opts.dim), prod(ods(opts.dim+1:end))]);
    opts.dim = 2;
    nd = ndims(d);
end
di = opts.dim;
ds = size(d);
if numel(ds) < 3
    ds(3) = 1;
end
dn = ds(di);

% if pre-filtering
kern = smoothkern(opts.kern, 1.42e-5 / opts.kern);
if numel(kern) < 3
    return;
end
kern = flexinterpn_method(kern, [Inf; 1; 1/4096; numel(kern)], 'cubic');
kern = 0.5 .* (kern + kern(end:-1:1));

% compute initial diff and std of diff
sdd = sqrt(varc(abs(diff(d, 1, opts.dim)), opts.dim));
sdd0 = 1000 * sdd;

% while smoothing necessary
smv = sdd0 ./ sdd;
smr = (smv > stdt);
while any(smr(:))

    % next round
    sdd0 = sdd;

    % third dim
    for d3 = 1:ds(3)

        % find indices
        si = smr(:, :, d3);
        if ds(3) > 1 && ...
            sum(si) == 0
            continue;
        end

        % data dim 1
        if di == 1
            for d2 = find(si(:, :, d3))
                d(:, d2, d3) = flexinterpn( ...
                    d(:, d2, d3), [Inf; 1; 1; dn], kern, 4096, 0);
            end
        else
            for d1 = find(si(:, :, d3))'
                d(d1, :, d3) = flexinterpn( ...
                    d(d1, :, d3)', [Inf; 1; 1; dn], kern, 4096, 0)';
            end
        end
    end

    % re-compute diff and std of diff
    sdd = sqrt(varc(abs(diff(d, 1, opts.dim)), opts.dim));

    % while smoothing necessary
    smv = sdd0 ./ sdd;
    smr = (smv > stdt);
end

% detect areas where diff itself is highly fluctuating
dd = diff(d, 1, di);
ddn = dn - 1;

% compute initial diff and std of diff
sdd = sqrt(varc(abs(diff(dd, 1, opts.dim)), opts.dim));
sdd0 = 1000 * sdd;

% while smoothing necessary
smv = sdd0 ./ sdd;
smr = (smv > stdt);
while any(smr(:))

    % next round
    sdd0 = sdd;

    % third dim
    for d3 = 1:ds(3)

        % find indices
        si = smr(:, :, d3);
        if ds(3) > 1 && ...
            sum(si) == 0
            continue;
        end

        % data dim 1
        if di == 1
            for d2 = find(si(:, :, d3))
                dd(:, d2, d3) = flexinterpn( ...
                    dd(:, d2, d3), [Inf; 1; 1; ddn], kern, 4096, 0);
            end
        else
            for d1 = find(si(:, :, d3))'
                dd(d1, :, d3) = flexinterpn( ...
                    dd(d1, :, d3)', [Inf; 1; 1; ddn], kern, 4096, 0)';
            end
        end
    end

    % re-compute diff and std of diff
    sdd = sqrt(varc(abs(diff(dd, 1, opts.dim)), opts.dim));

    % while smoothing necessary
    smv = sdd0 ./ sdd;
    smr = (smv > stdt);
end

% find areas where signal is highly noise (resonance)
addd = abs(diff(dd, 1, opts.dim));
dddn = dn - 2;

% compute initial diff and std of diff
sdd = sqrt(varc(abs(diff(addd, 1, opts.dim)), opts.dim));
sdd0 = 1000 * sdd;

% while smoothing necessary
smv = sdd0 ./ sdd;
smr = (smv > stdt);
while any(smr(:))

    % next round
    sdd0 = sdd;

    % third dim
    for d3 = 1:ds(3)

        % find indices
        si = smr(:, :, d3);
        if ds(3) > 1 && ...
            sum(si) == 0
            continue;
        end

        % data dim 1
        if di == 1
            for d2 = find(si(:, :, d3))
                addd(:, d2, d3) = flexinterpn( ...
                    addd(:, d2, d3), [Inf; 1; 1; dddn], kern, 4096, 0);
            end
        else
            for d1 = find(si(:, :, d3))'
                addd(d1, :, d3) = flexinterpn( ...
                    addd(d1, :, d3)', [Inf; 1; 1; dddn], kern, 4096, 0)';
            end
        end
    end

    % re-compute diff and std of diff
    sdd = sqrt(varc(abs(diff(addd, 1, opts.dim)), opts.dim));

    % while smoothing necessary
    smv = sdd0 ./ sdd;
    smr = (smv > stdt);
end

% rescale
maddd = median(addd, opts.dim) + sqrt(varc(addd, opts.dim));
w = ones(ds);
rm = repmat({':'}, 1, nd);
am = rm;
am{di} = ones(1, dn - 2);
rm{di} = 2:dn-1;
w(rm{:}) = (1 ./ maddd(am{:})) .* (maddd(am{:}) - addd);

% slide full window
wsz = round(opts.freq * opts.fwin);
hsz = 0.75 * wsz;
ssz = ceil(opts.freq * opts.fwslide);

% filtering matrix
cof = ceil(opts.fwin / opts.fwincut);
[wnull, f] = tempfilter(zeros(wsz, 1), struct('tempdt', false, 'tempsc', cof));
f = [f, ones(wsz, 1), ztrans((1:wsz)')];

% temporary arrays
dsum = zeros(ds(1:2));
ksum = zeros(ds(1:2));

% over 3rd dim
fcount = 0;
for d3 = 1:ds(3)

    % while window within data
    wf = 1;
    while wf < (dn - wsz)
        wt = wf + wsz - 1;
        wi = wf:wt;

        % dim = 1
        if di == 1

            % too bad overall weight
            outl = (sum(w(wi, :, d3)) < hsz);

            % for bad items, robustly estimate lower frequencies
            if any(outl)
                fcount = fcount + 1;
                dsum(wi, outl) = dsum(wi, outl) + f * fitrobustbisquare_img(f, d(wi, outl, d3)')';
            elseif ~all(outl)
                dsum(wi, ~outl) = dsum(wi, ~outl) + d(wi, ~outl, d3);
            end
            ksum(wi, :) = ksum(wi, :) + 1;

        % dim = 2
        else

            % too bad overall weight
            outl = (sum(w(:, wi, d3)) < hsz);

            % for bad items, robustly estimate lower frequencies
            if any(outl)
                dsum(outl, wi) = dsum(outl, wi) + f * fitrobustbisquare_img(f, d(outl, wi, d3));
            elseif ~all(outl)
                dsum(~outl, wi) = dsum(~outl, wi) + d(~outl, wi, d3);
            end
            ksum(:, wi) = ksum(:, wi) + 1;
        end

        % slide window
        wf = wf + ssz;
    end

    % replace where ksum >= 0.5 * fwin/fwslide
    kso = max(1, floor(0.5 * opts.fwin / opts.fwslide));
    if di == 1
        ksr = ksum(:, 1) >= kso;
        d(ksr, :, d3) = dsum(ksr, :) ./ ksum(ksr, :);
    else
        ksr = ksum(1, :)' >= kso;
        d(:, ksr, d3) = dsum(:, ksr) ./ ksum(:, ksr);
    end
end

% post filtering (recommended)
if opts.post > 0

    % create new kernel
    kern = smoothkern(opts.post, 1.42e-5 / opts.post);
    if numel(kern) < 3
        return;
    end
    kern = flexinterpn_method(kern, [Inf; 1; 1/4096; numel(kern)], 'cubic');
    kern = 0.5 .* (kern + kern(end:-1:1));

    % re-smooth data
    for d3 = 1:ds(3)
        if di == 1
            for d2 = 1:ds(2)
                d(:, d2, d3) = flexinterpn( ...
                    d(:, d2, d3), [Inf; 1; 1; dn], kern, 4096, 0);
            end
        else
            for d1 = 1:ds(1)
                d(d1, :, d3) = flexinterpn( ...
                    d(d1, :, d3)', [Inf; 1; 1; dn], kern, 4096, 0)';
            end
        end
    end
end

% reshape if necessary
if nd ~= numel(ods)
    d = reshape(d, ods);
end
