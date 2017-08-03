function k = smoothkern(f, t, ip, ipg)
% smoothkern  - build N-dim smoothing kernel from FWHM
%
% FORMAT:       k = smoothkern(f [, t [, ip [, method]])
%
% Input fields:
%
%       f           1xF FWHM values (e.g. [2, 2, 2])
%       t           optional weight threshold (default: eps ^ (1/numel(f)))
%       ip          interpolate values instead of N-d convolution
%       method      if given and one of the valid interpolation methods
%                   convolves with the interpolation kernel first
%                   if set to a positive, integer number compute moving
%                   average over nearest neighbor based kernel
%
% Output fields:
%
%       k           kernel for convolution
%
% Note: to use the method-based computation, ip must be set to false

% Version:  v0.9d
% Build:    14070711
% Date:     Jul-07 2014, 11:23 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

persistent i_kernel;
if isempty(i_kernel)
    i_kernel = struct;
    k21 = 1025:3073;
    k22 = [1:1024, 3074:4097];

    % cubic spline, see http://en.wikipedia.org/wiki/Bicubic_interpolation
    % create line from -2 ... 2 (with 1/4096 steps)
    k = abs(-2:1/1024:2)';

    % compute term for 0 <= |x| <= 1
    k(k21) =  1.5 * (k(k21) .^ 3) - 2.5 * (k(k21) .^ 2) + 1.0;

    % compute term for 1 < |x| <= 2
    k(k22) = -0.5 * (k(k22) .^ 3) + 2.5 * (k(k22) .^ 2) - 4.0 * k(k22) + 2.0;

    % set all "direct hit" nodes
    k(1:1024:end) = [0, 0, 1, 0, 0];

    % store in persistent struct
    i_kernel.cubic = k;


    % Lanczos kernels, see http://en.wikipedia.org/wiki/Lanczos_resampling
    % from kernel size 2:9
    for kc = 2:9

        % start with line again
        k = (-kc:1/1024:kc)';

        % set 0-value to 1 first
        k(kc * 1024 + 1) = 1;

        % compute sinc / (pi * k)
        ks = sin(pi * k) ./ ((pi * k) .^ 2);
        ks(1:1024:end) = 0;

        % Lanczos' addition
        ka = (kc * ks) .* sin((pi / kc) * k);

        % reset 0-value
        ka(kc * 1024 + 1) = 1;

        % and store it
        i_kernel.(sprintf('lanczos%d', kc)) = ka;
    end

    % store linear kernel as two empty arguments
    i_kernel.linear = ([0:(1/1024):1, (1 - 1/1024):-(1/1024):0])';

    % for nearest neighbor, use pseudo sampling
    i_kernel.nearest = [0.5; ones(1023, 1); 0.5];

    % poly3
    k = abs(-2:1/1024:2)';
    k(k21) = 1 + (-2.25 + 1.25 * k(k21)) .* (k(k21) .* k(k21));
    k(k22) = 3 + ((3.75 - 0.75 * k(k22)) .* (k(k22)) - 6) .* k(k22);
    k(1:1024:end) = [0,0,1,0,0]';
    i_kernel.poly3 = k;

    % spline2 (diff not contiguous, but acceptable approx.)
    k = abs(-2:1/1024:2)';
    k(k21) = 1 + ((-1.8 + k(k21)) .* k(k21) - 0.2) .* k(k21);
    k(k22) = (((-1/3) * (-1 + k(k22)) + 0.8) .* (-1 + k(k22)) - (7/15)) .* (-1 + k(k22));
    k(1:1024:end) = [0,0,1,0,0]';
    i_kernel.spline2 = k;

    % spline3 (diff not contiguous, but good approx.)
    k = abs(-3:1/1024:3)';
    k31 = 2049:4097;
    k32 = [1025:2048, 4098:5121];
    k33 = [1:1024, 5122:6145];
    k(k31) = 1 + (((13/11) * k(k31) - (453/209)) .* k(k31) - (3/209)) .* k(k31);
    k(k32) = (((-6/11) * (-1 + k(k32)) + (270/209)) .* (-1 + k(k32)) - (156/209)) .* (-1 + k(k32));
    k(k33) = (((1/11) * (-2 + k(k33)) - (45/209)) .* (-2 + k(k33)) + (26/209)) .* (-2 + k(k33));
    k(1:1024:end) = [0, 0, 0, 1, 0, 0, 0]';
    i_kernel.spline3 = k;
end

% argument check
if nargin < 1 || ...
   ~isa(f, 'double') || ...
    isempty(f) || ...
    numel(f) ~= max(size(f)) || ...
    ndims(f) > 2 || ...
    any(isinf(f) | isnan(f) | f < 0 | f > (35000 ^ (1 / numel(f))))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument.' ...
    );
end
f = f(:)';
if nargin < 2 || ...
   ~isa(t, 'double') || ...
    numel(t) ~= 1 || ...
    isinf(t) || ...
    isnan(t) || ...
    t < 0 || ...
    t > ((1/4) ^ numel(f))
    t = (1 / (numel(f) ^ numel(f))) * (eps ^ (1 / numel(f)));
end

% check for method argument
if nargin > 3 && ...
    ischar(ipg) && ...
    all(f <= 33) && ...
   ~isempty(ipg) && ...
   ~isempty(regexpi(ipg(:)', '^(cubic|lanczos[2-9]|linear|nearest|poly3|spline[23])$'))

    % get kernel
    k = i_kernel.(lower(ipg(:)'));

    % convolve kernel for each element
    ik = cell(1, numel(f));
    ikn = zeros(1, numel(f));
    for ic = 1:numel(f)
        if ic > 1 && ...
            any(f(ic) == f(1:ic-1))
            iki = findfirst(f(ic) == f(1:ic-1));
            ik{ic} = ik{iki};
            ikn(ic) = ikn(iki);
            continue;
        end
        ik{ic} = conv(k, smoothkern(1024 * f(ic), 1e-20));
        mp = 0.5 * (numel(ik{ic}) + 1);
        mp = mp - 1024 * floor((mp - 1) / 1024);
        ik{ic} = ik{ic}(mp:1024:end);
        ik{ic} = ik{ic} ./ sum(ik{ic});
        ikn(ic) = numel(ik{ic});
    end

    % create kernel
    k = ones([ikn, 1]);

    % then multiply in each dimension
    for ic = 1:numel(f)
        xd = ikn;
        rd = ikn;
        xd(ic) = 1;
        rd([1:ic-1, ic+1:end]) = 1;
        k = k .* repmat(reshape(ik{ic}, [rd, 1]), [xd, 1]);
    end

    % for single kernel and thresholding
    if numel(f) == 1 && ...
        t > 0

        % reduce kernel
        kred = findfirst(abs(k) >= t);
        if ~isempty(kred) && ...
            kred > 1
            k = k(kred:(numel(k)+1-kred));
        end
    end

    % return
    return;
end

% lower bound on FWHM is 0.29, for which the kernel function is 1
f = max(f, 0.29 * ones(size(f)));

% FWHM -> Gaussian parameter
f = f / sqrt(8 * log(2));

% max dist and dimension
md = round(6 .* max(1, log2(f)) .* f);
fd = 2 * md + 1;

% standard way
if nargin < 3 || ...
   ~islogical(ip) || ...
    numel(ip) ~= 1 || ...
   ~ip

    % kernel
    if nargin > 2 && ...
        isa(ip, 'double') && ...
        numel(ip) == 1 && ...
       ~isinf(ip) && ...
       ~isnan(ip) && ...
        ip == round(ip)
        iph = 0.5 * (ip - 1);
        sip = (1 / ip);
    else
        ip = 1;
        sip = 1;
    end
    k = ones([fd, 1]);
    mdb = md + 0.5 - 0.5 / ip;

    % non-nearest neighbor averaging
    if nargin < 4 || ...
       ~islogical(ipg) || ...
        numel(ipg) ~= 1
        ipg = false;
    end

    % for each dim multiply with linear kernel
    for nd = 1:numel(f)
        xd = fd;
        rd = fd;
        xd(nd) = 1;
        rd([1:nd-1, nd+1:end]) = 1;
        ed = exp(- (-mdb(nd):sip:mdb(nd)) .^ 2 ./ (2 * f(nd) .^ 2));
        if ip > 1
            ed = reshape(ed, ip, fd(nd));
            if ipg && ...
                ip > 3
                edw = exp(- (-iph:iph) .^ 2 ./ 2);
                edw = edw ./ sum(edw);
                ed = sum((edw(:) * ones(1, fd(nd))) .* ed, 1);
            else
                ed = sum(ed, 1);
            end
        end
        ed = ed ./ sum(ed);
        k = k .* repmat(reshape(ed, [rd, 1]), [xd, 1]);
    end

% alternatively, use flexinterpn
else

    % create large kernel for interpolation
    f0 = 64 / sqrt(8 * log(2));
    m0 = round(6 .* max(1, log2(f0)) .* f0);
    ed = exp(- (0:m0) .^ 2 ./ (2 * f0 .^ 2));
    ed = ed ./ (2 .* sum(ed) - ed(1));

    % create array for distances of kernel elements
    xd = md + 1;
    k = zeros([xd, 1]);

    % for each dim
    for nd = 1:numel(f)
        xd = md + 1;
        rd = xd;
        xd(nd) = 1;
        rd([1:nd-1, nd+1:end]) = 1;

        % add squared (relative to 64-size kernel)
        ik = (0:md(nd)) ./ (f(nd) / f0);
        k = k + repmat(reshape(ik .* ik, [rd, 1]), [xd, 1]);
    end

    % compute sqrt
    k = sqrt(k);

    % then interpolate from pre-made kernel
    k(:) = flexinterpn_method(ed(:), k(:) + 1, 'lanczos8');

    % then copy in each dimension
    for nd = 1:numel(f)
        rd = repmat({':'}, 1, numel(f));
        rd{nd} = size(k, nd):-1:2;
        k = cat(nd, subsref(k, struct('type', '()', 'subs', {rd})), k);
    end
end

% threshold
if t > 0
    k(k < t) = 0;

    % give back smallest possible kernel
    s = cell(1, numel(f));
    sk = s;
    md = ceil(0.5 .* size(k));
    for sc = 1:numel(s)
        sk{sc} = md(sc);
    end
    for sc = 1:numel(s)
        sks = sk;
        sks{sc} = ':';
        s{sc} = find(k(sks{:}) > 0);
    end
    k = subsref(k, struct('type', '()', 'subs', {s}));
end

% return early
if isa(ip, 'double') && ...
    ip > 1
    return;
end

% make sure all elements sum to 1
k = k ./ sum(k(:));
