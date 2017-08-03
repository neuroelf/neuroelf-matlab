function s = smoothdata3(d, f, t, skim)
% smoothdata3  - apply gaussian smoothing to 3D data
%
% FORMAT:       s = smoothdata3(d, f [, t [, skim]])
%
% Input fields:
%
%       d           3D double data (for N-dim, each 3D volume is smoothed)
%       f           FWHM kernel sizes
%       t           smoothing kernel threshold (default: 0)
%       skim        smoothing kernel interpolator (method, default: 'none')
%
% Output fields:
%
%       s           smoothed data
%
% Note: if skim is set to one of the interpolation methods supported by
%       flexinterpn_method, the smoothing kernel will be convolved with
%       the interpolation kernel (on a higher resolution) and then
%       sampled accordingly (as with spm_smoothkern)

% Version:  v0.9c
% Build:    13111413
% Date:     Nov-05 2011, 8:21 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2013, Jochen Weber
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
    ~isnumeric(d) || ...
    ndims(d) < 3 || ...
   ~isa(f, 'double') || ...
   ~any(numel(f) == [1, 3]) || ...
    any(isinf(f) | isnan(f) | f < 0 | f > 32)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if numel(f) == 1
    f = f([1, 1, 1]);
end
sd = size(d);
if nargin < 3 || ...
   ~isa(t, 'double') || ...
    numel(t) ~= 1 || ...
    isinf(t) || ...
    isnan(t) || ...
    t < 0
    t = 0;
elseif t > 0.25
    t = 0.25;
end
if nargin < 4 || ...
   ~ischar(skim)
    skim = '';
else
    skim = lower(skim(:)');
end
f = f(:)';

% for fft method
if strcmpi(skim, 'fft')
    kern = smoothkern(f, t, false, 'lanczos8');
    kdim = size(kern);
    rsd = 1 + round(0.5 .* (kdim - sd(1:3)));
    if kdim(1) >= sd(1)
        kern = kern(rsd(1):rsd(1)+2*(floor(0.5*(sd(1)-1))), :, :);
    end
    if kdim(2) >= sd(2)
        kern = kern(:, rsd(2):rsd(2)+2*(floor(0.5*(sd(2)-1))), :);
    end
    if kdim(3) >= sd(3)
        kern = kern(:, :, rsd(3):rsd(3)+2*(floor(0.5*(sd(3)-1))));
    end
    kdim = size(kern);
    kdimh = floor(kdim ./ 2);
    fftkern = zeros(sd(1:3));
    ddh = round((sd(1:3) + 1) / 2);
    fftkern(ddh(1)-kdimh(1):ddh(1)+kdimh(1), ...
            ddh(2)-kdimh(2):ddh(2)+kdimh(2), ...
            ddh(3)-kdimh(3):ddh(3)+kdimh(3)) = kern;
    fftkern = fftn(fftkern);

    % multiply
    if isa(d, 'double')
        s = d;
    else
        s = double(d);
    end
    if numel(sd) == 3
        s(:, :, :) = ifftshift(ifftn(fftn(s(:, :, :)) .* fftkern));
    else
        ns = prod(sd(4:end));
        for sc = 1:ns
            s(:, :, :, sc) = ifftshift(ifftn(fftn(s(:, :, :, sc)) .* fftkern));
        end
    end

    % return
    return;
end

% build 3 separate smoothing kernels for each dim
k1 = smoothkern(f(1), t, false, skim);
k2 = smoothkern(f(2), t, false, skim);
k3 = smoothkern(f(3), t, false, skim);

% for large datasets, we do not wish to create a copy in memory!
if prod(sd(1:3)) > 1e6 || ...
   (~isa(d, 'double') && ...
    ~isa(d, 'single'))

    % create flexinterpn arguments
    fis = [Inf; 1; 1; 1] * ones(1, 3);
    fis(end, :) = sd(1:3);
    skn = {1, 1, 1};
    sks = [0; 1; 0];

    % apply as three steps, but in one go
    s = d;
    if numel(sd) == 3
        s(:, :, :) = flexinterpn(flexinterpn(flexinterpn(d(:, :, :), ...
            fis, {k1, sks, sks}, skn, 0), fis, {sks, k2, sks}, skn, 0), ...
            fis, {sks, sks, k3}, skn, 0);
    else
        ns = prod(sd(4:end));
        for sc = 1:ns
            s(:, :, :, sc) = flexinterpn(flexinterpn(flexinterpn(s(:, :, :, sc), ...
                fis, {k1, sks, sks}, skn, 0), fis, {sks, k2, sks}, skn, 0), ...
                fis, {sks, sks, k3}, skn, 0);
        end
    end

% smaller datasets can be handles by conv3d
else

    % build 3D kernels
    kc1 = zeros(max(max(numel(k1), numel(k2)), numel(k3)) * ones(1, 3));
    kc2 = kc1;
    kc3 = kc1;

    % fill kernels
    nk = size(kc1, 1);
    mk = (nk + 1) / 2;
    kc1(1+(nk-numel(k1))/2:(nk+numel(k1))/2, mk, mk) = k1;
    kc2(mk, 1+(nk-numel(k2))/2:(nk+numel(k2))/2, mk) = k2;
    kc3(mk, mk, 1+(nk-numel(k3))/2:(nk+numel(k3))/2) = k3;

    % apply as three steps, but in one go
    if numel(sd) == 3
        s = conv3d(conv3d(conv3d(double(d), kc1), kc2), kc3);
    else
        s = double(d);
        ns = prod(sd(4:end));
        for sc = 1:ns
            s(:, :, :, sc) = ...
                conv3d(conv3d(conv3d(s(:, :, :, sc), kc1), kc2), kc3);
        end
    end
end
