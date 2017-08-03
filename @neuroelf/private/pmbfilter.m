function [out, bias, b, bo] = pmbfilter(in, order, mask, opts)
% pmfilter  - apply a poly-mask bias filter to an image
%
% FORMAT:       [out, bias, b, bo] = pmbfilter(in, order, mask [, opts])
%
% Input fields:
%
%       in          input image
%       order       polynomial order
%       mask        mask(s)
%       opts        optional fields
%        .bcutoff   bias cutoff value (lower values are removed, 0.2)
%        .cmask     flag, mask must be continuous (i.e. 1 cluster, false)
%        .debias    which de-biasing method ('biasonly', 'log', {'mult'})
%        .restrict  flag, if true restrict bias to (cutoff/bias^2)
%        .robust    perform robust regression (default: false)
%        .xmask     use given mask in conjunction with autodetected mask
%
% Output fields:
%
%       out         output image (after filtering)
%       bias        bias field in dims of input image
%       b           bias field weights
%       bo          bias field polynomial orders
%
% Note: if mask is empty, an auto-detection algorithm is employed

% Version:  v1.1
% Build:    16020515
% Date:     Feb-05 2016, 3:34 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2012, 2016, Jochen Weber
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
if nargin < 2 || ...
   ~isnumeric(in) || ...
    ndims(in) ~= 3 || ...
   ~isa(order, 'double') || ...
    numel(order) ~= 1 || ...
    isinf(order) || ...
    isnan(order) || ...
   ~any(1:7 == order)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument.' ...
    );
end
if ~isa(in, 'double')
    if istransio(in)
        in = double(resolve(in));
    else
        in = double(in);
    end
end
sz = size(in);
if nargin < 4 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'bcutoff') || ...
   ~isa(opts.bcutoff, 'double') || ...
    numel(opts.bcutoff) ~= 1 || ...
    isinf(opts.bcutoff) || ...
    isnan(opts.bcutoff) || ...
    opts.bcutoff < 0 || ...
    opts.bcutoff > 1
    opts.bcutoff = 0.2;
end
if ~isfield(opts, 'cmask') || ...
   ~islogical(opts.cmask) || ...
    numel(opts.cmask) ~= 1
    opts.cmask = false;
end
if ~isfield(opts, 'debias') || ...
   ~ischar(opts.debias) || ...
    isempty(opts.debias) || ...
   ~any(strcmpi(opts.debias, {'biasonly', 'l', 'log', 'm', 'mult', 'multi'}))
    opts.debias = 'm';
else
    opts.debias = lower(opts.debias(1));
end
if ~isfield(opts, 'robust') || ...
   ~islogical(opts.robust) || ...
    numel(opts.robust) ~= 1
    opts.robust = false;
end
if ~isfield(opts, 'xmask') || ...
   ~islogical(opts.xmask) || ...
    numel(opts.xmask) ~= 1
    opts.xmask = false;
end
if nargin < 3 || ...
   (~islogical(mask) && ...
    ~isa(mask, 'uint8')) || ...
    ndims(mask) < 3 || ...
    ndims(mask) > 4 || ...
    size(in, 1) ~= size(mask, 1) || ...
    size(in, 2) ~= size(mask, 2) || ...
    size(in, 3) ~= size(mask, 3) || ...
    opts.xmask
    if any(sz < 64)
        error( ...
            'neuroelf:BadArgument', ...
            'For small images no pre-mask can be estimated.' ...
        );
    end
    sh = round(0.5 * sz - 31);
    su = sh + 63;
    mx = smoothdata3(in, [3, 3, 3]);
    my = maxgrad(mx);
    mz = (1 / 3) * median(my(mx > mean(mx(:))));
    [dx, dy, dz] = ind2sub([64, 64, 64], maxpos(lsqueeze( ...
        in(sh(1):su(1), sh(2):su(2), sh(3):su(3)) .* ...
        (my(sh(1):su(1), sh(2):su(2), sh(3):su(3)) < mz))));
    dmask = floodfill3(my < mz, dx + sh(1) - 1, dy + sh(2) - 1, dz + sh(3) - 1);
    if sum(dmask(:)) < (0.001 * numel(dmask))
        [dx, dy, dz] = ind2sub([64, 64, 64], maxpos(lsqueeze( ...
            in(sh(1):su(1), sh(2):su(2), sh(3):su(3)) .* ...
            (my(sh(1):su(1), sh(2):su(2), sh(3):su(3)) < ((2 / 3) * mz)))));
        dmask = floodfill3(my < ((4 / 3) * mz), ...
            dx + sh(1) - 1, dy + sh(2) - 1, dz + sh(3) - 1);
    end
    my = mean(mx(:));
    dmask(mx < mean(mx(mx > my))) = 0;
    if opts.xmask && ...
        nargin > 2 && ...
       (islogical(mask) || ...
        isa(mask, 'uint8')) && ...
       ~isempty(mask) && ...
        size(mask, 1) == size(dmask, 1) && ...
        size(mask, 2) == size(dmask, 2) && ...
        size(mask, 3) == size(dmask, 3)
        if size(mask, 4) ~= 1
            mask = repmat(dmask, [1, 1, 1, size(mask, 4)]) & (uint8(mask) > 0);
        else
            mask = dmask & (uint8(mask) > 0);
        end
    else
        mask = dmask;
    end
end
if opts.cmask
    for nm = size(mask, 4):-1:1
        [mys, my] = clustercoordsc(mask(:, :, :, nm) ~= 0);
        if isempty(mys)
            warning( ...
                'neuroelf:BadArgument', ...
                'No voxels found in automatic mask generation.' ...
            );
            mask(:, :, :, nm) = [];
        else
            mask(:, :, :, nm) = (my == maxpos(mys));
        end
    end
end
nm = size(mask, 4);
if islogical(mask)
    mask = uint8(mask);
    for mc = 2:nm
        mask(:, :, :, 1) = uint8(mc) .* uint8(mask(:, :, :, mc) > 0) + ...
            uint8(mask(:, :, :, mc) == 0) .* mask(:, :, :, 1);
    end
else
    imask = mask;
    mask = mask(:, :, :, 1);
    mask(:) = 0;
    mskr = zeros(1, 256);
    nxtv = 1;
    for mc = 1:nm
        imaskv = imask(:, :, :, mc);
        mskv = unique(lsqueeze(imaskv));
        mskv(mskv == 0) = [];
        mskr(:) = 0;
        mskr(mskv) = nxtv:(nxtv + numel(mskv) - 1);
        mask(imaskv > 0) = mskr(imaskv(imaskv > 0));
        nxtv = nxtv + numel(mskv);
    end
    nm = nxtv - 1;
end
mask(:, :, :, 2:end) = [];

% design matrix sizes
bo = [4, 10, 20, 35, 56, 84, 120];

% get indices and voxel data (for each 3D mask)
mx = find(mask(:) > 0);
[mg, yd] = sort(mask(mx));
mx = mx(yd);
sv = [find(diff(mg(:))); numel(mg)];
nv = diff([0; sv]);
nn = sum(nv);
mv = [1; sv + 1];
yd = in(mx);
mg = maxgrad(in);
mg = mg(mx);
[mx, my, mz] = ind2sub(sz, mx);

% computer a center of gravity for each mask
cog = zeros(nm, 3);
for mc = 1:nm
    inval = yd(mv(mc):sv(mc));
    inval = inval - min(inval);
    inval = (eps + (1 / max(inval))) .* inval;
    inval = inval .* inval;
    snval = sum(inval);
    cog(mc, :) = [sum(mx(mv(mc):sv(mc)) .* inval) / snval, ...
        sum(my(mv(mc):sv(mc)) .* inval) / snval, ...
        sum(mz(mv(mc):sv(mc)) .* inval) / snval];
end
cog = (nv' * cog) ./ nn;

% get basis polynomials values for multiplication
szm = (sz - 1) ./ 2;
cog = cog - (szm + 1);
x = (0.5 / szm(1)) .* (-szm(1):szm(1)) - cog(1) / szm(1);
y = (0.5 / szm(2)) .* (-szm(2):szm(2)) - cog(2) / szm(2);
z = (0.5 / szm(3)) .* (-szm(3):szm(3)) - cog(3) / szm(3);
x = [ones(numel(x), 1), x(:), zeros(numel(x), order - 1)];
y = [ones(numel(y), 1), y(:), zeros(numel(y), order - 1)];
z = [ones(numel(z), 1), z(:), zeros(numel(z), order - 1)];
for oc = 2:order
    x(:, oc + 1) = x(:, 2) .^ oc;
    y(:, oc + 1) = y(:, 2) .^ oc;
    z(:, oc + 1) = z(:, 2) .^ oc;
end

% prepare design matrix
bo = bo(order);
X = zeros(nn, nm * bo);
XV = zeros(3, nm * bo);
xc = 0;

% for each mask
sf = zeros(1, nm);
for mc = 1:nm
    df = mv(mc);
    dt = sv(mc);

    % find peak
    mmm = minmaxmean(yd(df:dt));
    [hyd, hyp] = hist(yd(df:dt), round(1 + mmm(2) - mmm(1)));
    hyd = conv(hyd, [1/4, 1/2, 1/4]);
    sf(mc) = round(hyp(maxpos(hyd) - 1));

    % rescale for bias field
    yd(df:dt) = (1 / sf(mc)) .* yd(df:dt);

    % fill design matrix
    for kz = 0:order
        for ky = 0:order
            for kx = 0:order

                % if model order (in total) held
                if (kx + ky + kz) <= order
                    xc = xc + 1;
                    X(df:dt, xc) = ...
                        x(mx(df:dt), kx + 1) .* ...
                        y(my(df:dt), ky + 1) .* ...
                        z(mz(df:dt), kz + 1);
                    XV(:, xc) = [kx; ky; kz];
                end
            end
        end
    end
end
XV = XV(:, 1:bo);

% reject voxels with mg > diff
if nm > 1
    my = ceil(0.5 * min(diff(sort(sf))));
else
    my = mean(mg);
end
mx = false(nn, 1);
for mc = 1:nm
    df = mv(mc);
    dt = sv(mc);

    % either mean or diff
    mz = mean(mg(df:dt));
    if mz < my
        mz = my;
    end
    mx(df:dt) = (mg(df:dt) <= mz);
    nv(mc) = sum(mx(df:dt));
end
nn = sum(nv);

% remove from data and design
X = X(mx, :);
yd = yd(mx);

% reshape arguments
z = reshape(z, [1, 1, size(z)]);
mx = ones(1, sz(1));
my = ones(1, sz(2));
mz = ones(1, sz(3));

% which method
if opts.robust

    % perform robust regression
    b = fitrobustbisquare_img(X, yd');

% or
else

    b = pinv(X' * X) * X' * yd;
end

% weight betas
for mc = 1:nm
    b((mc - 1) * bo + 1:mc * bo) = b((mc - 1) * bo + 1:mc * bo) .* (nv(mc) / nn);
end
b = sum(reshape(b, bo, nm), 2);

% create bias field(s)
bias = zeros(size(in));
for xc = 1:numel(b)
    xy = x(:, XV(1, xc) + 1) * y(:, XV(2, xc) + 1)';
    bias = bias + b(xc) .* xy(:, :, mz) .* z(mx, my, :, XV(3, xc) + 1);
end

% compute output image
switch (opts.debias)
    case {'b'}
        out = in;
    case {'l'}
        out = in .* (1 + log2(2 - bias));
    case {'m'}
        out = in ./ max(bias, opts.bcutoff);
end
