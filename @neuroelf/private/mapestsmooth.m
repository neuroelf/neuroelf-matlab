function fwhm = mapestsmooth(map, mmpv, ml, ip2)
% mapestsmooth  - estimate smoothness of 3D map
%
% FORMAT:       fwhm = mapestsmooth(map, mmpv [, ml [, ip2]])
%
% Input fields:
%
%       map         map (image) to estimate smoothness of
%       mmpv        mm per voxel (voxel size)
%       ml          masking length (length of individual vectors, 16)
%       ip2         up-scaling flag, if true, rescale to half-voxels
%
% Output fields:
%
%       fwhm        [X, Y, Z] overall FWHM estimate
%
% Note: this function is experimental for the moment (although well
%       tested against smoothed versions of random data)

% Version:  v1.1
% Build:    16072114
% Date:     Jul-21 2016, 2:04 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2015, 2016, Jochen Weber
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
   ~isnumeric(map) || ...
    ndims(map) ~= 3 || ...
   ~isa(mmpv, 'double') || ...
   (numel(mmpv) ~= 1 && ...
    numel(mmpv) ~= 3) || ...
    any(isinf(mmpv) | isnan(mmpv) | mmpv <= 0)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if numel(mmpv) == 1
    mmpv = mmpv([1, 1, 1]);
end
if nargin < 3 || ...
   ~isa(ml, 'double') || ...
    numel(ml) ~= 1 || ...
    isinf(ml) || ...
    isnan(ml)
    ml = 16;
else
    ml = min(64, min(floor(0.5 * min(size(map))), max(6, ml)));
end

% make sure map is valid
map = double(map(:, :, :));

% increase size by 1 on each edge (to prevent bleeding over effects)
map = flexmask(map, 0, 2, 0, size(map) + 2, [1, 1, 1]);

% find valid voxels
vmap = (~isinf(map) & ~isnan(map) & map ~= 0);

% find central mass
[cs, cv] = clustercoordsc(vmap);
if isempty(cs)
    fwhm = 2 .* size(map) .* mmpv;
    return;
end
vmap = (cv == maxpos(cs));

% find holes in the central mass
[cs, cv] = clustercoordsc(~vmap);
if numel(cs) > 1
    hmap = (cv ~= 0 & cv ~= maxpos(cs));
else
    hmap = false(size(cv));
end

% now compute values for those
map(~vmap) = NaN;
smap = map;
fia = [Inf; 1; 1] * [1, 1, 1];
while any(hmap(:)) && ...
    any(isnan(smap(hmap)))
    smap = flexinterpn(smap, [fia; size(smap)], [0; 0.01; 0.98; 0.01; 0], 1);
end

% and replace the values in vmap
if any(hmap(:))
    map(hmap) = smap(hmap);
end

% then interpolate map
if nargin > 3 && ...
    islogical(ip2) && ...
    numel(ip2) == 1 && ...
    ip2
    vmap = isnan(map);
    vmap = vmap( ...
        round(0.75:0.5:size(vmap, 1) + 0.5), ...
        round(0.75:0.5:size(vmap, 2) + 0.5), ...
        round(0.75:0.5:size(vmap, 3) + 0.5));
    map = flexinterpn_method(map, [[Inf; 0.75; 0.5] * ones(1, 3); size(map) + 0.5], 'lanczos3');
    map(vmap) = NaN;
    ipf = 0.25;
else
    ipf = 0.5;
end

% estimate
fwhm = (ipf .* mmpv(:)') .* [ ...
    smoothest_1d(map, ml) + ...
    smoothest_1d(map(end:-1:1, :, :), ml), ...
    smoothest_1d(permute(map, [2, 3, 1]), ml) + ...
    smoothest_1d(permute(map(:, end:-1:1, :), [2, 3, 1]), ml), ...
    smoothest_1d(permute(map, [3, 1, 2]), ml) + ...
    smoothest_1d(permute(map(:, :, end:-1:1), [3, 1, 2]), ml)];



% subfunction



% perform the required computation on one dimension only
function fwhm = smoothest_1d(map, ml)

% create binary map (list)
bmap = (~isnan(map(:)));

% and make sure that only items of length at least ml are taken
bmap(:, 2:ml) = false;
for c = 2:ml
    bmap(1:end+1-c, c) = bmap(c:end, 1);
end
bmap = all(bmap, 2);

% find indices
fx = find(bmap);

% remove too much overlap
while any(diff(fx) == 1)
    dfx = fx(1 + [0; find(diff(fx) > 1)]);
    for c = 1:floor((ml - 1) / 2)
        bmap(dfx + c) = false;
    end
    fx = find(bmap);
end

% return HIGH for no values
if isempty(fx)
    fwhm = max(size(map));
    return;
end

% create oversampled z-transformed array of values
rmap = ztrans(reshape(...
    map(ones(ml, 1) * fx(:)' + (0:(ml-1))' * ones(1, numel(fx))), ml, numel(fx)));
dmap = flexinterpn(rmap, [inf,inf;1,1;1,1;size(rmap)], ...
        {smoothkern(0.3, 0, false, 'lanczos8'), [0;1;0]}, {1 ,1});

% compute periodograms
mmap = zeros(3200, 1);
c = 1;
while c <= numel(fx)
    li = min(numel(fx), c+999);
    mmap = mmap + sum(pwelch(rmap(:, c:li)-dmap(:, c:li), ...
        struct('nfft', 6400, 'window', ml, 'overlap', 0, 'detrend', 'n')), 2);
    c = c + 1000;
end

% average
pmap = maxpos(mmap(2:end));

% average point wasn't useful
if pmap > 2800

    % pre-smooth map then try again
    smap = smoothdata3(map, [2, 2, 2], 0, 'lanczos8');
    smap(isnan(map)) = NaN;
    fwhm = sqrt(max(0, smoothest_1d(smap, ml) .^ 2 - 4));

% determine if average is good
elseif pmap > 1600

    % if values before are all less
    if all(mmap(2:801) < 0.75 * mmap(pmap + 1))

        % accept value
        fwhm = 3200 / (pmap - 120);

    % potential problem
    else

        % see if there is a second peak
        tpmap = maxpos(mmap(2:1201));

        % not worth reporting
        if mmap(tpmap + 1) < 0.33 * mmap(pmap + 1)

            % compute as well
            fwhm = 3200 / (pmap - 120);

        % otherwise
        else

            % from alternative position
            fwhm = 3200 / (tpmap - 120);
        end

    end

% no problems
else

    % estimate directly
    fwhm = 3200 / (pmap - 120);
end
