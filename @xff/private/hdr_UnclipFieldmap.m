function nxo = hdr_UnclipFieldmap(xo, opts)
% HDR::UnclipFieldmap  - correct the clipped values in a B0 fieldmap
%
% FORMAT:       chdr = hdr.UnclipFieldmap([opts])
%
% Input fields:
%
%       opts        optional struct with settings
%        .imask     use intensity (second volume) as implicit mask (true)
%        .mask      explicit mask (must match in size, default: [])
%        .smk       smoothing kernel in voxels for detection (default: 2)
%
% Output fields:
%
%       chdr        corrected HDR/NII object
%
% Using: clustercoordsc.

% Version:  v1.1
% Build:    16102317
% Date:     Oct-23 2016, 5:40 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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

% neuroelf library
global ne_methods;
clustercoordsc = ne_methods.clustercoordsc;
dilate3d = ne_methods.dilate3d;
erode3d = ne_methods.erode3d;
indexarraynb = ne_methods.indexarraynb;
smoothdata3 = ne_methods.smoothdata3;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'hdr')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
if size(bc.VoxelData, 4) > 2
    error('neuroelf:xff:badArgument', 'HDR/NII object has more than 2 volumes.');
end
if nargin < 2 || numel(opts) ~= 1 || ~isstruct(opts)
    opts = struct;
end
if ~isfield(opts, 'imask') || ~islogical(opts.imask) || numel(opts.imask) ~= 1
    opts.imask = true;
end
if ~isfield(opts, 'smk') || ~isa(opts.smk, 'double') || numel(opts.smk) ~= 1 || ...
    isinf(opts.smk) || isnan(opts.smk) || opts.smk <= 0
    opts.smk = 2;
end

% get data
ovd = bc.VoxelData;
if istransio(ovd)
    ovd = resolve(ovd);
end
vd = double(ovd);
fm = vd(:, :, :, 1);

% image size and central coordinate
imsz = size(fm);
ccrd = 0.5 .* (imsz + 1);

% estimate range of field map
fmr = max(abs(fm(:)));

% implicit masking
if size(vd, 4) > 1 && opts.imask
    im = vd(:, :, :, 2);
    mim = mean(im(:));
    im = (im ./ median(im(im > mim))) >= 0.25;
    [tsizes, im] = clustercoordsc(im, 1);
    [tnull, tsize] = max(tsizes);
    im = (im == tsize);
    fm(~im) = 0;
else
    im = true(size(fm));
end

% explicit masking
if isfield(opts, 'mask') && (isnumeric(opts.mask) || islogical(opts.mask)) && ...
    isequal(size(opts.mask), size(fm))
    if islogical(opts.mask)
        fm(~opts.mask) = 0;
        im(~opts.mask) = false;
    else
        fm(opts.mask == 0) = 0;
        im(opts.mask == 0) = false;
    end
end

% get clipped tails
fmt = ceil(0.333 * fmr);
t1 = (fm >= fmt);
t2 = (fm <= -fmt);
if mod(fmr, 2) ~= 0
    fmc = 2 * (fmr + 1);
else
    fmc = 2 * fmr;
end

% cluster the parts we're interested in
[tnull, t1v] = clustercoordsc(t1, 1, 8);
[tnull, t2v] = clustercoordsc(t2, 1, 8);
t1vd = dilate3d(t1v > 0);
t2vd = dilate3d(t2v > 0);
tcv = t1vd & t2vd;

% cluster
[tnull, tcvv, cll, cllc] = clustercoordsc(tcv, 1, 8);

% sort by distance to center coordinate
cdist = zeros(numel(cllc), 1);
for cc = 1:numel(cllc)
    cds = sqrt(sum((cllc{cc} - repmat(ccrd, size(cllc{cc}, 1), 1)) .^ 2, 2));
    cdist(cc) = min(cds);
end
[cdist, cdo] = sort(cdist);

% replace in clustered volume
cdist(cdo) = 1:numel(cdist);
tcvv(tcvv > 0) = cdist(tcvv(tcvv > 0));

% for each cluster
for cc = 1:numel(cdo)

    % get cluster and check which part of the tail is more "central"
    tcvc = (tcvv == cc);
    t1c = tcvc & (t1v > 0);
    t2c = tcvc & (t2v > 0);
    if ~any(t1c(:)) || ~any(t2c(:))
        continue;
    end

    % construct rays from each voxel to center
    [t1cc, t1cy, t1cz] = ind2sub(size(t1c), find(t1c(:)));
    t1cc = [t1cc(:), t1cy(:), t1cz(:)];
    t1cy = ones(numel(t1cy), 1) * ccrd - t1cc;
    t1cz = sqrt(sum(t1cy .* t1cy, 2));
    t1cy = t1cy ./ (t1cz * ones(1, 3));
    [t2cc, t2cy, t2cz] = ind2sub(size(t2c), find(t2c(:)));
    t2cc = [t2cc(:), t2cy(:), t2cz(:)];
    t2cy = ones(numel(t2cy), 1) * ccrd - t2cc;
    t2cy = t2cy ./ (sqrt(sum(t2cy .* t2cy, 2)) * ones(1, 3));

    % sample toward the center
    t1cz = ceil(t1cz);
    for sc = 1:numel(t1cz)
        svals = diff(indexarraynb(fm, t1cc(sc .* ones(t1cz(sc), 1), :) + (0:(t1cz(sc)-1))' * t1cy(sc, :)));
        svals(abs(svals) >= fmr) = 0;
        t1cz(sc) = mean(svals(svals ~= 0));
    end
    t2cz = ceil(t2cz);
    for sc = 1:numel(t2cz)
        svals = diff(indexarraynb(fm, t2cc(sc .* ones(t2cz(sc), 1), :) + (0:(t2cz(sc)-1))' * t2cy(sc, :)));
        svals(abs(svals) >= fmr) = 0;
        t2cz(sc) = mean(svals(svals ~= 0));
    end

    % median diff
    mdiff = median([t1cz(:); t2cz(:)]);

    % get coordinates of different tails
    tcn1 = double(t1v(t1c));
    utcn1 = unique(tcn1);
    if numel(utcn1) ~= 1
        utcn1 = mode(tcn1);
    end
    tcx1 = (t1v == utcn1);
    tcn2 = double(t2v(t2c));
    utcn2 = unique(tcn2);
    if numel(utcn2) ~= 1
        utcn2 = mode(tcn2);
    end
    tcx2 = (t2v == utcn2);

    % patch
    if mdiff > 0
        fm(tcx1) = fm(tcx1) - fmc;
    elseif mdiff < 0
        fm(tcx2) = fm(tcx2) + fmc;
    end
    t1v(tcx1) = 0;
    t2v(tcx2) = 0;
end

% ultimately, replace data in voxel data and store in new object
nxo = aft_CopyObject(xo);
ovd(fm ~= 0) = fm(fm ~= 0);
nxo.C.VoxelData = ovd;
