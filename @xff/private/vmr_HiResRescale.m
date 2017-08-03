function xo2 = vmr_HiResRescale(xo, res, imeth, cutoff)
% VMR::HiResRescale  - upsample ISO-Voxel VMR to higher resolution
%
% FORMAT:       hresvmr = vmr.HiResRescale([res [, imeth [, cutoff]]])
%
% Input fields:
%
%       res         resolution, either {0.5} or 0.25 (in mm)
%       imeth       interpolation: {'linear'}, 'cubic', 'lanczos3'
%       cutoff      if given and false, do not cutoff all-zero planes
%
% Output fields:
%
%       hresvmr     HiRes iso-voxelated VMR
%
% Note: the framing cube will be chosen automatically:
%        512 for 0.5mm
%       1024 for 0.25mm
%
% Also: the function only supports VMR that are ISO-voxelated
%
% Using: flexinterpn_method.

% Version:  v1.1
% Build:    16021315
% Date:     Feb-13 2016, 3:53 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% argument check
if nargin < 1 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
if ~isfield(bc, 'VoxResX') || any([bc.VoxResX, bc.VoxResY, bc.VoxResZ] ~= 1) || ...
    bc.FramingCube ~= 256
    error('neuroelf:xff:invalidObject', ...
        'Only defined for 1mm ISOvoxel VMRs in a 256 FramingCube.');
end
if nargin < 2 || ~isa(res, 'double') || numel(res) ~= 1 || ...
    isinf(res) || isnan(res) || ~any(res == [0.25, 0.5])
    res = 0.5;
end
if nargin < 3 || ~ischar(imeth) || isempty(imeth) || ...
   ~any(strcmpi(imeth(:)', {'cubic', 'lanczos3', 'lanczos8', 'linear', 'nearest'}))
    imeth = 'linear';
else
    imeth = lower(imeth(:)');
end
if res < 0.5
    fcube = 1024;
else
    fcube = 512;
end
fch = fcube / 2;
if nargin < 4 || ~islogical(cutoff) || isempty(cutoff) || cutoff(1)
    cutoff = true;
else
    cutoff = false;
end
rtv = bc.RunTimeVars;
if isfield(rtv, 'SPMsn') && isstruct(rtv.SPMsn) && numel(rtv.SPMsn) == 1 && ...
    isfield(rtv.SPMsn, 'VG') && isfield(rtv.SPMsn, 'VF') && ...
    isfield(rtv.SPMsn, 'Tr') && isfield(rtv.SPMsn, 'Affine') && ...
    isstruct(rtv.SPMsn.VG) && isstruct(rtv.SPMsn.VF) && ...
    numel(rtv.SPMsn.VF) == 1 && ndims(rtv.SPMsn.Tr) == 4 && ...
    isequal(size(rtv.SPMsn.Affine), [4, 4])
    spmsn = rtv.SPMsn;
else
    spmsn = [];
end

% for cutoff
if cutoff

    % create temporary object
    xot = aft_CopyObject(xo);

    % and use MinBox
    vmr_MinBox(xo, 1, 1);

    % and get new data
    bc = xot.C;

    % and clear object
    delete(xot);
end

% which data
if ~isempty(bc.VMRData16) && isequal(size(bc.VMRData16), size(bc.VMRData))
    sdata = bc.VMRData16(:, :, :);
else
    sdata = bc.VMRData(:, :, :);
end
sds = size(sdata);

% find source offsets
if cutoff
    mxx = max(max(sdata, [], 2), [], 3);
    mxy = max(max(sdata, [], 1), [], 3);
    mxz = max(max(sdata, [], 1), [], 2);
    sox = find(mxx > 0);
    soy = find(mxy > 0);
    soz = find(mxz > 0);
else
    sox = (1:size(sdata, 1))';
    soy = (1:size(sdata, 2))';
    soz = (1:size(sdata, 3))';
end

% copy VMR, but set some fields to better values...
ndata = sdata(1);
ndata(1) = 0;
xo2 = aft_CopyObject(xo);
xo2.F = '';
xo2.C.FramingCube = fcube;
xo2.C.VoxResX = res;
xo2.C.VoxResY = res;
xo2.C.VoxResZ = res;
for fn = {'RenderBBox', 'RenderBBoxFull', 'SliceRanges', 'SPMsn', 'Trf', 'TrfPlus'}
    if isfield(xo2.C.RunTimeVars, fn{1})
        xo2.C.RunTimeVars = rmfield(xo2.C.RunTimeVars, fn{1});
    end
end

% empty VMR?
if cutoff
    xo2.C.VMRData = ndata;
    xo2.C.VMRData(1:4, 1:4, 1:4) = 0;
    xo2.C.OffsetX = fch - 2;
    xo2.C.OffsetY = fch - 2;
    xo2.C.OffsetZ = fch - 2;
    xo2.C.DimX = 4;
    xo2.C.DimY = 4;
    xo2.C.DimZ = 4;
else
    xo2.C.VMRData = ndata;
    xo2.C.VMRData(1:fcube, 1:fcube, 1:fcube) = 0;
    xo2.C.OffsetX = 0;
    xo2.C.OffsetY = 0;
    xo2.C.OffsetZ = 0;
    xo2.C.DimX = fcube;
    xo2.C.DimY = fcube;
    xo2.C.DimZ = fcube;
end
if isempty(sox)
    return;
end

% generate sampling grid
xsmp = res / bc.VoxResX;
ysmp = res / bc.VoxResY;
zsmp = res / bc.VoxResZ;
if fcube == 512
    ffac = 2;
else
    ffac = 4;
end
xs = 1 + xsmp .* (0:(ffac * sds(1) - 0.5));
ys = 1 + ysmp .* (0:(ffac * sds(2) - 0.5));
zs = 1 + zsmp .* (0:(ffac * sds(3) - 0.5));
xg = find(xs >= (sox(1) - 0.5) & xs <= (sox(end) + 0.5));
yg = find(ys >= (soy(1) - 0.5) & ys <= (soy(end) + 0.5));
zg = find(zs >= (soz(1) - 0.5) & zs <= (soz(end) + 0.5));

% generate output
ndata(1:numel(xs), 1:numel(ys), 1:numel(zs)) = ndata(1);
xs = [Inf; xs(xg(1)); xs(2) - xs(1); xs(xg(end))];
ys = [Inf; ys(yg(1)); ys(2) - ys(1); ys(yg(end))];
r = [xs, ys, ys];

% determine progress bar capabilities
try
    pbar = xprogress;
    xprogress(pbar, 'setposition', [80, 200, 640, 36]);
    xprogress(pbar, 'settitle', sprintf('%s interpolating VMR...', imeth));
    xprogress(pbar, 0, 'Building interpolation grid...', 'visible', 0, 1);
catch xfferror
    neuroelf_lasterr(xfferror);
    pbar = [];
end

% iterate over Z grid
pbv = 0;
if ~isempty(pbar)
    xprogress(pbar, 0, 'Interpolating...', 'visible', 0, numel(zg));
    pns = 1 / 86400;
    pnx = now + pns;
end
flexinterpn_method = ne_methods.flexinterpn_method;
for zc = zg(:)'
    r([2, 4], 3) = zs(zc);
    pbv = pbv + 1;
    if isempty(spmsn)
        ndata(xg, yg, zc) = flexinterpn_method(sdata, r, imeth);
    else
        rf = r;
        rf([2, 4], :) = (129 + repmat([bc.OffsetZ, bc.OffsetX, bc.OffsetY], 2, 1)) - r([2, 4], [3, 1, 2]);
        rf(3, :) = -r(3, [3, 1, 2]);
        ndata(xg, yg, zc) = aft_SampleData3D(xo, rf, struct('method', imeth, 'snmat', spmsn));
    end
    if ~isempty(pbar) && pnx < now
        xprogress(pbar, pbv);
        pnx = now + pns;
    end
end
if ~isempty(pbar)
    closebar(pbar);
end

% set into new VMR
xo2.C.OffsetX = ffac * bc.OffsetX;
xo2.C.OffsetY = ffac * bc.OffsetY;
xo2.C.OffsetZ = ffac * bc.OffsetZ;
xo2.C.DimX = size(ndata, 1);
xo2.C.DimY = size(ndata, 2);
xo2.C.DimZ = size(ndata, 3);

% simple if only 8-bit data
if ~xo2.C.VMR8bit || isa(ndata, 'uint8')
    xo2.C.VMRData = ndata;
    return;
else
    xo2.C.VMRData16 = ndata;
end

% otherwise VMR data has to be reconstructed
nh = cumsum(hist(single(ndata(:)), single(0:max(ndata(:)))));
ht = find(nh > 0.999 * numel(ndata));
if isempty(ht)
    ht = numel(nh);
end
ndata(ndata > ht(1)) = ht(1);
xo2.C.VMRData = uint8(225 * (single(ndata) ./ single(max(ndata(:)))));
