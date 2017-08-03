function dmr = hdr_Dyn3DToDMR(xo, vols)
% HDR::Dyn3DToDMR  - convert a dynamics 3D Analyze into DMR/DWI
%
% FORMAT:       [dmr] = hdr.Dyn3DToDMR([vols])
%
% Input fields:
%
%       vols        range of volumes in file (default: all)
%
% Output fields:
%
%       dmr         created DMR object (with data loaded)

% Version:  v1.1
% Build:    16020515
% Date:     Feb-05 2016, 3:52 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2011, 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'hdr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
hdrc = xo.C;
hdrf = hdr_CoordinateFrame(xo);
hdrt = hdrf.Trf;
[hdrpath, dwipref] = fileparts(xo.F);
if isempty(dwipref)
    dwipref = 'unsaved';
end
if nargin < 2 || ~isa(vols, 'double') || numel(vols) ~= max(size(vols)) || ...
    any(isinf(vols) | isnan(vols) | vols < 1 | vols > size(hdrc.VoxelData, 4) | vols ~= fix(vols))
    vols = 1:size(hdrc.VoxelData, 4);
else
    vols = unique(vols(:))';
end

% get settings from HDR and make sure the voxel data is loaded
iImg = hdrc.ImgDim;
iDim = size(hdrc.VoxelData);
iRes = iImg.PixSpacing(2:4);
dimX = iDim(1);
dimY = iDim(2);
dimZ = iDim(3);
resX = iRes(1);
resY = iRes(2);
resZ = iRes(3);
sl1c = hdrt * [0.5 .* (dimX + 1); 0.5 .* (dimY + 1); 1; 1];
slnc = hdrt * [0.5 .* (dimX + 1); 0.5 .* (dimY + 1); dimZ; 1];

% make some initial DMR project settings
dmr = xff('new:dmr');
dmrc = dmr.C;
dmrc.FileVersion = 3;
dmrc.DataStorageFormat = 3;
dmrc.DataType = 2;
dmrc.NrOfVolumes = numel(vols);
dmrc.NrOfSlices = dimZ;
dmrc.Prefix = dwipref;
dmrc.InterSliceTime = fix(dmrc.TR / dimZ);
dmrc.ResolutionX = dimX;
dmrc.ResolutionY = dimY;
dmrc.LayoutNColumns = ceil(sqrt(dimZ)) + 1;
dmrc.LayoutNRows = ceil(dimZ / dmrc.LayoutNColumns);
dmrc.InplaneResolutionX = resX;
dmrc.InplaneResolutionY = resY;
dmrc.SliceThickness = resZ;
dmrc.SliceGap = 0;
dmrc.Slice1CenterX = sl1c(1);
dmrc.Slice1CenterY = sl1c(2);
dmrc.Slice1CenterZ = sl1c(3);
dmrc.SliceNCenterX =  slnc(1);
dmrc.SliceNCenterY =  slnc(2);
dmrc.SliceNCenterZ =  slnc(3);
dmrc.RowDirX = hdrt(1, 1);
dmrc.RowDirY = hdrt(2, 1);
dmrc.RowDirZ = hdrt(3, 1);
dmrc.ColDirX = hdrt(1, 2);
dmrc.ColDirY = hdrt(2, 2);
dmrc.ColDirZ = hdrt(3, 2);
dmrc.NRows = dimX;
dmrc.NCols = dimY;
dmrc.FoVRows = dimX * resX;
dmrc.FoVCols = dimY * resY;
dmrc.GapThickness = 0;
dmrc.Convention = 'Neurological';
dmrc.FirstDataSourceFile = hdrs.F;
dmrc.GradientDirectionsVerified = 'NO';
dmrc.DWIData = single(0);
dmrc.DWIData(dimX, dimY, dimZ, numel(vols)) = 0;

% scaling
slope = iImg.ScalingSlope;
offset = iImg.ScalingIntercept;
if numel(slope) ~= 1 || isinf(slope) || isnan(slope) || slope == 0
    slope = 1;
end
if numel(offset) ~= 1 || isinf(offset) || isnan(offset)
    offset = 0;
end

% fill slice data
for vc = 1:numel(vols)
    if slope ~= 1 || offset ~= 0
        dmrc.DWIData(:, :, :, vc) = offset + slope .* double( ...
            hdrc.VoxelData(:, :, :, vols(vc)));
    else
        dmrc.DWIData(:, :, :, vc) = hdrc.VoxelData(:, :, :, vols(vc));
    end
end

% set content
dmr.C = dmrc;
