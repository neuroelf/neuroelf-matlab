function fmr = hdr_Dyn3DToFMR(xo, vols)
% HDR::Dyn3DToFMR  - convert a dynamics 3D Analyze into FMR/STC
%
% FORMAT:       [fmr] = hdr.Dyn3DToFMR([vols])
%
% Input fields:
%
%       vols        range of volumes in file (default: all)
%
% Output fields:
%
%       fmr         created FMR object (with slices loaded)

% Version:  v1.1
% Build:    16020515
% Date:     Feb-05 2016, 3:50 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'hdr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
hdrc = xo.C;
hdrf = hdr_CoordinateFrame(xo);
hdrt = hdrf.Trf;
[hdrpath, stcpref] = fileparts(xo.F);
if isempty(stcpref)
    stcpref = 'unsaved';
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

% make some initial FMR project settings
fmr = xff('new:fmr');
fmrc = fmr.C;
fmrc.FileVersion = 6;
fmrc.DataStorageFormat = 2;
fmrc.DataType = 2;
fmrc.NrOfVolumes = numel(vols);
fmrc.NrOfSlices = dimZ;
fmrc.Prefix = stcpref;
fmrc.TR = 2000;
fmrc.InterSliceTime = fix(fmrc.TR / dimZ);
fmrc.ResolutionX = dimX;
fmrc.ResolutionY = dimY;
fmrc.LayoutNColumns = ceil(sqrt(dimZ)) + 1;
fmrc.LayoutNRows = ceil(dimZ / fmrc.LayoutNColumns);
fmrc.InplaneResolutionX = resX;
fmrc.InplaneResolutionY = resY;
fmrc.SliceThickness = resZ;
fmrc.SliceGap = 0;
fmrc.Slice1CenterX = sl1c(1);
fmrc.Slice1CenterY = sl1c(2);
fmrc.Slice1CenterZ = sl1c(3);
fmrc.SliceNCenterX =  slnc(1);
fmrc.SliceNCenterY =  slnc(2);
fmrc.SliceNCenterZ =  slnc(3);
fmrc.RowDirX = hdrt(1, 1);
fmrc.RowDirY = hdrt(2, 1);
fmrc.RowDirZ = hdrt(3, 1);
fmrc.ColDirX = hdrt(1, 2);
fmrc.ColDirY = hdrt(2, 2);
fmrc.ColDirZ = hdrt(3, 2);
fmrc.NRows = dimX;
fmrc.NCols = dimY;
fmrc.FoVRows = dimX * resX;
fmrc.FoVCols = dimY * resY;
fmrc.GapThickness = 0;
fmrc.Convention = 'Neurological';
fmrc.FirstDataSourceFile = xo.F;
fmrc.Slice = struct('STCData', single(0));
fmrc.Slice.STCData(dimX, dimY, numel(vols), dimZ) = 0;

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
        fmrc.Slice.STCData(:, :, vc, :) = offset + slope .* double(reshape( ...
            hdrc.VoxelData(:, :, :, vols(vc)), [dimX, dimY, 1, dimZ]));
    else
        fmrc.Slice.STCData(:, :, vc, :) = reshape( ...
            hdrc.VoxelData(:, :, :, vols(vc)), [dimX, dimY, 1, dimZ]);
    end
end

% set content
fmr.C = fmrc;
