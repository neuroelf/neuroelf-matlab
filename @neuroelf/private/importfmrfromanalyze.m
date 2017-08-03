function fmr = importfmrfromanalyze(imgs, targetname, skipvols, tr, flipxy)
% importvtcfromanalyze  - import an FMR from Analzye files
%
% FORMAT:       fmr = importfmrfromanalyze(imgs, targetname [, skipvols, tr, flipxy])
%
% Input fields:
%
%       imgs        cell array with HDR filenames or pattern
%       targetname  target filename (possibly incl. path)
%       skipvols    number of skipped volumes, default: 0
%       tr          TR in milliseconds (default: 2000)
%       flipxy      string, if contains 'x', flipped in X, ..., default: ''
%
% Output fields:
%
%       fmr         created FMR object

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:20 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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
   (~iscell(imgs) && ...
     ~ischar(imgs)) || ...
   ~ischar(targetname) || ...
    isempty(targetname)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad input argument.' ...
    );
end
if ischar(imgs)
    imgs(imgs=='\') = '/';
    [cimgs{1:3}] = fileparts(imgs);
    if isempty(cimgs{1})
        cimgs{1} = pwd;
    end
    try
        imgs = findfiles(cimgs{1}, [cimgs{2} cimgs{3}], 'depth=1');
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        imgs = {};
    end
end
imgs = imgs(:);
nimg = numel(imgs);
if nimg < 3
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid image selection.' ...
    );
end
himg = cell(1, nimg);
vimg = zeros(20, nimg);
try
    for ic = 1:nimg
        if ~ischar(imgs{ic}) || ...
            isempty(imgs{ic})
            error('BAD_IMAGENAME');
        end
        himg{ic} = xff(imgs{ic});
        if numel(himg{ic}) ~= 1 ||...
           ~isxff(himg{ic}, 'hdr')
            error('BAD_IMAGECONT');
        end
        icf = himg{ic}.CoordinateFrame;
        vimg(:, ic) = [icf.Trf(:); icf.Dimensions(:)];
    end
    if any(any(diff(vimg, 1, 2)))
        error( ...
            'neuroelf:BadArgument', ...
            'Spatial orientation/dimensions mismatch between images.' ...
        );
    end
catch ne_eo;
    clearxffobjects(himg);
    error( ...
        'neuroelf:BadArgument', ...
        'Error loading image %d: %s.', ...
        ic, ne_eo.message ...
    );
end
bb = himg{1}.BoundingBox;
resx = bb.DimXYZ(1);
resy = bb.DimXYZ(2);
resxy = max(resx, resy);
sfn = himg{1}.FilenameOnDisk;
targetname = targetname(:)';
[targetpart{1:3}] = fileparts(targetname);
if nargin < 3 || ...
   ~isa(skipvols, 'double') || ...
    numel(skipvols) ~= 1 || ...
    isinf(skipvols) || ...
    isnan(skipvols) || ...
    skipvols < 0 || ...
    skipvols > nimg
    skipvols = 0;
else
    skipvols = round(skipvols);
end
if nargin < 4 || ...
   ~isa(tr, 'double') || ...
    numel(tr) ~= 1 || ...
    isinf(tr) || ...
    isnan(tr) || ...
    tr < 0 || ...
    tr > 30000
    tr = 2000;
else
    tr = round(tr);
end
nrvol = nimg - skipvols;
if nargin < 5 || ...
   ~ischar(flipxy)
    flipxy = ' ';
else
    flipxy = lower(flipxy(:)');
end
flipx = false;
if any(flipxy == 'x')
    flipx = true;
end
flipy = false;
if any(flipxy == 'y')
    flipy = true;
end

% create FMR object
fmr = xff('new:fmr');
fmr.FileVersion = 5;
fmr.NrOfVolumes = nrvol;
fmr.NrOfSlices = bb.DimXYZ(3);
fmr.NrOfSkippedVolumes = skipvols;
fmr.Prefix = targetpart{2};
fmr.DataStorageFormat = 2;
fmr.TR = tr;
fmr.InterSliceTime = floor(tr / fmr.NrOfSlices);
fmr.TimeResolutionVerified = 1;
fmr.TE = round(fmr.InterSliceTime / 2);
fmr.SliceAcquisitionOrder = 0;
fmr.SliceAcquisitionOrderVerified = 0;
fmr.ResolutionX = resxy;
fmr.ResolutionY = resxy;
fmr.LoadAMRFile = '';
fmr.ShowAMRFile = 0;
fmr.ImageIndex = 0;
fmr.LayoutNColumns = ceil(sqrt(fmr.NrOfSlices) + 1);
fmr.LayoutNRows = ceil(fmr.NrOfSlices / fmr.LayoutNColumns);
fmr.LayoutZoomLevel = 1;
fmr.SegmentSize = 10;
fmr.SegmentOffset = 0;
fmr.NrOfLinkedProtocols = 1;
fmr.ProtocolVersion = 2;
fmr.ProtocolFile = '<none>';
fmr.InplaneResolutionX = bb.ResXYZ(1);
fmr.InplaneResolutionY = bb.ResXYZ(2);
fmr.SliceThickness = bb.ResXYZ(3);
fmr.SliceGap = 0;
fmr.VoxelResolutionVerified = 1;
fmr.PosInfosVerified = 0;
fmr.CoordinateSystem = 0;
fmr.Slice1CenterX = 0;
fmr.Slice1CenterY = 0;
fmr.Slice1CenterZ = 0;
fmr.SliceNCenterX = 0;
fmr.SliceNCenterY = 0;
fmr.SliceNCenterZ = 0;
fmr.RowDirX = 0;
fmr.RowDirY = 0;
fmr.RowDirZ = 0;
fmr.ColDirX = 0;
fmr.ColDirY = 0;
fmr.ColDirZ = 0;
fmr.NRows = 0;
fmr.NCols = 0;
fmr.FoVRows = 0;
fmr.FoVCols = 0;
fmr.GapThickness = 0;
fmr.NrOfPastSpatialTransformations = 0;
fmr.Convention = 'Radiological';
fmr.FirstDataSourceFile = [sfn(1:end-3) 'img'];
fmr.Slice = struct;
stcdata = uint16(0);
fmr.Slice.STCData = stcdata;
stcdata(1:resxy, 1:resxy, 1:nrvol, 1:fmr.NrOfSlices) = stcdata(1);

% copy image data
if resx > resy
    fx = 1;
    fy = 1 + floor((resx - resy) / 2);
    tx = resx;
    ty = 1 + resx - fy;
elseif resy > resx
    fx = 1 + floor((resy - resx) / 2);
    fy = 1;
    tx = 1 + resy - fx;
    ty = resy;
else
    fx = 1;
    fy = 1;
    tx = resx;
    ty = resy;
end
if flipx
    ftx = tx:-1:fx;
else
    ftx = fx:tx;
end
if flipy
    fty = ty:-1:fy;
else
    fty = fy:ty;
end
rsh = [bb.DimXYZ(1:2), 1, bb.DimXYZ(3)];
for ic = (skipvols+1):nimg
    thi = himg{ic};
    thh = thi.ImgDim;
    thi = thi.VoxelData;
    if istransio(thi)
        thi = resolve(thi);
    end
    if thh.ScalingIntercept ~= 0 || ...
       ~any(thh.ScalingSlope == [0, 1])
        if thh.ScalingSlope == 0
            thh.ScalingSlope = 1;
        end
        thi = thh.ScalingIntercept + thh.ScalingSlope .* double(thi);
    end
    stcdata(ftx, fty, ic-skipvols, :) = reshape(thi, rsh);
end
clearxffobjects(himg);
fmr.Slice.STCData = stcdata;
try
    fmr.SaveAs(targetname);
    fmr.ClearObject;
    fmr = xff(targetname);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end
