function dmr = fmr_ConvertToDMR(xo, dmrfile)
% FMR::ConvertToDMR  - convert an FMR to DMR
%
% FORMAT:       dmr = fmr.ConvertToDMR(dmrfile)
%
% Input fields:
%
%       dmrfile     name of output DMR file
%
% Output fields:
%
%       dmr         created DMR object
%
% Note: the DWI data file will be created and then linked as transio.

% Version:  v1.1
% Build:    16020112
% Date:     Feb-01 2016, 12:11 PM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr') || ...
   ~ischar(dmrfile) || isempty(dmrfile)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
dmrfile = dmrfile(:)';

% try DMR save first
dmr = xff('new:dmr');
try
    aft_SaveAs(dmr, dmrfile);
catch xfferror
    delete(dmr);
    error('neuroelf:xff:saveFailed', 'Error saving DMR file: %s.', xfferror.message);
end

% get filename particles
[fnp{1:3}] = fileparts(dmrfile);
if isempty(fnp{1})
    fnp{1} = '.';
end
dwifile = [fnp{1} '/' fnp{2} '.dwi'];

% get FMR content
bc = xo.C;

% get resolution
res = [bc.ResolutionX, bc.ResolutionY, bc.NrOfSlices, bc.NrOfVolumes];

% get DWI data from Slice
if bc.DataStorageFormat == 1
    dwid = uint16(zeros(res));
    for slc = 1:res(3)
        dwid(:, :, slc, :) = reshape(bc.Slice(slc).STCData(:, :, :), [res(1:2), 1, res(4)]);
    end
else
    dwid = permute(bc.Slice.STCData(:, :, :, :), [1, 2, 4, 3]);
end

% create temporary DWI
dwi = xff('new:dwi');
dwic = dwi.C;
dwic.NrOfVolumes = res(4);
dwic.NrOfSlices = res(3);
dwic.ResolutionX = res(1);
dwic.ResolutionY = res(2);
dwic.DWIData = dwid;
dwi.C = dwic;
try
    aft_SaveAs(dwi, dwifile);
catch xfferror
    delete(dwi);
    error('neuroelf:xff:saveError', 'Error saving DWI data to file: %s', xfferror.message);
end
delete(dwi);

% copy fields to DMR
dmrc = dmr.C;
dmrc.NrOfVolumes = bc.NrOfVolumes;
dmrc.NrOfSlices = bc.NrOfSlices;
dmrc.NrOfSkippedVolumes = bc.NrOfSkippedVolumes;
dmrc.Prefix = fnp{2};
dmrc.DataType = bc.DataType;
dmrc.DataStorageFormat = 3;
dmrc.TR = bc.TR;
dmrc.InterSliceTime = bc.InterSliceTime;
dmrc.TimeResolutionVerified = bc.TimeResolutionVerified;
dmrc.TE = bc.TE;
dmrc.SliceAcquisitionOrder = bc.SliceAcquisitionOrder;
dmrc.SliceAcquisitionOrderVerified = bc.SliceAcquisitionOrderVerified;
dmrc.ResolutionX = bc.ResolutionX;
dmrc.ResolutionY = bc.ResolutionY;
dmrc.LoadAMRFile = bc.LoadAMRFile;
dmrc.ShowAMRFile = bc.ShowAMRFile;
dmrc.ImageIndex = bc.ImageIndex;
dmrc.LayoutNColumns = bc.LayoutNColumns;
dmrc.LayoutNRows = bc.LayoutNRows;
dmrc.LayoutZoomLevel = bc.LayoutZoomLevel;
dmrc.SegmentSize = bc.SegmentSize;
dmrc.SegmentOffset = bc.SegmentOffset;
dmrc.ProtocolFile = bc.ProtocolFile;
dmrc.InplaneResolutionX = bc.InplaneResolutionX;
dmrc.InplaneResolutionY = bc.InplaneResolutionY;
dmrc.SliceThickness = bc.SliceThickness;
dmrc.SliceGap = bc.SliceGap;
dmrc.VoxelResolutionVerified = bc.VoxelResolutionVerified;
dmrc.PosInfosVerified = bc.PosInfosVerified;
dmrc.CoordinateSystem = bc.CoordinateSystem;
dmrc.Slice1CenterX = bc.Slice1CenterX;
dmrc.Slice1CenterY = bc.Slice1CenterY;
dmrc.Slice1CenterZ = bc.Slice1CenterZ;
dmrc.SliceNCenterX = bc.SliceNCenterX;
dmrc.SliceNCenterY = bc.SliceNCenterY;
dmrc.SliceNCenterZ = bc.SliceNCenterZ;
dmrc.RowDirX = bc.RowDirX;
dmrc.RowDirY = bc.RowDirY;
dmrc.RowDirZ = bc.RowDirZ;
dmrc.ColDirX = bc.ColDirX;
dmrc.ColDirY = bc.ColDirY;
dmrc.ColDirZ = bc.ColDirZ;
dmrc.NRows = bc.NRows;
dmrc.NCols = bc.NCols;
dmrc.FoVRows = bc.FoVRows;
dmrc.FoVCols = bc.FoVCols;
dmrc.GapThickness = bc.GapThickness;
dmrc.NrOfPastSpatialTransformations = bc.NrOfPastSpatialTransformations;
dmrc.Trf = bc.Trf;
dmrc.Convention = bc.Convention;
dmrc.FirstDataSourceFile = bc.FirstDataSourceFile;
dmrc.GradientDirectionsVerified = 'NO';
dmrc.GradientInformationAvailable = 'NO';

% set content back
dmr.C = dmrc;

% save
aft_SaveAs(dmr, dmrfile);
try
    aft_ReloadFromDisk(dmr);
catch xfferror
    warning('neuroelf:xff:reloadError', ...
        'Error re-loading DMR from disk (for DWI access): %s', xfferror.message);
end
