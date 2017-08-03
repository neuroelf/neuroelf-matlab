function xo = hdr_UpdateHeader(xo)
% HDR::UpdateHeader  - update the header of an Analyze/NIftI file
%
% FORMAT:       hdr.UpdateHeader
%
% No input/output fields.
%
% Using: lsqueeze, unzerostring, writemat, zerodstring.

% Version:  v1.1
% Build:    16021216
% Date:     Feb-12 2016, 4:27 PM EST
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

% neuroelf library
global ne_methods;
zerodstring = ne_methods.zerodstring;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'hdr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get HDR file name -> IMG file
bc = xo.C;
fflp = xo.F;

% reject for empty files
if isempty(fflp)
    warning('neuroelf:xff:badArgument', 'Invalid for unsaved files.');
    return;
end

% and also for multi-file transio VoxelData objects
if istransio(bc.VoxelData) && iscell(filename(bc.VoxelData))
    error('neuroelf:xff:badArgument', 'Invalid for multi-file VoxelData objects.');
end

% some pre-write checks (from HDR.bff)
bc.HdrKey.DataType = zerodstring(bc.HdrKey.DataType, 10);
bc.HdrKey.DBName = zerodstring(bc.HdrKey.DBName, 18);
try
    bc.HdrKey.Regular = bc.HdrKey.Regular(1);
catch xfferror
    neuroelf_lasterr(xfferror);
    bc.HdrKey.Regular = 'r';
end
bc.ImgDim.Dim = [4, 1, 1, 1, 1, 0, 0, 0];
ndi = ndims(bc.VoxelData);
if ndi > 4
    bc.ImgDim.Dim(1) = ndi;
end
bc.ImgDim.Dim(2:1+ndims(bc.VoxelData)) = size(bc.VoxelData);
bc.DataHist.Description = zerodstring(bc.DataHist.Description, 80);
bc.DataHist.AuxFilename = zerodstring(bc.DataHist.AuxFilename, 24);
try
    bc.DataHist.Orientation = bc.DataHist.Orientation(1);
catch xfferror;
    neuroelf_lasterr(xfferror);
    bc.DataHist.Orientation = int8(0);
end
if numel(bc.DataHist.OriginString) > 5
    bc.DataHist.Originator = uint8(bc.DataHist.OriginString);
elseif numel(bc.DataHist.OriginSPM) == 5 && ~all(bc.DataHist.OriginSPM == 0)
    dh = min(32767, max(-32768, round(bc.DataHist.OriginSPM(:)')));
    dh(dh < 0) = dh(dh < 0) + 65536;
    bc.DataHist.Originator = [uint8(floor(dh / 256)); uint8(mod(dh, 256))];
    if ~isempty(strfind(lower(bc.Endian), 'le'))
        bc.DataHist.Originator = bc.DataHist.Originator([2, 1], :);
    end
    bc.DataHist.Originator = bc.DataHist.Originator(:)';
end
bc.DataHist.Generated = zerodstring(bc.DataHist.Generated, 10);
bc.DataHist.ScanNumber = zerodstring(bc.DataHist.ScanNumber, 10);
bc.DataHist.PatientID = zerodstring(bc.DataHist.PatientID, 10);
bc.DataHist.ExpiryDate = zerodstring(bc.DataHist.ExpiryDate, 10);
bc.DataHist.ExpiryTime = zerodstring(bc.DataHist.ExpiryTime, 10);
bc.DataHist.Hist_UN0 = zerodstring(bc.DataHist.Hist_UN0, 3);
bc.DataHist.NIftI1.IntentName = zerodstring(bc.DataHist.NIftI1.IntentName, 16);
bc.FileMagic = zerodstring(bc.FileMagic, 4);
bc.FileMagic(4) = 0;
bc.IntermedData = bc.IntermedData(:)';
nimed = numel(bc.IntermedData);

% depending on file magic
switch(ne_methods.unzerostring(bc.FileMagic))

    % ni1 / nii file
    case {'ni1', 'nii'}
        bc.NIIFileType = 1;

    % n+1 file
    case 'n+1'
        bc.NIIFileType = 2;

    % all other cases
    otherwise
        bc.NIIFileType = 0;
end

% check voxel offset
if bc.NIIFileType == 2
    newoffset = max(bc.ImgDim.VoxOffset, 348 + nimed);
    if newoffset > bc.ImgDim.VoxOffset
        error('neuroelf:xff:badArgument', ...
            'Cannot be used to add data to IntermedData field.');
    end
    bc.ImgDim.VoxOffset = newoffset;
end

% open temporary file
tfile = [tempname(tempdir) fflp(end-3:end)];
if isempty(strfind(lower(bc.Endian), 'le'))
    fid = fopen(tfile, 'w', 'ieee-be');
else
    fid = fopen(tfile, 'w', 'ieee-le');
end
if fid < 1
    error('neuroelf:xff:fileWriteError', 'Error opening temporary file.');
end

% open original header file
ofid = fopen(fflp, 'r+', 'ieee-le');
if ofid < 1
    fclose(fid);
    try
        delete(tfile);
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
    error('neuroelf:xff:fileWriteError', ...
        'Error opening target file in read/write mode.');
end

% start writing content
try
    fwrite(fid, bc.HdrKey.SizeOfHdr(1), 'uint32');
    fwrite(fid, bc.HdrKey.DataType(1:10), 'char');
    fwrite(fid, bc.HdrKey.DBName(1:18), 'char');
    fwrite(fid, bc.HdrKey.Extents(1), 'int32');
    fwrite(fid, bc.HdrKey.SessionError(1), 'int16');
    fwrite(fid, bc.HdrKey.Regular(1), 'char');
    fwrite(fid, bc.HdrKey.SliceOrdering(1), 'int8');
    fwrite(fid, bc.ImgDim.Dim(1:8), 'int16');
    fwrite(fid, bc.ImgDim.IntentParam1(1), 'single');
    fwrite(fid, bc.ImgDim.IntentParam2(1), 'single');
    fwrite(fid, bc.ImgDim.IntentParam3(1), 'single');
    fwrite(fid, bc.ImgDim.IntentCode(1), 'int16');
    fwrite(fid, bc.ImgDim.DataType(1), 'int16');
    fwrite(fid, bc.ImgDim.BitsPerPixel(1), 'uint16');
    fwrite(fid, bc.ImgDim.FirstSliceIndex(1) , 'int16');
    fwrite(fid, bc.ImgDim.PixSpacing(1:8), 'single');
    fwrite(fid, bc.ImgDim.VoxOffset(1), 'single');
    fwrite(fid, bc.ImgDim.ScalingSlope(1), 'single');
    fwrite(fid, bc.ImgDim.ScalingIntercept(1), 'single');
    fwrite(fid, bc.ImgDim.LastSliceIndex(1), 'int16');
    fwrite(fid, bc.ImgDim.SliceTimeOrder(1), 'int8');
    fwrite(fid, bc.ImgDim.NIftIUnits(1), 'int8');
    fwrite(fid, bc.ImgDim.CalMaxDisplay(1), 'single');
    fwrite(fid, bc.ImgDim.CalMinDisplay(1), 'single');
    fwrite(fid, bc.ImgDim.SliceDuration(1), 'single');
    fwrite(fid, bc.ImgDim.SliceTimeOffset(1), 'single');
    fwrite(fid, bc.ImgDim.GLMax(1), 'int32');
    fwrite(fid, bc.ImgDim.GLMin(1), 'int32');
    fwrite(fid, bc.DataHist.Description(1:80), 'char');
    fwrite(fid, bc.DataHist.AuxFilename(1:24), 'char');
    if bc.NIIFileType == 0
        fwrite(fid, bc.DataHist.Orientation(1), 'int8');
        fwrite(fid, bc.DataHist.Originator(1:10), 'uint8');
        fwrite(fid, bc.DataHist.Generated(1:10), 'char');
        fwrite(fid, bc.DataHist.ScanNumber(1:10), 'char');
        fwrite(fid, bc.DataHist.PatientID(1:10), 'char');
        fwrite(fid, bc.DataHist.ExpiryDate(1:10), 'char');
        fwrite(fid, bc.DataHist.ExpiryTime(1:10), 'char');
        fwrite(fid, bc.DataHist.Hist_UN0(1:3), 'char');
        fwrite(fid, bc.DataHist.Views(1), 'int32');
        fwrite(fid, bc.DataHist.VolumesAdded(1), 'int32');
        fwrite(fid, bc.DataHist.StartField(1), 'int32');
        fwrite(fid, bc.DataHist.FieldSkip(1), 'int32');
        fwrite(fid, bc.DataHist.OMax(1), 'int32');
        fwrite(fid, bc.DataHist.OMin(1), 'int32');
        fwrite(fid, bc.DataHist.SMax(1), 'int32');
        fwrite(fid, bc.DataHist.SMin(1), 'int32');
    else
        fwrite(fid, bc.DataHist.NIftI1.QFormCode(1), 'int16');
        fwrite(fid, bc.DataHist.NIftI1.SFormCode(1), 'int16');
        fwrite(fid, bc.DataHist.NIftI1.QuaternionB(1), 'single');
        fwrite(fid, bc.DataHist.NIftI1.QuaternionC(1), 'single');
        fwrite(fid, bc.DataHist.NIftI1.QuaternionD(1), 'single');
        fwrite(fid, bc.DataHist.NIftI1.QuatOffsetX(1), 'single');
        fwrite(fid, bc.DataHist.NIftI1.QuatOffsetY(1), 'single');
        fwrite(fid, bc.DataHist.NIftI1.QuatOffsetZ(1), 'single');
        fwrite(fid, bc.DataHist.NIftI1.AffineTransX(1:4), 'single');
        fwrite(fid, bc.DataHist.NIftI1.AffineTransY(1:4), 'single');
        fwrite(fid, bc.DataHist.NIftI1.AffineTransZ(1:4), 'single');
        fwrite(fid, bc.DataHist.NIftI1.IntentName(1:16), 'char');
        fwrite(fid, bc.FileMagic(1:4), 'char');
    end
catch xfferror
    oeo = xfferror;
    fclose(ofid);
    fclose(fid);
    try
        delete(tfile);
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
    rethrow(oeo);
end

% close file
fclose(fid);

% now, re-read this content (no error occurred so far!)
fid = fopen(tfile, 'r');
if fid < 1
    fclose(ofid);
    try
        delete(tfile);
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
    error('neuroelf:xff:fileReadError', 'Error reading back temporary file.');
end
fcont = fread(fid, [1, Inf], 'uint8=>uint8');
fclose(fid);
try
    delete(tfile);
catch xfferror
    neuroelf_lasterr(xfferror);
end

% check length
if numel(fcont) ~= 348
    fclose(ofid);
    error('neuroelf:xff:fileWriteError', ...
        'New header has wrong length; not overwriting.');
end

% write to original file, then close file (flush)
fwrite(ofid, fcont, 'uint8');
fclose(ofid);

% then work on any 4x4(xV) transformation matrix in a MAT file
if isfield(bc.RunTimeVars, 'Mat44') && isa(bc.RunTimeVars.Mat44, 'double') && ...
    size(bc.RunTimeVars.Mat44, 1) == 4 && size(bc.RunTimeVars.Mat44, 2) == 4 && ...
    size(bc.RunTimeVars.Mat44, 3) == size(bc.VoxelData, 4) && ...
   ~any(isinf(ne_methods.lsqueeze(bc.RunTimeVars.Mat44(4, 1:3, :))) | ...
        isnan(ne_methods.lsqueeze(bc.RunTimeVars.Mat44(4, 1:3, :))) | ...
        ne_methods.lsqueeze(bc.RunTimeVars.Mat44(4, 1:3, :)) ~= 0) && ...
    all(bc.RunTimeVars.Mat44(4, 4, :) == 1)

    % MAT file exist
    mflp = [fflp(1:end-4) '.mat'];
    if exist(mflp, 'file') == 2
        try
            mat = load(mflp);
            mat.mat = bc.RunTimeVars.Mat44;
        catch xfferror
            neuroelf_lasterr(xfferror);
            warning('neuroelf:xff:updateIncomplete', 'MAT-file (%s) not updated.', mflp);
            return;
        end
    else
        mat = struct('mat', bc.RunTimeVars.Mat44);
    end

    % save MAT file
    ne_methods.writemat(mflp, '-v7', fieldnames(mat), struct2cell(mat));
end
