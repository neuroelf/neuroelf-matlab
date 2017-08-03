function [dmr] = dicom4todmr(dcm4file, dmrfile, skipvols)
% dicom4todmr  - cheap try to convert a DICOM4 file into DMR/DWI
%
% FORMAT:       [dmr] = dicom4todmr(dcm4file, dmrfile [, skipvols])
%
% Input fields:
%
%       dcm4file    filename of DICOM4 image
%       dmrfile     filename of DMR/DWI to write (with .DMR extension)
%       skipvols    volumes to skip (default: empty [])
%
% Output fields:
%
%       dmr         DMR object with DWI loaded

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:19 AM EST
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
   ~ischar(dcm4file) || ...
    isempty(dcm4file) || ...
    exist(dcm4file(:)', 'file') ~= 2 || ...
   ~ischar(dmrfile) || ...
    isempty(dmrfile)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument for dicom4todmr.' ...
    );
end
dcm4file = dcm4file(:)';
dmrfile = dmrfile(:)';
if nargin < 3 || ...
   ~isa(skipvols, 'double') || ...
    isempty(skipvols) || ...
    any(isinf(skipvols(:)) | isnan(skipvols(:)) | skipvols(:) < 1 | skipvols(:) > 999) || ...
    any(skipvols(:) ~= round(skipvols(:))) || ...
    numel(skipvols) ~= numel(unique(skipvols(:)))
    skipvols = [];
else
    skipvols = skipvols(:)';
end

% try to read dicom file
dcm = [];
try
    dcm = xff(dcm4file);
    if ~isxff(dcm, 'dcm') || ...
        numel(dcm.Data) < 16 || ...
       ~all(dcm.Data(end).Key == [32736, 16])
        error('BAD_DICOM');
    end
    NrOfPixels = numel(dcm.Data(end).Value);
    NrOfRows = double(dcm.Value('Rows'));
    NrOfCols = double(dcm.Value('Columns'));
    PixelSpacing = dcm.Value('PixelSpacing');
    if ischar(PixelSpacing)
        PixelSpacing = splittocell(PixelSpacing, '\');
        PixSpcX = str2double(PixelSpacing{1});
        PixSpcY = str2double(PixelSpacing{2});
    else
        PixSpcX = double(PixelSpacing(1));
        PixSpcY = double(PixelSpacing(2));
    end
    SlcThick = dcm.Value('SliceThickness');
    if ischar(SlcThick)
        SlcThick = str2double(SlcThick);
    end
    GapThick = dcm.Value('SpacingBetweenSlices');
    if ischar(GapThick)
        GapThick = str2double(GapThick);
    end
    GapThick = GapThick - SlcThick;
    try
        NrOfSlices = double(dcm.Value('NrOfSlices'));
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        try
            NrOfSlices = double(dcm.Value(8193, 4120));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            error('MISSING_NROFSLICES');
        end
    end
    NrOfVolumes = NrOfPixels / (NrOfRows * NrOfCols * NrOfSlices);
    if NrOfVolumes ~= fix(NrOfVolumes)
        error('BAD_NROFVOLUMES');
    end
    if skipvols > NrOfVolumes
        skipvols = 0;
    end
    dcm.Data(end).Value = reshape(dcm.Data(end).Value, ...
        [NrOfRows, NrOfCols, NrOfVolumes, NrOfSlices]);
catch ne_eo;
    clearxffobjects({dcm});
    error( ...
        'neuroelf:BadArgument', ...
        'The given file couldn''t be detected as a DICOM4 file: %s.', ...
        ne_eo.message ...
    );
end

% check dmr filename
if numel(dmrfile) < 4 || ...
   ~strcmpi(dmrfile(end-3:end), '.dmr')
    dmrfile = [dmrfile '.dmr'];
end
[dmrp{1:3}] = fileparts(dmrfile);
sliceprefix = dmrp{2};

% create DMR
dmr = xff('new:dmr');

% fill structure
dmr.NrOfVolumes = NrOfVolumes - numel(skipvols);
dmr.NrOfSlices = NrOfSlices;
dmr.NrOfSkippedVolumes = numel(skipvols);
dmr.DataStorageFormat = 3;
dmr.Prefix = sliceprefix;
dmr.InterSliceTime = floor(dmr.TR / NrOfSlices);
dmr.ResolutionX = NrOfRows;
dmr.ResolutionY = NrOfCols;
dmr.LayoutNRows = ceil(NrOfSlices / dmr.LayoutNColumns);
dmr.InplaneResolutionX = PixSpcX;
dmr.InplaneResolutionY = PixSpcY;
dmr.SliceThickness = SlcThick;
dmr.SliceGap = GapThick;
dmr.Slice1CenterZ = -(NrOfSlices - 1) * (SlcThick + GapThick) / 2;
dmr.SliceNCenterZ = -dmr.Slice1CenterZ;
dmr.NRows = NrOfRows;
dmr.NCols = NrOfCols;
dmr.FoVRows = NrOfRows * PixSpcX;
dmr.FoVCols = NrOfCols * PixSpcY;
dmr.GapThickness = GapThick;
dmr.FirstDataSourceFile = dcm4file;

% try saving
try
    dmr.SaveAs(dmrfile);
catch ne_eo;
    clearxffobjects({dmr, dcm});
    error( ...
        'neuroelf:xffError', ...
        'Error saving DMR file to disk: %s.', ...
        ne_eo.message ...
    );
end

% creating slice files
try
    dwi = 0;
    dwi = fopen([dmrfile(1:end-4) '.dwi'], 'w', 'ieee-le');
    if dwi < 1
        error('ERROR_WRITING_DWI');
    end
    for vc = 1:NrOfVolumes
        if ~any(skipvols == vc)
            fwrite(dwi, dcm.Data(end).Value(:, :, vc, :), 'uint16');
        end
    end
    fclose(dwi);
catch ne_eo;
    clearxffobjects({dmr, dcm});
    if dwi > 0
        fclose(dwi);
    end
    error( ...
        'neuroelf:xffError', ...
        'Error saving DWI data: %s.', ...
        ne_eo.message ...
    );
end

% clear unrequired object
dcm.ClearObject;

% also clear DMR if no output
if nargout < 1
    dmr.ClearObject;
end
