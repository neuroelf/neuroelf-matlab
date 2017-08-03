function [fmr] = dicom4tofmr(dcm4file, fmrfile, skipvols)
% dicom4tofmr  - cheap try to convert a DICOM4 file into FMR/STC
%
% FORMAT:       [fmr] = dicom4tofmr(dcm4file, fmrfile [, skipvols])
%
% Input fields:
%
%       dcm4file    filename of DICOM4 image
%       fmrfile     filename of FMR/STC to write (with .FMR extension)
%       skipvols    number of volumes to skip (default: 0)
%
% Output fields:
%
%       fmr         FMR object with STCs loaded

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
   ~ischar(fmrfile) || ...
    isempty(fmrfile)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument for dicom4tofmr.' ...
    );
end
dcm4file = dcm4file(:)';
fmrfile = fmrfile(:)';
if nargin < 3 || ...
   ~isa(skipvols, 'double') || ...
    numel(skipvols) ~= 1 || ...
    isnan(skipvols) || ...
    skipvols < 0 || ...
    skipvols > 999
    skipvols = 0;
else
    skipvols = round(real(skipvols));
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
        'The given file couldn''t be detected as a DICOM4 file (%s).', ...
        ne_eo.message ...
    );
end

% check fmr filename
if numel(fmrfile) < 4 || ...
   ~strcmpi(fmrfile(end-3:end), '.fmr')
    fmrfile = [fmrfile '.fmr'];
end
[fmrp{1:3}] = fileparts(fmrfile);
sliceprefix = fmrp{2};

% create FMR
fmr = xff('new:fmr');

% fill structure
fmr.NrOfVolumes = NrOfVolumes - skipvols;
fmr.NrOfSlices = NrOfSlices;
fmr.NrOfSkippedVolumes = skipvols;
fmr.DataStorageFormat = 2;
fmr.Prefix = sliceprefix;
fmr.InterSliceTime = floor(fmr.TR / NrOfSlices);
fmr.ResolutionX = NrOfRows;
fmr.ResolutionY = NrOfCols;
fmr.LayoutNRows = ceil(NrOfSlices / fmr.LayoutNColumns);
fmr.InplaneResolutionX = PixSpcX;
fmr.InplaneResolutionY = PixSpcY;
fmr.SliceThickness = SlcThick;
fmr.SliceGap = GapThick;
fmr.Slice1CenterZ = -(NrOfSlices - 1) * (SlcThick + GapThick) / 2;
fmr.SliceNCenterZ = -fmr.Slice1CenterZ;
fmr.NRows = NrOfRows;
fmr.NCols = NrOfCols;
fmr.FoVRows = NrOfRows * PixSpcX;
fmr.FoVCols = NrOfCols * PixSpcY;
fmr.GapThickness = GapThick;
fmr.FirstDataSourceFile = dcm4file;

% try saving
try
    fmr.SaveAs(fmrfile);
catch ne_eo;
    clearxffobjects({fmr, dcm});
    error( ...
        'neuroelf:xffError', ...
        'Error saving FMR file to disk: %s.', ...
        ne_eo.message ...
    );
end

% creating slice files
try
    stc = 0;
    stc = fopen([fmrfile(1:end-4) '.stc'], 'w', 'ieee-le');
    if stc < 1
        error('ERROR_WRITING_STC');
    end
    fwrite(stc, dcm.Data(end).Value, 'uint16');
    fclose(stc);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    clearxffobjects({fmr, dcm});
    if stc > 0
        fclose(stc);
    end
    error( ...
        'neuroelf:xffError', ...
        'Error saving STC data.' ...
    );
end

% create a firstvol fmr
firstvol = fmr.CopyObject;

% reload slices
try
    fmr = fmr.LoadSTC;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    clearxffobjects({firstvol, fmr, dcm});
    error( ...
        'neuroelf:xffError', ...
        'Error test re-loading slice data.' ...
    );
end

% edit firstvol
firstvol.NrOfVolumes = 1;
firstvol.NrOfSkippedVolumes = 0;
firstvol.Prefix = [sliceprefix '_firstvol'];
try
    firstvol.SaveAs([fmrfile(1:end-4) '_firstvol.fmr']);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    clearxffobjects({firstvol, fmr, dcm});
    warning( ...
        'neuroelf:xffError', ...
        'Error saving _firstvol.fmr' ...
    );
    return;
end

% create firstvol slice files
try
    stc = 0;
    stc = fopen([fmrfile(1:end-4) '_firstvol.stc'], 'w', 'ieee-le');
    if stc < 1
        error('ERROR_WRITING_STC');
    end
    fwrite(stc, dcm.Data(end).Value(:, :, 1, :), 'uint16');
    fclose(stc);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    clearxffobjects({firstvol, fmr, dcm});
    if stc > 0
        fclose(stc);
    end
    error( ...
        'neuroelf:xffError', ...
        'Error saving _firstvol slice %d.', ...
        sc ...
    );
end

% clear unrequired object
clearxffobjects({firstvol, dcm});

% also clear FMR if no output
if nargout < 1
    fmr.ClearObject;
end
