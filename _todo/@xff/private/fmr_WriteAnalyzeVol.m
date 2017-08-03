function rvalue = fmr_WriteAnalyzeVol(xo, volnum, filename, flip)
% FMR::WriteAnalyzeVol  - write an Analyze image from one volume
%
% FORMAT:       [success] = fmr.WriteAnalyzeVol(volume, filename [, flip]);
%
% Input fields:
%
%       volume      volume number to write
%       filename    analyze filename
%       flip        char string for flipping (e.g. 'xy', default: '')
%
% Output fields:
%
%       success     true if write was successful

% Version:  v1.1
% Build:    16020311
% Date:     Feb-03 2016, 11:08 AM EST
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
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr')
    error('xff:BadArgument', 'Invalid call to %s.', mfilename);
end

% check further arguments
bc = xo.C;
if nargin < 2 || ~isa(volnum, 'double') || numel(volnum) ~= 1 || isnan(volnum) || ...
    isinf(volnum) || volnum ~= fix(volnum) || volnum < 1 || volnum > bc.NrOfVolumes
    error('neuroelf:xff:badArgument', ...
        'Invalid or missing volnum argument for call to %s.', mfilename);
end
if nargin < 3 || ~ischar(filename) || isempty(filename) || ...
    numel(filename) ~= length(filename) || length(filename) < 5 || ...
  (~strcmpi(filename(end-3:end), '.hdr') && ~strcmpi(filename(end-3:end), '.img'))
    error('neuroelf:xff:badArgument', ...
        'Invalid or missing filename argument for call to %s.', mfilename);
end

% flipping
if nargin < 4 || ~ischar(flip) || numel(flip) > 3
    flip = '';
else
    flip = lower(flip(:)');
end

% default rvalue
rvalue = false;

% get reference
ofile = xo.F;
[opath{1:2}] = fileparts(ofile);
ofile = opath{2};

% slice data loaded at all?
if isempty(bc.Slice) || ~isstruct(bc.Slice) || ~isfield(bc.Slice, 'STCData')
    try
        fmr_LoadSTC(xo);
    catch xfferror
        rethrow(xfferror);
    end
    bc = xo.C;
    if isempty(bc.Slice) || ~isstruct(bc.Slice) || ~isfield(bc.Slice, 'STCData')
        error('neuroelf:xff:invalidObject', 'Slice data cannot be loaded.');
    end
end

% filename
filename = filename(:)';
filename = [filename(1:end-3) 'img'];
hdrfname = [filename(1:end-3) 'hdr'];
matfname = [filename(1:end-3) 'mat'];

% dimension and datatype
xpix = bc.ResolutionX;
ypix = bc.ResolutionY;
zpix = bc.NrOfSlices;
if (bc.DataStorageFormat == 1 && zpix ~= length(bc.Slice)) || ...
   (bc.DataStorageFormat == 2 && zpix ~= size(bc.Slice.STCData, 4))
    error('neuroelf:xff:invalidArraySize', ...
        'Invalid number of slices in header, cannot save data.');
end
isiz = [xpix, ypix, zpix, 1];

% X/Y resolution
xres = bc.InplaneResolutionX;
yres = bc.InplaneResolutionY;
zres = bc.SliceThickness + bc.SliceGap;
tres = bc.TR / 1000;

% positional information
xvec = xres .* [bc.RowDirX; bc.RowDirY; bc.RowDirZ];
yvec = yres .* [bc.ColDirX; bc.ColDirY; bc.ColDirZ];
ovec = [bc.Slice1CenterX; bc.Slice1CenterY; bc.Slice1CenterZ];
lvec = [bc.SliceNCenterX; bc.SliceNCenterY; bc.SliceNCenterZ];
zspan = bc.NrOfSlices - 1;
zvec = (1 / zspan) .* (lvec - ovec);
svec = ovec - ((xpix + 1) / 2) * xvec - ((ypix + 1) / 2) * yvec;

% build mat
tmat = [[xvec, yvec, zvec], svec;  0, 0, 0, 1];

% radiological convention ?
if isfield(bc, 'Convention') && lower(bc.Convention(1)) == 'r'
    tfmat = eye(4);
    tfmat(2, 2) = -1;
    tmat = tfmat * tmat;
end

% volume creation
hdr = xff('new:hdr');
hdrc = hdr.C;

% set dims and data
hdrc.ImgDim.Dim(1:5) = [4, isiz];
if bc.DataType == 1
    hdrc.ImgDim.DataType = 4;
    hdrc.ImgDim.BitsPerPixel = 16;
    hdrc.VoxelData = int16(zeros(isiz));
else
    hdrc.ImgDim.DataType = 16;
    hdrc.ImgDim.BitsPerPixel = 32;
    hdrc.VoxelData = single(zeros(isiz));
end
hdrc.ImgDim.PixSpacing(2:5) = [xres, yres, zres, tres];
hdrc.ImgDim.CalMaxDisplay = 32767;
hdrc.ImgDim.CalMinDisplay = 0;
hdrc.ImgDim.GLMax = 32767;
hdrc.ImgDim.GLMin = 0;
hdrc.DataHist.Description = sprintf('Volume %d of FMR %s', volnum, ofile);
hdrc.DataHist.ScanNumber = sprintf('%d', volnum);

% which format
switch bc.DataStorageFormat

    % old format (BrainVoyager QX <= 1.8.7)
    case 1
        for sc = 1:zpix
            hdrc.VoxelData(:,:,sc,1) = bc.Slice(sc).STCData(:, :, volnum);
        end

    % new format (BrainVoyager QX >= 1.9.9)
    case 2
        hdrc.VoxelData(:, :, :) = squeeze(bc.Slice.STCData(:, :, volnum, :));

    % unsupported
    otherwise
        delete(hdr);
        error('neuroelf:xff:invalidObject', 'Unsupported DataStorageFormat of FMR.');
end
hdr.C = hdrc;

% save hdr/img
try
    rvalue = false;
    aft_SaveAs(hdr, hdrfname);
    rvalue = hdr_SaveVoxelData(hdr, flip);
catch xfferror
    warning('neuroelf:xff:internalError', ...
        'Error writing Analyze header/image ''%s'': %s.', hdrfname, xfferror.message);
end
delete(hdr);

% leave early ?
if ~rvalue
    return;
end

% save mat
M = tmat;
mat = tmat;
try
    save(matfname, 'M', 'mat', '-v6');
catch xfferror
    neuroelf_lasterr(xfferror);
    rvalue = false;
end
if exist(matfname, 'file') ~= 2
    rvalue = false;
end
