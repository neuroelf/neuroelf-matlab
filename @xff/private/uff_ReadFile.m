function imdata = uff_ReadFile(xo, filename, numimages)
% UFF::ReadFile  - read a UFF specified file
%
% FORMAT:       imdata = uff.ReadFile(filename [, numimages])
%
% Input fields:
%
%       filename    filename of image data file
%       numimages   number of images to read from file
%
% Output fields:
%
%       imdata      image data (N-D array)

% Version:  v1.1
% Build:    16020314
% Date:     Feb-03 2016, 2:46 PM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'uff') || ...
   ~ischar(filename) || isempty(filename) || exist(filename(:)', 'file') ~= 2
    error('neuroelf:xff:badArgument', 'Invalid call to %s or file not found.', mfilename);
end
ui = xo.C;
if nargin < 3 || ~isa(numimages, 'double') || numel(numimages) ~= 1 || ...
    isinf(numimages) || isnan(numimages) || numimages < 1
    numimages = Inf;
else
    numimages = fix(real(numimages));
end

% get shortcuts
sh = ui.SubHeaderSize;

% reject DICOM files
if ui.IsDICOM
    error('neuroelf:xff:fileTypeNotSupported', 'Please use a real DICOM reader.');
end

% open file according to syntax
filename = filename(:)';
try
    if ui.IsBigEndian
        fid = fopen(filename, 'rb', 'ieee-be');
    else
        fid = fopen(filename, 'rb', 'ieee-le');
    end
    if fid < 1
        error('INVALID_FID');
    end
catch xfferror
    neuroelf_lasterr(xfferror);
    error('neuroelf:xff:fileOpenError', 'Error opening specified file.');
end

% get file size
fseek(fid, 0, 1);
flen = ftell(fid);

% position IO pointer after header
fseek(fid, ui.HeaderSize, -1);

% element size
switch (ui.PixelFormat(1))
    case 1
        rclass = '*uint8';
        rcsize = 1;
        imdata = uint8(0);
    case 2
        rclass = '*uint16';
        rcsize = 2;
        imdata = uint16(0);
    case 3
        rclass = '*uint32';
        rcsize = 4;
        imdata = uint32(0);
    case 4
        rclass = '*single';
        rcsize = 4;
        imdata = single(0);
    otherwise
        error('neuroelf:xff:invalidValue', 'Invalid PixelFormat value.');
end

% slice size
noc = ui.NrOfCols;
nor = ui.NrOfRows;
sls = noc * nor * rcsize;
fdt = ui.FirstDimIsTime;

% if time dim is first, calculate size
if any([1, 2] == ui.SingleFuncType(1))
    tdim = floor((flen - ui.HeaderSize) / sls);
else
    tdim = 1;
end

% skip images ?
if ~isempty(ui.ImageIndex) && fix(ui.ImageIndex(1)) > 1 && fdt == 0
    for sc = 2:fix(ui.ImageIndex)
        try
            fseek(fid, sh + sls, 0);
        catch xfferror
            warning('neuroelf:xff:fileReadError', ...
                'Error skipping slice images UFF file: %s.', xfferror.message);
            return;
        end
    end
end

% read slices (time is not first)
if fdt == 0

    % estimate number of slices if needed
    if isinf(numimages)
        sm = 1 + ceil((flen - ui.HeaderSize) / (sh + sls));
    else
        sm = numimages;
    end

    % initialize imdata
    imdata(1:nor, 1:noc, 1:sm) = imdata(1);

    sc = 1;
    while numimages > 0
        try
            fseek(fid, sh, 0);
            imdata(:, :, sc) = fread(fid, [nor, noc], rclass);
            sc = sc + 1;
        catch xfferror
            warning('neuroelf:xff:fileReadError', ...
                'Error reading slice %d in file: %s.', sc, xfferror.message);
            return;
        end
        if ~isinf(numimages)
            numimages = numimages - 1;
        end
        if (ftell(fid) + sh + sls) > flen
            break;
        end
    end

    % remove unneeded slices
    imdata(:, :, sc:end) = [];

% time is first
else

    % ignore sub header size!
    imdata = reshape(fread(fid, [tdim, nor * noc], rclass), ...
        [tdim, nor, noc]);
end
