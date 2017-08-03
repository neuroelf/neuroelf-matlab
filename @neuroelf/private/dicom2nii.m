function nii = dicom2nii(imafiles, flags)
% dicom2nii  - create a NIftI object from a set of IMA/DCM files
%
% FORMAT:       nii = dicom2nii(imafiles [, flags])
%
% Input fields:
%
%       imafiles    FxC char or Fx1/1xF cell array with IMA filenames
%                   multiple-slice files must be sorted correctly
%       flags       optional settings
%        .disdaqs   if given, remove first N volumes of 4D data (0)
%        .dtype     datatype, either of 'int16' or {'single'}
%        .flip      additionally flip data (e.g. 'xy')
%        .mosaic    flag to force mosaic processing
%        .mosdimord mosaic dimension order (default: [1, 2])
%        .nslices   number of slices for 1-slice DCM files
%        .xyres     functional x/y resolution (default [64, 64])
%
% Output fields:
%
%       nii         NIftI object

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:18 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2016, Jochen Weber
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
if nargin < 1 || ...
   (~ischar(imafiles) && ...
    ~iscell(imafiles))
    error( ...
        'neuroelf:BadArguments', ...
        'Bad or missing argument.' ...
    );
end
if ischar(imafiles)
    if size(imafiles, 1) == 1 && ...
        any(imafiles == '*') && ...
        any(imafiles == '/' | imafiles == '\')
        [imapath, imapatt, imaext] = fileparts(imafiles);
        if isempty(imaext) || ...
            strcmp(imaext, '.')
            imaext = '.*';
        end
        imafiles = findfiles(imapath, [imapatt imaext], 'depth=1');
        % temp workaround
        rimafiles = findfiles(imapath, [imapatt imaext], 'depth=1');
        if ~isequal(imafiles, rimafiles)
            imafiles = union(imafiles(:), rimafiles(:));
        end
        imafiles = imafiles(:);
        disp(sprintf('%d files found for conversion.', numel(imafiles)));
    else
        imafiles = cellstr(imafiles);
    end
end
if numel(imafiles{1}) < 5
    error( ...
        'neuroelf:BadArguments', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
    numel(flags) ~= 1 || ...
   ~isstruct(flags)
    flags = struct;
end
if ~isfield(flags, 'disdaqs') || ...
   ~isa(flags.disdaqs, 'double') || ...
    numel(flags.disdaqs) ~= 1 || ...
    isinf(flags.disdaqs) || ...
    isnan(flags.disdaqs) || ...
    flags.disdaqs <= 0
    flags.disdaqs = 0;
else
    flags.disdaqs = floor(flags.disdaqs);
end
if ~isfield(flags, 'dtype') || ...
   ~ischar(flags.dtype) || ...
   ~any(strcmpi(flags.dtype(:)', {'int16', 'single'}))
    flags.dtype = 'single';
else
    flags.dtype = lower(flags.dtype(:)');
end
if ~isfield(flags, 'flip') || ...
    isempty(flags.flip) || ...
   ~ischar(flags.flip)
    flip = '';
else
    flip = lower(flags.flip(:)');
end
if ~isfield(flags, 'mosaic') || ...
   ~islogical(flags.mosaic) || ...
    numel(flags.mosaic) ~= 1
    flags.mosaic = false;
end
if ~isfield(flags, 'mosdimord') || ...
   ~isa(flags.mosdimord, 'double') || ...
    numel(flags.mosdimord) ~= 2 || ...
   (~all(flags.mosdimord(:) == [1; 2]) && ...
    ~all(flags.mosdimord(:) == [2; 1]))
    flags.mosdimord = [1, 2];
else
    flags.mosdimord = flags.mosdimord(:)';
end
if ~isfield(flags, 'nslices') || ...
   ~isa(flags.nslices, 'double') || ...
    numel(flags.nslices) ~= 1 || ...
    isinf(flags.nslices) || ...
    isnan(flags.nslices) || ...
    flags.nslices < 2
    flags.nslices = [];
end
if ~isfield(flags, 'xyres') || ...
   ~isa(flags.xyres, 'double') || ...
    numel(flags.xyres) ~= 2 || ...
    any(isinf(flags.xyres) | isnan(flags.xyres) | flags.xyres < 16 | flags.xyres > 512)
    flags.xyres = [64, 64];
else
    flags.xyres = round(flags.xyres(:))';
end

% make sure files exist
for sc = 1:numel(imafiles)
    if ~ischar(imafiles{sc}) || ...
        exist(imafiles{sc}(:)', 'file') ~= 2
        error( ...
            'neuroelf:FileNotFound', ...
            'File for slice %d not found.', ...
            sc ...
        );
    end
    imafiles{sc} = imafiles{sc}(:)';
end

% datatype
if flags.dtype ~= 's'
    vdt = int16(0);
else
    vdt = single(0);
end

% the next big scope decides on the type of file, reads in the
% file(s) and determines a few settings.
% after the scope, the following variables must be defined and
% set correctly to procede:
%
% DimX, DimY, DimZ   : dimensions of dataset
% gaps :             : gap size
% iop : double(4, 3) : [rowdir, coldir, slcdir, slcdircomp]
% pxs : double(1, 2) : pixel spacing within slice
% pxz :              : pixel spacing between slices
% sl1c / slnc        : slice 1 and N center coordinates (in world space)
% vdt                : data array with read voxel data

% what kind of files
fid = 0;
cll = cell(1, 3);
exten = regexprep(imafiles{1}, '^.*(\.[^\.]+)$', '$1');
if numel(exten) > 1 && ...
    all(exten(2:end) >= '0' & exten(2:end) <= '9')
    exten = '.dcm';
end
if numel(imafiles) == 1

    % check extension
    switch (lower(exten))

        % PAR/REC
        case {'.par', '.rec'}
            try

                % read file
                p = readpar(imafiles{1});

                % get resolution
                DimX = p.SliceResolution(1);
                DimY = p.SliceResolution(2);
                DimZ = numel(p.RECData.Dyn(1).Slice);
                DimT = numel(p.RECData.Dyn);

                % get some settings
                try
                    angi = find(strcmp(p.MatrixHeaders, 'image_angulation_1'));
                    rotV = p.MatrixValues(1, [angi + 2, angi, angi + 1]);
                    rotX = [ ...
                        1, 0, 0; ...
                        0,  cos(rotV(1) * pi / 180), -sin(rotV(1) * pi / 180); ...
                        0,  sin(rotV(1) * pi / 180),  cos(rotV(1) * pi / 180)];
                    rotY = [ ...
                         cos(rotV(2) * pi / 180), 0,  sin(rotV(2) * pi / 180); ...
                        0, 1, 0;  ...
                        -sin(rotV(2) * pi / 180), 0,  cos(rotV(2) * pi / 180)]; ...
                    rotZ = [ ...
                         cos(rotV(3) * pi / 180), -sin(rotV(3) * pi / 180), 0; ...
                         sin(rotV(3) * pi / 180),  cos(rotV(3) * pi / 180), 0; ...
                         0, 0, 1]; ...
                    iop = (rotX * rotY * rotZ)';
                    iop(4, :) = -cross(iop(1, :), iop(2, :));
                    iop = iop ./ repmat(sqrt(sum(iop .* iop, 2)), [1, 3]);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    iop = [0, 1, 0; 0, 0, -1; 1 0 0; 1 0 0];
                end

                % and some more settings
                resi = find(strcmp(p.MatrixHeaders, 'pixel_spacing_1'));
                if isempty(resi)
                    error('UNSUPPORTED_PARREC');
                end
                pxs  = p.MatrixValues(1, [resi, resi + 1]);
                slci = find(strcmp(p.MatrixHeaders, 'image_offcentre_1'));
                if isempty(slci)
                    error('UNSUPPORTED_PARREC');
                end
                sl1c = p.MatrixValues( 1 , [slci + 2, slci, slci + 1]);
                slnc = p.MatrixValues(end, [slci + 2, slci, slci + 1]);
                spbs = sqrt(sum((slnc - sl1c) .^ 2)) / (DimZ - 1);

                % initialize array size
                vdt(DimX, DimY, DimZ, DimT) = vdt(1);

                % get data access struct
                for tc = 1:DimT
                    slio = p.RECData.Dyn(tc).Slice;
                    for sc = 1:DimZ
                        vdt(:, :, sc, tc) = slio(sc).IO(:, :);
                    end
                end
            catch ne_eo;
                error( ...
                    'neuroelf:ErrorReadingFile', ...
                    'Error reading specified PAR/REC file (%s).', ...
                    ne_eo.message ...
                );
            end

        % unsupported
        otherwise
            error( ...
                'neuroelf:BadArgument', ...
                'Single file of this type not supported.' ...
            );
    end
else

    % also check extension
    switch lower(exten)

        % DICOM / IMA
        case {'.dcm', '.ima'}
            try

                % assume number of files is number of time points
                if isempty(flags.nslices) || ...
                    flags.mosaic
                    DimT = numel(imafiles);
                    DimZ = 1;
                else
                    DimT = floor(numel(imafiles) / flags.nslices);
                    DimZ = flags.nslices;
                end

                % read in first and last file
                im1 = xff(imafiles{1});
                cll{1} = im1;
                iml = xff(imafiles{end});
                cll{2} = iml;

                % check files
                if ~isxff(im1, 'dcm') || ...
                   ~isxff(iml, 'dcm') || ...
                    im1.DataLittleEndian ~= iml.DataLittleEndian
                    error('BAD_FILES');
                end

                % check some crucial settings
                dle = im1.DataLittleEndian;
                tse = im1.DataTSExplicit;
                DimX = im1.Value('Columns');
                DimY = im1.Value('Rows');
                if ischar(DimX) || ...
                    ischar(DimY)
                    error('BAD_VALUES');
                end

                % check pixel data of first image
                DimS = DimX * DimY;
                IInt = im1.PixelData;
                if (~isa(IInt, 'uint16') && ...
                    ~isa(IInt, 'uint8')) || ...
                   (numel(IInt) ~= DimS && ...
                    numel(IInt) ~= (2 * DimS))
                    error('WRONG_DIM_OR_TYPE');
                end

                % patch for mosaic
                mosxy = [1, 1];
                if (DimX > flags.xyres(1) && ...
                    DimY > flags.xyres(2)) || ...
                    flags.mosaic
                    flags.mosaic = true;
                    mosxy = [DimX, DimY] ./ flags.xyres;
                    if all(mod(mosxy, 1) == 0) && ...
                        mosxy(1) == mosxy(2) && ...
                        DimZ <= prod(mosxy)
                        if DimZ > 1
                            DimT = DimT * DimZ;
                        else
                            testi = unpackmosaic(reshape(IInt, [DimX, DimY]), ...
                                flags.mosdimord, 3, mosxy(1));
                            if ~isempty(flags.nslices)
                                DimZ = flags.nslices;
                            else
                                DimZ = size(testi, 3);
                            end
                        end
                        DimX = DimX / mosxy(1);
                        DimY = DimY / mosxy(1);
                    end
                end

                % try to get some crucial parameters from the images
                try
                    iop = reshape(im1.Value('ImageOrientationPatient'), [3, 2])';
                    ipp = im1.Value('ImagePositionPatient');
                    pxs = im1.Value('PixelSpacing');
                    try
                        spbs = im1.Value('SpacingBetweenSlices');
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                        spbs = im1.Value('SliceThickness');
                    end
                catch ne_eo;
                    error( ...
                        'NeuroElf:DICOMError', ...
                        'Missing fields in DICOM file (%s).', ...
                        ne_eo.message ...
                    );
                end

                % compute orthogonal 3rd dimension vector
                iop(3, :) = cross(iop(1, :), iop(2, :));

                % and re-orthogonalize second vector to full precision
                iop(2, :) = -cross(iop(1,:), iop(3, :));

                % scale each direction to 1
                iop = iop ./ repmat(sqrt(sum(iop .* iop, 2)), [1, 3]);

                % compute BrainVoyager QX settings
                sl1c = ipp + pxs(1) * (DimX / 2) * iop(1, :) + pxs(2) * (DimY / 2) * iop(2, :);

                % initialize array size
                vdt(DimX, DimY, DimZ, DimT) = vdt(1);

                % iterate over time points, slices
                fc = 1;
                for tc = 1:DimT
                    if all(mosxy == 1)
                        for sc = 1:DimZ
                            fid = 0;
                            if dle
                                fid = fopen(imafiles{fc}, 'rb', 'ieee-le');
                            else
                                fid = fopen(imafiles{fc}, 'rb', 'ieee-be');
                            end
                            if fid < 3
                                error('FILE_OPEN_ERROR');
                            end
                            fc = fc + 1;
                            if tse
                                fseek(fid, -(12 + 2 * DimS), 1);
                                dtag = fread(fid, [1, 2], 'uint16=>double');
                                fread(fid, [1, 2], '*uint8');
                                dtls = fread(fid, [1, 1], 'uint16=>double');
                                dtll = fread(fid, [1, 1], 'uint32=>double');
                            else
                                fseek(fid, -(8 + 2 * DimS), 1);
                                dtag = fread(fid, [1, 2], 'uint16=>double');
                                dtls = 0;
                                dtll = fread(fid, [1, 1], 'uint32=>double');
                            end
                            if ~all(dtag == [32736, 16]) || ...
                                dtls ~= 0 || ...
                                dtll ~= (2 * DimS)
                                error('BAD_TAG_OR_LENGTH');
                            end
                            fcnt = fread(fid, [1, Inf], '*uint16');
                            if numel(fcnt) ~= DimS
                                error('FILE_TOO_SHORT');
                            end
                            vdt(:, :, sc, tc) = reshape(fcnt, [DimX, DimY]);
                            fclose(fid);
                        end
                    else
                        fid = 0;
                        if dle
                            fid = fopen(imafiles{fc}, 'rb', 'ieee-le');
                        else
                            fid = fopen(imafiles{fc}, 'rb', 'ieee-be');
                        end
                        if fid < 3
                            error('FILE_OPEN_ERROR');
                        end
                        fc = fc + 1;
                        if tse
                            fseek(fid, -(12 + 2 * DimS), 1);
                            dtag = fread(fid, [1, 2], 'uint16=>double');
                            fread(fid, [1, 2], '*uint8');
                            dtls = fread(fid, [1, 1], 'uint16=>double');
                            dtll = fread(fid, [1, 1], 'uint32=>double');
                        else
                            fseek(fid, -(8 + 2 * DimS), 1);
                            dtag = fread(fid, [1, 2], 'uint16=>double');
                            dtls = 0;
                            dtll = fread(fid, [1, 1], 'uint32=>double');
                        end
                        if ~all(dtag == [32736, 16]) || ...
                            dtls ~= 0 || ...
                            dtll ~= (2 * DimS)
                            error('BAD_TAG_OR_LENGTH');
                        end
                        fcnt = fread(fid, [1, Inf], '*uint16');
                        if numel(fcnt) ~= DimS
                            error('FILE_TOO_SHORT');
                        end
                        testi = unpackmosaic(reshape(fcnt, mosxy(1) * [DimX, DimY]), ...
                            flags.mosdimord, 3, mosxy(1));
                        vdt(:, :, 1:min(DimZ, size(testi, 3)), tc) = ...
                            testi(:, :, 1:min(DimZ, size(testi, 3)));
                        fclose(fid);
                    end
                end

            catch ne_eo;
                clearxffobjects(cll);
                if fid > 0
                    fclose(fid);
                end
                error( ...
                    'neuroelf:xffError', ...
                    'Error creating NII from DICOM files (%s).', ...
                    ne_eo.message ...
                );
            end

        % unsupported
        otherwise
            error( ...
                'neuroelf:BadArgument', ...
                'Single file of this type not supported.' ...
            );
    end
end

% remove disdaqs
if flags.disdaqs > 0 && ...
    size(vdt, 4) > 1
    vdt(:, :, :, 1:min(size(vdt, 4) - 1, flags.disdaqs)) = [];
end

% squeeze for one-slice data
if size(vdt, 3) == 1
    vdt = squeeze(vdt);
end

% clear possible objects
clearxffobjects(cll);

% reject invalid centers
sl1c(abs(sl1c) > 128) = 0;

% compute transformation matrix offset and matrix
off = sl1c - 0.5 * (pxs(1) * (DimX + 1) * iop(1, :) + pxs(2) * (DimY + 1) * iop(2, :));
pxs(3) = spbs;
mat = [repmat(pxs(:)', [3, 1]) .* (iop(1:3, 1:3)'), off(:)];

% create object
nii = xff('new:hdr');

% set object properties
nii.FileMagic = 'n+1';
nii.NIIFileType = 2;
nii.ImgDim.Dim(2:1+ndims(vdt)) = size(vdt);
if flags.dtype(1) ~= 's'
    nii.ImgDim.DataType = 4;
    nii.ImgDim.BitsPerPixel = 16;
else
    nii.ImgDim.DataType = 16;
    nii.ImgDim.BitsPerPixel = 32;
end
nii.ImgDim.PixSpacing(2:8) = [pxs(:)', 2, 1, 1, 1];
nii.DataHist.Description = 'NeuroElf-imported NIftI image';
nii.DataHist.NIftI1.QFormCode = 2;
nii.DataHist.NIftI1.SFormCode = 2;
nii.DataHist.NIftI1.AffineTransX = mat(1, :);
nii.DataHist.NIftI1.AffineTransY = mat(2, :);
nii.DataHist.NIftI1.AffineTransZ = mat(3, :);

% extend RunTimeVars.Map
nii.RunTimeVars.Map = repmat(nii.RunTimeVars.Map(1), [1, size(vdt, 4)]);

% data flips
if any(flip == 'x')
    vdt = vdt(end:-1:1, :, :, :);
end
if any(flip == 'y')
    vdt = vdt(:, end:-1:1, :, :);
end
if any(flip == 'z')
    vdt = vdt(:, :, end:-1:1, :);
end
nii.VoxelData = vdt;
