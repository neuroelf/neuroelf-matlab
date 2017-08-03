function [fmr, fvol] = createfmr(imafiles, flags)
% createfmr  - create a FMR from a set of IMA/DCM files
%
% FORMAT:       [fmr, fvol] = createfmr(imafiles [, flags])
%
% Input fields:
%
%       imafiles    FxC char or Fx1/1xF cell array with IMA filenames
%                   multiple-slice files must be sorted correctly
%       flags       optional settings
%        .fileorder slice-wise DCM file order, either 'tz' or {'zt'}
%        .flip      additionally flip data (e.g. 'xy')
%        .mosaic    flag to force mosaic processing
%        .mosdimord mosaic dimension order (default: [1, 2])
%        .nslices   number of slices for 1-slice DCM files
%        .skipvols  number of skipped volumes (default: 0)
%        .xyres     functional x/y resolution (default [64, 64])
%
% Output fields:
%
%       fmr         FMR object

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:15 AM EST
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

% xff factory
global xffsngl;

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
    imafiles = cellstr(imafiles);
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
if ~isfield(flags, 'fileorder') || ...
   ~ischar(flags.fileorder) || ...
    isempty(flags.fileorder) || ...
    lower(flags.fileorder(1)) ~= 't'
    forderz = true;
else
    forderz = false;
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
if ~isfield(flags, 'skipvols') || ...
   ~isa(flags.skipvols, 'double') || ...
    numel(flags.skipvols) ~= 1 || ...
    isinf(flags.skipvols) || ...
    isnan(flags.skipvols) || ...
    flags.skipvols < 0
    flags.skipvols = 0;
else
    flags.skipvols = round(flags.skipvols);
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

% assume TR/TE of 2000/50ms
atr = 2000;
ate = 50;

% prepare array and set un-used file ID
% intal = false; (not used yet)
xff;
dtype = xffsngl.CONF.settings.DataTypes.STC;
if dtype == 1
    vdt = uint16(0);
else
    vdt = single(0);
end
fid = 0;

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
cll = cell(1, 3);
exten = imafiles{1}(end-3:end);
if numel(imafiles) == 1

    % check extension
    switch (lower(exten))

        % Analyze
        case {'.hdr', '.img', '.nii'}

            % for old analyze check HDR whether HDR exists
            gaps = 0;

            % read file
            afile = imafiles{1};
            if any(exten == upper(exten))
                hfile = [afile(1:end-4) '.HDR'];
            else
                hfile = [afile(1:end-4) '.hdr'];
            end
            if strcmpi(exten, '.img')
                afile = hfile;
            end
            if exist(afile, 'file') ~= 2
                error( ...
                    'neuroelf:FileNotFound', ...
                    'Required Analyze file not found: %s.', ...
                    afile ...
                );
            end

            % get coordinate frame
            try
                afileo = xff(afile);
                cll{1} = afileo;
                aframe = afileo.CoordinateFrame;
                vdt = resolve(afileo.VoxelData);

                % get dims, resolution, and settings
                DimX = aframe.DimX;
                DimY = aframe.DimY;
                DimZ = aframe.DimZ;
                pxs = [aframe.ResX, aframe.ResY];
                pxz = aframe.ResZ;
                iop = [aframe.RowDir; aframe.ColDir; aframe.SlcDir];
                iop = iop ./ repmat(sqrt(sum(iop .* iop, 2)), [1, 3]);
                iop(4, :) = -cross(iop(1, :), iop(2, :));
                sl1c = aframe.Slice1Center;
                slnc = aframe.SliceNCenter;

            catch ne_eo;
                clearxffobjects(cll);
                error( ...
                    'neuroelf:ErrorReadingFile', ...
                    'Error reading specified Analyze file: %s.', ...
                    ne_eo.message ...
                );
            end

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
                slti = find(strcmp(p.MatrixHeaders, 'slice_thickness'));
                if isempty(slti)
                    error('UNSUPPORTED_PARREC');
                end
                pxz  = p.MatrixValues(1, slti(1));
                slgi = find(strcmp(p.MatrixHeaders, 'slice_gap'));
                if isempty(slgi)
                    error('UNSUPPORTED_PARREC');
                end
                gaps = p.MatrixValues(1, slgi(1));

                % possibly recalc gap?
                if abs(spbs - (pxz + gaps)) > 5e-4
                    gaps = spbs - pxz;
                end

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
                    'Error reading specified PAR/REC file: %s.', ...
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

        % Analyze
        case {'.hdr', '.img', '.nii'}

            % for old analyze check HDR whether HDR exists
            gaps = 0;

            % read file
            afile = imafiles{1};
            if any(exten == upper(exten))
                hfile = [afile(1:end-4) '.HDR'];
            else
                hfile = [afile(1:end-4) '.hdr'];
            end
            if strcmpi(exten, '.img')
                afile = hfile;
            end
            if exist(afile, 'file') ~= 2
                error( ...
                    'neuroelf:FileNotFound', ...
                    'Required Analyze file not found: %s.', ...
                    afile ...
                );
            end

            % get coordinate frame
            try
                afileo = xff(afile);
                cll{1} = afileo;
                aframe = afileo.CoordinateFrame;
                vdt = resolve(afileo.VoxelData);
                vds = size(vdt);

                % get dims, resolution, and settings
                DimX = aframe.DimX;
                DimY = aframe.DimY;
                DimZ = aframe.DimZ;
                pxs = [aframe.ResX, aframe.ResY];
                pxz = aframe.ResZ;
                iop = [aframe.RowDir; aframe.ColDir; aframe.SlcDir];
                iop = iop ./ repmat(sqrt(sum(iop .* iop, 2)), [1, 3]);
                iop(4, :) = -cross(iop(1, :), iop(2, :));
                sl1c = aframe.Slice1Center;
                slnc = aframe.SliceNCenter;

            catch ne_eo;
                clearxffobjects(cll);
                error( ...
                    'neuroelf:ErrorReadingFile', ...
                    'Error reading specified Analyze file: %s.', ...
                    ne_eo.message ...
                );
            end

            % load further analyze files into vdt
            vdt(end, end, end, numel(imafiles)) = 0;
            for ic = 2:numel(imafiles)
                try
                    afile = imafiles{ic};
                    exten = afile(end-3:end);
                    if any(exten == upper(exten))
                        hfile = [afile(1:end-4) '.HDR'];
                    else
                        hfile = [afile(1:end-4) '.hdr'];
                    end
                    if strcmpi(exten, '.img')
                        afile = hfile;
                    end
                    if exist(afile, 'file') ~= 2
                        error( ...
                            'neuroelf:FileNotFound', ...
                            'Required Analyze file not found: %s.', ...
                            afile ...
                        );
                    end
                    cll{2} = xff(afile);
                    vdr = resolve(cll{2}.VoxelData);
                    if ~isequal(size(vdr), vds)
                        error( ...
                            'neuroelf:SizeMismatch', ...
                            'Size of Analyze series image mismatch.' ...
                        );
                    end
                    vdt(:, :, :, ic) = vdr;
                    cll{2}.ClearObject;
                    cll{2} = [];
                catch ne_eo;
                    clearxffobjects(cll);
                    rethrow(ne_eo);
                end
            end


        % DICOM / IMA
        case {'.dcm', '.ima'}
            try

                % assume number of files is number of time points
                if isempty(flags.nslices) || ...
                    flags.mosaic
                    DimT = numel(imafiles);
                    DimZ = 1;
                else
                    DimT = numel(imafiles) / flags.nslices;
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
                DimX = im1.Value('Columns');
                DimY = im1.Value('Rows');
                if ischar(DimX) || ...
                    ischar(DimY)
                    error('BAD_VALUES');
                end

                % check pixel data of first image
                IInt = im1.PixelData;
                if flags.mosaic && ...
                    DimX == 64 && ...
                    DimY == 64
                    DimS = numel(IInt);
                    nummos = ceil(sqrt(flags.nslices - sqrt(eps)));
                    DimXt = sqrt(DimS / (nummos * nummos));
                    if DimXt ~= 64
                        error('WRONG_DIM_OR_TYPE');
                    end
                    DimX = sqrt(DimS);
                    DimY = DimX;
                else
                    DimS = DimX * DimY;
                end
                if (~isa(IInt, 'uint16') && ...
                    ~isa(IInt, 'uint8')) || ...
                    numel(IInt) ~= DimS
                    error('WRONG_DIM_OR_TYPE');
                end

                % also override TR/TE if available
                try
                    atr = im1.Value('RepetitionTime');
                    ate = im1.Value('EchoTime');
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end

                % patch for mosaic
                oDimX = DimX;
                oDimY = DimY;
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
                            DimZ = size(testi, 3);
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
                    pxz = im1.Value('SliceThickness');
                    if ~flags.mosaic
                        ipl = iml.Value('ImagePositionPatient');
                    else
                        iop(3, :) = cross(iop(1, :), iop(2, :));
                        iop(3, :) = iop(3, :) ./ sqrt(sum(iop(3, :) .^ 2));
                        try
                            spbs = im1.Value('SpacingBetweenSlices');
                            iop(3, :) = spbs .* iop(3, :);
                            ipl = ipp + (DimZ - 1) .* iop(3, :);
                        catch ne_eo;
                            ipl = ipp + (DimZ - 1) .* (pxz .* iop(3, :));
                            rethrow(ne_eo);
                        end
                    end
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    iop = [0, 1, 0; 0, 0, -1];
                    ipp = [-(DimZ - 1) / 2, -DimX / 2, DimY / 2];
                end
                iop(3, :) = -cross(iop(1, :), iop(2, :));
                iop = iop ./ repmat(sqrt(sum(iop .* iop, 2)), [1, 3]);

                % compute BrainVoyager QX settings
                sl1c = ipp + pxs(1) * (oDimX / 2) * iop(1, :) + pxs(2) * (oDimY / 2) * iop(2, :);
                slnc = ipl + pxs(1) * (oDimX / 2) * iop(1, :) + pxs(2) * (oDimY / 2) * iop(2, :);
                spbs = sqrt(sum((slnc - sl1c) .^ 2)) / (DimZ - 1);
                gaps = spbs - pxz;
                if abs(gaps) < 0.005
                    gaps = 0;
                end
                iop(4, :) = (slnc - sl1c) / (DimZ - 1);

                % initialize array size
                vdt(DimX, DimY, DimZ, DimT) = vdt(1);

                % initialize fast read offset
                dtse = im1.DataTSExplicit;
                if dtse
                    dcmdatao = 12;
                else
                    dtls = 0;
                    dcmdatao = 8;
                end

                % iterate over time points, slices
                for tc = 1:DimT
                    if all(mosxy == 1)
                        for sc = 1:DimZ
                            fid = 0;
                            if forderz
                                fc = (tc - 1) * DimZ + sc;
                            else
                                fc = (sc - 1) * DimT + tc;
                            end
                            if dle
                                fid = fopen(imafiles{fc}, 'rb', 'ieee-le');
                            else
                                fid = fopen(imafiles{fc}, 'rb', 'ieee-be');
                            end
                            if fid < 3
                                error('FILE_OPEN_ERROR');
                            end
                            fseek(fid, -(dcmdatao + 2 * DimS), 1);
                            dtag = fread(fid, [1, 2], 'uint16=>double');
                            if dtse
                                fread(fid, [1, 2], '*uint8');
                                dtls = fread(fid, [1, 1], 'uint16=>double');
                            end
                            dtll = fread(fid, [1, 1], 'uint32=>double');
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
                            fid = fopen(imafiles{tc}, 'rb', 'ieee-le');
                        else
                            fid = fopen(imafiles{tc}, 'rb', 'ieee-be');
                        end
                        if fid < 3
                            error('FILE_OPEN_ERROR');
                        end
                        fseek(fid, -(dcmdatao + 2 * DimS), 1);
                        dtag = fread(fid, [1, 2], 'uint16=>double');
                        if dtse
                            fread(fid, [1, 2], '*uint8');
                            dtls = fread(fid, [1, 1], 'uint16=>double');
                        end
                        dtll = fread(fid, [1, 1], 'uint32=>double');
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
                    'Error creating FMR from DICOM files: %s.', ...
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

% clear possible objects
clearxffobjects(cll);

% check orientation
iopr = corrcoef(iop(3, :), iop(4, :));
iopr = iopr(1, 2);
if abs(iopr) < 0.75
    error( ...
        'neuroelf:InternalError', ...
        'Error getting image orientation.' ...
    );
end
if iopr < 0
    cnvt = 0;
else
    cnvt = 1;
end

% first-vol FMR?
if nargout > 1
    fvold = vdt(:, :, :, 1);
end

% skip volumes
if flags.skipvols > 0
    novols = size(vdt, 4);
    vdt(:, :, :, 1:min(novols, flags.skipvols)) = [];
    noskipped = novols - size(vdt, 4);
else
    noskipped = 0;
end

% get things into a square matrix with a power-2
psdim = [64, 96, 128, 192, 256];
mxdim = max(DimX, DimY);
mndim = find(mxdim <= psdim);
if isempty(mndim)
    error( ...
        'neuroelf:DataTooLarge', ...
        'Data too large to fit into square matrix.' ...
    );
end
acdim = psdim(mndim(1));

% create object
fmr = xff('new:fmr');

% set object properties
xffroot = xff();
updstate = xffroot.UpdateState('fmr', false);
fmr.NrOfVolumes = size(vdt, 4);
fmr.NrOfSlices = DimZ;
fmr.NrOfSkippedVolumes = noskipped;
fmr.TR = atr;
fmr.InterSliceTime = floor(0.25 + atr / DimZ);
fmr.TE = ate;
fmr.ResolutionX = acdim;
fmr.ResolutionY = acdim;
fmr.LayoutNColumns = 1 + ceil(sqrt(fmr.NrOfSlices));
fmr.LayoutNRows = ceil(fmr.NrOfSlices / fmr.LayoutNColumns);
fmr.InplaneResolutionX = pxs(1);
fmr.InplaneResolutionY = pxs(2);
fmr.SliceThickness = pxz;
fmr.SliceGap = gaps;
fmr.VoxelResolutionVerified = 1;
fmr.PosInfosVerified = 1;
% fmr.CoordinateSystem = 1;
fmr.Slice1CenterX = sl1c(1);
fmr.Slice1CenterY = sl1c(2);
fmr.Slice1CenterZ = sl1c(3);
fmr.SliceNCenterX = slnc(1);
fmr.SliceNCenterY = slnc(2);
fmr.SliceNCenterZ = slnc(3);
fmr.RowDirX = iop(1, 1);
fmr.RowDirY = iop(1, 2);
fmr.RowDirZ = iop(1, 3);
fmr.ColDirX = iop(2, 1);
fmr.ColDirY = iop(2, 2);
fmr.ColDirZ = iop(2, 3);
fmr.NRows = acdim;
fmr.NCols = acdim;
fmr.FoVRows = acdim * pxs(1);
fmr.FoVCols = acdim * pxs(2);
fmr.GapThickness = gaps;
if cnvt
    fmr.Convention = 'Radiological';
else
    fmr.Convention = 'Neurological';
end
if dtype == 1
    stcd = uint16(0);
else
    stcd = single(0);
end
stcd1 = stcd;
stcd(acdim, acdim, fmr.NrOfVolumes, fmr.NrOfSlices) = 0;
stcd1(acdim, acdim, 1, fmr.NrOfSlices) = 0;
if acdim > DimX
    dxf = ceil((acdim - DimX) / 2);
    dxt = dxf + DimX - 1;
else
    dxf = 1;
    dxt = acdim;
end
if acdim > DimY
    dyf = ceil((acdim - DimY) / 2);
    dyt = dyf + DimY - 1;
else
    dyf = 1;
    dyt = acdim;
end
if any(flip == 'x')
    dxr = dxt:-1:dxf;
else
    dxr = dxf:dxt;
end
if any(flip == 'y')
    dyr = dyt:-1:dyf;
else
    dyr = dyf:dyt;
end
if any(flip == 'z')
    dzr = DimZ:-1:1;
else
    dzr = 1:DimZ;
end

% set data
if dtype < 2
    mmm = minmaxmean(vdt);
    if mmm(1) < 0
        vdt = single(vdt);
        vdt = vdt - mmm(1);
        if nargout > 1
            fvold = single(fvold) - mmm(1);
        end
        mmm(2) = mmm(2) - mmm(1);
    end
    if mmm(2) > 32767
        vdt = single(vdt);
        vdt = (32767 / mmm(2)) .* vdt;
        if nargout > 1
            fvold = (32767 / mmm(2)) .* single(fvold);
        end
    end
    fmr.FileVersion = 5;
end
stcd(dxr, dyr, :, dzr) = permute(vdt, [1, 2, 4, 3]);
fmr.DataType = dtype;
if nargout > 1
    stcd1(dxr, dyr, :, dzr) = permute(fvold, [1, 2, 4, 3]);
    fvol = fmr.CopyObject;
    fvol.NrOfVolumes = 1;
    fvol.Slice.STCData = stcd1;
end
fmr.Slice.STCData = stcd;

% try to perform good limitting
xffroot.UpdateState('fmr', updstate);
