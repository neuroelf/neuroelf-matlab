function vmr = createvmr(imafiles, flags)
% createvmr  - create a VMR from a set of IMA/DCM files
%
% FORMAT:       vmr = createvmr(imafiles [, flags])
%
% Input fields:
%
%       imafiles    FxC char or Fx1/1xF cell array with IMA filenames
%                   multiple-slice files must be sorted correctly
%       flags       optional settings
%        .flipped   override assume flipped flag of xff
%
% Output fields:
%
%       vmr         VMR object

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:16 AM EST
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
xffroot = xff();
if nargin < 2 || ...
    numel(flags) ~= 1 || ...
   ~isstruct(flags)
    flags = struct;
end
if ~isfield(flags, 'flipped') || ...
    isempty(flags.flipped) || ...
   (~isnumeric(flags.flipped) && ...
    ~islogical(flags.flipped))
    flags.flipped = xffroot.Config('hdr', 'assumeflipped');
else
    flags.flipped = flags.flipped(1);
end
if ~isfield(flags, 'interp') || ...
   ~ischar(flags.interp) || ...
   ~any(strcmpi(flags.interp, {'cubic', 'lanczos3', 'linear', 'nearest'}))
    flags.interp = 'linear';
else
    flags.interp = lower(flags.interp(:)');
end
if ~isfield(flags, 'ismni') || ...
    isempty(flags.ismni) || ...
   (~isnumeric(flags.ismni) && ...
    ~islogical(flags.ismni))
    flags.ismni = false;
else
    flags.ismni = false || flags.ismni(1);
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

% prepare array and set un-used file ID
intal = false;
vdt = uint16(0);
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
if numel(imafiles) == 1

    % check extension
    exten = imafiles{1}(end-3:end);
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
            try
                oassfl = xffroot.Config('hdr', 'assumeflipped', flags.flipped);
                afileo = xff(afile);
                cll{1} = afileo;
                aframe = afileo.CoordinateFrame;
                vdt = afileo.VoxelData(:, :, :);

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
                xffroot.Config('hdr', 'assumeflipped', oassfl);

            catch ne_eo;
                clearxffobjects(cll);
                xffroot.Config('hdr', 'assumeflipped', oassfl);
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
                DimZ = size(p.MatrixValues, 1);

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
                    iop = [0, 1, 0; 0, 0, -1; 1, 0, 0; 1, 0, 0];
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
                vdt(DimX, DimY, DimZ) = vdt(1);

                % get data access struct
                switch lower(p.RECDataType)
                    case {'dyn.slice'}
                        slio = p.RECData.Dyn(1).Slice;
                    case {'slice'}
                        slio = p.RECData.Slice;
                    otherwise
                        error('UNSUPPORTED_PARREC');
                end

                % iterate over slices
                for sc = 1:DimZ
                    vdt(:, :, sc) = slio(sc).IO(:, :);
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
    switch lower(imafiles{1}(end-3:end))

        % DICOM / IMA
        case {'.dcm', '.ima'}
            try

                % assume number of files is number of slices
                DimZ = numel(imafiles);

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
                DimS = DimX * DimY;
                IInt = im1.PixelData;
                if (~isa(IInt, 'uint16') && ...
                    ~isa(IInt, 'uint8')) || ...
                    numel(IInt) ~= DimS
                    error('WRONG_DIM_OR_TYPE');
                end

                % try to get some crucial parameters from the images
                try
                    iop = reshape(im1.Value('ImageOrientationPatient'), [3, 2])';
                    ipp = im1.Value('ImagePositionPatient');
                    ipl = iml.Value('ImagePositionPatient');
                    pxs = im1.Value('PixelSpacing');
                    pxz = im1.Value('SliceThickness');
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    iop = [0, 1, 0; 0, 0, -1];
                    ipp = [-(DimZ - 1) / 2, -DimX / 2, DimY / 2];
                end
                iop(3, :) = -cross(iop(1, :), iop(2, :));
                iop = iop ./ repmat(sqrt(sum(iop .* iop, 2)), [1, 3]);

                % compute BrainVoyager QX settings
                sl1c = ipp + pxs(1) * (DimX / 2) * iop(1, :) + pxs(2) * (DimY / 2) * iop(2, :);
                slnc = ipl + pxs(1) * (DimX / 2) * iop(1, :) + pxs(2) * (DimY / 2) * iop(2, :);
                spbs = sqrt(sum((slnc - sl1c) .^ 2)) / (DimZ - 1);
                gaps = pxz - spbs;
                if abs(gaps) < 0.005
                    gaps = 0;
                end
                iop(4, :) = (slnc - sl1c) / (DimZ - 1);

                % initialize array size
                vdt(DimX, DimY, DimZ) = vdt(1);

                % iterate over slices
                for sc = 1:DimZ
                    fid = 0;
                    if dle
                        fid = fopen(imafiles{sc}, 'rb', 'ieee-le');
                    else
                        fid = fopen(imafiles{sc}, 'rb', 'ieee-be');
                    end
                    if fid < 3
                        error('FILE_OPEN_ERROR');
                    end
                    fseek(fid, -(12 + 2 * DimS), 1);
                    dtag = fread(fid, [1, 2], 'uint16=>double');
                    fread(fid, [1, 2], '*uint8');
                    dtls = fread(fid, [1, 1], 'uint16=>double');
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
                    vdt(:, :, sc) = reshape(fcnt, [DimX, DimY]);
                    fclose(fid);
                end

            catch ne_eo;
                clearxffobjects(cll);
                if fid > 0
                    fclose(fid);
                end
                error( ...
                    'neuroelf:xffError', ...
                    'Error creating VMR from DICOM files: %s.', ...
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

% create object
vmr = xff('new:vmr');

% set object properties
xffroot = xff();
updstate = xffroot.UpdateState('vmr', false);
vmr.DimX = DimX;
vmr.DimY = DimY;
vmr.DimZ = DimZ;
vmr.VMRData = [];
vmr.VMRData16 = [];
vmr.FramingCube = max([DimX, DimY, DimZ]);
fc = vmr.FramingCube / 2;
vmr.PosInfoVerified = 1;
vmr.OffsetX = max(0, ceil(fc - ceil(DimX / 2)));
vmr.OffsetY = max(0, ceil(fc - ceil(DimY / 2)));
vmr.OffsetZ = max(0, ceil(fc - ceil(DimZ / 2)));
vmr.Slice1CenterX = sl1c(1);
vmr.Slice1CenterY = sl1c(2);
vmr.Slice1CenterZ = sl1c(3);
vmr.SliceNCenterX = slnc(1);
vmr.SliceNCenterY = slnc(2);
vmr.SliceNCenterZ = slnc(3);
vmr.RowDirX = iop(1, 1);
vmr.RowDirY = iop(1, 2);
vmr.RowDirZ = iop(1, 3);
vmr.ColDirX = iop(2, 1);
vmr.ColDirY = iop(2, 2);
vmr.ColDirZ = iop(2, 3);
vmr.NRows = DimX;
vmr.NCols = DimY;
vmr.FoVRows = DimX * pxs(1);
vmr.FoVCols = DimY * pxs(2);
vmr.SliceThickness = pxz;
vmr.GapThickness = gaps;
vmr.Convention = cnvt;
vmr.VoxResX = pxs(1);
vmr.VoxResY = pxs(2);
vmr.VoxResZ = pxz + gaps;
vmr.VoxResInTalairach = intal;
vmr.VoxResVerified = true;

% set data
vdt = double(vdt);
mmm = minmaxmean(vdt);
if mmm(1) ~= 0
    vdt = vdt - mmm(1);
    mmm = mmm - mmm(1);
end
if mmm(2) > 32767
    vdt = floor((32767 / mmm(2)) .* vdt);
else
    vdt = round(vdt);
end
vmr.VMRData = uint8([]);
vmr.VMRData(size(vdt, 1), size(vdt, 2), size(vdt, 3)) = 0;
vmr.VMRData16 = uint16(vdt);

% try to perform good limitting
vmr.LimitVMR(struct('recalc8b', true));

% reset updatestate flag
xffroot.UpdateState('vmr', updstate);
