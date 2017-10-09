function xo = aft_ExportNifti(xo, niifile, talorder, range)
% AnyFileType::ExportNifti  - save object as NII object
%
% FORMAT:       obj.ExportNifti(niifile [, talorder [, range]]);
%
% Input fields:
%
%       niffile     output filename (must end in .nii)
%       talorder    if given and 1x1 boolean true, permute BV data dims
%       range       range of volumes (default: [1:NrOfVolumes])
%
% No output fields.
%
% TYPES: AVA, CMP, DDT, DMR, FMR, GLM, HEAD, MGH, MSK, VDW, VMP, VMR, VTC
%
% Note: this method does *not* resample data, and map/volume order
%       is as in the source object (see file format documentation)
%
% Using: multimatch, spmm2q.

% Version:  v1.1
% Build:    16031615
% Date:     Mar-16 2016, 3:34 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011 - 2016, Jochen Weber
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

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'ava', 'cmp', 'ddt', 'dmr', ...
    'fmr', 'glm', 'head', 'mgh', 'msk', 'vdw', 'vmp', 'vmr', 'vtc'}) || ~ischar(niifile) || ...
    isempty(niifile) || isempty(regexpi(niifile(:)', '\.(hdr|img|nii)$'))
    error('neuroelf:xff:badArguments', 'Invalid call to %s.', mfilename );
end
niifile = niifile(:)';
ftype = lower(xo.S.Extensions{1});
if nargin < 3 || ~islogical(talorder) || numel(talorder) ~= 1 || ...
    any(strcmpi(ftype, {'dmr', 'fmr', 'head'}))
    talorder = false;
end

% what NIftI type
if strcmpi(niifile(end-2:end), 'nii')
    niitype = 2;
    niimagic = 'n+1';
    niioff = 352;
else
    niitype = 1;
    niimagic = 'ni1';
    niifile = regexprep(niifile, '\.img$', '.hdr', 'preservecase');
    niioff = 0;
end

% get objects super-struct, content, and filetype
bc = xo.C;
fname = xo.F;
if isempty(fname)
    fname = sprintf('%s object', upper(ftype));
elseif numel(fname) > 70
    fname = [fname(1:34) '...' fname(end-32:end)];
end
fname = ['exported ' fname];

% get bounding box
if ~any(strcmp(ftype, {'dmr', 'fmr'}))
    bbox = aft_BoundingBox(xo);
end

% for functional files
if any(strcmp(ftype, {'dmr', 'fmr', 'vtc'}))

    numvol = bc.NrOfVolumes;

% for mask/structural file
elseif any(strcmp(ftype, {'mgh', 'msk', 'vmr'}))

    % only one volume
    numvol = 1;

% for DDT
elseif strcmp(ftype, 'ddt')

    % only one volume
    numvol = 12;

% for stats files
else

    % get map names
    mnames = aft_MapNames(xo);
    numvol = numel(mnames);
end
if nargin < 4 || ~isa(range, 'double') || isempty(range) || ...
    any(isinf(range(:)) | isnan(range(:)) | range(:) < 1 | range(:) > numvol)
    range = 1:numvol;
else
    range = sort(unique(round(range(:))))';
end
numvol = numel(range);

% some more settings
tr = 1;
switch (ftype)
    case {'ava', 'cmp', 'ddt', 'glm', 'vdw', 'vmp'}
        dclass = 'single';
        dtype = 16;
        bpv = 32;
        res = bbox.ResXYZ;
    case {'dmr', 'fmr', 'vtc'}
        if bc.DataType == 1
            dclass = 'int16';
            dtype = 4;
            bpv = 16;
        else
            dclass = 'single';
            dtype = 16;
            bpv = 32;
        end
        if strcmp(ftype, 'vtc')
            res = bbox.ResXYZ;
        else
            res = [bc.InplaneResolutionX, bc.InplaneResolutionY, ...
                bc.SliceThickness + bc.SliceGap];
        end
        tr = 0.001 * bc.TR;
    case {'head', 'mgh'}
        dclass = 'single';
        dtype = 16;
        bpv = 32;
        if ftype(1) == 'h'
            cfr = head_CoordinateFrame(xo);
        else
            cfr = mgh_CoordinateFrame(xo);
        end
        res = cfr.Resolution;
    case 'msk'
        dclass = 'uint8';
        dtype = 2;
        bpv = 8;
        res = bbox.ResXYZ;
    case 'vmr'
        if isempty(bc.VMRData16) && isa(bc.VMRData, 'uint8')
            dclass = 'uint8';
            dtype = 2;
            bpv = 8;
        else
            dclass = 'int16';
            dtype = 4;
            bpv = 16;
        end
        res = bbox.ResXYZ;
end

% transformation matrix for BV voxel-space objects
if any(strcmp(ftype, {'ava', 'cmp', 'ddt', 'glm', 'msk', 'vdw', 'vmp', 'vmr', 'vtc'}))

    % from bounding box
    trf = bbox.QuatB2T;

% for HEAD objects
elseif strcmp(ftype, 'head') || strcmp(ftype, 'mgh')

    % get from coordinate frame
    trf = cfr.Trf;

% for DMR/FMR objects
else

    % get from CoordinateFrame (also works for DMRs!)
    cfr = fmr_CoordinateFrame(xo);
    trf = cfr.Trf;
end

% get volume size
vsz = aft_GetVolumeSize(xo);

% correct stuff
if talorder
    newoffset = trf * [vsz(:) + 1; 1];
	trf = [0, -1, 0, 0; 0, 0, -1, 0; -1, 0, 0, 0; 0, 0, 0, 1] * trf;
    trf(1:3, 4) = newoffset(1:3);
    res = res(1, [3, 1, 2]);
    vsz = vsz(1, [3, 1, 2]);
end

% compute offset
offset = trf(1:3, 4) + trf(1:3, 1:3) * [1; 1; 1];
trf(1:3, 4) = offset;
Q = ne_methods.spmm2q(trf);

% create temporary nii object
tnii = xff('new:nii');

% fill NII
tnic = tnii.C;
tnic.FileMagic = niimagic;
tnic.NIIFileType = niitype;
tnic.ImgDim.Dim(1:5) = [3 + double(numvol > 1), vsz, numvol];
tnic.ImgDim.DataType = dtype;
tnic.ImgDim.BitsPerPixel = bpv;
if talorder
    tnic.ImgDim.PixSpacing(1:5) = [1, res, tr];
else
    tnic.ImgDim.PixSpacing(1:5) = [-1, res, tr];
end
tnic.ImgDim.VoxOffset = niioff;
tnic.DataHist.Description = fname;
tnic.DataHist.NIftI1.QFormCode = 2;
tnic.DataHist.NIftI1.SFormCode = 3;
tnic.DataHist.NIftI1.QuaternionB = Q(1);
tnic.DataHist.NIftI1.QuaternionC = Q(2);
tnic.DataHist.NIftI1.QuaternionD = Q(3);
tnic.DataHist.NIftI1.QuatOffsetX = offset(1);
tnic.DataHist.NIftI1.QuatOffsetY = offset(2);
tnic.DataHist.NIftI1.QuatOffsetZ = offset(3);
tnic.DataHist.NIftI1.AffineTransX = trf(1, :);
tnic.DataHist.NIftI1.AffineTransY = trf(2, :);
tnic.DataHist.NIftI1.AffineTransZ = trf(3, :);
tnic.DataHist.NIftI1.IntentName = [upper(ftype) '2NII export'];

% also store RunTimeVars (for CMP/VMP type)
if any(strcmp(ftype, {'cmp', 'glm', 'vmp'}))
    oofields = struct;
    for cf = {'ProjectType', 'ProjectTypeRFX', 'NrOfSubjects', ...
            'NrOfSubjectPredictors', 'NrOfTimePoints', 'NrOfConfounds', ...
            'NrOfStudiesWithConfounds', 'NrOfConfoundsPerStudy', ...
            'SeparatePredictors', 'TransformationType', 'SerialCorrelation', ...
            'Resolution', 'XStart', 'XEnd', 'YStart', 'YEnd', ...
            'ZStart', 'ZEnd', 'NrOfVoxelsForBonfCorrection', ...
            'Study', 'Predictor'}
        if isfield(bc, cf{1})
            oofields.(cf{1}) = bc.(cf{1});
        end
    end
    tnic.RunTimeVars.OriginalObjectFields = oofields;
    mf = fieldnames(tnic.RunTimeVars.Map);
    tnic.RunTimeVars.Map = tnic.RunTimeVars.Map(ones(1, numvol));
    if ~strcmp(ftype, 'glm')
        srcmap = bc.Map;
    else
        srcmap = bc.RunTimeVars.Map;
    end
    mf(ne_methods.multimatch(mf, fieldnames(srcmap)) < 1) = [];
    for vc = 1:numvol
        for fc = 1:numel(mf)
            tnic.RunTimeVars.Map(vc).(mf{fc}) = srcmap(vc).(mf{fc});
        end
        tnic.RunTimeVars.Map(vc).Name = mnames{vc};
    end
    tnic.RunTimeVars.StatsObject = true;
    tnic.RunTimeVars.StatsObjectType = upper(ftype);
    mf = fieldnames(bc.RunTimeVars);
    for fc = 1:numel(mf)
        if ~any(strcmp(mf{fc}, {'Map'}))
            tnic.RunTimeVars.(mf{fc}) = bc.RunTimeVars.(mf{fc});
        end
    end
end

% save NII file
tnii.C = tnic;
try
    aft_SaveAs(tnii, niifile);
    if any(strcmp(ftype, {'cmp', 'glm', 'vmp'}))
        aft_SaveRunTimeVars(tnii);
    end
catch xfferror
    delete(tnii);
    rethrow(xfferror);
end

% clear object
delete(tnii);

% update .Dim field
niifid = fopen(niifile, 'r+', 'ieee-le');
fseek(niifid, 40, -1);
fwrite(niifid, tnic.ImgDim.Dim, 'int16');
fclose(niifid);

% perform as transio
if niitype < 2
    niifile = regexprep(niifile, '\.hdr$', '.img', 'preservecase');
end

% create transio object (extend/create if necessary)
tio = transio(niifile, 'ieee-le', dclass, niioff, [vsz, numvol], true);

% write each volume
for vc = 1:numvol
    try
        vcont = double(aft_GetVolume(xo, range(vc)));
        if talorder
            tio(:, :, :, vc) = permute(vcont(end:-1:1, end:-1:1, end:-1:1), [3, 1, 2]);
        else
            tio(:, :, :, vc) = vcont;
        end
    catch xfferror
        rethrow(xfferror);
    end
end
