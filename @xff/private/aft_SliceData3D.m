function [sag, cor, tra, sc, odata, oalpha] = aft_SliceData3D(xo, crd, opts)
% AFT::SliceData3D  - slice spatial data
%
% FORMAT:       [sag, cor, tra, sc, data] = obj.SliceData3D(crd [, opts]);
%
% Input fields:
%
%       crd         1x3 coordinate to create oblique slices at
%       opts        optional settings
%        .cctype    cross-correlation map data, either of 'lag' or {'r'}
%        .conv      convention, either 'neurological' or {'radiological'}
%        .dir       direction, either of {'all'}, 'cor', 'sag', 'tra'
%        .draw      either of {'0'}, 'c', 's', 't', or '3'
%        .drawcode  1x4 cell array with code, factor, smoothness, and undo
%                   values, default: no drawing (setting .draw to '0')
%        .drawrad   drawing radius (only required once, default: 3mm)
%        .fmr       required for AMR/MAP objects
%        .frame     min/max values (default: 256 box)
%        .gradient  boolean flag, if true return 2D gradient (false)
%        .mapvol    required for CMP/VMP objects, defaults to 1
%        .method    resampling, one of {'cubic'}, 'lanczos3', 'linear', 'nearest'
%        .snmat     SPM-compatible normalization structure
%        .space     either of 'bvc', 'bvi', 'bvs', {'tal'}
%        .step      step size (default: 1)
%        .trans     4x4 transformation matrix, used for volume datatypes
%        .trf       TRF object(s) passed on to samplefmrspace
%        .v16       for VMR objects, if true access available V16 data
%
% Output fields:
%
%       sag         SAG slice data
%       cor         COR slice data
%       tra         TRA slice data
%       sc          internal sampling coordinate
%       data        entire data slab (used for sampling)
%
% TYPES: AMR, AVA, CMP, DMR, DDT, FMR, GLM, HDR, HEAD, MAP, MGH, MSK, NLF, TVL, VDW, VMP, VMR, VTC
%
% Note: for AMR/DMR/FMR/MAP objects, the trf option must be set with
%       at least the IA/FA (or FA alone) TRF object(s) to work
%
% Using: applyspmsnc, bvcoordconv, findfirst, flexinterpn_method,
%        indexarraynb, limitrangec, lsqueeze, samplefmrspace.

% Version:  v1.1
% Build:    16071915
% Date:     Jul-19 2016, 3:46 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% global neuroelf methods
global ne_methods;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true) || ~isa(crd, 'double') || ...
    numel(crd) ~= 3 || any(isinf(crd) | isnan(crd) | crd < -256 | crd > 256)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'cctype') || ~ischar(opts.cctype) || isempty(opts.cctype) || ...
   ~any(lower(opts.cctype(1)) == 'lr')
    opts.cctype = 'r';
else
    opts.cctype = lower(opts.cctype(1));
end
if ~isfield(opts, 'conv') || ~ischar(opts.conv) || isempty(opts.conv) || ...
   ~any(lower(opts.conv(1)) == 'nr')
    opts.conv = 'r';
else
    opts.conv = lower(opts.conv(1));
end
if ~isfield(opts, 'dir') || ~ischar(opts.dir) || ...
   ~any(strcmpi(opts.dir(:)', {'a', 'all', 'c', 'cor', 's', 'sag', 't', 'tra'}))
    opts.dir = 'a';
else
    opts.dir = lower(opts.dir(1));
end
if ~isfield(opts, 'draw') || ~ischar(opts.draw) || numel(opts.draw) ~= 1 || ...
   ~any(opts.draw == '03cstCST')
    opts.draw = '0';
else
    opts.draw = lower(opts.draw);
end
if ~isfield(opts, 'drawcode') || ~iscell(opts.drawcode) || numel(opts.drawcode) ~= 4 || ...
   ~isa(opts.drawcode{1}, 'double') || isempty(opts.drawcode{1}) || ...
    any(isinf(opts.drawcode{1}(:)) | isnan(opts.drawcode{1}(:))) || ...
   ~isa(opts.drawcode{2}, 'double') || numel(opts.drawcode{2}) ~= numel(opts.drawcode{1}) || ...
    any(isinf(opts.drawcode{2}(:)) | isnan(opts.drawcode{2}(:))) || ...
   ~isa(opts.drawcode{3}, 'double') || numel(opts.drawcode{3}) ~= 1 || ...
    isinf(opts.drawcode{3}) || isnan(opts.drawcode{3}) || opts.drawcode{3} < 0 || ...
    numel(opts.drawcode{4}) ~= 1 || ~islogical(opts.drawcode{4})
    opts.draw = '0';
end
if ~isfield(opts, 'drawrad') || ~isa(opts.drawrad, 'double') || numel(opts.drawrad) ~= 1 || ...
    isinf(opts.drawrad) || isnan(opts.drawrad)
    opts.drawrad = 0;
else
    opts.drawrad = min(128, opts.drawrad);
end
if ~isfield(opts, 'fmr') || numel(opts.fmr) ~= 1 || ...
   (~xffisobject(opts.fmr, true, 'dmr') && ~xffisobject(opts.fmr, true, 'fmr'))
    opts.fmr = [];
end
if ~isfield(opts, 'frame') || ~isa(opts.frame, 'double') || ~isequal(size(opts.frame), [2, 3]) || ...
    any(isinf(opts.frame(:)) | isnan(opts.frame(:)) | opts.frame(:) < -512 | opts.frame(:) > 512)
    opts.frame = [];
end
if ~isfield(opts, 'gradient') || numel(opts.gradient) ~= 1 || ~islogical(opts.gradient)
    opts.gradient = false;
end
if ~isfield(opts, 'gradtype') || ~ischar(opts.gradtype) || ...
   ~any(strcmpi(opts.gradtype(:)', {'1', '1s', '1l', '2s', '2l', '3', '3s', '3l'}))
    opts.gradtype = '1s';
else
    opts.gradtype = lower(opts.gradtype(:)');
    if numel(opts.gradtype) == 1
        opts.gradtype(2) = 's';
    end
end
if ~isfield(opts, 'mapvol') || numel(opts.mapvol) ~= 1 || ~isa(opts.mapvol, 'double') || ...
    isinf(opts.mapvol) || isnan(opts.mapvol) || opts.mapvol < 1
    opts.mapvol = 1;
end
if ~isfield(opts, 'method') || ~ischar(opts.method) || ...
   ~any(strcmpi(opts.method, {'cubic', 'lanczos3', 'linear', 'nearest'}))
    opts.method = 'cubic';
else
    opts.method = lower(opts.method(:)');
end
rtv = xo.C.RunTimeVars;
if isfield(opts, 'snmat') && islogical(opts.snmat) && numel(opts.snmat) == 1 && ...
    opts.snmat && isfield(rtv, 'SPMsn') && numel(rtv.SPMsn) == 1 && isstruct(rtv.SPMsn)
    opts.snmat = rtv.SPMsn;
end
if ~isfield(opts, 'snmat') || numel(opts.snmat) ~= 1 || ~isstruct(opts.snmat) || ...
   ~isfield(opts.snmat, 'VF') || ~isfield(opts.snmat, 'VG') || ...
   ~isfield(opts.snmat, 'Tr') || ~isfield(opts.snmat, 'Affine') || ...
   ~isstruct(opts.snmat.VF) || numel(opts.snmat.VF) ~= 1 || ...
   ~isfield(opts.snmat.VF, 'mat') || ~isa(opts.snmat.VF.mat, 'double') || ...
   ~isequal(size(opts.snmat.VF.mat), [4, 4]) || ~isstruct(opts.snmat.VG) || ...
    isempty(opts.snmat.VG) || ~isfield(opts.snmat.VG, 'dim') || ~isfield(opts.snmat.VG, 'mat') || ...
   ~isa(opts.snmat.VG(1).dim, 'double') || ~isequal(size(opts.snmat.VG(1).dim), [1, 3]) || ...
   ~isa(opts.snmat.VG(1).mat, 'double') || ~isequal(size(opts.snmat.VG(1).mat), [4, 4]) || ...
   ~isa(opts.snmat.Tr, 'double') || ndims(opts.snmat.Tr) ~= 4 || ...
   ~isa(opts.snmat.Affine, 'double') || ~isequal(size(opts.snmat.Affine), [4, 4])
    opts.snmat = [];
end
snm = opts.snmat;
if ~isfield(opts, 'space') || ~ischar(opts.space) || ...
   ~any(strcmpi(opts.space(:)', {'bvc', 'bvi', 'bvs', 'tal'}))
    opts.space = 'tal';
else
    opts.space = lower(opts.space(:)');
end
if ~isfield(opts, 'step') || ~isa(opts.step, 'double') || ...
    isinf(opts.step) || isnan(opts.step) || opts.step <= 0 || opts.step >= 16
    opts.step = 1;
end
if ~isfield(opts, 'tal') || numel(opts.tal) ~= 1 || ~xffisobject(opts.tal, true, 'tal')
    opts.tal = {};
end
if ~isfield(opts, 'trans') || ~isa(opts.trans, 'double') || ~isequal(size(opts.trans), [4, 4]) || ...
    any(isinf(opts.trans(:)) | isnan(opts.trans(:))) || any(opts.trans(4, :) ~= [0, 0, 0, 1])
    opts.trans = [];
end
if ~isfield(opts, 'trf') || (~iscell(opts.trf) && (numel(opts.trf) ~= 1 || ...
    ~xffisobject(opts.trf, true, 'trf')))
    opts.trf = {};
elseif ~iscell(opts.trf)
    opts.trf = {opts.trf};
end
for tc = numel(opts.trf):-1:1
    if numel(opts.trf{tc}) ~= 1 || (~xffisobject(opts.trf{tc}, true, 'tal') && ...
        ~xffisobject(opts.trf{tc}, true, 'trf'))
        opts.trf(tc) = [];
    end
end
opts.trf = opts.trf(:)';
if ~isfield(opts, 'v16') || ~islogical(opts.v16) || numel(opts.v16) ~= 1
    opts.v16 = false;
end
epst = eps ^ 0.75;

% alpha defaults to 1
oalpha = 1;

% switch on filetype
bc = xo.C;
rtv = bc.RunTimeVars;
ft = lower(xo.S.Extensions{1});

% check TrfPlus being different from eye(4)
trfplus = [];
if isfield(rtv, 'TrfPlus') && isequal([4, 4], size(rtv.TrfPlus)) && ...
    any(any(rtv.TrfPlus ~= eye(4)))
    trfplus = inv(rtv.TrfPlus);
end

% for DMR/FMR/VTC objects, resolve transio, otherwise accept map
resetc = false;
sci = 0;
scs = 1;
switch (ft)
    case 'amr'
        if isempty(opts.fmr)
            error('neuroelf:xff:missingOption', ...
                'Option field .fmr must be set to access AMR in 3D space.');
        end
        odata = bc.Slice(1).AMRData(1);
        odata(1) = 0;
        odata = repmat(odata, [size(bc.Slice(1).AMRData), numel(bc.Slice)]);
        for slc = 1:numel(bc.Slice)
            if istransio(bc.Slice(slc).AMRData)
                bc.Slice(slc).AMRData = resolve(bc.Slice(slc).AMRData);
                resetc = true;
            end
            odata(:, :, slc) = bc.Slice(slc).AMRData;
        end
        if opts.draw ~= '0' && (~isfield(rtv, 'UndoBuffer') || ...
            ~isequal(size(rtv.UndoBuffer), size(odata)))
            bc.RunTimeVars.UndoBuffer = odata;
            resetc = true;
        end
    case 'ava'
        if bc.ProjectType > 2
            error('neuroelf:xff:notImplemented', ...
                'Slicing of SRF-based component Maps not supported.');
        end
        if bc.ProjectType == 1 && isempty(opts.fmr)
            error('neuroelf:xff:missingOption', ...
                'Option field .fmr must be set to access AVA in 3D space.');
        end
        fvol = fieldnames(bc.Maps);
        nvol = numel(fvol) + (size(bc.Maps.CellMeans, ndims(bc.Maps.CellMeans)) - 1);
        opts.mapvol = round(min(nvol, opts.mapvol));
        if opts.mapvol < numel(fvol)
            if istransio(bc.Maps.(fvol{opts.mapvol}))
                bc.Maps.(fvol{opts.mapvol}) = resolve(bc.Maps.(fvol{opts.mapvol}));
                xo.C = bc;
            end
            odata = bc.Maps.(fvol{opts.mapvol});
        else
            if istransio(bc.Maps.CellMeans)
                bc.Maps.CellMeans = resolve(bc.Maps.CellMeans);
            end
            odata = bc.Maps.CellMeans;
            odatar = repmat({':'}, 1, ndims(odata));
            odatar{end} = opts.mapvol + 1 - numel(fvol);
            odata = odata(odatar{:});
        end
    case 'cmp'
        if bc.DocumentType > 1
            error('neuroelf:xff:notImplemented', ...
                'Slicing of SRF-based component Maps not supported.');
        end
        if bc.DocumentType == 0 && isempty(opts.fmr)
            error('neuroelf:xff:missingOption', ...
                'Option field .fmr must be set to access fCMP in 3D space.');
        end
        nvol = numel(bc.Map);
        opts.mapvol = round(min(nvol, opts.mapvol));
        if istransio(bc.Map(opts.mapvol).CMPData)
            bc.Map(opts.mapvol).CMPData = resolve(bc.Map(opts.mapvol).CMPData);
            xo.C = bc;
        end
        odata = bc.Map(opts.mapvol).CMPData;
        if opts.draw ~= '0' && ...
           (~isfield(bc.Map(opts.mapvol).RunTimeVars, 'UndoBuffer') || ...
            ~isequal(size(bc.Map(opts.mapvol).RunTimeVars.UndoBuffer), size(odata)))
            bc.Map(opts.mapvol).RunTimeVars.UndoBuffer = odata;
            resetc = true;
        end
        if bc.Map(opts.mapvol).EnableClusterCheck && ...
            isequal(size(odata), size(bc.Map(opts.mapvol).CMPDataCT))
            odata = odata .* bc.Map(opts.mapvol).CMPDataCT;
        end
    case 'dmr'
        if istransio(bc.DWIData)
            bc.DWIData = resolve(bc.DWIData);
            xo.C = bc;
        end
        switch bc.DataStorageFormat
            case 2
                svol = 3;
            case 3
                svol = 4;
            case 4
                svol = 1;
            otherwise
                error('neuroelf:xff:illegalSetting', ...
                    'Illegal DataStorageFormat setting in DMR object.');
        end
        nvol = size(bc.DWIData, svol);
        opts.mapvol = round(min(nvol, opts.mapvol));
        switch (svol)
            case 1
                odata = squeeze(bc.DWIData(opts.mapvol, :, :, :));
            case 3
                odata = squeeze(bc.DWIData(:, :, opts.mapvol, :));
            case 4
                odata = squeeze(bc.DWIData(:, :, :, opts.mapvol));
        end
    case 'fmr'
        for slc = 1:numel(bc.Slice)
            if istransio(bc.Slice(slc).STCData)
                bc.Slice(slc).STCData = resolve(bc.Slice(slc).STCData);
                resetc = true;
            end
        end
        szdt = size(bc.Slice(1).STCData);
        nvol = szdt(3);
        opts.mapvol = round(min(nvol, opts.mapvol));
        switch bc.DataStorageFormat
            case 1
                odata = zeros([szdt(1:2), numel(bc.Slice)]);
                for slc = 1:numel(bc.Slice)
                    odata(:, :, slc) = bc.Slice(slc).STCData(:, :, opts.mapvol);
                end
            case 2
                odata = squeeze(bc.Slice.STCData(:, :, opts.mapvol, :));
            otherwise
                error('neuroelf:xff:illegalSetting', ...
                    'Illegal DataStorageFormat setting in DMR object.');
        end
    case 'glm'
        nvol = numel(bc.Predictor);
        opts.mapvol = round(min(nvol, opts.mapvol));
        [odata, trfplusi, snmat] = aft_GetVolume(xo, opts.mapvol);
        if isempty(trfplus)
            trfplus = inv(trfplusi);
        else
            trfplus = inv(trfplusi) * trfplus;
        end
        if ~isempty(snmat)
            snm = snmat;
        end
    case 'hdr'
        if any(bc.ImgDim.DataType == [32, 128, 1536, 1792, 2048, 2304])
            switch (bc.ImgDim.DataType)
                case {32, 1792}
                    if istransio(bc.VoxelData)
                        bc.VoxelData = resolve(bc.VoxelData);
                        resetc = true;
                    end
                    if istransio(bc.VoxelDataComplex)
                        bc.VoxelDataComplex = resolve(bc.VoxelDataComplex);
                        resetc = true;
                    end
                    odata = complex(bc.VoxelData, bc.VoxelDataComplex);
                case {128, 2304}
                    if istransio(bc.VoxelDataRGBA)
                        bc.VoxelDataRGBA = resolve(bc.VoxelDataRGBA);
                        resetc = true;
                    end
                    odata = bc.VoxelDataRGBA;
                otherwise
                    error('neuroelf:xff:unsupported', ...
                        'Slicing of this datatype not yet supported.');
            end
        else
            if istransio(bc.VoxelData)
                bc.VoxelData = resolve(bc.VoxelData);
                resetc = true;
            end
            odata = bc.VoxelData;
        end
        if opts.draw ~= '0' && (~isfield(rtv, 'UndoBuffer') || ...
            ~isequal(size(rtv.UndoBuffer), size(odata)))
            bc.RunTimeVars.UndoBuffer = odata;
            resetc = true;
        end
        nvol = size(odata, 4);
        if nvol > 1
            opts.mapvol = round(min(nvol, opts.mapvol));
            odata = odata(:, :, :, opts.mapvol, :);
        else
            opts.mapvol = 1;
        end

        % scaling
        if any([2, 4, 8, 130, 132, 136, 256, 512, 768] == bc.ImgDim.DataType)
            sci = bc.ImgDim.ScalingIntercept;
            scs = bc.ImgDim.ScalingSlope;
        end

        % clustered ?
        if rtv.Map(opts.mapvol).EnableClusterCheck && ...
            numel(bc.VoxelDataCT) >= opts.mapvol && ...
            isequal(size(odata), size(bc.VoxelDataCT{opts.mapvol}))
            if isa(odata, 'double') || isa(odata, 'single')
                odata = odata .* bc.VoxelDataCT{opts.mapvol};
            else
                odata(~bc.VoxelDataCT{opts.mapvol}) = 0;
            end
        end

        % also get transformation matrix and set space to bvc!
        cfr = hdr_CoordinateFrame(xo, opts.mapvol);
        if isempty(opts.trans)
            opts.trans = inv(cfr.Trf);
        else
            opts.trans = inv(cfr.Trf) * opts.trans;
        end
        if isempty(opts.frame)
            opts.frame = [128, 128, 128; -127.99, -127.99, -127.99];
        end
        opts.space = 'bvc';
    case 'head'
        nvol = numel(bc.Brick);
        opts.mapvol = round(min(nvol, opts.mapvol));
        if istransio(bc.Brick(opts.mapvol).Data)
            bc.Brick(opts.mapvol).Data = resolve(bc.Brick(opts.mapvol).Data);
            resetc = true;
        end
        odata = bc.Brick(opts.mapvol).Data;
        if opts.draw ~= '0' && (~isfield(rtv, 'UndoBuffer') || ...
            ~isequal(size(rtv.UndoBuffer), size(odata)))
            bc.Brick(opts.mapvol).RunTimeVars.UndoBuffer = odata;
        end

        % scaling
        if ~isempty(regexpi(class(odata), 'int')) && bc.Brick(opts.mapvol).ScalingFactor ~= 1
            scs = bc.Brick(opts.mapvol).ScalingFactor;
        end

        % clustered ?
        if rtv.Map(opts.mapvol).EnableClusterCheck && ...
            isequal(size(odata), size(bc.Brick(opts.mapvol).DataCT))
            odata = odata .* bc.Brick(opts.mapvol).DataCT;
        end

        % also get transformation matrix and set space to bvc!
        cfr = head_CoordinateFrame(xo);
        if isempty(opts.trans)
            opts.trans = inv(cfr.Trf);
        else
            opts.trans = inv(cfr.Trf) * opts.trans;
        end
        opts.trans(abs(opts.trans) < epst) = 0;
        if isempty(opts.frame)
            opts.frame = [128, 128, 128; -127.99, -127.99, -127.99];
        end
        opts.space = 'bvc';
    case 'map'
        if isempty(opts.fmr)
            error('neuroelf:xff:missingOption', ...
                'Option field .fmr must be set to access MAP in 3D space.');
        end
        odata = bc.Map(1).Data(1);
        odata(1) = 0;
        odata = repmat(odata, [size(bc.Map(1).Data), numel(bc.Map)]);
        for slc = 1:numel(bc.Map)
            if istransio(bc.Map(slc).Data)
                bc.Map(slc).Data = resolve(bc.Map(slc).Data);
                resetc = true;
            end
            odata(:, :, slc) = bc.Map(slc).Data;
        end
        if opts.draw ~= '0' && (~isfield(rtv, 'UndoBuffer') || ...
            ~isequal(size(rtv.UndoBuffer), size(odata)))
            bc.RunTimeVars.UndoBuffer = odata;
            resetc = true;
        end
    case 'mgh'
        opts.mapvol = 1;
        if istransio(bc.MGHData)
            bc.MGHData = resolve(bc.MGHData);
            resetc = true;
        end
        odata = bc.MGHData;
        if opts.draw ~= '0' && (~isfield(rtv, 'UndoBuffer') || ...
            ~isequal(size(rtv.UndoBuffer), size(odata)))
            bc.RunTimeVars.UndoBuffer = odata;
        end

        % also get transformation matrix and set space to bvc!
        cfr = mgh_CoordinateFrame(xo);
        if isempty(opts.trans)
            opts.trans = inv(cfr.Trf);
        else
            opts.trans = inv(cfr.Trf) * opts.trans;
        end
        opts.trans(abs(opts.trans) < epst) = 0;
        if isempty(opts.frame)
            opts.frame = [128, 128, 128; -127.99, -127.99, -127.99];
        end
        opts.space = 'bvc';
    case 'msk'
        if istransio(bc.Mask)
            bc.Mask = resolve(bc.Mask);
            resetc = true;
        end
        odata = bc.Mask;
        if opts.draw ~= '0' && (~isfield(rtv, 'UndoBuffer') || ...
            ~isequal(size(rtv.UndoBuffer), size(odata)))
            bc.RunTimeVars.UndoBuffer = odata;
            resetc = true;
        end
    case 'nlf'
        if isempty(regexpi(bc.DimMeaning, 'xyz'))
            error('neuroelf:xff:badSubType', 'Invalid subtype for NLF::SliceData3D.');
        end
        xd = find(bc.DimMeaning == 'x');
        td = find(bc.DimMeaning ~= 'x' & bc.DimMeaning ~= 'y' & bc.DimMeaning ~= 'z');
        if istransio(bc.Data)
            bc.Data = resolve(bc.Data);
            resetc = true;
        end
        odata = bc.Data;
        if opts.draw ~= '0' && (~isfield(rtv, 'UndoBuffer') || ...
            ~isequal(size(rtv.UndoBuffer), size(odata)))
            bc.RunTimeVars.UndoBuffer = odata;
            resetc = true;
        end
        if ~isempty(td)
            nvol = size(odata, td(1));
        else
            nvol = 1;
        end
        if nvol > 1
            opts.mapvol = round(min(nvol, opts.mapvol));
            sra = {1, 1, 1, 1, 1, 1, 1, 1};
            sra{xd} = ':';
            sra{xd+1} = ':';
            sra{xd+2} = ':';
            sra{td} = opts.mapvol;
            odata = odata(sra{:});
        else
            opts.mapvol = 1;
        end

        % scaling
        if (bc.ScalingIntercept ~= 0 || all([0, 1] ~= bc.ScalingSlope))
            sci = bc.ScalingIntercept;
            scs = bc.ScalingSlope;
        end

        % also get transformation matrix and set space to bvc!
        cfr = nlf_CoordinateFrame(xo);
        if isempty(opts.trans)
            opts.trans = inv(cfr.Trf);
        else
            opts.trans = inv(cfr.Trf) * opts.trans;
        end
        opts.trans(abs(opts.trans) < epst) = 0;
        if isempty(opts.frame)
            opts.frame = [128, 128, 128; -127.99, -127.99, -127.99];
        end
        opts.space = 'bvc';
    case 'srf'
        if isempty(bc.VertexVMRData)
            srf_BackToVMR(xo);
            bc = xo.C;
        end
        odata = bc.VertexVMRData;
    case 'vdw'
        if istransio(bc.VDWData)
            bc.VDWData = resolve(bc.VDWData);
            xo.C = bc;
        end
        nvol = size(bc.VDWData, 1);
        opts.mapvol = round(min(nvol, opts.mapvol));
        odata = squeeze(bc.VDWData(opts.mapvol, :, :, :));
    case 'vmp'
        nvol = numel(bc.Map);
        opts.mapvol = round(min(nvol, opts.mapvol));
        if istransio(bc.Map(opts.mapvol).VMPData)
            bc.Map(opts.mapvol).VMPData = resolve(bc.Map(opts.mapvol).VMPData);
            xo.C = bc;
        end
        odata = bc.Map(opts.mapvol).VMPData;
        if isfield(bc.Map(opts.mapvol), 'RunTimeVars') && ...
            isstruct(bc.Map(opts.mapvol).RunTimeVars) && ...
            isfield(bc.Map(opts.mapvol).RunTimeVars, 'VMPDataError') && ...
            isfield(bc.Map(opts.mapvol).RunTimeVars, 'LowerErrorThresh') && ...
            isfield(bc.Map(opts.mapvol).RunTimeVars, 'UpperErrorThresh')
            lerror = bc.Map(opts.mapvol).RunTimeVars.LowerErrorThresh;
            uerror = bc.Map(opts.mapvol).RunTimeVars.UpperErrorThresh;
            escale = 1 ./ max(sqrt(eps), uerror - lerror);
            oalpha = ne_methods.limitrangec(escale .* (abs(odata) ./ ...
                bc.Map(opts.mapvol).RunTimeVars.VMPDataError - lerror), 0, 1, 0);
        end
        if opts.draw ~= '0' && ...
           (~isfield(bc.Map(opts.mapvol).RunTimeVars, 'UndoBuffer') || ...
            ~isequal(size(bc.Map(opts.mapvol).RunTimeVars.UndoBuffer), size(odata)))
            bc.Map(opts.mapvol).RunTimeVars.UndoBuffer = odata;
            resetc = true;
        end
        if bc.Map(opts.mapvol).EnableClusterCheck && ...
            isequal(size(odata), size(bc.Map(opts.mapvol).VMPDataCT))
            odata = odata .* bc.Map(opts.mapvol).VMPDataCT;
        end
        if bc.Map(opts.mapvol).Type == 3
            if opts.cctype == 'r'
                odata = odata - floor(odata);
            else
                odata = floor(odata);
            end
        end
    case 'vmr'
        if opts.draw ~= '0'
            if ~isempty(bc.VMRData16) && (~isfield(rtv, 'UndoBuffer16') || ...
                ~isequal(size(rtv.UndoBuffer16), size(bc.VMRData16)))
                bc.RunTimeVars.UndoBuffer16 = bc.VMRData16(:, :, :);
                resetc = true;
            end
            if ~isfield(rtv, 'UndoBuffer') || ~isequal(size(rtv.UndoBuffer), size(bc.VMRData))
                bc.RunTimeVars.UndoBuffer = bc.VMRData(:, :, :);
                resetc = true;
            end
        end
        if ~isempty(bc.VMRData16) && opts.v16
            if istransio(bc.VMRData16)
                bc.VMRData16 = resolve(bc.VMRData16);
                resetc = true;
            end
            odata = bc.VMRData16;
        else
            if istransio(bc.VMRData)
                bc.VMRData = resolve(bc.VMRData);
                resetc = true;
            end
            odata = bc.VMRData;
        end
    case 'vtc'
        if istransio(bc.VTCData)
            bc.VTCData = resolve(bc.VTCData);
            xo.C = bc;
        end
        if ~isfield(rtv, 'AvgVTC') || ~rtv.AvgVTC
            nvol = size(bc.VTCData, 1);
            opts.mapvol = round(min(nvol, opts.mapvol));
            avgvtc = false;
            condi = 1;
            tcpc = 0;
            vptc = 0;
        else
            avgvtc = true;
            nvol = numel(rtv.ConditionNames);
            if ~isfield(rtv, 'SubMapVol') || ~isa(rtv.SubMapVol, 'double') || ...
                numel(rtv.SubMapVol) ~= 1
                rtv.SubMapVol = 1;
            end
            vptc = rtv.NrOfVolumesPerTC;
            smv = max(min(rtv.SubMapVol, vptc), 1);
            condi = round(min(nvol, opts.mapvol));
            tcpc = rtv.NrOfTCsPerCondition;
            cthr = rtv.ConditionThresholds(condi, 2, :);
            opts.mapvol = smv;
            if opts.mapvol == round(opts.mapvol)
                oalpha = double(squeeze(bc.VTCData((condi - 1) * tcpc * vptc + opts.mapvol + vptc, :, :, :)));
            else
                smv = round(opts.mapvol);
                smd = opts.mapvol - smv;
                rweights = ne_methods.flexinterpn_method([0; 0; 0; 1; 0; 0; 0], ...
                    [Inf; 1 - smd;1;7], 'cubic');
                rweight1 = ne_methods.findfirst(rweights ~= 0);
                rweights = rweights(rweight1:ne_methods.findfirst(rweights ~= 0, -1));
                smv = smv - 4 + rweight1;
                if smv < 1
                    rweights = rweights(2 - smv:end);
                    smv = 1;
                end
                while (smv + numel(rweights) - 1) > rtv.NrOfVolumesPerTC
                    rweights(end) = [];
                end
                smv = smv + (condi - 1) * tcpc * vptc;
                if isempty(rweights)
                    rweights = 1;
                end
                odatasz = size(bc.VTCData);
                odatasz(1) = [];
                oalpha = zeros(odatasz);
            end
        end
        if opts.mapvol == round(opts.mapvol)
            odata = squeeze(bc.VTCData(opts.mapvol + (condi - 1) * tcpc * vptc, :, :, :));
            if avgvtc
                if strcmpi(rtv.TCNames{2}, 'sd')
                    oalpha = (1 / sqrt(rtv.NrOfConditionOnsets(condi))) .* oalpha;
                end
                if lower(rtv.TCNames{2}(1)) == 's'
                    oalpha = abs(odata) ./ abs(oalpha);
                    oalpha(isnan(oalpha)) = 0;
                end
                oalpha = ne_methods.limitrangec((1 / max(cthr(2) - cthr(1), sqrt(eps))) .* (oalpha - cthr(1)), 0, 1, 0);
            end
        else
            odata = zeros(odatasz);
            for rwc = 1:numel(rweights)
                odatap = squeeze(bc.VTCData(smv+rwc-1, :, :, :));
                oalphap = squeeze(bc.VTCData(vptc+smv+rwc-1, :, :, :));
                if strcmpi(rtv.TCNames{2}, 'sd')
                    oalphap = (1 / sqrt(rtv.NrOfConditionOnsets(condi))) .* oalphap;
                end
                if lower(rtv.TCNames{2}(1)) == 's'
                    oalphap = abs(odatap) ./ abs(oalphap);
                    oalphap(isinf(oalphap) | isnan(oalphap)) = 0;
                end
                odata = odata + rweights(rwc) .* odatap;
                oalpha = oalpha + rweights(rwc) .* ...
                    ne_methods.limitrangec((1 / max(cthr(2) - cthr(1), sqrt(eps))) .* (oalphap - cthr(1)), 0, 1, 0);
            end
        end
    otherwise
        error('neuroelf:xff:objectTypeUnsupported', ...
            'Unsupported object type for ::SliceData3D method.');
end

% udpate contents?
if resetc
    xo.C = bc;
end

% check frame/stepsize
opts.step = opts.step([1, 1, 1]);
if isempty(opts.frame)
    switch (opts.space)
        case 'tal'
            opts.frame = [128, 128, 128; -127.9999, -127.9999, -127.9999];
        case {'bvi', 'bvs'}
            opts.frame = [0, 0, 0; 255.9999, 255.9999, 255.9999];
        case 'bvc'
            opts.frame = 1 + [0, 0, 0; (size(odata) - 0.01)];
    end
end
if opts.conv ~= 'r'
    opts.frame(:, 1) = opts.frame([2, 1], 1);
end
off = diff(opts.frame) < 0;
opts.step(off) = -opts.step(off);

% sample true BV spatial objects for nearest/linear in a fast way
if isempty(snm) && isempty(opts.trans) && isempty(trfplus) && ...
   (any(strcmp(ft, {'glm', 'msk', 'vdw', 'vmp', 'vmr', 'vtc'})) || ...
    (strcmp(ft, 'cmp') && bc.DocumentType == 1))

    % create three coord lists
    xc = opts.frame(1, 1):opts.step(1):opts.frame(2, 1);
    yc = opts.frame(1, 2):opts.step(2):opts.frame(2, 2);
    zc = opts.frame(1, 3):opts.step(3):opts.frame(2, 3);

    % get object layout
    if ~strcmp(ft, 'vmr')
        ol = [size(odata), 1, bc.XStart, bc.YStart, bc.ZStart, ...
            bc.XEnd, bc.YEnd, bc.ZEnd, bc.Resolution([1, 1, 1])];
    else
        ol = aft_Layout(xo);
        if any(ol(11:13) < 1) && all(ol(11:13) <= 1)
            ol(5:7) = ol(5:7) .* ol(11:13);
        end
    end

    % depending on space
    switch (opts.space)
        case 'tal'
            crd = (1 ./ ol([13, 11, 12])) .* ((128 - crd) - ol([7, 5, 6]));
            xc = (1 / ol(13)) .* ((128 - xc) - ol(7));
            yc = (1 / ol(11)) .* ((128 - yc) - ol(5));
            zc = (1 / ol(12)) .* ((128 - zc) - ol(6));
            nstep = -opts.step(:)' ./ ol([13, 11, 12]);
        case 'bvs'
            crd = (1 ./ ol([13, 11, 12])) .* (crd - ol([7, 5, 6]));
            xc = (1 / ol(13)) .* (xc - ol(7));
            yc = (1 / ol(11)) .* (yc - ol(5));
            zc = (1 / ol(12)) .* (zc - ol(6));
            nstep = opts.step(:)' ./ ol([13, 11, 12]);
        case 'bvi'
            crd = (1 ./ ol([13, 11, 12])) .* crd([3, 1, 2]);
            xs = xc;
            xc = (1 / ol(13)) .* (zc - ol(7));
            zc = (1 / ol(12)) .* (yc - ol(6));
            yc = (1 / ol(11)) .* (xs - ol(5));
            nstep = opts.step([3, 1, 2]) ./ ol([13, 11, 12]);
        case 'bvc'
            crd = crd([3, 1, 2]);
            xs = xc;
            xc = zc;
            zc = yc;
            yc = xs;
            nstep = opts.step([3, 1, 2]);
    end
    sc = crd([2, 3, 1]) + 1;

    % override method if possible
    if all(abs(nstep - round(nstep)) <= 0.00001) && ...
        all(abs([xc(1), yc(1), zc(1)] - round([xc(1), yc(1), zc(1)])) <= 0.002)
        opts.method = 'nearest';
    end

    % prepare outputs
    if opts.dir == 'a'
        sag = zeros(numel(zc), numel(yc));
        cor = zeros(numel(zc), numel(xc));
        tra = zeros(numel(yc), numel(xc));
        if numel(oalpha) > 1
            asag = sag;
            acor = cor;
            atra = tra;
        end
    elseif opts.dir == 's'
        sag = zeros(numel(zc), numel(yc));
        if numel(oalpha) > 1
            asag = sag;
        end
    elseif opts.dir == 'c'
        sag = zeros(numel(zc), numel(xc));
        if numel(oalpha) > 1
            asag = sag;
        end
    else
        sag = zeros(numel(yc), numel(xc));
        if numel(oalpha) > 1
            asag = sag;
        end
    end

    % sample only within valid coordinates
    xv = (xc > -0.5) & (xc < (size(odata, 3)) - 0.5);
    yv = (yc > -0.5) & (yc < (size(odata, 1)) - 0.5);
    zv = (zc > -0.5) & (zc < (size(odata, 2)) - 0.5);
    xf = xc(ne_methods.findfirst(xv));
    yf = yc(ne_methods.findfirst(yv));
    zf = zc(ne_methods.findfirst(zv));
    if isempty(xf) || isempty(yf) || isempty(zf)
        if scs ~= 1
            if sci ~= 0
                sag = scs .* sag + sci;
                if opts.dir == 'a'
                    cor = scs .* cor + sci;
                    tra = scs .* tra + sci;
                end
            else
                sag = scs .* sag;
                if opts.dir == 'a'
                    cor = scs .* cor;
                    tra = scs .* tra;
                end
            end
        elseif sci ~= 0
            sag = sag + sci;
            if opts.dir == 'a'
                cor = cor + sci;
                tra = tra + sci;
            end
        end
        if opts.dir ~= 'a'
            cor = sc;
            tra = odata;
            sc = 0;
        end
        return;
    end

    % now sample only required parts, nearest
    if strcmpi(opts.method, 'nearest')
        crd = 1 + round(crd);
        xc = 1 + round(xc(xv));
        yc = 1 + round(yc(yv));
        zc = 1 + round(zc(zv));
        if crd(1) > 0 && crd(1) <= size(odata, 3) && any('as' == opts.dir)
            sag(zv, yv) = odata(yc, zc, crd(1))';
            if numel(oalpha) > 1
                asag(zv, yv) = oalpha(yc, zc, crd(1))';
            end
        end
        if crd(2) > 0 && crd(2) <= size(odata, 1)
            if opts.dir == 'a'
                cor(zv, xv) = odata(crd(2), zc, xc);
                if numel(oalpha) > 1
                    acor(zv, xv) = oalpha(crd(2), zc, xc);
                end
            elseif opts.dir == 'c'
                sag(zv, xv) = odata(crd(2), zc, xc);
                if numel(oalpha) > 1
                    asag(zv, xv) = oalpha(crd(2), zc, xc);
                end
            end
        end
        if crd(3) > 0 && crd(3) <= size(odata, 2)
            if opts.dir == 'a'
                tra(yv, xv) = odata(yc, crd(3), xc);
                if numel(oalpha) > 1
                    atra(yv, xv) = oalpha(yc, crd(3), xc);
                end
            elseif opts.dir == 't'
                sag(yv, xv) = odata(yc, crd(3), xc);
                if numel(oalpha) > 1
                    asag(yv, xv) = oalpha(yc, crd(3), xc);
                end
            end
        end

    % or higher interp
    else
        crd = crd + 1;
        xt = xf + (sum(xv) - 0.25) * nstep(1);
        yt = yf + (sum(yv) - 0.25) * nstep(2);
        zt = zf + (sum(zv) - 0.25) * nstep(3);

        crdf = [Inf, Inf, Inf; yf, zf, xf; nstep([2, 3, 1]); yt, zt, xt];
        crdf([2, 4], :) = 1 + crdf([2, 4], :);
        if crd(1) > -0.5 && crd(1) < size(odata, 3) && any('as' == opts.dir)
            crds = crdf;
            crds([2, 4], 3) = crd(1);
            sag(zv, yv) = ne_methods.flexinterpn_method(odata, crds, 0, opts.method)';
            if numel(oalpha) > 1
                asag(zv, yv) = ne_methods.flexinterpn_method(oalpha, crds, 0, opts.method)';
            end
        end
        if crd(2) > -0.5 && crd(2) < size(odata, 1)
            crds = crdf;
            crds([2, 4], 1) = crd(2);
            if opts.dir == 'a'
                cor(zv, xv) = ne_methods.flexinterpn_method(odata, crds, 0, opts.method);
                if numel(oalpha) > 1
                    acor(zv, xv) = ne_methods.flexinterpn_method(oalpha, crds, 0, opts.method);
                end
            elseif opts.dir == 'c'
                sag(zv, xv) = ne_methods.flexinterpn_method(odata, crds, 0, opts.method);
                if numel(oalpha) > 1
                    asag(zv, xv) = ne_methods.flexinterpn_method(oalpha, crds, 0, opts.method);
                end
            end
        end
        if crd(3) > -0.5 && crd(3) < size(odata, 2)
            crds = crdf;
            crds([2, 4], 2) = crd(3);
            if opts.dir == 'a'
                tra(yv, xv) = ne_methods.flexinterpn_method(odata, crds, 0, opts.method);
                if numel(oalpha) > 1
                    atra(yv, xv) = ne_methods.flexinterpn_method(oalpha, crds, 0, opts.method);
                end
            elseif opts.dir == 't'
                sag(yv, xv) = ne_methods.flexinterpn_method(odata, crds, 0, opts.method);
                if numel(oalpha) > 1
                    asag(yv, xv) = ne_methods.flexinterpn_method(oalpha, crds, 0, opts.method);
                end
            end
        end
    end

    % stack alpha values
    if numel(oalpha) > 1
        if opts.dir == 'a'
            oalpha = {asag, acor, atra, oalpha};
        else
            oalpha = {asag, [], [], oalpha};
        end
    end

% sample FMR based data differently
elseif any(strcmp(ft, {'amr', 'dmr', 'fmr', 'map'})) || ...
   (strcmp(ft, 'cmp') && bc.DocumentType == 0)

    % build 1-d vectors, then pairs of coords
    xv = opts.frame(1, 1):opts.step(1):opts.frame(2, 1);
    yv = opts.frame(1, 2):opts.step(2):opts.frame(2, 2);
    zv = opts.frame(1, 3):opts.step(3):opts.frame(2, 3);
    [ys, xs] = ndgrid(zv, yv);
    [yc, xc] = ndgrid(zv, xv);
    [yt, xt] = ndgrid(yv, xv);
    if opts.dir == 'a'
        crd = [[crd(1) .* ones(numel(xs), 1), xs(:), ys(:)]; ...
               [xc(:), crd(2) .* ones(numel(xc), 1), yc(:)]; ...
               [xt(:), yt(:), crd(3) .* ones(numel(xt), 1)]; crd];
    elseif opts.dir == 's'
        crd = [[crd(1) .* ones(numel(xs), 1), xs(:), ys(:)]; crd];
    elseif opts.dir == 'c'
        crd = [[xc(:), crd(2) .* ones(numel(xc), 1), yc(:)]; crd];
    else
        crd = [[xt(:), yt(:), crd(3) .* ones(numel(xt), 1)]; crd];
    end

    % convert coordinates to BV system
    if strcmp(opts.space, 'tal')
        if ~isempty(snm)
            crd = ne_methods.applyspmsnc(crd, snm.Tr, ...
                snm.VG(1).dim, inv(snm.VG(1).mat), snm.VF.mat * snm.Affine);
        end
        crd = 128 - crd;
    elseif strcmp(opts.space, 'bvi')
        crd = crd(:, [3, 1, 2]);
        if ~isempty(snm)
            crd = 128 - ne_methods.applyspmsnc(128 - crd, snm.Tr, ...
                snm.VG(1).dim, inv(snm.VG(1).mat), snm.VF.mat * snm.Affine);
        end
    elseif ~strcmp(opts.space, 'bvs')
        error('neuroelf:xff:invalidOption', 'Sampling of DMR/FMR based data must not have bvc space.');
    elseif ~isempty(snm)
        crd = 128 - ne_methods.applyspmsnc(128 - crd, snm.Tr, ...
            snm.VG(1).dim, inv(snm.VG(1).mat), snm.VF.mat * snm.Affine);
    end

    % use AMR/MAP option DMR/FMR
    if ~isempty(opts.fmr)
        xo = opts.fmr;
    end

    % sample data with samplefmrspace function
    [data, sc] = ne_methods.samplefmrspace(odata, crd, xo, [opts.trf, opts.tal], opts.method);
    sc = sc(end, :);

    % split into sag, cor, tra
    if any('as' == opts.dir)
    end
    if opts.dir == 'a'
        sag = reshape(data(1:numel(xs)), size(xs));
        cor = reshape(data(numel(xs)+1:numel(xs)+numel(xc)), size(xc));
        tra = reshape(data(numel(xs)+numel(xc)+1:numel(xs)+numel(xc)+numel(xt)), size(xt));
    elseif opts.dir == 's'
        sag = reshape(data(1:numel(xs)), size(xs));
    elseif opts.dir == 'c'
        sag = reshape(data(1:numel(xc)), size(xc));
    elseif opts.dir == 't'
        sag = reshape(data(1:numel(xt)), size(xt));
    end

% sample "already 3D" formats
else

    % get bounding box for valid objects and get transformation
    if ~strcmp(opts.space, 'bvc')
        ftrf = {ne_methods.bvcoordconv(zeros(0, 3), [opts.space '2bvc'], aft_BoundingBox(xo))};
    else
        ftrf = {};
    end

    % any additional transformation?
    if ~isempty(opts.trans)
        if isempty(ftrf)
            ftrf = {opts.trans};
        else
            ftrf = {ftrf{1} * opts.trans};
            ftrf{1}(abs(ftrf{1}) < epst) = 0;
        end
    end

    % and yet additional from TrfPlus
    if ~isempty(trfplus)
        if isempty(ftrf)
            ftrf = {trfplus};
        else
            ftrf = {ftrf{1} * trfplus};
            ftrf{1}(abs(ftrf{1}) < epst) = 0;
        end
    end

    % make sure ones are ones!
    if ~isempty(ftrf)
        ftrf{1}(abs(ftrf{1} - 1) < epst) = 1;

        % and that only to apply if necessary
        if all(all(ftrf{1} == diag(ones(4, 1))))
            ftrf = {};
        end
    end

    % sample data (without normalization)
    if isempty(snm)

        % regular sampling (full coordinates)
        if ~isempty(ftrf)
            ftrf1611 = ne_methods.lsqueeze(ftrf{1}([1, 6, 11]))';
        else
            ftrf1611 = [1, 1, 1];
        end
        if all(abs(opts.frame(:) - round(opts.frame(:))) < 0.001) && ...
           all(abs(ftrf1611 .* opts.step - round(ftrf1611 .* opts.step)) < 0.0001) && ...
           (isempty(ftrf) || (all(ftrf{1}([2, 3, 4, 5, 7, 8, 9, 10, 12]) == 0) && ...
             all(ftrf{1}([1, 6, 11, 13, 14, 15]) == round(ftrf{1}([1, 6, 11, 13, 14, 15])))))
            if isempty(ftrf)
                ftrf = {eye(4)};
            end
            opts.frame = repmat(ftrf1611, 2, 1) .* opts.frame + ftrf{1}(1:3, [4, 4])';
            opts.step = opts.step .* ftrf1611;
            s5 = size(odata, 5);
            if any('as' == opts.dir)
                sf = opts.frame;
                sf(:, 1) = ftrf{1}(1) * crd(1) + ftrf{1}(13);
                if s5 > 1
                    sag = double(permute(squeeze(ne_methods.indexarraynb(squeeze(odata), ...
                        [Inf, Inf, Inf, Inf; sf(1, :), 1; opts.step, 1; sf(2, :), s5])), ...
                        [2, 1, 3]));
                else
                    sag = double(squeeze(ne_methods.indexarraynb(squeeze(odata), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)]))');
                end
            end
            sf = opts.frame;
            sf(:, 2) = ftrf{1}(6) * crd(2) + ftrf{1}(14);
            if opts.dir == 'a'
                if s5 > 1
                    cor = double(permute(squeeze(ne_methods.indexarraynb(squeeze(odata), ...
                        [Inf, Inf, Inf, Inf; sf(1, :), 1; opts.step, 1; sf(2, :), s5])), ...
                        [2, 1, 3]));
                else
                    cor = double(squeeze(ne_methods.indexarraynb(squeeze(odata), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)]))');
                end
            elseif opts.dir == 'c'
                if s5 > 1
                    sag = double(permute(squeeze(ne_methods.indexarraynb(squeeze(odata), ...
                        [Inf, Inf, Inf, Inf; sf(1, :), 1; opts.step, 1; sf(2, :), s5])), ...
                        [2, 1, 3]));
                else
                    sag = double(squeeze(ne_methods.indexarraynb(squeeze(odata), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)]))');
                end
            end
            sf = opts.frame;
            sf(:, 3) = ftrf{1}(11) * crd(3) + ftrf{1}(15);
            if opts.dir == 'a'
                if s5 > 1
                    tra = double(permute(squeeze(ne_methods.indexarraynb(squeeze(odata), ...
                        [Inf, Inf, Inf, Inf; sf(1, :), 1; opts.step, 1; sf(2, :), s5])), ...
                        [2, 1, 3]));
                else
                    tra = double(squeeze(ne_methods.indexarraynb(squeeze(odata), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)]))');
                end
            elseif opts.dir == 't'
                if s5 > 1
                    sag = double(permute(squeeze(ne_methods.indexarraynb(squeeze(odata), ...
                        [Inf, Inf, Inf, Inf; sf(1, :), 1; opts.step, 1; sf(2, :), s5])), ...
                        [2, 1, 3]));
                else
                    sag = double(squeeze(ne_methods.indexarraynb(squeeze(odata), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)]))');
                end
            end

        % interpolation anyway
        else
            if any('as' == opts.dir)
                sf = opts.frame;
                sf(:, 1) = crd(1);
                if ndims(odata) == 3
                    sag = squeeze(ne_methods.flexinterpn_method(odata, ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                else
                    sag = squeeze(ne_methods.flexinterpn_method(odata(:, :, :, 1), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                end
                for d5c = 2:size(odata, 5)
                    sag(:, :, d5c) = ...
                        squeeze(ne_methods.flexinterpn_method(odata(:, :, :, 1, d5c), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                end
            end
            sf = opts.frame;
            sf(:, 2) = crd(2);
            if opts.dir == 'a'
                if ndims(odata) == 3
                    cor = squeeze(ne_methods.flexinterpn_method(odata, ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                else
                    cor = squeeze(ne_methods.flexinterpn_method(odata(:, :, :, 1), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                end
                for d5c = 2:size(odata, 5)
                    cor(:, :, d5c) = ...
                        squeeze(ne_methods.flexinterpn_method(odata(:, :, :, 1, d5c), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                end
            elseif opts.dir == 'c'
                if ndims(odata) == 3
                    sag = squeeze(ne_methods.flexinterpn_method(odata, ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                else
                    sag = squeeze(ne_methods.flexinterpn_method(odata(:, :, :, 1), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                end
                for d5c = 2:size(odata, 5)
                    sag(:, :, d5c) = ...
                        squeeze(ne_methods.flexinterpn_method(odata(:, :, :, 1, d5c), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                end
            end
            sf = opts.frame;
            sf(:, 3) = crd(3);
            if opts.dir == 'a'
                if ndims(odata) == 3
                    tra = squeeze(ne_methods.flexinterpn_method(odata, ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                else
                    tra = squeeze(ne_methods.flexinterpn_method(odata(:, :, :, 1), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                end
                for d5c = 2:size(odata, 5)
                    tra(:, :, d5c) = ...
                        squeeze(ne_methods.flexinterpn_method(odata(:, :, :, 1, d5c), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                end
            elseif opts.dir == 't'
                if ndims(odata) == 3
                    sag = squeeze(ne_methods.flexinterpn_method(odata, ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                else
                    sag = squeeze(ne_methods.flexinterpn_method(odata(:, :, :, 1), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                end
                for d5c = 2:size(odata, 5)
                    sag(:, :, d5c) = ...
                        squeeze(ne_methods.flexinterpn_method(odata(:, :, :, 1, d5c), ...
                        [Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                        0, ftrf{:}, opts.method))';
                end
            end
        end
        if ~isempty(ftrf)
            crd(4) = 1;
            crd = ftrf{1} * crd(:);
            sc = crd(1:3);
            sc = sc(:)';
        end

    % with normalization
    else
        if isempty(ftrf)
            ftrf = eye(4);
        else
            ftrf = ftrf{1};
        end
        ivgm = inv(snm.VG(1).mat);
        sf = opts.frame;
        nx = numel(sf(1, 1):opts.step(1):sf(2, 1));
        ny = numel(sf(1, 2):opts.step(2):sf(2, 2));
        nz = numel(sf(1, 3):opts.step(3):sf(2, 3));
        if any('as' == opts.dir)
            sf(:, 1) = crd(1);
            sag = reshape(ne_methods.flexinterpn_method(odata(:, :, :, 1), ...
                ne_methods.applyspmsnc([Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                snm.Tr, snm.VG(1).dim, ivgm, ...
                ftrf * snm.VF.mat * snm.Affine), 0, opts.method), ny, nz)';
            for d5c = 2:size(odata, 5)
                sag(:, :, d5c) = ...
                    reshape(ne_methods.flexinterpn_method(odata(:, :, :, 1, d5c), ...
                    ne_methods.applyspmsnc([Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                    snm.Tr, snm.VG(1).dim, ivgm, ...
                    ftrf * snm.VF.mat * snm.Affine), 0, opts.method), ny, nz)';
            end
            sf = opts.frame;
        end
        if any(opts.dir == 'ac')
            sf(:, 2) = crd(2);
            cor = reshape(ne_methods.flexinterpn_method(odata(:, :, :, 1), ...
                ne_methods.applyspmsnc([Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                snm.Tr, snm.VG(1).dim, ivgm, ...
                ftrf * snm.VF.mat * snm.Affine), 0, opts.method), nx, nz)';
            for d5c = 2:size(odata, 5)
                cor(:, :, d5c) = ...
                    reshape(ne_methods.flexinterpn_method(odata(:, :, :, 1, d5c), ...
                    ne_methods.applyspmsnc([Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                    snm.Tr, snm.VG(1).dim, ivgm, ...
                    ftrf * snm.VF.mat * snm.Affine), 0, opts.method), nx, nz)';
            end
            if opts.dir == 'c'
                sag = cor;
            end
            sf = opts.frame;
        end
        if any(opts.dir == 'at')
            sf(:, 3) = crd(3);
            tra = reshape(ne_methods.flexinterpn_method(odata(:, :, :, 1), ...
                ne_methods.applyspmsnc([Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                snm.Tr, snm.VG(1).dim, ivgm, ...
                ftrf * snm.VF.mat * snm.Affine), 0, opts.method), nx, ny)';
            for d5c = 2:size(odata, 5)
                tra(:, :, d5c) = ...
                    reshape(ne_methods.flexinterpn_method(odata(:, :, :, 1, d5c), ...
                    ne_methods.applyspmsnc([Inf, Inf, Inf; sf(1, :); opts.step; sf(2, :)], ...
                    snm.Tr, snm.VG(1).dim, ivgm, ...
                    ftrf * snm.VF.mat * snm.Affine), 0, opts.method), nx, ny)';
            end
            if opts.dir == 't'
                sag = tra;
            end
        end
        sc = ne_methods.applyspmsnc(crd(:)', snm.Tr, snm.VG(1).dim, ivgm, ...
            ftrf * snm.VF.mat * snm.Affine);
    end
end
if opts.dir ~= 'a'
    cor = sc;
    tra = odata;
    sc = oalpha;
end

% scale/intercept data
if scs ~= 1
    if sci ~= 0
        sag = scs .* sag + sci;
        if opts.dir == 'a'
            cor = scs .* cor + sci;
            tra = scs .* tra + sci;
        end
    else
        sag = scs .* sag;
        if opts.dir == 'a'
            cor = scs .* cor;
            tra = scs .* tra;
        end
    end
elseif sci ~= 0
    sag = sag + sci;
    if opts.dir == 'a'
        cor = cor + sci;
        tra = tra + sci;
    end
end

% compute gradients over data
if opts.gradient

    % override method
    if opts.method(1) == 'n'
        opts.method = 'linear';
    end

    % sampling grid
    smpgx = [[inf; 1; 1] * ones(1, ndims(sag)); size(sag) + 0.4];
    smpgy = smpgx;
    smpgx(2:3, 1) = [0.625; 0.125];
    smpgy(2:3, 2) = [0.625; 0.125];
    gdx = ne_methods.flexinterpn_method(sag, smpgx, opts.method);
    gdy = ne_methods.flexinterpn_method(sag, smpgy, opts.method);
    if opts.gradtype(1) == '1'
        grx = gdx(5:8:end, :, :) - gdx(3:8:end, :, :);
        gry = gdy(:, 5:8:end, :) - gdy(:, 3:8:end, :);
        sag = 4 .* sqrt(grx .* grx + gry .* gry);
    elseif opts.gradtype(1) == '2'
    else
    end
    if opts.dir == 'a'
        smpgx(end, :) = size(cor) + 0.4;
        smpgy(end, :) = size(cor) + 0.4;
        gdx = ne_methods.flexinterpn_method(cor, smpgx, opts.method);
        gdy = ne_methods.flexinterpn_method(cor, smpgy, opts.method);
        if opts.gradtype(1) == '1'
            grx = gdx(5:8:end, :, :) - gdx(3:8:end, :, :);
            gry = gdy(:, 5:8:end, :) - gdy(:, 3:8:end, :);
            cor = 4 .* sqrt(grx .* grx + gry .* gry);
        elseif opts.gradtype(1) == '2'
        else
        end
        smpgx(end, :) = size(tra) + 0.4;
        smpgy(end, :) = size(tra) + 0.4;
        gdx = ne_methods.flexinterpn_method(tra, smpgx, opts.method);
        gdy = ne_methods.flexinterpn_method(tra, smpgy, opts.method);
        if opts.gradtype(1) == '1'
            grx = gdx(5:8:end, :, :) - gdx(3:8:end, :, :);
            gry = gdy(:, 5:8:end, :) - gdy(:, 3:8:end, :);
            tra = 4 .* sqrt(grx .* grx + gry .* gry);
        elseif opts.gradtype(1) == '2'
        else
        end
    end

    % log
    if opts.gradtype(2) == 'l'
        sag = log(1 + sag);
        if opts.dir == 'a'
            cor = log(1 + cor);
            tra = log(1 + tra);
        end
    end
end
