function [data, crd, ftrf] = aft_SampleData3D(xo, crd, opts)
% AFT::SampleData3D  - sample spatial data
%
% FORMAT:       [data, crd, ftrf] = obj.SampleData3D(crd [, opts]);
%
% Input fields:
%
%       crd         either a 4x3 range definition or Cx3 coordinates
%       opts        optional settings
%        .fmr       required for AMR/MAP objects
%        .mapvol    required for CMP/VMP objects, defaults to 1
%        .method    resampling, one of {'cubic'}, 'lanczos3', 'linear', 'nearest'
%        .snmat     SPM-compatible normalization structure
%        .space     either of 'bvc', 'bvi', 'bvs', {'tal'}
%        .trans     4x4 transformation matrix, used for volume datatypes
%        .trf       TRF object(s) passed on to samplefmrspace
%        .v16       for VMR objects, if true access available V16 data
%
% Output fields:
%
%       data        either XxYxZ or Cx1 data
%
% TYPES: AMR, AVA, CMP, DDT, DMR, FMR, GLM, HDR, HEAD, MAP, MGH, MSK, NLF, SRF, TVL, VDW, VMP, VMR, VTC
%
% Note: for AMR/DMR/FMR/MAP objects, the trf option must be set with
%       at least the IA/FA (or FA alone) TRF object(s) to work
%
% Using: applyspmsnc, bvcoordconv, findfirst, flexinterpn_method, samplefmrspace.

% Version:  v1.1
% Build:    16031615
% Date:     Mar-16 2016, 3:45 PM EST
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

% neuroelf library
global ne_methods;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true) || ...
   ~isa(crd, 'double') || isempty(crd) || (~isequal(size(crd), [4, 3]) && ...
    (ndims(crd) ~= 2 || size(crd, 2) ~= 3))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'fmr') || numel(opts.fmr) ~= 1 || ...
   (~xffisobject(opts.fmr, true, 'dmr') && ~xffisobject(opts.fmr, true, 'fmr'))
    opts.fmr = [];
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
   ~isfield(opts.snmat, 'VF') || ~isfield(opts.snmat, 'VG') || ~isfield(opts.snmat, 'Tr') || ...
   ~isfield(opts.snmat, 'Affine') || ...
   ~isstruct(opts.snmat.VF) || numel(opts.snmat.VF) ~= 1 || ~isfield(opts.snmat.VF, 'mat') || ...
   ~isa(opts.snmat.VF.mat, 'double') || ~isequal(size(opts.snmat.VF.mat), [4, 4]) || ...
   ~isstruct(opts.snmat.VG) || isempty(opts.snmat.VG) || ~isfield(opts.snmat.VG, 'dim') || ...
   ~isfield(opts.snmat.VG, 'mat') || ...
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
if ~isfield(opts, 'tal') || numel(opts.tal) ~= 1 || ...
   ~xffisobject(opts.tal, true, 'tal')
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
    if numel(opts.trf{tc}) ~= 1 || ...
       (~xffisobject(opts.trf{tc}, true, 'tal') && ~xffisobject(opts.trf{tc}, true, 'trf'))
        opts.trf(tc) = [];
    end
end
opts.trf = opts.trf(:)';
if ~isfield(opts, 'v16') || ~islogical(opts.v16) || numel(opts.v16) ~= 1
    opts.v16 = false;
end
epst = eps ^ 0.75;

% switch on filetype
bc = xo.C;
ft = lower(xo.S.Extensions{1});

% check TrfPlus being different from eye(4)
trfplus = [];
if isfield(bc.RunTimeVars, 'TrfPlus') && isequal([4, 4], size(bc.RunTimeVars.TrfPlus)) && ...
    any(any(bc.RunTimeVars.TrfPlus ~= eye(4)))
    trfplus = inv(bc.RunTimeVars.TrfPlus);
end

% for DMR/FMR/VTC objects, resolve transio, otherwise accept map
switch (ft)
    case 'amr'
        if isempty(opts.fmr)
            error('neuroelf:xff:missingOption', ...
                'Option field .fmr must be set to access AMR in 3D space.');
        end
        resetc = false;
        odata = zeros([size(bc.Slice(1).AMRData), numel(bc.Slice)]);
        for slc = 1:numel(bc.Slice)
            if istransio(bc.Slice(slc).AMRData)
                bc.Slice(slc).AMRData = resolve(bc.Slice(slc).AMRData);
                resetc = true;
            end
            odata(:, :, slc) = bc.Slice(slc).AMRData;
        end
        if resetc
            xo.C = bc;
        end
    case 'ava'
        if bc.ProjectType > 2
            error('neuroelf:xff:notImplemented', ...
                '3D data sampling of SRF-based component Maps not supported.');
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
        resetc = false;
        for slc = 1:numel(bc.Slice)
            if istransio(bc.Slice(slc).STCData)
                bc.Slice(slc).STCData = resolve(bc.Slice(slc).STCData);
                resetc = true;
            end
        end
        if resetc
            xo.C = bc;
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
    case 'hdr'
        if any(bc.ImgDim.DataType == [32, 128, 1536, 1792, 2048, 2304])
            switch (bc.ImgDim.DataType)
                case {32, 1792}
                    nvol = size(bc.VoxelData, 4);
                    opts.mapvol = round(min(nvol, opts.mapvol));
                    odata = complex(bc.VoxelData(:, :, :, opts.mapvol, :), ...
                        bc.VoxelDataComplex(:, :, :, opts.mapvol, :));
                case {128, 2304}
                    nvol = size(bc.VoxelDataRGBA, 4);
                    opts.mapvol = round(min(nvol, opts.mapvol));
                    odata = bc.VoxelDataRGBA(:, :, :, opts.mapvol, :);
                otherwise
                    error('neuroelf:xff:unsupported', ...
                        'Slicing of this datatype not yet supported.');
            end
        else
            nvol = size(bc.VoxelData, 4);
            opts.mapvol = round(min(nvol, opts.mapvol));
            odata = bc.VoxelData(:, :, :, opts.mapvol, :);
        end

        % scaling
        if any([2, 4, 8, 130, 132, 136, 256, 512, 768] == bc.ImgDim.DataType) && ...
           (bc.ImgDim.ScalingIntercept ~= 0 || ...
            all([0, 1] ~= bc.ImgDim.ScalingSlope))
            if bc.ImgDim.ScalingSlope ~= 0
                odata = bc.ImgDim.ScalingIntercept + bc.ImgDim.ScalingSlope .* double(odata);
            else
                odata = bc.ImgDim.ScalingIntercept + double(odata);
            end
        end

        % clustered ?
        if bc.RunTimeVars.Map(opts.mapvol).EnableClusterCheck && ...
            numel(bc.VoxelDataCT) >= mapvol && ...
            isequal(size(odata), size(bc.VoxelDataCT{mapvol}))
            odata = odata .* bc.VoxelDataCT{mapvol};
        end

        % also get transformation matrix and set space to bvc!
        cfr = hdr_CoordinateFrame(xo, opts.mapvol);
        if isempty(opts.trans)
            opts.trans = inv(cfr.Trf);
        else
            opts.trans = inv(cfr.Trf) * opts.trans;
        end
        opts.trans(abs(opts.trans) < epst) = 0;
        opts.space = 'bvc';
    case 'head'
        nvol = numel(bc.Brick);
        opts.mapvol = round(min(nvol, opts.mapvol));
        odata = bc.Brick(opts.mapvol).Data;
        if istransio(odata)
            odata = resolve(odata);
        end

        % clustered ?
        if bc.RunTimeVars.Map(opts.mapvol).EnableClusterCheck && ...
            isequal(size(odata), size(bc.Brick(mapvol).DataCT))
            odata = odata .* bc.Brick(mapvol).DataCT;
        end

        % also get transformation matrix and set space to bvc!
        cfr = head_CoordinateFrame(xo);
        if isempty(opts.trans)
            opts.trans = inv(cfr.Trf);
        else
            opts.trans = inv(cfr.Trf) * opts.trans;
        end
        opts.trans(abs(opts.trans) < epst) = 0;
        opts.space = 'bvc';
    case 'map'
        if isempty(opts.fmr)
            error('neuroelf:xff:missingOption', ...
                'Option field .fmr must be set to access MAP in 3D space.');
        end
        resetc = false;
        odata = zeros([size(bc.Map(1).Data), numel(bc.Map)]);
        for slc = 1:numel(bc.Map)
            if istransio(bc.Map(slc).Data)
                bc.Map(slc).Data = resolve(bc.Map(slc).Data);
                resetc = true;
            end
            odata(:, :, slc) = bc.Map(slc).Data;
        end
        if resetc
            xo.C = bc;
        end
    case 'mgh'
        opts.mapvol = 1;
        odata = bc.MGHData;
        if istransio(odata)
            odata = resolve(odata);
        end

        % also get transformation matrix and set space to bvc!
        cfr = mgh_CoordinateFrame(xo);
        if isempty(opts.trans)
            opts.trans = inv(cfr.Trf);
        else
            opts.trans = inv(cfr.Trf) * opts.trans;
        end
        opts.trans(abs(opts.trans) < epst) = 0;
        opts.space = 'bvc';
    case 'msk'
        if istransio(bc.Mask)
            bc.Mask = resolve(bc.Mask);
            xo.C = bc;
        end
        odata = bc.Mask;
    case 'nlf'
        odata = bc.Data;
        if istransio(odata)
            odata = resolve(odata);
        end

        % scaling
        if bc.ScalingIntercept ~= 0 || all([0, 1] ~= bc.ScalingSlope)
            if bc.ScalingSlope ~= 0
                odata = bc.ScalingIntercept + bc.ScalingSlope .* double(odata);
            else
                odata = bc.ScalingIntercept + double(odata);
            end
        end

        % also get transformation matrix and set space to bvc!
        cfr = nlf_CoordinateFrame(xo);
        if isempty(opts.trans)
            opts.trans = inv(cfr.Trf);
        else
            opts.trans = inv(cfr.Trf) * opts.trans;
        end
        opts.trans(abs(opts.trans) < epst) = 0;
        opts.space = 'bvc';
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
        if bc.Map(opts.mapvol).EnableClusterCheck && ...
            isequal(size(odata), size(bc.Map(opts.mapvol).VMPDataCT))
            odata = odata .* bc.Map(opts.mapvol).VMPDataCT;
        end
    case 'vmr'
        if ~isempty(bc.VMRData16) && opts.v16
            if istransio(bc.VMRData16)
                bc.VMRData16 = resolve(bc.VMRData16);
                xo.C = bc;
            end
            odata = bc.VMRData16;
        else
            if istransio(bc.VMRData)
                bc.VMRData = resolve(bc.VMRData);
                xo.C = bc;
            end
            odata = bc.VMRData;
        end
    case {'vtc'}
        if istransio(bc.VTCData)
            bc.VTCData = resolve(bc.VTCData);
            xo.C = bc;
        end
        nvol = size(bc.VTCData, 1);
        opts.mapvol = round(min(nvol, opts.mapvol));
        odata = squeeze(bc.VTCData(opts.mapvol, :, :, :));
    otherwise
        error('neuroelf:xff:objectTypeUnsupported', ...
            'Unsupported object type for ::SampleData3D method.');
end

% sample FMR based data differently
if any(strcmp(ft, {'amr', 'dmr', 'fmr', 'map'})) || (strcmp(ft, 'cmp') && bc.DocumentType == 0)

    % take care of special coordinates
    if size(crd, 1) == 4 && all(isinf(crd(1, :)))
        crdx = cell(1, 3);
        crds = zeros(1, 3);
        crdx{1} = crd(2, 1):crd(3, 1):crd(4, 1);
        crdx{2} = crd(2, 2):crd(3, 2):crd(4, 2);
        crdx{3} = crd(2, 3):crd(3, 3):crd(4, 3);
        crds(1) = numel(crdx{1});
        crds(2) = numel(crdx{2});
        crds(3) = numel(crdx{3});
        [crdx{1:3}] = ndgrid(crdx{1}, crdx{2}, crdx{3});
        crd = [crdx{1}(:), crdx{2}(:), crdx{3}(:)];
    else
        crds = [size(crd, 1), 1];
    end

    % convert coordinates to BV system (incl. SPM normalization)
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
        error('neuroelf:xff:invalidOption', ...
            'Sampling of DMR/FMR based data must not have bvc space.');
    elseif ~isempty(snm)
        crd = 128 - ne_methods.applyspmsnc(128 - crd, snm.Tr, ...
            snm.VG(1).dim, inv(snm.VG(1).mat), snm.VF.mat * snm.Affine);
    end

    % use AMR/MAP option DMR/FMR
    if ~isempty(opts.fmr)
        xo = opts.fmr;
    end

    % sample data with samplefmrspace function
    ftrf = [opts.trf, opts.tal];
    data = reshape(ne_methods.samplefmrspace(odata, crd, xo, ftrf, opts.method), crds);

% sample "already 3D" formats
else

    % get bounding box for valid objects and get transformation
    if ~strcmp(opts.space, 'bvc')
        bbox = aft_BoundingBox(xo);
        ftrf = {ne_methods.bvcoordconv(zeros(0, 3), [opts.space '2bvc'], bbox)};
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

    % no normalization
    if isempty(snm)

        % manually apply transformation?
        if ~isempty(ftrf) && (size(crd, 1) ~= 4 || ~all(isinf(crd(1, :))))
            crd(:, 4) = 1;
            crd = crd * ftrf{1}';
            crd(:, 4) = [];
            ftrf = {};
        end

        % sample data
        data = ne_methods.flexinterpn_method(odata(:, :, :, 1), crd, 0, ftrf{:}, opts.method);
        if size(odata, 5) > 1
            if numel(data) < 2
                dataci = {};
            else
                dataci = repmat({':'}, 1, ne_methods.findfirst(size(data) > 1, -1));
            end
        end
        for d5c = 2:size(odata, 5)
            data(dataci{:}, d5c) = ne_methods.flexinterpn_method( ...
                odata(:, :, :, 1, d5c), crd, 0, ftrf{:}, opts.method);
        end

    % 3D data with normalization
    else

        % apply transformation to coordinates
        if isempty(ftrf)
            ftrf = eye(4);
        else
            ftrf = ftrf{1};
        end
        ivgm = inv(snm.VG(1).mat);
        if all(isinf(crd(1, :)))
            data = ne_methods.flexinterpn_method(odata(:, :, :, 1), ...
                ne_methods.applyspmsnc(crd, snm.Tr, snm.VG(1).dim, ivgm, ...
                ftrf * snm.VF.mat * snm.Affine), 0, opts.method);
            for d5c = 2:size(odata, 5)
                data(:, d5c) = ne_methods.flexinterpn_method(odata(:, :, :, 1, d5c), ...
                    ne_methods.applyspmsnc(crd, snm.Tr, snm.VG(1).dim, ivgm, ...
                    ftrf * snm.VF.mat * snm.Affine), 0, opts.method);
            end
        else
            crd = ne_methods.applyspmsnc(crd, snm.Tr, snm.VG(1).dim, ivgm, ...
                ftrf * snm.VF.mat * snm.Affine);
            data = ne_methods.flexinterpn_method(odata(:, :, :, 1), ...
                crd, 0, opts.method);
            for d5c = 2:size(odata, 5)
                data(:, d5c) = ne_methods.flexinterpn_method(odata(:, :, :, 1, d5c), ...
                    crd, 0, opts.method);
            end
        end

        % for coordinate range
        if size(crd, 1) == 4 && all(isinf(crd(1, :)))

            % reshape output from vector to range size
            lx = numel(crd(2, 1):crd(3, 1):crd(4, 1));
            ly = numel(crd(2, 2):crd(3, 2):crd(4, 2));
            lz = numel(crd(2, 3):crd(3, 3):crd(4, 3));
            data = reshape(data, [lx, ly, lz, size(data, 2)]);
        end

        % add to ftrf
        if nargout > 1
            if nargout > 2
                ftrf = {ftrf, snm};
            end
        end
    end
end
