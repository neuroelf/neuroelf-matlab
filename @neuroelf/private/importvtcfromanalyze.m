function vtc = importvtcfromanalyze(imgs, bbox, res, imeth, trans, opts)
% importvtcfromanalyze  - import a VTC from Analzye files
%
% FORMAT:       vtc = importvtcfromanalyze(imgs [, bb [, res [, im [, t [, o]]]]])
%
% Input fields:
%
%       imgs        cell array with HDR filenames for xff
%       bb          optional 2x3 bounding box (default: MNI covering box)
%                   must be given in BrainVoyager's axes order!!
%       res         optional resolution (default: 3)
%       imeth       interpolation 'cubic', 'lanczos3', {'linear'}, 'nearest'
%       t           4x4 transformation matrix (also stored in RunTimeVars)
%       o           additional options
%        .gsmasks   files used to extract subject-specific global signals
%        .gsnames   global signal mask names (default: filenames w/o ext.)
%        .gssnmat   spatial-normalization MAT content for GS masks
%        .mppca     compute motion-parameter based PCA components (12)
%        .raw       store raw matrix and apply transformation
%        .snmat     spatial-normalization MAT content
%
% Output fields:
%
%       vtc         created VTC object
%
% Note: this function requires the MEX file flexinterpn.

% Version:  v1.1
% Build:    16061121
% Date:     Jun-11 2016, 9:34 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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
if nargin < 1 || ~iscell(imgs) || numel(imgs) < 1 || isempty(imgs{1}) || ...
   (~ischar(imgs{1}) && (numel(imgs{1}) ~= 1 || ~isxff(imgs{1}, {'hdr', 'head'})))
    error('neuroelf:general:badArgument', 'Bad input argument.');
end

% go on with loading images and sampling VTC data
imgs = imgs(:);
nimg = numel(imgs);
hclr = true(1, nimg);
himg = cell(1, nimg);
hrtv = repmat({struct}, 1, nimg);
htyp = himg;
vimg = zeros(20, nimg);
nvol = ones(1, nimg);
try
    pbar = xprogress;
    xprogress(pbar, 'settitle', 'Converting Analyze to VTC...');
    xprogress(pbar, 0, 'Checking HDRs...', 'visible', 0, 8 * nimg);
    pgstep = 1 / 150000;
    nextpg = now + pgstep;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    pbar = [];
end

% load data
rps = [];
try
    for ic = 1:nimg
        if isempty(imgs{ic}) || (~ischar(imgs{ic}) && (numel(imgs{ic}) ~= 1 || ...
            ~isxff(imgs{ic}, {'hdr', 'head'})))
            error('BAD_IMAGENAME');
        end
        if ischar(imgs{ic})
            imgs{ic} = strrep(strrep(imgs{ic}, '.img', '.hdr'), '.IMG', '.HDR');
            himg{ic} = xff(imgs{ic});
        else
            himg(ic) = imgs(ic);
            hclr(ic) = false;
            hrtv{ic} = himg{ic}.RunTimeVars;
        end
        if numel(himg{ic}) ~= 1 || ~isxff(himg{ic}, {'hdr', 'head'})
            error('BAD_IMAGECONT');
        end
        htyp{ic} = himg{ic}.Filetype;
        icf = himg{ic}.CoordinateFrame;
        vimg(:, ic) = [icf.Trf(:); icf.Dimensions(:)];
        if strcmpi(htyp{ic}, 'hdr')
            nvol(ic) = size(himg{ic}.VoxelData, 4);
        else
            nvol(ic) = numel(himg{ic}.Brick);
        end
        if ic == 1
            [hpath, hfile] = fileparts(himg{1}.FilenameOnDisk);
            frfile = '';
            if numel(hfile) > 2 && (hfile(1) == 'w' || strcmp(hfile(1:2), 'sw'))
                rfile = regexprep(hfile, '^s?wr?', '');
                if exist([hpath '/rp_' rfile '.txt'], 'file')
                    frfile = [hpath '/rp_' rfile '.txt'];
                elseif exist([hpath '/rp-' rfile '.txt'], 'file')
                    frfile = [hpath '/rp-' rfile '.txt'];
                end
            elseif numel(hfile) > 1 && hfile(1) == 'r'
                if exist([hpath '/rp_' hfile(2:end) '.txt'], 'file')
                    frfile = [hpath '/rp_' hfile(2:end) '.txt'];
                elseif exist([hpath '/rp-' hfile(2:end) '.txt'], 'file')
                    frfile = [hpath '/rp-' hfile(2:end) '.txt'];
                end
            end
            if ~isempty(frfile)
                try
                    rps = load(frfile);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            end
        end
        if ~isempty(pbar) && now >= nextpg
            xprogress(pbar, ic);
            nextpg = now + pgstep;
        end
    end
    if any(any(diff(vimg, 1, 2)))
        warning('neuroelf:general:badArgument', ...
            'Spatial orientation/dimensions mismatch between images.');
    end
catch ne_eo;
    clearxffobjects(himg(hclr));
    if ~isempty(pbar)
        closebar(pbar);
    end
    error('neuroelf:general:badArgument', ...
        'Error loading image %d (%s).', ic, ne_eo.message);
end

% checking other arguments
sfn = himg{1}.FilenameOnDisk;
tbb = false;
if nargin < 2 || ~isa(bbox, 'double') || ~isequal(size(bbox), [2, 3]) || ...
    any(isnan(bbox(:)) | bbox(:) < -128 | bbox(:) > 255)
    bbox = [];
else
    bbox = round(bbox);
    if any(bbox(:) < 0)
        bbox = 128 - bbox([2, 1], [2, 3, 1]);
        tbb = true;
    end
end
if nargin < 3 || ~isa(res, 'double') || numel(res) ~= 1 || ~any((1:15) == res)
    res = 3;
end
if tbb
    bboxe = (1 + bbox(2, :) - bbox(1, :)) ./ res;
    if any(bboxe ~= round(bboxe))
        bbox(2, :) = bbox(1, :) + res .* round(bboxe) - 1;
    end
end
if nargin < 4 || ~ischar(imeth) || ...
   ~any(strcmpi(imeth(:)', {'cubic', 'lanczos3', 'linear', 'nearest'}))
    imeth = 'linear';
else
    imeth = lower(imeth(:)');
end
if nargin < 5 || ~isa(trans, 'double') || ~isequal(size(trans), [4, 4]) || ...
    any(isinf(trans(:)) | isnan(trans(:))) || any(trans(4, :) ~= [0, 0, 0, 1])
    trans = [];
end
if nargin < 6 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'gsmasks') || ~iscell(opts.gsmasks) || ...
   ~all(cellfun(@ischar, opts.gsmasks(:))) || any(cellfun('isempty', opts.gsmasks(:)))
    opts.gsmasks = {};
    gsmi = {};
else
    gsmi = opts.gsmasks;
    try
        for ic = 1:numel(opts.gsmasks)
            gsmi{ic} = xff(gsmi{ic}(:)');
            if ~isxff(gsmi{ic}, 'hdr')
                error('neuroelf:xff:badObjectType', 'Not a HDR/NII file.');
            end
            gsmisz = size(gsmi{ic}.VoxelData);
            if numel(gsmisz) < 3
                error('neuroelf:xff:badObjectType', 'Not a 3D/4D HDR/NII file.');
            end
        end
    catch ne_eo;
        if ic > 1
            clearxffobjects(gsmi(1:ic-1));
        end
        rethrow(ne_eo);
    end
end
if ~isfield(opts, 'gsnames') || ~iscell(opts.gsnames) || ...
   ~all(cellfun(@ischar, opts.gsnames(:))) || any(cellfun('isempty', opts.gsnames(:)))
    if ~isempty(opts.gsmasks)
        [gspaths, opts.gsnames] = mfileparts(opts.gsmasks(:));
    else
        opts.gsnames = {};
    end
end
if ~isfield(opts, 'gssnmat') || isempty(opts.gssnmat) || ...
   (~ischar(opts.gssnmat) && ~isstruct(opts.gssnmat))
    opts.gssnmat = [];
end
if ischar(opts.gssnmat)
    try
        opts.gssnmat = load(opts.gssnmat);
    catch ne_eo;
        clearxffobjects(himg(hclr));
        if ~isempy(gsmi)
            clearxffobjects(gsmi);
        end
        if ~isempty(pbar)
            closebar(pbar);
        end
        rethrow(ne_eo);
    end
end
if ~isfield(opts, 'mppca') || ~isa(opts.mppca, 'double') || numel(opts.mppca) ~= 1
    opts.mppca = 12;
end
if ~isfield(opts, 'raw') || ~islogical(opts.raw) || numel(opts.raw) ~= 1
    opts.raw = false;
end
if ~isfield(opts, 'snmat') || isempty(opts.snmat) || ...
   (~ischar(opts.snmat) && ~isstruct(opts.snmat))
    opts.snmat = [];
end
if ischar(opts.snmat)
    try
        opts.snmat = load(opts.snmat);
    catch ne_eo;
        clearxffobjects(himg(hclr));
        if ~isempy(gsmi)
            clearxffobjects(gsmi);
        end
        if ~isempty(pbar)
            closebar(pbar);
        end
        rethrow(ne_eo);
    end
end
if isstruct(opts.snmat)
    if numel(opts.snmat) ~= 1 || ~isfield(opts.snmat, 'VG') || ~isstruct(opts.snmat.VG) || ...
        isempty(opts.snmat.VG) || ~isfield(opts.snmat.VG, 'dim') || ~isfield(opts.snmat.VG, 'mat') || ...
       ~isfield(opts.snmat, 'VF') || ~isstruct(opts.snmat.VF) || numel(opts.snmat.VF) ~= 1 || ...
       ~isfield(opts.snmat.VF, 'dim') || ~isfield(opts.snmat.VF, 'mat') || ...
       ~isfield(opts.snmat, 'Tr') || ~isa(opts.snmat.Tr, 'double') || ndims(opts.snmat.Tr) ~= 4 || ...
        any(isinf(opts.snmat.Tr(:)) | isnan(opts.snmat.Tr(:))) || ...
       ~isfield(opts.snmat, 'Affine') || ~isa(opts.snmat.Affine, 'double') || ...
       ~isequal(size(opts.snmat.Affine), [4, 4]) || ...
        any(isinf(opts.snmat.Affine(:)) | isnan(opts.snmat.Affine(:))) || ...
        any(opts.snmat.Affine(4, :) ~= [0, 0, 0, 1])
        clearxffobjects(himg(hclr));
        if ~isempy(gsmi)
            clearxffobjects(gsmi);
        end
        if ~isempty(pbar)
            closebar(pbar);
        end
        error('neuroelf:general:badArgument', ...
            'SPM-based SN-mat structure not correctly specified.');
    end
end

% get global setting to figure out DataType/FileVersion
global xffsngl;
dtype = xffsngl.CONF.settings.DataTypes.VTC;

% guess scaling
if dtype == 1
    if strcmpi(htyp{1}, 'hdr')
        if istransio(himg{1}.VoxelData)
            vdm = minmaxmean(himg{1}.VoxelData(:, :, :, :));
        else
            vdm = minmaxmean(himg{1}.VoxelData);
        end
    else
        vdm = minmaxmean(himg{1}.Brick(1).Data(:, :, :));
    end
    vdm = vdm(2);
    if vdm > 16384
        vdf = 16384 / double(vdm);
    else
        vdf = [];
    end
else
    vdf = [];
end

% create VTC
vtc = xff('new:vtc');
vtc.DataType = dtype;
if dtype > 1
    vtc.FileVersion = 3;
end
vtc.NameOfSourceFMR = sfn;
vtc.Resolution = res;

% raw storage
if opts.raw

    % get first volume
    if strcmpi(htyp{1}, 'hdr')
        vtd = shiftdim(single(himg{1}.VoxelData(:, :, :, 1)), -1);
    else
        vtd = shiftdim(single(himg{1}.Brick(1).Data(:, :, :)), -1);
    end
    if dtype == 1
        vtd = uint16(round(vtd));
    end

    % create VTC data
    vtd(sum(nvol), 1, 1, 1) = 0;
    vts = size(vtd);
    vtf = cell(vts(1), 1);
    vtcf = zeros(4, 4, vts(1));

    % storing
    if ~isempty(pbar)
        xprogress(pbar, nimg, 'Importing images...');
        nextpg = now + pgstep;
    end
    icc = 1;
    dsc = [];
    for ic = 1:nimg
        himgf = himg{ic}.FilenameOnDisk;
        if isfield(hrtv{ic}, 'Discard')
            dsc = [dsc(:); hrtv{ic}.Discard(:) + sum(nvol(1:ic-1))];
        end
        for vc = 1:nvol(ic)
            try
                crdf = himg{ic}.CoordinateFrame(vc);
                vtcf(:, :, icc) = crdf.Trf;
                if strcmpi(htyp{ic}, 'hdr')
                    hy = shiftdim(single(himg{ic}.VoxelData(:, :, :, vc)), -1);
                else
                    hy = shiftdim(single(himg{ic}.Brick(vc).Data(:, :, :)), -1);
                end
                if ~isempty(vdf)
                    hy = vdf * hy;
                end
                hy(isinf(hy) | isnan(hy)) = 0;
                if ~isempty(pbar) && now >= nextpg
                    xprogress(pbar, nimg + 7 * ((ic - 1) + vc / nvol(ic)));
                    nextpg = now + pgstep;
                end
                vtd(icc, :, :, :) = hy;
                if nvol(ic) > 1
                    vtf{icc} = sprintf('%s,%d', himgf, vc);
                else
                    vtf{icc} = himgf;
                end
                icc = icc + 1;
            catch ne_eo;
                clearxffobjects(himg(hclr));
                if ~isempy(gsmi)
                    clearxffobjects(gsmi);
                end
                if ~isempty(pbar)
                    closebar(pbar);
                end
                vtc.ClearObject;
                error('neuroelf:internalError:storingFailed', ...
                    'Error storing data of volume %d (%s).', ic, ne_eo.message);
            end
        end
        if hclr(ic)
            himg{ic}.ClearObject;
        end
        himg{ic} = [];
    end

    % store simplified frame
    if all(all(abs(diff(vtcf, 1, 3)) <= sqrt(eps)))
        vtcf = vtcf(:, :, 1);
    elseif isstruct(opts.snmat)
        clearxffobjects(himg(hclr));
        if ~isempty(pbar)
            closebar(pbar);
        end
        vtc.ClearObject;
        if ~isempy(gsmi)
            clearxffobjects(gsmi);
        end
        error('neuroelf:general:invalidCombination', ...
            'SPM-based SN-mat requires unique spatial orientation.');
    end

    % VTC settings
    xyzstart = max(0, 128 - (res / 2) .* vts(2:4));
    vtc.XStart = xyzstart(1);
    vtc.XEnd = xyzstart(1) + res * vts(2);
    vtc.YStart = xyzstart(2);
    vtc.YEnd = xyzstart(2) + res * vts(3);
    vtc.ZStart = xyzstart(3);
    vtc.ZEnd = xyzstart(3) + res * vts(4);

% re-sampling
else

    % try to sample first vol (check)
    try
        [vtd, obox] = himg{1}.SampleBVBox( ...
            struct('BBox', bbox, 'ResXYZ', res), 1, imeth, trans);
        if dtype == 1
            vtd = shiftdim(uint16(round(vtd)), -1);
        else
            vtd = shiftdim(single(vtd), -1);
        end
    catch ne_eo;
        clearxffobjects(himg(hclr));
        if ~isempy(gsmi)
            clearxffobjects(gsmi);
        end
        rethrow(ne_eo);
    end

    % create VTC data
    vtd(sum(nvol), 1, 1, 1) = 0;
    vts = size(vtd);
    vtf = cell(vts(1), 1);

    % sampling
    if ~isempty(pbar)
        xprogress(pbar, nimg, 'Sampling images...');
        nextpg = now + pgstep;
    end
    icc = 1;
    dsc = [];
    for ic = 1:nimg
        himgf = himg{ic}.FilenameOnDisk;
        if isfield(hrtv{ic}, 'Discard')
            dsc = [dsc(:); hrtv{ic}.Discard(:) + sum(nvol(1:ic-1))];
        end
        for vc = 1:nvol(ic)
            try
                hy = shiftdim(himg{ic}.SampleBVBox( ...
                    struct('BBox', bbox, 'ResXYZ', res), vc, imeth, trans), -1);
                if ~isempty(vdf)
                    hy = vdf * hy;
                end
                if ~isempty(pbar) && now >= nextpg
                    xprogress(pbar, nimg + 7 * ((ic - 1) + vc / nvol(ic)));
                    nextpg = now + pgstep;
                end
                vtd(icc, :, :, :) = hy;
                if nvol(ic) > 1
                    vtf{icc} = sprintf('%s,%d', himgf, vc);
                else
                    vtf{icc} = himgf;
                end
                icc = icc + 1;
            catch ne_eo;
                clearxffobjects(himg(hclr));
                if ~isempy(gsmi)
                    clearxffobjects(gsmi);
                end
                if ~isempty(pbar)
                    closebar(pbar);
                end
                vtc.ClearObject;
                error('neuroelf:internalError:samplingFailed', ...
                    'Error sampling data of volume %d (%s).', ic, ne_eo.message);
            end
        end
        if hclr(ic)
            himg{ic}.ClearObject;
        end
        himg{ic} = [];
    end

    % VTC settings
    vtc.XStart = obox.BBox(1, 1);
    vtc.XEnd = obox.BBox(1, 1) + res * vts(2);
    vtc.YStart = obox.BBox(1, 2);
    vtc.YEnd = obox.BBox(1, 2) + res * vts(3);
    vtc.ZStart = obox.BBox(1, 3);
    vtc.ZEnd = obox.BBox(1, 3) + res * vts(4);
end

% store data
vtc.VTCData = vtd;
vtc.NrOfVolumes = size(vtc.VTCData, 1);

% auto-mask
mdata = lsqueeze(mean(vtc.VTCData));
mmask = mdata > (0.5 .* mean(mdata));

% subject specific masks
if ~isempty(gsmi)
    globsigs = repmat(opts.gsnames(:), 1, 3);
    for mc = 1:size(globsigs, 1)
        if opts.raw
            gsmaxv = 0.5 * max(gsmi{mc}.VoxelData(:));
            globsigs{mc, 1} = find(gsmi{mc}.VoxelData(:) >= gsmaxv);
        else
            globsigs{mc, 1} = find(lsqueeze(gsmi{mc}.SampleBVBox(vtc.BoundingBox, ...
                1, 'linear', [], opts.snmat)) >= 0.5);
        end
        globsigs{mc, 2} = ztrans(meannoinfnan(vtc.VTCData(:, globsigs{mc, 1}), 2));
        globsigs{mc, 3} = [globsigs{mc, 3} ',1'];
    end
    vtc.RunTimeVars.GlobSigs = [vtc.RunTimeVars.GlobSigs; globsigs];
    clearxffobjects(gsmi);
end

% default masks
if all([44, 242, 38, 194, 44, 212] == ...
        [vtc.XStart, vtc.XEnd, vtc.YStart, vtc.YEnd, vtc.ZStart, vtc.ZEnd]) && ...
    any([2, 3] == res)
    mskfiles = findfiles(neuroelf_path('masks'), sprintf('*%dmm.msk', res));
    mskfiles = mskfiles(:);
    for mc = 1:size(mskfiles, 1)
        msko = xff(mskfiles{mc, 1});
        [mskpath, mskfname] = fileparts(mskfiles{mc, 1});
        mskfiles{mc, 1} = find(msko.Mask(:));
        mskfiles{mc, 2} = ztrans(meannoinfnan(vtc.VTCData(:, mskfiles{mc, 1}), 2));
        mskfiles{mc, 3} = sprintf('%s,%d', mskfname, 1);
        msko.ClearObject;
    end
else
    mskfiles = cell(0, 3);
end

% add RunTimeVars
vtc.RunTimeVars.AutoSave = true;
vtc.RunTimeVars.DVARS = {find(mmask(:)), ...
    sqrt(mean(diff(psctrans(vtc.VTCData(:, mmask))) .^ 2, 2))};
vtc.RunTimeVars.Discard = dsc(:);
if ~isempty(mskfiles)
    vtc.RunTimeVars.GlobSigs = [vtc.RunTimeVars.GlobSigs; mskfiles];
end
if size(rps, 1) == vtc.NrOfVolumes
    vtc.RunTimeVars.MotionParameters = rps;
    if opts.mppca > 0
        if ~isempty(pbar)
            xprogress(pbar, 7.9 * nimg, 'Computing high-MP-based PCA...');
        end
        drps = diff(rps, 1);
        xrps = [rps, [zeros(1, size(rps, 2)); drps]];
        frps = tempfilter(ztrans(xrps), struct('tempdct', ceil(100000 / vtc.TR)));
        frps = orthvecs(frps);
        wrps = min(1, 4.685 ./ sum(frps .* frps, 2));
        frps = ztrans(frps);
        swrps = (sqrt(wrps) * ones(1, 12)) .* frps;
        rmaps = ((swrps' * swrps) \ swrps') * (diag(sqrt(wrps)) * ztrans(tempfilter( ...
            reshape(vtc.VTCData, vtc.NrOfVolumes, prod(vts(2:end))), ...
            struct('tempdct', ceil(100000 / vtc.TR)))));
        rmaps = sum(rmaps .* rmaps);
        rmskthresh = 0.5;
        rmsk = (rmaps >= rmskthresh);
        while sum(rmsk(:)) < (2 * opts.mppca)
            rmskthresh = 0.75 * rmskthresh;
            if rmskthresh < 0.1
                opts.mppca = 0;
                break;
            end
            rmsk = (rmaps >= rmskthresh);
        end
        if opts.mppca > 0
            nmsk = sum(rmsk);
            rmaps = rmaps(rmsk);
            spdiag = (1:nmsk)';
            pcad = ztrans(vtc.VTCData(:, rmsk)) * sparse(spdiag, spdiag, rmaps, nmsk, nmsk, nmsk);
            pcad = ne_fastica(double(pcad), struct('step', 'pca'));
            vtc.RunTimeVars.MotionParametersPCA = pcad(:, end:-1:(end+1-opts.mppca));
            vtc.RunTimeVars.MotionParametersPCAThresh = rmskthresh;
        end
    end
    vtc.RunTimeVars.MPFD = sum(abs(diff([rps(:, 1:3), 50 .* rps(:, 4:6)])), 2);
end
vtc.RunTimeVars.SourceFiles = vtf;
if isstruct(opts.snmat)
    vtc.RunTimeVars.SPMsn = opts.snmat;
end
if opts.raw
    vtcbb = vtc.BoundingBox;
    for d3c = 1:size(vtcf, 3)
        vtcf(:, :, d3c) = vtcf(:, :, d3c) * vtcbb.QuatT2B;
    end
    vtc.RunTimeVars.TrfPlus = vtcf;
end

% close bar
if ~isempty(pbar)
    closebar(pbar);
end
