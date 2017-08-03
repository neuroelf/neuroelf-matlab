function [astruct, adata] = averagestruct(v, opts)
% averagestruct  - average structural images
%
% FORMAT:       astruct = averagestruct(v [, opts])
%
% Input fields:
%
%       v           cell array with list of structural image filenames
%       opts        optional settings
%        .avgtype   averaging type, {'mean'}, 'mean50', 'median', 'robmean'
%        .bbox      bounding box (default: [-127 .. 128] for X/Y/Z)
%        .clean     clean data (default: true)
%        .cnorm     normalize contrast on modes of histogram (default: true)
%        .ihc       perform inhomogeneity correction on data (default: true)
%        .imeth     interpolation (passed to OBJ::SampleTalBox, 'cubic')
%        .indivvmr  pattern, write individual datasets as VMRs (default: '')
%        .maskc1c2  GM+WM only, valid for HDR objects (default: false)
%        .maskc3    also mask with CSF (requires maskc1c2 to be true, false)
%        .median    median value for voxels (default: 125)
%        .outtype   output object type, either of 'HDR' or {'VMR'}
%        .pbar      progress bar object (if not given, created)
%        .res       resolution, either 0.5 or {1}
%        .snorm     load normalization file for each struct (default: true)
%
% Output fields:
%
%       astruct     average VMR with maximum bounding box
%
% Note: the .indivvmr field is used to replace the original filename
%       extension (e.g. '_norm.vmr')

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:13 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, 2016, Jochen Weber
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
   ~iscell(v) || ...
    numel(v) < 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
v = v(:);
nv = numel(v);
for vc = 1:nv
    if ~ischar(v{vc}) || ...
        isempty(v{vc}) || ...
        exist(v{vc}(:)', 'file') ~= 2
        clearxffobjects(v);
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid list of files.' ...
        );
    end
    try
        v{vc} = xff(v{vc}(:)', 't');
    catch ne_eo;
        clearxffobjects(v);
        rethrow(ne_eo);
    end
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'avgtype') || ...
   ~ischar(opts.avgtype) || ...
    isempty(opts.avgtype) || ...
   ~any(strcmpi(opts.avgtype(:)', {'mean', 'mean50', 'median', 'robmean'}))
    opts.avgtype = 'mean';
else
    opts.avgtype = lower(opts.avgtype(:)');
end
if ~isfield(opts, 'bbox') || ...
   ~isa(opts.bbox, 'double') || ...
   ~isequal(size(opts.bbox), [2, 3]) || ...
    any(isinf(opts.bbox(:)) | isnan(opts.bbox(:)) | opts.bbox(:) < -256 | opts.bbox(:) > 256) || ...
    any(diff(opts.bbox) <= 0)
    opts.bbox = [-127, -127, -127; 128, 128, 128];
else
    opts.bbox(1, :) = floor(opts.bbox(1, :));
    opts.bbox(2, :) = ceil(opts.bbox(2, :));
end
if ~isfield(opts, 'clean') || ...
   ~islogical(opts.clean) || ...
    numel(opts.clean) ~= 1
    opts.clean = true;
end
if ~isfield(opts, 'cnorm') || ...
   ~islogical(opts.cnorm) || ...
    numel(opts.cnorm) ~= 1
    opts.cnorm = true;
end
if ~isfield(opts, 'ihc') || ...
   ~islogical(opts.ihc) || ...
    numel(opts.ihc) ~= 1
    opts.ihc = true;
end
if ~isfield(opts, 'imeth') || ...
   ~ischar(opts.imeth) || ...
    isempty(opts.imeth)
    opts.imeth = 'cubic';
else
    opts.imeth = lower(opts.imeth(:)');
end
if ~isfield(opts, 'indivvmr') || ...
   ~ischar(opts.indivvmr) || ...
    numel(opts.indivvmr) < 4 || ...
    numel(opts.indivvmr) ~= size(opts.indivvmr, 2) || ...
   ~strcmpi(opts.indivvmr(end-3:end), '.vmr')
    opts.indivvmr = '';
end
if ~isfield(opts, 'maskc1c2') || ...
   ~islogical(opts.maskc1c2) || ...
    numel(opts.maskc1c2) ~= 1
    opts.maskc1c2 = false;
end
if ~isfield(opts, 'median') || ...
   ~isa(opts.median, 'double') || ...
    numel(opts.median) > 1 || ...
    any(isinf(opts.median) | isnan(opts.median) | opts.median < 1)
    opts.median = 125;
end
if ~isfield(opts, 'outtype') || ...
   ~ischar(opts.outtype) || ...
    isempty(opts.outtype) || ...
   ~any(strcmpi(opts.outtype(:)', {'hdr', 'vmr'}))
    opts.outtype = 'vmr';
else
    opts.outtype = lower(opts.outtype(:)');
end
if ~isfield(opts, 'pbar') || ...
    numel(opts.pbar) ~= 1 || ...
   (~isxfigure(opts.pbar, true) && ...
    ~isa(opts.pbar, 'xprogress'))
    opts.pbar = [];
end
if ~isfield(opts, 'res') || ...
   ~isa(opts.res, 'double') || ...
    numel(opts.res) ~= 1 || ...
    isinf(opts.res) || ...
    isnan(opts.res) || ...
   ~any(opts.res == [0.5, 1])
    opts.res = 1;
end
if ~isfield(opts, 'snorm') || ...
   ~islogical(opts.snorm) || ...
    numel(opts.snorm) ~= 1
    opts.snorm = true;
end

% check additional files
if opts.maskc1c2
    c1c2 = cell(nv, 2);
    for vc = 1:nv
        if ~isxff(v{vc}, 'hdr')
            clearxffobjects(v);
            error( ...
                'neuroelf:MissingFile', ...
                'GM/WM masking requires all files to be HDR/NII.' ...
            );
        end
        c1 = fnpp(v{vc}.FilenameOnDisk, 'c1');
        c2 = fnpp(v{vc}.FilenameOnDisk, 'c2');
        if exist(c1, 'file') ~= 2 || ...
            exist(c2, 'file') ~= 2
            clearxffobjects(v);
            error( ...
                'neuroelf:MissingFile', ...
                'Missing GM/WM file for file %d.', ...
                vc ...
            );
        end
        c1c2(vc, :) = {c1, c2};
    end
    try
        c1c2 = reshape(xff(c1c2, 't'), nv, 2);
    catch ne_eo;
        clearxffobjects(v);
        rethrow(ne_eo);
    end
else
    c1c2 = {};
end
if opts.maskc1c2 && ...
    opts.maskc3
    c3 = cell(nv, 1);
    for vc = 1:nv
        c3f = fnpp(v{vc}.FilenameOnDisk, 'c3');
        if exist(c3f, 'file') ~= 2
            clearxffobjects(v);
            clearxffobjects(c1c2);
            error( ...
                'neuroelf:MissingFile', ...
                'Missing CSF file for file %d.', ...
                vc ...
            );
        end
        c3{vc} = c3f;
    end
    try
        c3 = reshape(xff(c3, 't'), nv, 1);
    catch ne_eo;
        clearxffobjects(v);
        clearxffobjects(c1c2);
        rethrow(ne_eo);
    end
else
    c3 = {};
end
if opts.snorm
    snorm = cell(nv, 1);
    for vc = 1:nv
        if ~isxff(v{vc}, 'hdr')
            snorm = [];
            break;
        end
        hf = v{vc}.FilenameOnDisk;
        s1 = [hf(1:end-4) '_seg_sn.mat'];
        s2 = [hf(1:end-4) '_sn.mat'];
        if exist(s1, 'file') == 2
            snorm{vc} = load(s1);
        elseif exist(s2, 'file') == 2
            snorm{vc} = load(s2);
        else
            clearxffobjects(v);
            clearxffobjects(c1c2);
            clearxffobjects(c3);
            error( ...
                'neuroelf:MissingFile', ...
                'Normalization file for file %d is missing.', ...
                vc ...
            );
        end
    end
else
    snorm = {};
end

% test xprogress
try
    pbarvis = '';
    if isempty(opts.pbar)
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 200, 640, 36]);
        xprogress(pbar, 'settitle', 'Averaging structurals...');
        xprogress(pbar, 0, 'Averaging structurals...', 'visible', 0, 1);
    else
        pbar = opts.pbar;
        pbarvis = pbar.Visible;
    end
    pbar.Progress(0, 'Averaging structurals...');
    pbar.Visible = 'on';
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    pbar = [];
end

% create target data
if opts.res == 1
    adata = zeros(1 + diff(opts.bbox));
    res = [1, 1, 1];
else
    opts.bbox(1, :) = opts.bbox(1, :) - 0.5;
    adata = zeros(1 + round(2 .* diff(opts.bbox)));
    res = [0.5, 0.5, 0.5];
end

% transio needed for all but straight mean
tiofile = '';
if ~strcmp(opts.avgtype, 'mean')

    % create transio
    tiofile = [tempname '.tio'];
    try
        if ~isempty(pbar)
            pbar.Progress(0, 'Preparing array for non-mean averaging...');
        end
        tioobj = transio(tiofile, 'ieee-le', 'single', 0, [size(adata), nv], true);
    catch ne_eo;
        clearxffobjects(v);
        clearxffobjects(c1c2);
        clearxffobjects(c3);
        error( ...
            'neuroelf:TransIOError', ...
            'Couldn''t create transio for averaging: %s', ...
            ne_eo.message ...
        );
    end
end

% iterate over data
vmr = [];
cvmr = [];
vstruct = [];
try
    for vc = 1:nv

        % filename
        [fp, fn] = fileparts(v{vc}.FilenameOnDisk);

        % normalize data?
        sn = cell(1, 2);
        if ~isempty(snorm)
            sn{2} = snorm{vc};
        end

        % sample data
        if ~isempty(pbar)
            pbar.Progress((vc - 1) / nv, sprintf('Sampling from %s...', fn));
        end
        vdata = v{vc}.SampleTalBox( ...
            struct('BBox', opts.bbox, 'ResXYZ', res), 1, opts.imeth, sn{:});
        v{vc}.ClearObject;
        v{vc} = [];

        % restrict to brain
        if ~isempty(c1c2)
            if ~isempty(pbar)
                pbar.Progress((vc - 1) / nv, sprintf('Brain restricting %s...', fn));
            end
            cdata = single(c1c2{vc, 1}.SampleTalBox( ...
                struct('BBox', opts.bbox, 'ResXYZ', res), 1, 'linear', sn{:})) + ...
                    single(c1c2{vc, 2}.SampleTalBox( ...
                struct('BBox', opts.bbox, 'ResXYZ', res), 1, 'linear', sn{:}));
            c1c2{vc, 1}.ClearObject;
            c1c2{vc, 2}.ClearObject;
            c1c2(vc, :) = {[], []};
            if ~isempty(c3)
                cdata = cdata + single(c3{vc}.SampleTalBox( ...
                    struct('BBox', opts.bbox, 'ResXYZ', res), 1, 'linear', sn{:}));
                c3{vc}.ClearObject;
                c3{vc} = [];
            end
            cdata = (smoothdata3(cdata, [2, 2, 2] ./ opts.res) > 0.5);
            [cls, cdata] = clustercoordsc(~cdata, 1);
            vdata(cdata == maxpos(cls)) = 0;
        end

        % compute preliminary min/max/mean
        v0 = (vdata == 0);
        if opts.clean || ...
            opts.cnorm || ...
            opts.ihc || ...
           ~isempty(opts.indivvmr)
            mm = minmaxmean(vdata(~v0));
        end

        % clean data
        if opts.clean
            if ~isempty(pbar)
                pbar.Progress((vc - 1) / nv, sprintf('Cleaning %s...', fn));
            end
            vmr = xff('new:vmr');
            vmr.VMRData = uint8(min(224, round((35 / mm(3)) .* vdata)));
            vmr.VMRData(v0) = 0;
            cvmr = vmr.CleanVMR;
            vmr.ClearObject;
            v0 = (cvmr.VMRData == 0);
            cvmr.ClearObject;
            cvmr = 0;
            vdata(v0) = 0;
            v0 = (vdata == 0);
        end

        % inhomogeneity correction
        if opts.ihc
            if ~isempty(pbar)
                pbar.Progress((vc - 1) / nv, sprintf('IHC for %s...', fn));
            end
            vmr = xff('new:vmr');
            vmr.VMRData = uint8(round((224 / (mm(2) - mm(1))) .* (vdata - mm(1))));
            vmr.VMRData(v0) = 0;
            vdd = vmr.VMRData;
            vmr.InhomogeneityCorrect;
            vdd = double(vmr.VMRData) ./ double(vdd);
            vmr.ClearObject;
            vmr = 0;
            vdd(isinf(vdd) | isnan(vdd)) = 1;
            vdd(v0 | vdd == 0) = 1;
            vdd = smoothdata3(vdd, [3, 3, 3]);
            vdata = vdata .* vdd;
            vdata(v0) = 0;
            vdd = [];
            v0 = (vdata == 0);

            % recompute min/max/mean
            mm = minmaxmean(vdata(~v0));
        end

        % median
        if ~isempty(opts.median)
            vdata = (opts.median / median(vdata(vdata > 0))) .* vdata;
            mm = minmaxmean(vdata(~v0));
        end
        vdata(vdata < 0) = 0;

        % contrast normalization
        if opts.cnorm || ...
           ~isempty(opts.indivvmr)

            % rescale mean to 100
            if ~isempty(pbar)
                pbar.Progress((vc - 1) / nv, sprintf('Contrast normalizing %s...', fn));
            end
            vdata = (100 / (mm(3) - mm(1))) .* (vdata - mm(1));
            mm(1:3) = (100 / (mm(3) - mm(1))) .* (mm(1:3) - mm(1));
            vdata(v0) = 0;

            % compute histogram (up to 5 * mean)
            mh = histcount(vdata(~v0), 0, 500, 0.5);
            mh(1:175) = 0;
            mh(end-99:end) = 0;
            if mm(2) < 500
                mm(round(2 * mm(2))-10:end) = 0;
            end

            % estimate modes
            mh = conv(mh(:), smoothkern(8, 0.01));
            md = diff(mh);
            m1 = maxpos(mh);
            mh(m1:m1+findfirst(md(m1:end) >= 0)) = 0;
            mh(m1:-1:m1-findfirst(md(m1-1:-1:1) <= 0)) = 0;
            m2 = maxpos(mh);
            if m1 > m2
                m21 = m2;
                m2 = m1;
                m1 = m21;
            end

            % renormalize to 125 and 175
            if opts.cnorm
                if m2 > m1
                    vdata = 125 + (100 / (m2 - m1)) .* (vdata - 0.5 * m1);
                else
                    vdata = 125 + (vdata - 0.5 * m1);
                end
                vdata(v0) = 0;
                vdata(vdata < 0) = 0;
                v0 = (vdata == 0);
            end

            % saving
            if ~isempty(opts.indivvmr)

                % generate object
                vstruct = xff('new:vmr');

                % set data
                if opts.cnorm
                    vstruct.VMRData = uint8(round(min(224, ...
                        permute(vdata(end:-1:1, end:-1:1, end:-1:1), [2, 3, 1]))));
                    vstruct.VMRData16 = uint16(round(min(4095, ...
                        10 .* permute(vdata(end:-1:1, end:-1:1, end:-1:1), [2, 3, 1]))));
                else
                    vstruct.VMRData = uint8(max(0, min(225, ...
                        round(125 + (100 / (m2 - m1)) .* ...
                        (permute(vdata(end:-1:1, end:-1:1, end:-1:1), [2, 3, 1]) - 0.5 * m1)))));
                    vstruct.VMRData(v0) = 0;
                    vstruct.VMRData16 = uint16(max(0, min(4095, ...
                        round((1000 / mm(3)) .* ...
                        permute(vdata(end:-1:1, end:-1:1, end:-1:1), [2, 3, 1])))));
                    vstruct.VMRData16(v0) = 0;
                end

                % set bounding box
                vstruct.OffsetX = round(1 / opts.res) * (128 - opts.bbox(2, 2));
                vstruct.OffsetY = round(1 / opts.res) * (128 - opts.bbox(2, 3));
                vstruct.OffsetZ = round(1 / opts.res) * (128 - opts.bbox(2, 1));
                vstruct.FramingCube = round(256 / opts.res);
                vstruct.RowDirY = opts.res;
                vstruct.ColDirZ = -opts.res;
                vstruct.NRows = vstruct.FramingCube;
                vstruct.NCols = vstruct.FramingCube;
                vstruct.SliceThickness = opts.res;
                vstruct.VoxResX = opts.res;
                vstruct.VoxResY = opts.res;
                vstruct.VoxResZ = opts.res;
                vstruct.VoxResVerified = true;

                % save object
                if any(opts.indivvmr == '/' | opts.indivvmr == '\')
                    [pp, pf, pe] = fileparts(opts.indivvmr);
                    vstruct.SaveAs([pp, '/', fn, pf, pe]);
                else
                    vstruct.SaveAs([fp, '/', fn, opts.indivvmr]);
                end

                % clear object
                vstruct.ClearObject;
                vstruct = [];
            end
        end

        % add data to average
        if ~isempty(pbar)
            pbar.Progress(vc / nv, sprintf('Averaging %s...', fn));
        end
        if strcmp(opts.avgtype, 'mean')
            adata = adata + vdata;
        else
            tioobj(:, :, :, vc) = single(vdata);
        end
    end
catch ne_eo;
    clearxffobjects(v);
    clearxffobjects(c1c2);
    clearxffobjects(c3);
    clearxffobjects({vmr, cvmr, vstruct});
    if ~isempty(pbar) && ...
        isempty(opts.pbar)
        closebar(pbar);
    end
    rethrow(ne_eo);
end

% averaging
if ~isempty(pbar)
    pbar.Progress(1, 'Averaging...');
end
switch (opts.avgtype)

    % plain average
    case {'mean'}
        adata = (1 / nv) .* adata;

    % average from middle 50 percentile
    case {'mean50'}

        % compute 50-percentile indices
        i1 = 1 + floor(0.25 * nv);
        i2 = nv + 1 - ceil(0.25 * nv);

        % iterate over Z slices
        for zc = 1:size(adata, 3)

            % access data
            zdata = squeeze(tioobj(:, :, zc, :));

            % sort in 3rd dim
            zdata = sort(zdata, 3, 'ascend');

            % average
            adata(:, :, zc) = mean(zdata(:, :, i1:i2), 3);
        end

    % median
    case {'median'}

        % iterate over Z slices
        for zc = 1:size(adata, 3)

            % access data
            zdata = squeeze(tioobj(:, :, zc, :));

            % compute median
            adata(:, :, zc) = median(zdata, 3);
        end

    % robust mean
    case {'robmean'}

        % iterate over Z slices
        for zc = 1:size(adata, 3)

            % access data
            zdata = squeeze(tioobj(:, :, zc, :));

            % compute
            adata(:, :, zc) = robustmean(zdata, 3);
        end
end

% remote transio file
if ~isempty(tiofile)
    try
        delete(tiofile);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

% close bar
if ~isempty(pbar) && ...
   ~isempty(pbarvis)
    pbar.Visible = pbarvis;
end
if ~isempty(pbar) && ...
    isempty(opts.pbar)
    closebar(pbar);
end

% VMR
if opts.outtype(1) == 'v'

    % generate object
    astruct = xff('new:vmr');

    % requires adaptation?
    if ~opts.cnorm
        c1 = ceil(0.25 .* size(adata));
        c2 = ceil(0.75 .* size(adata));
        adata = (150 / mean(lsqueeze(adata(c1(1):c2(1), c1(2):c2(2), c1(3):c2(3))))) .* adata;
    end

    % set data
    adata = permute(adata(end:-1:1, end:-1:1, end:-1:1), [2, 3, 1]);
    astruct.VMRData = uint8(round(min(225, adata)));
    astruct.VMRData16 = uint16(round(min(32767, 50 .* adata)));

    % set bounding box
    astruct.OffsetX = round(1 / opts.res) * (128 - opts.bbox(2, 2));
    astruct.OffsetY = round(1 / opts.res) * (128 - opts.bbox(2, 3));
    astruct.OffsetZ = round(1 / opts.res) * (128 - opts.bbox(2, 1));
    astruct.FramingCube = round(256 / opts.res);
    astruct.RowDirY = opts.res;
    astruct.ColDirZ = -opts.res;
    astruct.NRows = astruct.FramingCube;
    astruct.NCols = astruct.FramingCube;
    astruct.SliceThickness = opts.res;
    astruct.VoxResX = opts.res;
    astruct.VoxResY = opts.res;
    astruct.VoxResZ = opts.res;
    astruct.VoxResVerified = true;

% NII
else

    % generate object
    astruct = xff('new:hdr');
    astruct.ImgDim.Dim(2:5) = [size(adata), 1];
    astruct.ImgDim.DataType = 16;
    astruct.ImgDim.BitsPerPixel = 32;
    astruct.ImgDim.PixSpacing(1:5) = [1, res, 1];
    astruct.DataHist.Description = sprintf('Averaged struct of %d datasets', nv);
    astruct.DataHist.NIftI1.QFormCode = 2;
    astruct.DataHist.NIftI1.SFormCode = 2;
    astruct.DataHist.NIftI1.QuaternionB = 0;
    astruct.DataHist.NIftI1.QuaternionC = 0;
    astruct.DataHist.NIftI1.QuaternionD = 0;
    astruct.DataHist.NIftI1.QuatOffsetX = opts.bbox(1, 1) - opts.res;
    astruct.DataHist.NIftI1.QuatOffsetY = opts.bbox(1, 2) - opts.res;
    astruct.DataHist.NIftI1.QuatOffsetZ = opts.bbox(1, 3) - opts.res;
    astruct.DataHist.NIftI1.AffineTransX = [opts.res, 0, 0, opts.bbox(1, 1) - opts.res];
    astruct.DataHist.NIftI1.AffineTransY = [0, opts.res, 0, opts.bbox(1, 2) - opts.res];
    astruct.DataHist.NIftI1.AffineTransZ = [0, 0, opts.res, opts.bbox(1, 3) - opts.res];
    astruct.VoxelData = single(adata);
end


% sub function
function f = fnpp(f, p)
[fp{1:3}] = fileparts(f);
if isempty(fp{1})
    fp{1} = pwd;
end
f = [fp{1}, filesep, p, fp{2}, fp{3}];
