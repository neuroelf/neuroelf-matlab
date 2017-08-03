function [cs, ctab, vmp, voi] = aft_ClusterTable(xo, mapno, threshold, opts)
% VMP::ClusterTable  - generate a table with clusters
%
% FORMAT:       [c, t, v, vo] = obj.ClusterTable(mapno [, thresh [, opts]])
%
% Input fields:
%
%       mapno       map number (1 .. NrOfMaps)
%       thresh      either p/r values (0 .. 1) or t/F value (1 .. Inf)
%                   if not given or 0, uses the LowerThreshold of map
%       opts        optional settings
%        .altmaps   alternative maps to extract values from (default: [])
%        .altstat   either of 'mean' or {'peak'}
%        .cclag     flag, interpret the threshold as lag number (false)
%        .clconn    cluster connectivity ('face', {'edge'}, 'vertex')
%        .icbm2tal  flag, VMP coords are passed to icbm2tal (default: false)
%        .localmax  break down larger clusters threshold (default: Inf)
%        .localmaxi iterate on sub-clusters (default: false)
%        .localmin  minimum size for sub-clusters (default: 2)
%        .localmsz  print sub-cluster sizes (default: false)
%        .lupcrd    look-up coordinate, either or 'center', 'cog', {'peak'}
%        .minsize   minimum cluster size (map resolution, default by map)
%        .mni2tal   flag, VMP coords are passed to mni2tal (default: false)
%        .showneg   flag, negative values are considered (default: false)
%        .showpos   flag, positive values are considered (default: true)
%        .sorting   either of 'maxstat', {'maxstats'}, 'size', 'x', 'y', 'z'
%        .svc       small-volume correction (either VOI, voxels, or image)
%        .svcdilate dilate SVC mask for visual interpolation (default true)
%        .svcthresh threshold for SVC (default: 0.05)
%        .tdclient  flag, lookup closest talairach label (default false)
%
% Output fields:
%
%       c           1xC struct with properties of clusters
%       t           text table output
%       v           thresholded map (in VMP resolution!)
%       vo          if requested, VOI structure with (TAL) coords
%
% TYPES: CMP, HDR, HEAD, MSK, VMP
%
% Note: if only one output is requested, the table is text table is
%       returned!!
%
% Note: icbm2tal overrides mni2tal!
%
% Note: svc only valid for VMP (given FWHM estimate in RunTimeVars!)
%
% Using: bvcoordconv, clustervol, correlinvtstat, correlpvalue, dilate3d,
%        flexinterpn_method, limitrangec, lsqueeze, minmaxmean, sdist,
%        smoothdata3, smoothkern.

% Version:  v1.1
% Build:    16021413
% Date:     Feb-14 2016, 1:09 PM EST
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

% use neuroelf library
global ne_methods;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'cmp', 'hdr', 'head', 'msk', 'vmp'}) || ...
   ~isa(mapno, 'double') || numel(mapno) ~= 1 || isinf(mapno) || isnan(mapno) || mapno < 1
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
stype = lower(xo.S.Extensions{1});
if strcmp(stype, 'msk')
    mapno = 1;
end
if (any(strcmp(stype, {'cmp', 'vmp'})) && mapno > numel(bc.Map)) || ...
   (any(strcmp(stype, {'hdr', 'head'})) && mapno > numel(bc.RunTimeVars.Map))
    error('neuroelf:xff:badArgument', 'Map number out of bounds.');
end

% get bounding box
if any(strcmp(stype, {'cmp', 'msk', 'vmp'}))
    bb = aft_BoundingBox(xo);
else
    if strcmp(stype, 'hdr')
        cfr = hdr_CoordinateFrame(xo);
    else
        cfr = head_CoordinateFrame(xo);
    end
    bb = struct('QuatB2T', cfr.Trf, 'ResXYZ', cfr.Resolution);
end

% get map structure and data
mapno = fix(mapno);
switch (stype)
    case 'cmp'
        nmap = numel(bc.Map);
        map = bc.Map(mapno);
        vmp = double(map.CMPData);
    case 'hdr'
        if any(bc.ImgDim.DataType == [32, 128, 1792, 2048, 2304])
            error('neuroelf:xff:unsupported', 'Cluster tables not supported for complex datatypes.');
        end
        nmap = numel(bc.RunTimeVars.Map);
        map = bc.RunTimeVars.Map(mapno);
        vmp = double(bc.VoxelData(:, :, :, mapno));
        if bc.ImgDim.ScalingSlope ~= 1
            vmp = vmp .* bc.ImgDim.ScalingSlope;
        end
        if bc.ImgDim.ScalingIntercept ~= 0
            vmp = vmp + bc.ImgDim.ScalingIntercept;
        end
    case 'head'
        nmap = numel(bc.RunTimeVars.Map);
        map = bc.RunTimeVars.Map(mapno);
        vmp = double(bc.Brick(mapno).Data);
    case 'msk'
        [mskpath, mskname] = fileparts(xo.F);
        if isempty(mskname)
            mskname = '<unsaved.msk>';
        else
            mskname = [mskname '.msk'];
        end
        nmap = 1;
        map = struct('Name', mskname, 'Type', 15, 'EnableClusterCheck', 0, ...
            'ClusterSize', 1, 'ShowPositiveNegativeFlag', 1, 'LowerThreshold', 0.01);
        vmp = double(bc.Mask(:, :, :) ~= 0);
        vmp = ne_methods.smoothdata3(vmp, 3.5) .* vmp;
    case 'vmp'
        nmap = numel(bc.Map);
        map = bc.Map(mapno);
        vmp = double(map.VMPData);
end

% no RunTimeVars.FWHMResEst/FWHMResImg in map -> no SVC!
if ~isfield(map, 'RunTimeVars') || ~isstruct(map.RunTimeVars) || numel(map.RunTimeVars) ~= 1 || ...
   ~isfield(map.RunTimeVars, 'FWHMResEst') || numel(map.RunTimeVars.FWHMResEst) ~= 3 || ...
   ~isfield(map.RunTimeVars, 'FWHMResImg') || ...
   (~isa(map.RunTimeVars.FWHMResImg, 'double') && ~isa(map.RunTimeVars.FWHMResImg, 'single')) || ...
    ndims(map.RunTimeVars.FWHMResImg) ~= 3
    if nargin > 3 && isstruct(opts) && numel(opts) == 1 && isfield(opts, 'svc') && ~isempty(opts.svc)
        disp('NeuroElf''s SVC requires a smoothness estimate in vmp.Map(M).RunTimeVars.FWHMResImg.');
        opts.svc = [];
    end
end

% check threshold argument
if nargin < 3 || ~isa(threshold, 'double') || numel(threshold) ~= 1 || ...
    isinf(threshold) || isnan(threshold)
    threshold = 0;
end

% options check
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'altmaps') || ~isa(opts.altmaps, 'double') || isempty(opts.altmaps) || ...
    any(isinf(opts.altmaps(:)) | isnan(opts.altmaps(:)) | opts.altmaps(:) < 1 | opts.altmaps(:) > nmap)
    opts.altmaps = [];
else
    altmaps = unique(round(opts.altmaps(:)));
    opts.altmapsp = true(1, numel(altmaps));
    opts.altmaps = double(zeros([size(vmp), numel(altmaps)]));
    for amc = 1:numel(altmaps)
        switch (stype)
            case 'cmp'
                vmpalt = double(bc.Map(altmaps(amc)).CMPData);
            case 'hdr'
                vmpalt = double(bc.VoxelData(:, :, :, altmaps(amc)));
            case 'head'
                vmpalt = double(bc.Brick(altmaps(amc)).Data);
            case 'vmp'
                vmpalt = double(bc.Map(altmaps(amc)).VMPData);
        end
        opts.altmaps(:, :, :, amc) = vmpalt;
    end
end
if ~isfield(opts, 'altstat') || ~ischar(opts.altstat) || isempty(opts.altstat) || ...
   ~any(lower(opts.altstat(1)) == 'mp')
    opts.altstat = 'peak';
else
    opts.altstat = lower(opts.altstat(1));
end
if ~isfield(opts, 'cclag') || ~islogical(opts.cclag) || numel(opts.cclag) ~=1
    opts.cclag = false;
end
if ~isfield(opts, 'clconn') || ~ischar(opts.clconn) || isempty(opts.clconn) || ...
   ~any(strcmpi(opts.clconn(:)', {'edge', 'face', 'vertex'}))
    opts.clconn = 'edge';
else
    opts.clconn = lower(opts.clconn(:)');
end
if ~isfield(opts, 'icbm2tal') || numel(opts.icbm2tal) ~= 1 || ...
   (~isa(opts.icbm2tal, 'double') && ~islogical(opts.icbm2tal))
    opts.icbm2tal = false;
else
    opts.icbm2tal = true && opts.icbm2tal;
end
if ~isfield(opts, 'localmax') || ~isa(opts.localmax, 'double') || numel(opts.localmax) ~= 1 || ...
    isnan(opts.localmax) || opts.localmax < 3
    opts.localmax = Inf;
else
    opts.localmax = round(opts.localmax);
end
if ~isfield(opts, 'localmaxi') || ~islogical(opts.localmaxi) || numel(opts.localmaxi) ~= 1
    opts.localmaxi = false;
end
if ~isfield(opts, 'localmin') || ~isa(opts.localmin, 'double') || numel(opts.localmin) ~= 1 || ...
    isinf(opts.localmin) || isnan(opts.localmin) || opts.localmin < 1
    opts.localmin = 2;
end
if ~isfield(opts, 'localmsz') || ~islogical(opts.localmsz) || numel(opts.localmsz) ~= 1
    opts.localmsz = false;
end
if ~isfield(opts, 'lupcrd') || ~ischar(opts.lupcrd) || ~any(strcmpi(opts.lupcrd(:)', {'center', 'cog', 'peak'}))
    opts.lupcrd = 'peak';
else
    opts.lupcrd = lower(opts.lupcrd(:)');
end
if ~isfield(opts, 'minsize') || ~isa(opts.minsize, 'double') || numel(opts.minsize) ~= 1 || ...
    isinf(opts.minsize) || isnan(opts.minsize) || opts.minsize < 1 || opts.minsize > (numel(vmp) / 2)
    if map.EnableClusterCheck
        opts.minsize = map.ClusterSize;
    else
        opts.minsize = 1;
    end
else
    opts.minsize = fix(opts.minsize);
end
if ~isfield(opts, 'mni2tal') || numel(opts.mni2tal) ~= 1 || ...
   (~isa(opts.mni2tal, 'double') && ~islogical(opts.mni2tal)) || opts.icbm2tal
    opts.mni2tal = false;
else
    opts.mni2tal = true && opts.mni2tal;
end
if ~isfield(opts, 'showneg') || ~islogical(opts.showneg) || numel(opts.showneg) ~= 1
    opts.negative = (map.ShowPositiveNegativeFlag > 1);
else
    opts.negative = opts.showneg;
end
if ~isfield(opts, 'showpos') || ~islogical(opts.showpos) || numel(opts.showpos) ~= 1
    opts.positive = (mod(map.ShowPositiveNegativeFlag, 2) == 1);
else
    opts.positive = opts.showpos;
end
if ~isfield(opts, 'sorting') || ~ischar(opts.sorting) || isempty(opts.sorting) || ...
   ~any(strcmpi(opts.sorting, {'maxstat', 'maxstats', 'size', 'x', 'y', 'z'}))
    opts.sorting = 'maxstats';
else
    opts.sorting = lower(opts.sorting(:)');
end
if numel(opts.sorting) == 1
    switch (opts.sorting)
        case 'x'
            opts.sorting = 'z';
        case 'y'
            opts.sorting = 'x';
        case 'z'
            opts.sorting = 'y';
    end
end
if ~isfield(opts, 'svc') || ((numel(opts.svc) ~= 1 || ...
     ~xffisobject(opts.svc, true, {'hdr', 'head', 'msk', 'vmr', 'voi', 'vtc'})) && ...
    (~isa(opts.svc, 'double') || isempty(opts.svc) || ndims(opts.svc) > 2 || size(opts.svc, 2) ~= 3 || ...
     any(isinf(opts.svc(:)) | isnan(opts.svc(:)) | opts.svc(:) < -128 | opts.svc(:) > 256)))
    opts.svc = [];
end
if ~isfield(opts, 'svcdilate') || numel(opts.svcdilate) ~= 1 || ~islogical(opts.svcdilate)
    opts.svcdilate = true;
end
if ~isfield(opts, 'svcthresh') || numel(opts.svcthresh) ~= 1 || ~isa(opts.svcthresh, 'double') || ...
    isinf(opts.svcthresh) || isnan(opts.svcthresh) || opts.svcthresh <= 0
    opts.svcthresh = 0.05;
else
    opts.svcthresh = min(0.1, opts.svcthresh);
end
if ~isfield(opts, 'tdclient') || ~islogical(opts.tdclient) || numel(opts.tdclient) ~= 1
    opts.tdclient = false;
end
if ~isfield(opts, 'tpnoconv') || ~islogical(opts.tpnoconv) || numel(opts.tpnoconv) ~= 1
    opts.tpnoconv = false;
end
if threshold == 0
    threshold = map.LowerThreshold;
    opts.tpnoconv = true;
end

% some initial checks on threshold
if threshold < 0
    if nargin < 4
        opts.negative = true;
        opts.positive = false;
    end
    threshold = -threshold;
end

% remove invalid entries first
vmp(isinf(vmp) | isnan(vmp)) = 0;

% put some additional fields in opts
if ~isfield(bc.RunTimeVars, 'TrfPlus') || ~isequal(size(bc.RunTimeVars.TrfPlus), [4, 4]) || ...
    all(all(bc.RunTimeVars.TrfPlus == eye(4)))
    opts.mat = bb.QuatB2T;
else
    opts.mat = bc.RunTimeVars.TrfPlus * bb.QuatB2T;
end
thdegf = 'n/a';
switch (map.Type)
    case 1
        thtype = 't-Map';
        thdegf = sprintf('%d', map.DF1);
        opts.tptype = 't';
        opts.tptypedf = map.DF1;
    case 2
        thtype = 'correlation Map';
        opts.tptype = 'r';
        if stype(1) == 'v' && bc.FileVersion < 5
            thdegf = sprintf('%d (:= t[%d])', map.DF1, map.DF1 - 2);
            opts.tptypedf = map.DF1 - 2;
        else
            thdegf = sprintf('%d', map.DF1);
            opts.tptypedf = map.DF1;
        end
    case 3
        thtype = 'CC Map';
        thdegf = sprintf('%d', map.DF1);

        % in this case also patch the VMP data, and set fields
        if ~opts.cclag
            vmp = mod(vmp + 2, 1000) - 2;
            opts.tptype = 'r';
            if stype(1) == 'v' && bc.FileVersion < 5
                thdegf = sprintf('%d (:= t[%d])', map.DF1, map.DF1 - 2);
                opts.tptypedf = map.DF1 - 2;
            else
                thdegf = sprintf('%d', map.DF1);
                opts.tptypedf = map.DF1;
            end
        else
            opts.tptype = 'l';
            vmp = floor(0.001 * (vmp + 2));
        end

        % and disallow SVC
        opts.svc = [];
    case 4
        thtype = 'F-Map';
        thdegf = sprintf('%d, %d', map.DF1, map.DF1);
        opts.tptype = 'f';
        opts.tptypedf = [map.DF1, map.DF2];
    case 11
        thtype = 'PSC Map';
        opts.svc = [];
    case 12
        thtype = 'z(ica) Map';
        thdegf = sprintf('%d', map.DF1);
    case 13
        thtype = 'Thickness Map';
        opts.svc = [];
    case 15
        thtype = 'beta Map';
        opts.svc = [];
    case 16
        thtype = 'Probability Map';
        thdegf = sprintf('%d', map.DF1);
    case 20
        thtype = 'Mean Diffusivity Map';
        opts.svc = [];
    case 21
        thtype = 'Fractional Anisotropy Map';
        opts.svc = [];
    case 99
        thtype = 'Mask';
        opts.svc = [];
    otherwise
        thtype = 'unknown';
        opts.svc = [];
end
mthreshp = nan;
if any([1, 2, 3, 4] == map.Type)
    if threshold <= 0.1 && ~opts.tpnoconv
        mthreshp = threshold;
        switch (opts.tptype)
            case 'f'
                threshold = ne_methods.sdist('finv', threshold, ...
                    opts.tptypedf(1), opts.tptypedf(2), true);
            case 'r'
                if opts.negative && opts.positive
                    threshold = ne_methods.correlinvtstat(-ne_methods.sdist('tinv', 0.5 * threshold, ...
                        opts.tptypedf), opts.tptypedf + 2);
                else
                    threshold = ne_methods.correlinvtstat(-ne_methods.sdist('tinv', threshold, ...
                        opts.tptypedf), opts.tptypedf + 2);
                end
            case 't'
                if opts.negative && opts.positive
                    threshold = -ne_methods.sdist('tinv', 0.5 * threshold, opts.tptypedf);
                else
                    threshold = -ne_methods.sdist('tinv', threshold, opts.tptypedf);
                end
        end
    else
        switch (opts.tptype)
            case 'f'
                mthreshp = 1 - ne_methods.sdist('fcdf', threshold, ...
                    opts.tptypedf(1), opts.tptypedf(2));
            case 'r'
                if opts.negative && opts.positive
                    mthreshp = ne_methods.correlpvalue(threshold, opts.tptypedf + 2);
                else
                    mthreshp = 0.5 * ne_methods.correlpvalue(threshold, opts.tptypedf + 2);
                end
            case 't'
                if opts.negative && opts.positive
                    mthreshp = 2 - 2 * ne_methods.sdist('tcdf', threshold, opts.tptypedf);
                else
                    mthreshp = 1 - ne_methods.sdist('tcdf', threshold, opts.tptypedf);
                end
        end
    end
end

% SVC
if ~isempty(opts.svc)

    % object
    if numel(opts.svc) == 1

        % voi
        if xffisobject(opts.svc, true, 'voi')

            % collect all coordinates
            voi = opts.svc.C;
            opts.svc = voi.VOI;
            opts.svc = {opts.svc.Voxels};
            opts.svc = cat(1, opts.svc{:});
            if strcmpi(voi.ReferenceSpace, 'bvs')
                opts.svc = 128 - opts.svc;
            elseif strcmpi(voi.ReferenceSpace, 'bvi')
                opts.svc = 128 - opts.svc(:, [3, 1, 2]);
            end

        % BV image
        elseif xffisobject(opts.svc, true, {'msk', 'vmr', 'vtc'})

            % sample BVBox (first volume)
            opts.svc = ne_methods.bvcoordconv(find(ne_methods.lsqueeze( ...
                aft_SampleBVBox(opts.svc, bb)) > 0), 'bvx2tal', bb);

        % HDR/HEAD image
        else

            % get volume
            if xffisobject(opts.svc, true, 'hdr')
                svctrf = hdr_CoordinateFrame(opts.svc, 1);
            else
                svctrf = head_CoordinateFrame(opts.svc, 1);
            end
            mskcoord = (aft_GetVolume(opts.svc, 1) > 0);
            [mskcoord, mskcrdy, mskcrdz] = ind2sub(size(mskcoord), ...
                find(ne_methods.lsqueeze(mskcoord)));
            opts.svc = (svctrf.Trf * ...
                [mskcoord(:)'; mskcrdy(:)'; mskcrdz(:)'; ones(1, numel(mskcrdy))])';
            opts.svc(:, 4) = [];
        end
    end

    % coordinate list by now, if empty
    if isempty(opts.svc)

        % show error
        error('neuroelf:xff:badArgument', 'Empty Small-Volume, no correction possible.');
    end

    % augment with estimate
    svccoord = (bb.QuatB2T \ ([opts.svc, ones(size(opts.svc, 1), 1)])')';
    svcvals = 0.5 .* ne_methods.flexinterpn_method( ...
        map.RunTimeVars.FWHMResImg, svccoord(:, 1:3), 'linear');
    svcvals = ne_methods.limitrangec(svcvals, 0, 3 * median(svcvals(svcvals > 0)), 0);

    % mask results to SVC
    svcnvx = size(svccoord, 1);
    svcvox = zeros(27 * svcnvx, 3);
    svcvtc = 1;
    for svcx = -.3:.3:.3
        for svcy = -.3:.3:.3
            for svcz = -.3:.3:.3
                svcvox(svcvtc:svcvtc+svcnvx-1, :) = svccoord(:, 1:3) + ...
                    ones(svcnvx, 1) * [svcx, svcy, svcz];
                svcvtc = svcvtc + svcnvx;
            end
        end
    end
    szrmap = size(map.RunTimeVars.FWHMResImg);
    svcvox = round(svcvox);
    svcvox(any(svcvox < 1, 2) | svcvox(:, 1) > szrmap(1) | svcvox(:, 2) > szrmap(2) | svcvox(:, 3) > szrmap(3), :) = [];
    svcvox = sub2ind(szrmap, svcvox(:, 1), svcvox(:, 2), svcvox(:, 3));
    svcvox = unique(svcvox);
    svcvmsk = true(size(map.RunTimeVars.FWHMResImg));
    svcvmsk(svcvox) = false;
    if opts.svcdilate
        svcvmsk(ne_methods.dilate3d(~svcvmsk)) = false;
    end
    switch (stype)
        case 'cmp'
            map.CMPData(svcvmsk) = 0;
            vmp = double(map.CMPData);
        case 'hdr'
            if istransio(bc.VoxelData)
                bc.VoxelData = resolve(bc.VoxelData);
            end
            svcvoxd = bc.VoxelData(:, :, :, mapno);
            svcvoxd(svcvmsk) = -(bc.ImgDim.ScalingIntercept / bc.ImgDim.ScalingSlope);
            bc.VoxelData(:, :, :, mapno) = svcvoxd;
            vmp = double(bc.VoxelData(:, :, :, mapno));
            if bc.ImgDim.ScalingSlope ~= 1
                vmp = vmp .* bc.ImgDim.ScalingSlope;
            end
            if bc.ImgDim.ScalingIntercept ~= 0
                vmp = vmp + bc.ImgDim.ScalingIntercept;
            end
        case 'head'
            svcvoxd = bc.Brick(mapno).Data(:, :, :);
            svcvoxd(svcvmsk) = 0;
            bc.Brick(mapno).Data = svcvoxd;
            vmp = double(bc.Brick(mapno).Data);
        case 'vmp'
            map.VMPData(svcvmsk) = 0;
            vmp = double(map.VMPData);
    end

    % max dist
    svcmd = ceil(max(svcvals));

    % coords
    vcf = 1 + svcmd - min(ne_methods.lsqueeze(opts.svc(:, 1:3)));
    svccoord = round(opts.svc(:, 1:3) + vcf);

    % generate volume
    szvol = max(svccoord) + svcmd;
    res2m = false;
    if prod(szvol) > 2e6
        svcmd = ceil(max(svcvals));
        opts.svc = 0.5 .* (opts.svc + 1);
        svcvals = 0.5 .* svcvals;
        vcf = 1 + svcmd - min(ne_methods.lsqueeze(opts.svc(:, 1:3)));
        svccoord = round(opts.svc(:, 1:3) + vcf);
        szvol = max(svccoord) + svcmd;
        res2m = true;
    end
    volume = false(szvol);

    % and coordinate list
    [scx, scy, scz] = ndgrid(-svcmd:svcmd, -svcmd:svcmd, -svcmd:svcmd);
    scxyz = [scx(:), scy(:), scz(:)];
    scd = sqrt(sum(scxyz .* scxyz, 2));
    scxyz = scx(:) + szvol(1) .* scy(:) + (szvol(1) * szvol(2)) .* scz(:);
    svccoord = unique(sub2ind(szvol, svccoord(:, 1), svccoord(:, 2), svccoord(:, 3)));
    scxyz(scd > svcmd) = [];
    scd(scd > svcmd) = [];

    % for every coordinate, set in volume
    for cc = 1:numel(svccoord)
        volume(svccoord(cc) + scxyz(scd <= svcvals(cc))) = true;
    end

    % sample volume
    [scx, scy, scz] = ind2sub(szvol, find(volume(:)));
    if res2m
        scx = 2 .* (scx(:)' - vcf) + 1;
        scy = 2 .* (scy(:)' - vcf) + 1;
        scz = 2 .* (scz(:)' - vcf) + 1;
    else
        scx = scx(:)' - vcf;
        scy = scy(:)' - vcf;
        scz = scz(:)' - vcf;
    end
    svccoord = (bb.QuatB2T \ [scx; scy; scz; ones(1, numel(scx))])';

    % resample values
    svcvals = ne_methods.flexinterpn_method(map.RunTimeVars.FWHMResImg, svccoord(:, 1:3), 'linear');

    % remove bad entries
    svcvals(svcvals <= 1) = [];
    svcvals = ne_methods.limitrangec(svcvals, 0, 5 * median(svcvals), 0);

    % compute the central kernel weights for range
    svmmm = ne_methods.minmaxmean(svcvals);
    svmmf = (1 / 100) * (svmmm(2) - svmmm(1));
    svr = svmmm(1):svmmf:svmmm(2);
    svcvals = 1 + (svcvals - svmmm(1)) ./ svmmf;
    svrm = zeros(101, 1);
    for cc = 1:101
        svrm(cc) = max(ne_methods.smoothkern(svr(cc), 0, true, 'lanczos3'));
    end
    svrm = svrm .^ 3;
    svcvals = ne_methods.flexinterpn_method(svrm, svcvals(:), 'linear');

    % store as upper limit
    if res2m
        opts.svc = ceil(8 .* sum(svcvals));
    else
        opts.svc = ceil(sum(svcvals));
    end
    map.RunTimeVars.SVCResels = opts.svc;
    map.Name = sprintf('SVC (%d resels, p<%g): %s', opts.svc, opts.svcthresh, map.Name);

    % update threshold
    switch (map.Type)

        % t map (and r via t)
        case {1, 2}

            % two-tailed
            if opts.negative && opts.positive

                % compute threshold
                threshold = abs(ne_methods.sdist('tinv', 0.5 * opts.svcthresh / opts.svc, map.DF1));

            % one-tailed
            else
                threshold = abs(ne_methods.sdist('tinv', opts.svcthresh / opts.svc, map.DF1));
            end
            if map.Type == 2
                threshold = ne_methods.correlinvtstat(threshold, map.DF1 + 2);
            end

        % F map
        case 4
            threshold = ne_methods.sdist('finv', opts.svcthresh / opts.svc, map.DF1, map.DF2, true);

        % z map
        case 12
            if opts.negative && opts.positive
                threshold = -ne_methods.sdist('norminv', 0.5 * opts.svcthresh / opts.svc, 0, 1);
            else
                threshold = -ne_methods.sdist('norminv', opts.svcthresh / opts.svc, 0, 1);
            end

        % prob map
        case 16
            threshold = opts.svcthresh / opts.svc;
    end

    % remove invalid entries first
    vmp(isinf(vmp) | isnan(vmp)) = 0;

    % set in data
    map.EnableClusterCheck = 0;
    map.LowerThreshold = threshold;
    map.UpperThreshold = 2 * threshold;
    map.RunTimeVars.FWHMResImg(svcvmsk) = NaN;
    if any(strcmp(stype, {'cmp', 'vmp'}))
        bc.Map(mapno) = map;
    elseif any(strcmp(stype, {'hdr', 'head'}))
        bc.RunTimeVars.Map(mapno) = map;
    end
    xo.C = bc;
end

% ensure threshold is > 0
threshold = max(threshold, sqrt(eps));

% output header
thead = cell(1, 6);
thead{1} = sprintf('   Clustertable of map:   "%s"', map.Name);
thead{2} = sprintf('           Type of map:   %s', thtype);
thead{3} = sprintf('    Degrees of freedom:   %s', thdegf);
thead{4} = sprintf('   Cluster k-threshold:   %d mm^3 (%d voxel)', opts.minsize * prod(bb.ResXYZ), opts.minsize);
thead{5} = sprintf(' Applied map threshold:   %.5f (p < %.5f)', threshold, mthreshp);
thead{6} = sprintf('   Type of coordinates:   %s', opts.lupcrd);

% pass on to clustervol
[ctab, cs, vmp] = ne_methods.clustervol(vmp, threshold, opts.minsize, opts);

% store back in struct if cluster-threshold enabled
if opts.minsize > 1
    switch(stype)
        case 'cmp'
            bc.Map(mapno).CMPDataCT = ne_methods.dilate3d(vmp ~= 0);
        case 'hdr'
            bc.VoxelDataCT{mapno} = ne_methods.dilate3d(vmp ~= 0);
        case 'head'
            bc.Brick(mapno).DataCT = ne_methods.dilate3d(vmp ~= 0);
        case 'vmp'
            bc.Map(mapno).VMPDataCT = ne_methods.dilate3d(vmp ~= 0);
    end
    xo.C = bc;
end

% extend ctab
ctab = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n\n%s', thead{:}, ctab);

% create VOI structure ?
if nargout > 3
    voi = xff('new:voi');
    voic = voi.C;
    voic.RunTimeVars.ClusterTableResXYZ = bb.ResXYZ;
    if ~isempty(cs)
        voic.VOI(numel(cs)).Name = '';
    else
        voic.VOI(:) = [];
    end
    cnum = 1;
    scnum = 1;
    for cc = 1:numel(cs)
        if cs(cc).localmax ~= 'L'
            voic.VOI(cc).Name = sprintf('Cluster%04d_%g_%g_%g_%s', ...
                cnum, 0.1 .* round(10 .* cs(cc).rwcoords(1, :)), ...
                strrep(regexprep(char(cs(cc).talout), '\(.*$', ''), ' ', '_'));
            cnum = cnum + 1;
            scnum = 1;
        else
            voic.VOI(cc).Name = sprintf('SC%04d_%04d_%g_%g_%g_%s', ...
                cnum - 1, scnum, 0.1 .* round(10 .* cs(cc).rwcoords(1, :)), ...
                strrep(regexprep(char(cs(cc).talout), '\(.*$', ''), ' ', '_'));
            scnum = scnum + 1;
        end
        if cs(cc).values(1) >= 0
            voic.VOI(cc).Color = floor([255, 127.999 * rand(1, 2)]);
        else
            voic.VOI(cc).Color = floor([127.999 * rand(1, 2), 255]);
        end
        voic.VOI(cc).NrOfVoxels = size(cs(cc).rwcoords, 1);
        voic.VOI(cc).Voxels = cs(cc).rwcoords;
        voic.VOI(cc).VoxelValues = cs(cc).values;
        voic.VOI(cc).IsLocalMax = double(cs(cc).localmax == 'L');
    end
    voic.NrOfVOIs = numel(voic.VOI);
    voi.C = voic;
end

% return table if only one output is requested
if nargout < 2
    cs = ctab;
end
