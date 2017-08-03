function xo = aft_RenderToTransimg(xo, mview, ti, opts)
% AFT::RenderToTransimg  - render spatial data into transimg object
%
% FORMAT:       [obj = ] obj.SliceToTransimg(mview, ti [, opts]);
%
% Input fields:
%
%       mview       4x4 viewing matrix (see renderv3d: opts.mview)
%       ti          1x1 transimg object(s)
%       opts        optional settings
%        .alpha     alpha blending vector between min/max
%        .avol      alpha volume (alternative to alpha setting)
%        .backcolor 1x3 background color
%        .bbox      3x2 voxel-index based bounding box
%        .collut    Cx3 color palette for grayscale data (default: [])
%        .layer     target layer in transimg (not checked, errors thrown)
%        .mapvol    required for CMP/VMP objects (default: 1)
%        .msk       masking object (HDR, HEAD, MSK, VMR)
%        .mskrange  1x2 range to scale msk to [0, 1] range
%        .mskthresh 1x2 with minimum and maximum alpha for [0, 1] range
%        .max       max intensity (default: from RunTimeVars)
%        .min       min intensity (default: from RunTimeVars)
%        .preview   perform preview (default: false)
%        .prevres   preview resolution (in voxels per dim, default: 128)
%        .sranges   slicing ranges (see 'help sliceranges')
%        .svinterp  stats volume interpolation method (default: 'linear')
%        .svol      stats volume (Vx3 cell array with xff, mapvol, and opts)
%        .svolafac  stats volume alpha factor (default: 1)
%        .update    update transio (default: false, optional give filename)
%        .v16       use VMRData16 on VMRs if available (default: false)
%        .warpip    warping interpolation (see renderv3d, default: 'cubic')
%
% Output fields:
%
%       sval        sampled values at coordinate
%       scrd        sampling coordinate
%
% TYPES: HDR, HEAD, MGH, VMR, VTC
%
% Using: bvcoordconv, flexinterpn, limitrangec, minmaxmean, renderv3d,
%        smoothkern, spmtrf.

% Version:  v1.1
% Build:    16031616
% Date:     Mar-16 2016, 4:17 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, Jochen Weber
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
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'hdr', 'head', 'mgh', 'vmr', 'vtc'}) || ...
   ~isa(mview, 'double') || ~isequal(size(mview), [4, 4]) || ...
    any(isinf(mview(:)) | isnan(mview(:))) || any(mview(4, :) ~= [0, 0, 0, 1]) || ...
   ~isa(ti, 'transimg') || numel(ti) ~= 1
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
se = lower(xo.S.Extensions{1});
bc = xo.C;
h = ti(1).Height;
w = ti(1).Width;
l = numel(ti(1).Layer);
for tc = 2:numel(ti)
    if ti(tc).Height ~= h || ti(tc).Width ~= w || numel(ti(tc).Layer) ~= l
        error('neuroelf:xff:badArgument', 'transimg objects must match in size.');
    end
end
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'alpha') || ~isa(opts.alpha, 'double')
    opts.alpha = [0; 1];
end
if ~isfield(opts, 'avol') || (~isnumeric(opts.avol) && ...
    (numel(opts.avol) ~= 1 || ~xffisobject(opts.avol, true)))
    opts.avol = [];
end
if ~isfield(opts, 'backcolor') || ~isa(opts.backcolor, 'double') || numel(opts.backcolor) ~= 3 || ...
    any(isinf(opts.backcolor) | isnan(opts.backcolor) | opts.backcolor < 0 | opts.backcolor >= 256)
    opts.backcolor = [0, 0, 0];
else
    opts.backcolor = floor(opts.backcolor(:)');
end
if ~isfield(opts, 'bbox') || ~isa(opts.bbox, 'double') || ~isequal(size(opts.bbox), [3, 2]) || ...
    any(isinf(opts.bbox(:)) | isnan(opts.bbox(:)) | opts.bbox(:) < 1)
    opts.bbox = [];
else
    opts.bbox = round(opts.bbox);
end
if ~isfield(opts, 'collut') || ~isa(opts.collut, 'double') || size(opts.collut, 2) ~= 3 || ...
    any(isinf(opts.collut(:)) | isnan(opts.collut(:)) | opts.collut(:) < 0 | opts.collut(:) > 255)
    opts.collut = [0, 0, 0; 255, 255, 255];
elseif all(opts.collut(:) <= 1)
    opts.collut = round(255 .* double(opts.collut(:, :)));
end
if ~isfield(opts, 'layer') || ~isa(opts.layer, 'double') || numel(opts.layer) ~= 1 || ...
    isinf(opts.layer) || isnan(opts.layer) || opts.layer < 1
    opts.layer = 1;
else
    opts.layer = round(opts.layer);
end
if ~isfield(opts, 'mapvol') || ~isa(opts.mapvol, 'double') || numel(opts.mapvol) ~= 1 || ...
    isinf(opts.mapvol) || isnan(opts.mapvol) || opts.mapvol < 1
    opts.mapvol = 1;
end
if ~isfield(opts, 'max') || ~isa(opts.max, 'double') || numel(opts.max) ~= 1 || ...
    isinf(opts.max) || isnan(opts.max)
    opts.max = [];
end
if ~isfield(opts, 'min') || ~isa(opts.min, 'double') || numel(opts.min) ~= 1 || ...
    isinf(opts.min) || isnan(opts.min)
    opts.min = [];
end
if ~isfield(opts, 'v16') || ~islogical(opts.v16) || numel(opts.v16) ~= 1 || ...
   ~strcmp(se, 'vmr') || ~isequal(size(bc.VMRData), size(bc.VMRData16))
    opts.v16 = false;
end
if isempty(opts.max) || isempty(opts.min)
    if ~isfield(bc.RunTimeVars, 'ScalingWindow') || numel(bc.RunTimeVars.ScalingWindow) ~= 2
        aft_SetScalingWindow(xo);
        bc = xo.C;
    end
    if opts.v16
        if isempty(opts.max)
            opts.max = bc.RunTimeVars.ScalingWindow16(2);
        end
        if isempty(opts.min)
            opts.min = bc.RunTimeVars.ScalingWindow16(1);
        end
    else
        if isempty(opts.max)
            opts.max = bc.RunTimeVars.ScalingWindow(2);
        end
        if isempty(opts.min)
            opts.min = bc.RunTimeVars.ScalingWindow(1);
        end
    end
end
if ~isfield(opts, 'preview') || ~islogical(opts.preview) || numel(opts.preview) ~= 1
    opts.preview = false;
end
if ~isfield(opts, 'prevres') || ~isa(opts.prevres, 'double') || numel(opts.prevres) ~= 1 || ...
    isinf(opts.prevres) || isnan(opts.prevres) || opts.prevres <= 0
    opts.prevres = 128;
else
    opts.prevres = round(min(256, max(16, opts.prevres)));
end
if ~isfield(opts, 'sranges') || ~isstruct(opts.sranges) || numel(opts.sranges) ~= 1
    opts.sranges = struct;
end
if ~isfield(opts, 'svinterp') || ~ischar(opts.svinterp) || isempty(opts.svinterp) || ...
   ~any(strcmpi(opts.svinterp(:)', {'cubic', 'lanczos3', 'linear', 'nearest'}))
    opts.svinterp = 'linear';
else
    opts.svinterp = opts.svinterp(:)';
end
if ~isfield(opts, 'svol') || ~iscell(opts.svol) || size(opts.svol, 2) ~= 3
    opts.svol = cell(0, 3);
end
svol = opts.svol;
for svc = size(svol, 1):-1:1
    if numel(svol{svc, 1}) ~= 1 || ...
       ~xffisobject(svol{svc, 1}, true, {'glm', 'hdr', 'head', 'vmp', 'vmr', 'vtc'}) || ...
        numel(svol{svc, 2}) ~= 1 || ~isa(svol{svc, 2}, 'double') || ...
        isinf(svol{svc, 2}) || isnan(svol{svc, 2}) || svol{svc, 2} < 1
        svol(svc, :) = [];
        continue;
    end
    svol{svc, 1} = svol{svc, 1}.L;
end
if ~isfield(opts, 'svolafac') || ~isa(opts.svolafac, 'double') || numel(opts.svolafac) ~= 1 || ...
    isinf(opts.svolafac) || isnan(opts.svolafac) || opts.svolafac <= 0
    opts.svolafac = 1;
else
    opts.svolafac = min(1, opts.svolafac);
end
if ~isfield(opts, 'update') || ~islogical(opts.update) || numel(opts.update) ~= 1
    opts.update = false;
end
if ~isfield(opts, 'warpip') || ~ischar(opts.warpip) || isempty(opts.warpip) || ...
   ~any(strcmpi(opts.warpip(:)', {'linear', 'cubic', 'lanczos3'}))
    opts.warpip = '';
else
    opts.warpip = lower(opts.warpip(:)');
end

% alpha volume
if ~isfield(xo.H, 'RenderAVol') || ~iscell(xo.H.RenderAVol)
    xo.H.RenderAVol = cell(1, 2);
end
bbx = aft_BoundingBox(xo);
bbxbx = bbx.BBox;
if strcmp(se, 'vmr')
    res = [bc.VoxResX, bc.VoxResY, bc.VoxResZ];
    bbxbx(2, :) = bbxbx(2, :) + res;
    trans = bbx.QuatB2T;
elseif ~any(strcmp(se, {'hdr', 'head', 'mgh'}))
    res = bc.Resolution .* ones(1, 3);
    trans = bbx.QuatB2T;
else
    if strcmp(se, 'hdr')
        vsz = size(bc.VoxelData);
        trans = hdr_CoordinateFrame(xo, opts.mapvol);
        if any(vsz == 0)
            vsz = trans.Dimensions(1:3);
        end
    elseif strcmp(se, 'head')
        vsz = size(bc.Brick(min(numel(bc.Brick), round(opts.mapvol))).Data);
        trans = head_CoordinateFrame(xo, opts.mapvol);
    else
        vsz = size(bc.MGHData);
        trans = mgh_CoordinateFrame(xo, 1);
    end
    if numel(vsz) < 3
        vsz(3) = 1;
    end
    trans = trans.Trf;
end
if numel(opts.avol) == 1 && xffisobject(opts.avol, true, {'hdr', 'head', 'msk', 'vmr', 'vtc'})
    if ~isequal(xo.H.RenderAVol{1}, opts.avol)

        % store for future reference
        xo.H.RenderAVol{1} = opts.avol;

        % for all but HDR/HEAD
        if ~any(strcmp(se, {'hdr', 'head', 'mgh'}))

            % resample data
            xo.H.RenderAVol{2} = aft_SampleBVBox(opts.avol, struct('BBox', bbxbx, 'ResXYZ', res), 1);

        % sample data in TAL space with inverse matrix
        else
            xo.H.RenderAVol{2} = aft_SampleTalBox(opts.avol, ...
                struct('BBox', [ones(1, 3); vsz], 'ResXYZ', [1, 1, 1]), 1, 'linear', trans);
        end
        avolbc = opts.avol.C;
        avolscl = avolbc.RunTimeVars.ScalingWindow;
        xo.H.RenderAVol{2} = single(ne_methods.limitrangec((1 / (avolscl(2) - avolscl(1))) .* ...
            (xo.H.RenderAVol{2} - avolscl(1)), 0, 1, 0));
    end
    opts.avol = xo.H.RenderAVol(2);
else
    opts.avol = {opts.avol};
end

% compare to previously selected volumes
if ~isfield(xo.H, 'RenderSVol') || ~iscell(xo.H.RenderSVol)
    xo.H.RenderSVol = cell(0, 4);
end
if ~isequal(xo.H.RenderSVol(:, 1:2), svol(:, 1:2))

    % initialize volumes
    xo.H.RenderSVol = [svol(:, 1:2), cell(size(svol, 1), 3)];

    % resample data
    for svc = 1:size(svol, 1)
        if svc > 1 && isequal(svol{svc-1, 1}, svol{svc, 1}) && ...
            isequal(svol{svc-1, 2}, svol{svc, 2})
            if isstruct(opts.svol{svc, 3}) && numel(opts.svol{svc, 3}) == 1 && ...
                isfield(opts.svol{svc, 3}, 'max') && isfield(opts.svol{svc, 3}, 'min') && ...
                isnumeric(opts.svol{svc, 3}.max) && isnumeric(opts.svol{svc, 3}.min) && ...
                numel(opts.svol{svc, 3}.max) == 1 && numel(opts.svol{svc, 3}.min) == 1 && ...
                opts.svol{svc, 3}.max < 0 && opts.svol{svc, 3}.min < opts.svol{svc, 3}.max
                svolo = opts.svol{svc, 3};
                xo.H.RenderSVol(svc, 3:5) = xo.H.RenderSVol(svc-1, 3:5);
                xo.H.RenderSVol{svc, 4} = ne_methods.limitrangec(xo.H.RenderSVol{svc, 4}, svolo.min + double(32 * eps('single')), 0, 0);
            else
                xo.H.RenderSVol(svc, 3:5) = xo.H.RenderSVol(svc-1, 3:5);
            end
        else
            if ~any(strcmp(se, {'hdr', 'head', 'mgh'}))
                if ~xffisobject(opts.svol{svc, 1}, true, 'vtc')
                    xo.H.RenderSVol{svc, 4} = single(aft_SampleBVBox(opts.svol{svc, 1}, ...
                        struct('BBox', bbxbx, 'ResXYZ', res), svol{svc, 2}, opts.svinterp));
                else
                    [xo.H.RenderSVol{svc, 4}, svolbb, xo.H.RenderSVol{svc, 5}] = ...
                        aft_SampleBVBox(opts.svol{svc, 1}, ...
                        struct('BBox', bbxbx, 'ResXYZ', res), svol{svc, 2}, ...
                        opts.svinterp, [], [], 6);
                    xo.H.RenderSVol{svc, 4} = single(xo.H.RenderSVol{svc, 4});
                    if opts.svolafac ~= 1
                        xo.H.RenderSVol{svc, 5} = opts.svolafac .* xo.H.RenderSVol{svc, 5};
                    end
                    xo.H.RenderSVol{svc, 5} = single(xo.H.RenderSVol{svc, 5});
                end
            else
                if ~xffisobject(opts.svol{svc, 1}, true, 'vtc')
                    xo.H.RenderSVol{svc, 4} = single(aft_SampleTalBox(opts.svol{svc, 1}, ...
                        struct('BBox', [ones(1, 3); vsz(1:3)], 'ResXYZ', [1, 1, 1]), svol{svc, 2}, ...
                        opts.svinterp, trans));
                else
                    [xo.H.RenderSVol{svc, 4}, svolbb, xo.H.RenderSVol{svc, 5}] = ...
                        aft_SampleTalBox(opts.svol{svc, 1}, ...
                        struct('BBox', [ones(1, 3); vsz(1:3)], 'ResXYZ', [1, 1, 1]), svol{svc, 2}, ...
                        opts.svinterp, trans, [], 6);
                    xo.H.RenderSVol{svc, 4} = single(xo.H.RenderSVol{svc, 4});
                    if opts.svolafac ~= 1
                        xo.H.RenderSVol{svc, 5} = opts.svolafac .* xo.H.RenderSVol{svc, 5};
                    end
                    xo.H.RenderSVol{svc, 5} = single(xo.H.RenderSVol{svc, 5});
                end
            end
            if size(xo.H.RenderSVol{svc, 4}, 4) > 1
                rsvolsz = size(xo.H.RenderSVol{svc, 4});
                xo.H.RenderSVol{svc, 4} = reshape( ...
                    xo.H.RenderSVol{svc, 4}, [rsvolsz(1:3), 1, rsvolsz(4:end)]);
            end

            % limit if necessary
            if isstruct(opts.svol{svc, 3}) && numel(opts.svol{svc, 3}) == 1 && ...
                isfield(opts.svol{svc, 3}, 'max') && isfield(opts.svol{svc, 3}, 'min') && ...
                isnumeric(opts.svol{svc, 3}.max) && isnumeric(opts.svol{svc, 3}.min) && ...
                numel(opts.svol{svc, 3}.max) == 1 && numel(opts.svol{svc, 3}.min) == 1 && ...
                opts.svol{svc, 3}.max <= 0 && opts.svol{svc, 3}.min < opts.svol{svc, 3}.max
                xo.H.RenderSVol{svc, 4} = ne_methods.limitrangec(xo.H.RenderSVol{svc, 4}, opts.svol{svc, 3}.min + double(32 * eps('single')), 0, 0);
            end
            mmm = ne_methods.minmaxmean(xo.H.RenderSVol{svc, 4});
            xo.H.RenderSVol{svc, 3} = mmm(1:2);
        end
    end
    for svc = 1:size(svol, 1)
        if isempty(xo.H.RenderSVol{svc, 5})
            continue;
        end
        if isstruct(opts.svol{svc, 3}) && numel(opts.svol{svc, 3}) == 1 && ...
            isfield(opts.svol{svc, 3}, 'max') && isfield(opts.svol{svc, 3}, 'min') && ...
            isnumeric(opts.svol{svc, 3}.max) && isnumeric(opts.svol{svc, 3}.min) && ...
            numel(opts.svol{svc, 3}.max) == 1 && numel(opts.svol{svc, 3}.min) == 1 && ...
            opts.svol{svc, 3}.max <= 0 && opts.svol{svc, 3}.min < opts.svol{svc, 3}.max
            xo.H.RenderSVol{svc, 5}(xo.H.RenderSVol{svc, 4} >= 0) = 0;
        else
            xo.H.RenderSVol{svc, 5}(xo.H.RenderSVol{svc, 4} <= 0) = 0;
        end
    end
end

% for switch on filetype
bc = xo.C;
ft = lower(xo.S.Extensions{1});

% get number of volumes
switch (ft)
    case 'hdr'
        nvol = size(bc.VoxelData, 4);
    case 'head'
        nvol = numel(bc.Brick);
    case 'mgh'
        nvol = 1;
    case 'vmr'
        if opts.v16 && isequal(size(bc.VMRData), size(bc.VMRData16))
            nvol = 2;
        else
            nvol = 1;
        end
        opts.mapvol = nvol;
    case 'vtc'
        nvol = size(bc.VTCData, 1);
    otherwise
        error('neuroelf:xff:objectTypeUnsupported', 'Unsupported object type for ::SliceToTransimg method.');
end
if opts.mapvol > nvol
    error('neuroelf:xff:badArgument', 'Requested mapvol option out of bounds.');
end

% get volume data
if ~opts.preview
    vols = [{aft_GetVolume(xo, opts.mapvol)}; xo.H.RenderSVol(:, 4)];
else
    vols = {aft_GetVolume(xo, opts.mapvol)};
end

% preview for failing to check size
if isempty(opts.avol)
    opts.avol = {[]};
elseif (size(vols{1}, 1) ~= size(opts.avol{1}, 1) || ...
    size(vols{1}, 2) ~= size(opts.avol{1}, 2) || size(vols{1}, 3) ~= size(opts.avol{1}, 3))

    % no alpha volume
    opts.avol = {[]};
end
if ~isempty(xo.H.RenderSVol)
    opts.avol(2:(size(xo.H.RenderSVol, 1)+1)) = xo.H.RenderSVol(:, 5)';
end

% preview
osize = size(vols{1});
pstep = (osize - 1) ./ (opts.prevres - 1);
pfrom = 0.5 .* (pstep + 1);
szvol = size(vols{1});
cr1 = [Inf, Inf, Inf; ones(2, 3); szvol(1:3)];
cr2 = [Inf, Inf, Inf; ones(2, 3); szvol(1:3)];
cr3 = [Inf, Inf, Inf; ones(2, 3); szvol(1:3)];
cr1(2:3, 1) = [pfrom(1); pstep(1)];
cr2(2:3, 2) = [pfrom(2); pstep(2)];
cr2(4, 1) = opts.prevres;
cr3(2:3, 3) = [pfrom(3); pstep(3)];
cr3(4, 1:2) = opts.prevres;
if ~isfield(xo.H, 'RenderPreview') || isempty(xo.H.RenderPreview)
    gkern1 = ne_methods.smoothkern(pstep(1), 0.001);
    gkern2 = ne_methods.smoothkern(pstep(2), 0.001);
    gkern3 = ne_methods.smoothkern(pstep(3), 0.001);
    if numel(szvol) < 4
        xo.H.RenderPreview = {single(ne_methods.flexinterpn(ne_methods.flexinterpn(ne_methods.flexinterpn(vols{1}, ...
            cr1, {gkern1, [0; 1; 0], [0; 1; 0]}, {1, 1, 1}), ...
            cr2, {[0; 1; 0], gkern2, [0; 1; 0]}, {1, 1, 1}), ...
            cr3, {[0; 1; 0], [0; 1; 0], gkern3}, {1, 1, 1}))};
    else
        fivol = single(zeros([opts.prevres, opts.prevres, opts.prevres]));
        ficnt = prod(szvol(4:end));
        fisrc = reshape(vols{1}, [szvol(1:3), ficnt]);
        for fic = 1:ficnt
            fivol = fivol + single(ne_methods.flexinterpn(ne_methods.flexinterpn(ne_methods.flexinterpn(fisrc(:, :, :, fic), ...
            cr1, {gkern1, [0; 1; 0], [0; 1; 0]}, {1, 1, 1}), ...
            cr2, {[0; 1; 0], gkern2, [0; 1; 0]}, {1, 1, 1}), ...
            cr3, {[0; 1; 0], [0; 1; 0], gkern3}, {1, 1, 1}));
        end
        fivol = single(1 / ficnt) .* fivol;
        xo.H.RenderPreview = {fivol};
    end
    xo.H.RenderPView = ne_methods.spmtrf([0, 0, 0], [0, 0, 0], pstep(1:3));
end
if ~isempty(opts.avol) && ~isempty(opts.avol{1}) && (numel(xo.H.RenderPreview) < 2 || ...
    ~isequal(size(xo.H.RenderPreview{1}), size(xo.H.RenderPreview{2})))
    gkern1 = ne_methods.smoothkern(pstep(1), 0.001);
    gkern2 = ne_methods.smoothkern(pstep(2), 0.001);
    gkern3 = ne_methods.smoothkern(pstep(3), 0.001);
    xo.H.RenderPreview{2} = single(ne_methods.flexinterpn(ne_methods.flexinterpn(ne_methods.flexinterpn(opts.avol{1}, ...
        cr1, {gkern1, [0; 1; 0], [0; 1; 0]}, {1, 1, 1}), ...
        cr2, {[0; 1; 0], gkern2, [0; 1; 0]}, {1, 1, 1}), ...
        cr3, {[0; 1; 0], [0; 1; 0], gkern3}, {1, 1, 1}));
end
if ~isempty(opts.warpip)
    warpip = opts.warpip;
else
    warpip = '';
end
if opts.preview;
    vols = xo.H.RenderPreview(1);
    pview = xo.H.RenderPView;
    if ~isempty(opts.avol) && ~isempty(opts.avol{1})
        opts.avol = xo.H.RenderPreview(2);
    end
    if isempty(warpip)
        warpip = 'linear';
    end
else
    pview = eye(4);
    if isempty(warpip)
        warpip = 'cubic';
    end
end

% determine view matrix
if ~isfield(xo.H, 'RenderMView') || ~isa(xo.H.RenderMView, 'double') || ...
   ~isequal(size(xo.H.RenderMView), [4, 4])
    if ~any(strcmp(se, {'hdr', 'head', 'mgh'}))
        xo.H.RenderMView = ne_methods.bvcoordconv(zeros(0, 3), 'bvc2tal', bbx);
    elseif strcmp(se, 'hdr')
        xo.H.RenderMView = hdr_CoordinateFrame(xo, opts.mapvol);
        xo.H.RenderMView = xo.H.RenderMView.Trf;
    elseif strcmp(se, 'head')
        xo.H.RenderMView = head_CoordinateFrame(xo, opts.mapvol);
        xo.H.RenderMView = xo.H.RenderMView.Trf;
    else
        xo.H.RenderMView = mgh_CoordinateFrame(xo, opts.mapvol);
        xo.H.RenderMView = xo.H.RenderMView.Trf;
    end
    szvol = size(vols{1});
    xo.H.RenderMView(1:3, 4) = xo.H.RenderMView(1:3, 4) + ...
        xo.H.RenderMView(1:3, 1:3) * (0.5 .* (szvol(1:3)' + 1));
end

% general options
ropts = struct('atable', opts.alpha, 'avolume', {opts.avol}, 'backcolor', (1 / 255) .* opts.backcolor, ...
    'ctable', opts.collut, 'imagesize', [h, w], 'imax', opts.max, 'imin', opts.min, ...
    'mview', ne_methods.spmtrf([0, 0, 0], [0, -pi/2, pi]) * mview * xo.H.RenderMView * pview, ...
    'rendtype', 'color', 'warpip', warpip);

% update
if islogical(opts.update) && opts.update
    ropts.update = {@setlayer, ti, opts.layer, '%p', 1; ...
        @render, ti, [], [], []; @display, ti, [], [], []; @drawnow, [], [], [], []};
end

% extend options
if numel(vols) > 1
    ropts.atable = {ropts.atable};
    ropts.ctable = {ropts.ctable};
    ropts.imax = [ropts.imax; zeros(size(svol, 1), 1)];
    ropts.imin = [ropts.imin; zeros(size(svol, 1), 1)];
    for svc = 1:size(svol, 1)
        svt = svc + 1;
        ropts.atable{svt} = (0:0.25:1)';
        if numel(ropts.avolume) < svt
            ropts.avolume{svt} = [];
        end
        ropts.ctable{svt} = (0:0.25:1)' * ones(1, 3);
        ropts.imax(svt) = xo.H.RenderSVol{svc, 3}(2);
        ropts.imin(svt) = xo.H.RenderSVol{svc, 3}(1);
        if ~isstruct(svol{svc, 3}) || numel(svol{svc, 3}) ~= 1
            continue;
        end
        if isfield(svol{svc, 3}, 'alpha') && isa(svol{svc, 3}.alpha, 'double') && ...
           ~isempty(svol{svc, 3}.alpha)
            ropts.atable{svt} = svol{svc, 3}.alpha(:);
        end
        if isfield(svol{svc, 3}, 'collut') && isa(svol{svc, 3}.collut, 'double') && ...
           ~isempty(svol{svc, 3}.collut)
            ropts.ctable{svt} = svol{svc, 3}.collut(:, :, 1);
        end
        if isfield(svol{svc, 3}, 'max') && isa(svol{svc, 3}.max, 'double') && ...
           ~isempty(svol{svc, 3}.max)
            ropts.imax(svt) = svol{svc, 3}.max(1);
        end
        if isfield(svol{svc, 3}, 'min') && isa(svol{svc, 3}.min, 'double') && ...
           ~isempty(svol{svc, 3}.min)
            ropts.imin(svt) = svol{svc, 3}.min(1);
        end
    end
end

% add sliceranges
if ~opts.preview && isfield(opts.sranges, 'xfromy') && isfield(opts.sranges, 'xtoy') && ...
    isfield(opts.sranges, 'xfromz') && isfield(opts.sranges, 'xtoz') && ...
    isfield(opts.sranges, 'yfromx') && isfield(opts.sranges, 'ytox') && ...
    isfield(opts.sranges, 'yfromz') && isfield(opts.sranges, 'ytoz') && ...
    isfield(opts.sranges, 'zfromx') && isfield(opts.sranges, 'ztox') && ...
    isfield(opts.sranges, 'zfromy') && isfield(opts.sranges, 'ztoy')
    sopts = opts.sranges;
    if isempty(opts.bbox)
        ropts.xfromy = sopts.xfromy;
        ropts.xtoy   = sopts.xtoy;
        ropts.xfromz = sopts.xfromz;
        ropts.xtoz   = sopts.xtoz;
        ropts.yfromx = sopts.yfromx;
        ropts.ytox   = sopts.ytox;
        ropts.yfromz = sopts.yfromz;
        ropts.ytoz   = sopts.ytoz;
        ropts.zfromx = sopts.zfromx;
        ropts.ztox   = sopts.ztox;
        ropts.zfromy = sopts.zfromy;
        ropts.ztoy   = sopts.ztoy;
    else
        ropts.xfromy = max(opts.bbox(2, 1), sopts.xfromy) + 0 .* sopts.xfromy;
        ropts.xtoy   = min(opts.bbox(2, 2), sopts.xtoy) + 0 .* sopts.xtoy;
        ropts.xfromz = max(opts.bbox(3, 1), sopts.xfromz) + 0 .* sopts.xfromz;
        ropts.xtoz   = min(opts.bbox(3, 2), sopts.xtoz) + 0 .* sopts.xtoz;
        ropts.yfromx = max(opts.bbox(1, 1), sopts.yfromx) + 0 .* sopts.yfromx;
        ropts.ytox   = min(opts.bbox(1, 2), sopts.ytox) + 0 .* sopts.ytox;
        ropts.yfromz = max(opts.bbox(3, 1), sopts.yfromz) + 0 .* sopts.yfromz;
        ropts.ytoz   = min(opts.bbox(3, 2), sopts.ytoz) + 0 .* sopts.ytoz;
        ropts.zfromx = max(opts.bbox(1, 1), sopts.zfromx) + 0 .* sopts.zfromx;
        ropts.ztox   = min(opts.bbox(1, 2), sopts.ztox) + 0 .* sopts.ztox;
        ropts.zfromy = max(opts.bbox(2, 1), sopts.zfromy) + 0 .* sopts.zfromy;
        ropts.ztoy   = min(opts.bbox(2, 2), sopts.ztoy) + 0 .* sopts.ztoy;
        if opts.bbox(1, 1) > 1
            ropts.xfromy(1:opts.bbox(1, 1)-1) = NaN;
            ropts.xtoy(1:opts.bbox(1, 1)-1) = NaN;
            ropts.xfromz(1:opts.bbox(1, 1)-1) = NaN;
            ropts.xtoz(1:opts.bbox(1, 1)-1) = NaN;
        end
        if opts.bbox(1, 2) < numel(ropts.xfromy)
            ropts.xfromy(opts.bbox(1, 2)+1:end) = NaN;
            ropts.xtoy(opts.bbox(1, 2)+1:end) = NaN;
            ropts.xfromz(opts.bbox(1, 2)+1:end) = NaN;
            ropts.xtoz(opts.bbox(1, 2)+1:end) = NaN;
        end
        if opts.bbox(2, 1) > 1
            ropts.yfromx(1:opts.bbox(2, 1)-1) = NaN;
            ropts.ytox(1:opts.bbox(2, 1)-1) = NaN;
            ropts.yfromz(1:opts.bbox(2, 1)-1) = NaN;
            ropts.ytoz(1:opts.bbox(2, 1)-1) = NaN;
        end
        if opts.bbox(2, 2) < numel(ropts.yfromx)
            ropts.yfromx(opts.bbox(2, 2)+1:end) = NaN;
            ropts.ytox(opts.bbox(2, 2)+1:end) = NaN;
            ropts.yfromz(opts.bbox(2, 2)+1:end) = NaN;
            ropts.ytoz(opts.bbox(2, 2)+1:end) = NaN;
        end
        if opts.bbox(3, 1) > 1
            ropts.zfromx(1:opts.bbox(3, 1)-1) = NaN;
            ropts.ztox(1:opts.bbox(3, 1)-1) = NaN;
            ropts.zfromy(1:opts.bbox(3, 1)-1) = NaN;
            ropts.ztoy(1:opts.bbox(3, 1)-1) = NaN;
        end
        if opts.bbox(3, 2) < numel(ropts.zfromx)
            ropts.zfromx(opts.bbox(3, 2)+1:end) = NaN;
            ropts.ztox(opts.bbox(3, 2)+1:end) = NaN;
            ropts.zfromy(opts.bbox(3, 2)+1:end) = NaN;
            ropts.ztoy(opts.bbox(3, 2)+1:end) = NaN;
        end
    end
end

% render
rimg = ne_methods.renderv3d(vols, ropts);

% set in transimg
setlayer(ti, opts.layer, rimg, 1);
