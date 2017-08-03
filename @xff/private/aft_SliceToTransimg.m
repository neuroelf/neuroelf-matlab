function [sval, scrd, oval] = aft_SliceToTransimg(xo, crd, ti, opts)
% AFT::SliceToTransimg  - slice spatial data into transimg object(s)
%
% FORMAT:       [sval, scrd = ] obj.SliceToTransimg(crd, ti [, opts]);
%
% Input fields:
%
%       crd         1x3 coordinate to create oblique slices at
%       ti          1x1 or 1x3 transimg object(s)
%       opts        optional settings
%        .cctype    cross-correlation map data, either of 'lag' or {'r'}
%        .conv      convention, either 'neurological' or {'radiological'}
%        .dir       direction, either of {'all'}, 'cor', 'sag', 'tra'
%        .fmr       required for AMR/MAP/fCMP objects
%        .frame     min/max values (default: 256 box)
%        .gcolblend gray color blend mode, one of 'hsv', 'lut', or {'rgb'}
%        .gcolhigh  low-intensity color (default: [255, 255, 255])
%        .gcollow   low-intensity color (default: [0, 0, 0])
%        .gcollut   256x3 color palette for grayscale data (default: [])
%        .grayalpha gray alpha (multiplied by 1/255, default: Inf = opaque)
%        .layers    target layers in transimg (not checked, errors thrown)
%        .mapvol    required for CMP/VMP objects, defaults to all
%        .mapvola   optional alpha volume (will be sliced as well)
%        .method    resampling, one of 'cubic', 'lanczos3', {'linear'}, 'nearest'
%        .orient    orientation of cor and tra slice, either of {'n'} or 'r'
%        .rgbalpha  alpha level (-Inf...1) on top of Map.TransColorFactor
%        .rgbbars   1x4 double position for color bars (last layer, opaque)
%        .rgbcol    Cx3 RGB colors (double/uint8), C must be even
%        .rgbctails coloring the tails
%        .rgblthr   lower threshold for scaling to RGB
%        .rgbuthr   upper threshold for scaling to RGB
%        .space     either of 'bvs' or {'tal'}
%        .tensorc   tensor computation, either of
%        .trans     4x4 transformation matrix, used for volume datatypes
%        .trf       TRF object(s) passed on to samplefmrspace
%        .type      either of 'dsq', {'gray'}, or 'rgb'
%        .v16       for VMR objects, if true access available V16 data
%        .xscale    scaling for type-x objects (e.g. DARTEL flowfields)
%
% Output fields:
%
%       sval        sampled values at coordinate
%       scrd        sampling coordinate
%
% TYPES: AMR, AVA, CMP, DMR, DDT, FMR, GLM, HDR, HEAD, MAP, MGH, MSK, NLF, SRF, TVL, VDW, VMP, VMR, VTC
%
% Note: for AMR/DMR/FMR/MAP objects, the trf option must be set with
%       at least the IA/FA (or FA alone) TRF object(s) to work
%
% Note: other than for the OBJ::SliceData3D method, the space can not be
%       bvi or bvc!
%
% Note: the rgbbars defaults to empty, coordinates are *relative* positions
%
% Using: findfirst, flexinterpn, hsvconv, limitrangec, minmaxmean, threshlutc,
%        threshmapc, winsorize.

% Version:  v1.1
% Build:    16061321
% Date:     Jun-13 2016, 9:35 PM EST
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

% persistent VMR coloring
persistent stt_vmrlut;
if isempty(stt_vmrlut)
    stt_vmrlut = xffnewcont('olt');
    stt_vmrlut = stt_vmrlut.Colors;
    stt_vmrlut(end+1:30, :) = 255;
end

% neuroelf library
global ne_methods;

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true) || ~isa(ti, 'transimg') || ...
   ~any([1, 3] == numel(ti)) || ~isa(crd, 'double') || numel(crd) ~= 3 || ...
    any(isinf(crd) | isnan(crd) | crd < -256 | crd > 256)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
h = ti(1).Height;
w = ti(1).Width;
if h ~= w
    error('neuroelf:xff:badArgument', 'The transimg object(s) must be square.');
end
l = numel(ti(1).Layer);
for tc = 2:numel(ti)
    if ti(tc).Height ~= h || ti(tc).Width ~= w || numel(ti(tc).Layer) ~= l
        error('neuroelf:xff:badArgument', 'transimg objects must match in size.');
    end
end
bc = xo.C;
rtv = bc.RunTimeVars;
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dir') || ~ischar(opts.dir) || ...
   ~any(strcmpi(opts.dir(:)', {'a', 'all', 'c', 'cor', 's', 'sag', 't', 'tra'}))
    opts.dir = 'a';
else
    opts.dir = lower(opts.dir(1));
end
if ~isfield(opts, 'frame') || ~isa(opts.frame, 'double') || ~isequal(size(opts.frame), [2, 3]) || ...
    any(isinf(opts.frame(:)) | isnan(opts.frame(:)) | opts.frame(:) < -512 | opts.frame(:) > 512)
    opts.frame = [];
end
if ~isfield(opts, 'gcolblend') || ~ischar(opts.gcolblend) || ...
   ~any(strcmpi(opts.gcolblend(:)', {'hsv', 'lut', 'rgb'}))
    opts.gcolblend = 'r';
    opts.gcollut = [];
else
    opts.gcolblend = lower(opts.gcolblend(1));
    if opts.gcolblend == 'l'
        if ~isfield(opts, 'gcollut') || ~isa(opts.gcollut, 'double') || ...
           ~isequal(size(opts.gcollut), [256, 3]) || ...
            any(isinf(opts.gcollut(:)) | isnan(opts.gcollut(:)))
            opts.gcolblend = 'r';
            opts.gcollut = [];
        else
            opts.gcollut = min(255, max(0, opts.gcollut));
            if all(opts.gcollut(:) <= 1)
                opts.gcollut = 255 .* opts.gcollut;
            end
            opts.gcollut = uint8(round(opts.gcollut));
            opts.gcolhigh = [225, 225, 225];
        end
    else
        opts.gcollut = [];
    end
end
if ~isfield(opts, 'gcolhigh') || ~isa(opts.gcolhigh, 'double') || numel(opts.gcolhigh) ~= 3 || ...
    any(isinf(opts.gcolhigh) | isnan(opts.gcolhigh) | opts.gcolhigh < 0 | opts.gcolhigh > 255)
    opts.gcolhigh = [255, 255, 255];
elseif all(opts.gcolhigh <= 1)
    opts.gcolhigh = round(255 .* opts.gcolhigh(:)');
else
    opts.gcolhigh = round(opts.gcolhigh(:)');
end
if ~isfield(opts, 'gcollow') || ~isa(opts.gcollow, 'double') || numel(opts.gcollow) ~= 3 || ...
    any(isinf(opts.gcollow) | isnan(opts.gcollow) | opts.gcollow < 0 | opts.gcollow > 255)
    opts.gcollow = [0, 0, 0];
elseif all(opts.gcollow <= 1)
    opts.gcollow = round(255 .* opts.gcollow(:)');
else
    opts.gcollow = round(opts.gcollow(:)');
end
if ~isfield(opts, 'grayalpha') || ~isa(opts.grayalpha, 'double') || numel(opts.grayalpha) ~= 1 || ...
    isnan(opts.grayalpha) || opts.grayalpha <= 0
    opts.grayalpha = Inf;
end
if ~isfield(opts, 'layers') || ~isa(opts.layers, 'double') || isempty(opts.layers) || ...
    any(isinf(opts.layers(:)) | isnan(opts.layers(:)) | opts.layers(:) < 1)
    opts.layers = 1;
else
    opts.layers = round(opts.layers(:));
    if numel(opts.layers) ~= numel(unique(opts.layers))
        error('neuroelf:xff:badArgument', 'Option layers must contain unique numbers.');
    end
end
if ~isfield(opts, 'mapvol') || ~isa(opts.mapvol, 'double') || ...
    any(isinf(opts.mapvol(:)) | isnan(opts.mapvol(:)) | opts.mapvol(:) < 1)
    opts.mapvol = [];
end
if ~isfield(opts, 'mapvola') || ~isa(opts.mapvola, 'double') || ...
    any(isinf(opts.mapvola(:)) | isnan(opts.mapvola(:)) | opts.mapvola(:) < 1)
    opts.mapvola = [];
else
    opts.mapvola = round(opts.mapvola(:)');
end
if ~isfield(opts, 'method') || ~ischar(opts.method) || ...
   ~any(strcmpi(opts.method, {'cubic', 'lanczos3', 'linear', 'nearest'}))
    opts.method = 'linear';
else
    opts.method = lower(opts.method(:)');
end
if ~isfield(opts, 'orient') || ~ischar(opts.orient) || isempty(opts.orient) || ...
   ~any(lower(opts.orient(1)) == 'nr')
    opts.orient = 'n';
else
    opts.orient = lower(opts.orient(1));
end
if ~isfield(opts, 'rgbalpha') || numel(opts.rgbalpha) ~= 1 || ~isa(opts.rgbalpha, 'double') || ...
    isnan(opts.rgbalpha) || opts.rgbalpha > 1
    opts.rgbalpha = 1;
end
if ~isfield(opts, 'rgbbars') || ~isa(opts.rgbbars, 'double') || ~isequal(size(opts.rgbbars), [1, 4]) || ...
    any(isinf(opts.rgbbars) | isnan(opts.rgbbars) | opts.rgbbars < 0 | opts.rgbbars > 1) || ...
    any(opts.rgbbars(3:4) <= opts.rgbbars(1:2))
    opts.rgbbars = [];
else
    opts.rgbbars = [1 + floor(opts.rgbbars(1:2) .* [h, w]), ceil(opts.rgbbars(3:4) .* [h, w])];
end
if ~isfield(opts, 'rgbcol') || isempty(opts.rgbcol) || ~isa(opts.rgbcol, 'double') || ...
    size(opts.rgbcol, 2) ~= 3 || mod(size(opts.rgbcol, 1), 2) ~= 0
    opts.rgbcol = [];
end
if ~isfield(opts, 'rgbctails') || numel(opts.rgbctails) ~= 1 || ...
   ~isa(opts.rgbctails, 'double') || ~any(0:3 == opts.rgbctails)
    opts.rgbctails = [];
end
if ~isfield(opts, 'rgblthr') || ~isa(opts.rgblthr, 'double') || numel(opts.rgblthr) ~= 1 || ...
    isinf(opts.rgblthr) || isnan(opts.rgblthr)
    opts.rgblthr = [];
elseif opts.rgblthr == 0
    opts.rgblthr = sqrt(eps);
else
    opts.rgblthr = abs(opts.rgblthr);
end
if ~isfield(opts, 'rgbuthr') || ~isa(opts.rgbuthr, 'double') || numel(opts.rgbuthr) ~= 1 || ...
    isinf(opts.rgbuthr) || isnan(opts.rgbuthr)
    opts.rgbuthr = [];
elseif opts.rgbuthr == 0
    opts.rgbuthr = sqrt(eps);
else
    opts.rgbuthr = abs(opts.rgbuthr);
end
if ~isfield(opts, 'space') || ~ischar(opts.space) || ~any(strcmpi(opts.space(:)', {'bvs', 'tal'}))
    opts.space = 'tal';
else
    opts.space = lower(opts.space(:)');
end
opts.step = [];
if isempty(opts.frame)
    if opts.space(1) == 't'
        opts.frame = [128, 128, 128; -127.9999, -127.9999, -127.9999];
    else
        opts.frame = [0, 0, 0; 255.9999, 255.9999, 255.9999];
    end
    opts.step = 256 / w;
end
if any(diff(abs(diff(opts.frame))))
    error('neuroelf:xff:badArgument', 'The frame must be cubic in size.');
end
if opts.orient == 'n'
    if opts.space(1) == 'b'
        opts.frame(:, 3) = opts.frame([2, 1], 3);
    else
        opts.frame(:, 1) = opts.frame([2, 1], 1);
    end
end
if isempty(opts.step)
    stepi = abs(diff(opts.frame(:, 1)));
    if all(abs(stepi - round(stepi)) < 0.001)
        opts.step = round(abs(diff(opts.frame(:, 1)))) / w;
    else
        opts.step = abs(diff(opts.frame(:, 1))) / (w - 0.0001);
    end
    d = sort(opts.frame(:, 1));
    if numel(d(1):opts.step:d(2)) > w
        opts.frame(2, 1) = opts.frame(2, 1) - 0.1 * opts.step * sign(diff(opts.frame(:, 1)));
    end
    d = sort(opts.frame(:, 2));
    if numel(d(1):opts.step:d(2)) > w
        opts.frame(2, 2) = opts.frame(2, 2) - 0.1 * opts.step * sign(diff(opts.frame(:, 2)));
    end
    d = sort(opts.frame(:, 3));
    if numel(d(1):opts.step:d(2)) > w
        opts.frame(2, 3) = opts.frame(2, 3) - 0.1 * opts.step * sign(diff(opts.frame(:, 3)));
    end
end
if ~isfield(opts, 'type') || ~ischar(opts.type) || ...
   ~any(strcmpi(opts.type(:)', {'d', 'dsq', 'g', 'gray', 'grey', 'l', 'ldsq', 'lsq', 'r', 'rgb'}))
    if isfield(rtv, 'DisplayType') && ischar(rtv.DisplayType) && ...
        any(strcmpi(rtv.DisplayType(:)', {'d', 'dsq', 'g', 'gray', 'grey', 'l', 'ldsq', 'lsq', 'r', 'rgb'}))
        opts.type = lower(rtv.DisplayType(1));
    else
        opts.type = 'g';
    end
else
    opts.type = lower(opts.type(1));
end
if ~isfield(opts, 'xscale') || ~isa(opts.xscale, 'double') || ~any(numel(opts.xscale) == [1, 2]) || ...
    any(isinf(opts.xscale) | isnan(opts.xscale)) || ...
   (numel(opts.xscale) == 1 && opts.xscale <= 0) || ...
   (numel(opts.xscale) == 2 && opts.xscale(1) == opts.xscale(2))
    opts.xscale = [];
end

% for switch on filetype
ft = lower(xo.S.Extensions{1});

% special HDRs
if strcmp(ft, 'hdr')

    % RGB/RGBA
    if bc.ImgDim.DataType == 128
        opts.type = '3';
    elseif bc.ImgDim.DataType == 2304
        opts.type = '4';

    elseif bc.ImgDim.Dim(1) > 4 && bc.ImgDim.Dim(5) == 1 && bc.ImgDim.Dim(6) > 1
        opts.type = 'x';
    end
end

% coloring of gray-scale values
if any(opts.type == 'dgl') && (~isequal(opts.gcolhigh, [255, 255, 255]) || ...
    ~isequal(opts.gcollow, [0, 0, 0]))

    % LUT coloring
    if ~isempty(opts.gcollut)

        gcol = opts.gcollut;

    % RGB coloring
    elseif opts.gcolblend == 'r'

        % create 256 colors
        hcol = opts.gcolhigh;
        lcol = opts.gcollow;
        dcol = (hcol - lcol) ./ 255;
        gcol = uint8(round([(lcol(1):dcol(1):hcol(1))', ...
            (lcol(2):dcol(2):hcol(2))', (lcol(3):dcol(3):hcol(3))']));

    % HSV coloring
    else

        % convert first
        hcol = ne_methods.hsvconv(opts.gcolhigh, 2);
        lcol = ne_methods.hsvconv(opts.gcollow, 2);
        if hcol(3) == 0
            hcol(1) = lcol(1);
        end
        if abs(hcol(1) - lcol(1)) > 0.5
            if hcol(1) > lcol(1)
                lcol(1) = lcol(1) + 1;
            else
                hcol(1) = hcol(1) + 1;
            end
        end
        dcol = (hcol - lcol) ./ 255;
        gcol = ne_methods.hsvconv(mod([(lcol(1):dcol(1):hcol(1))', ...
            (lcol(2):dcol(2):hcol(2))', (lcol(3):dcol(3):hcol(3))'], 1));
    end
else
    gcol = [];
end

% ensure that the scaling window is known
if any(strcmp(ft, {'amr', 'dmr', 'fmr', 'hdr', 'head', 'mgh', 'vmr', 'vtc'})) && ...
   (~isfield(rtv, 'ScalingWindow') || numel(rtv.ScalingWindow) ~= 2) && ...
   (~strcmp(ft, 'hdr') || ~isfield(rtv, 'StatsObject') || ~islogical(rtv.StatsObject) || ...
    numel(rtv.StatsObject) ~= 1 || ~rtv.StatsObject)
    aft_SetScalingWindow(xo);
    bc = xo.C;
    rtv = bc.RunTimeVars;
end

% defaults for non-VMPs
alpha = opts.rgbalpha;
lut = opts.rgbcol;
tls = opts.rgbctails;
lthr = double(opts.rgblthr);
uthr = double(opts.rgbuthr);

% get number of volumes
dsc = false;
switch (ft)
    case {'amr', 'map', 'mgh', 'msk', 'srf', 'vmr'}
        nvol = 1;
    case 'ava'
        nvol = numel(fieldnames(bc.Maps)) + ...
            (size(bc.Maps.CellMeans, ndims(bc.Maps.CellMeans)) - 1);
    case {'cmp', 'vmp'}
        nvol = numel(bc.Map);
    case 'hdr'
        nvol = size(bc.VoxelData, 4);
        if any([2, 4, 8, 130, 132, 136, 256, 512, 768] == bc.ImgDim.DataType)
            dscf = bc.ImgDim.ScalingSlope;
            if isinf(dscf) || isnan(dscf) || dscf == 0
                dscf = 1;
            end
            dsci = bc.ImgDim.ScalingIntercept;
            if isinf(dsci) || isnan(dsci)
                dsci = 0;
            end
            if dscf ~= 1 || dsci ~= 0
                dsc = true;
            end
        end
    case 'nlf'
        xd = bc.DimMeaning;
        xd = ne_methods.findfirst(xd ~= 'x' & xd ~= 'y' & xd ~= 'z');
        if isempty(xd)
            nvol = 1;
        else
            nvol = size(bc.Data, xd(1));
        end
    case 'head'
        nvol = numel(bc.Brick);
        if ~isempty(bc.Brick) && ~isempty(regexpi(class(bc.Brick(1).Data), 'int')) && ...
            bc.Brick(1).ScalingFactor ~= 1
            dscf = {bc.Brick.ScalingFactor};
            dsci = repmat({0}, size(dscf));
            dsc = true;
        end
    case 'glm'
        nvol = numel(bc.Predictor);
    case 'vtc'
        if ~isfield(bc.RunTimeVars, 'AvgVTC') || ~bc.RunTimeVars.AvgVTC
            nvol = size(bc.VTCData, 1);
        else
            nvol = numel(bc.RunTimeVars.ConditionNames);
            if ~isfield(bc.RunTimeVars, 'SubMapVol') || ...
               ~isa(bc.RunTimeVars.SubMapVol, 'double') || ...
                numel(bc.RunTimeVars.SubMapVol) ~= 1
                bc.RunTimeVars.SubMapVol = 1;
            end
        end
    case 'fmr'
        nvol = size(bc.Slice(1).STCData, 3);
    case 'vdw'
        nvol = size(bc.VDWData, 1);
    case 'dmr'
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
    otherwise
        error('neuroelf:xff:objectTypeUnsupported', ...
            'Unsupported object type for ::SliceToTransimg method.');
end

% thresholds update flag
thrc = false;

% for RGB mode and non-VMPs, determine thresholds
if opts.type == 'r' && ~any(strcmp(ft, {'ava', 'cmp', 'glm', 'hdr', 'head', 'vmp'}))

    % no thresholds handle yet
    if ~isfield(bc.RunTimeVars, 'Thresholds') || ...
       ~isequal(size(bc.RunTimeVars.Thresholds), [nvol, 5])

        % then create and store
        bc.RunTimeVars.Thresholds = nan(nvol, 5);
        xo.C = bc;
    end

    % get thresholds
    thrh = bc.RunTimeVars.Thresholds;
end

% default options
if isempty(opts.mapvol)
    mapvol = 1:nvol;
elseif any(opts.mapvol(:) > nvol)
    error('neuroelf:xff:badArgument', 'Requested mapvol option out of bounds.');
else
    mapvol = round(opts.mapvol(:));
end
if ~isempty(opts.mapvola)
    opts.mapvola(opts.mapvola > nvol) = Inf;
elseif any(strcmp(ft, {'cmp', 'vmp'})) && isfield(bc.Map, 'RunTimeVars') && ...
    isstruct(bc.Map(opts.mapvol(1)).RunTimeVars) && ...
    isfield(bc.Map(opts.mapvol(1)).RunTimeVars, 'AlphaMap')
    mapvolai = {bc.Map(opts.mapvol).RunTimeVars};
    for sc = 1:numel(mapvolai)
        if isstruct(mapvolai{sc}) && isfield(mapvolai{sc}, 'AlphaMap') && ...
            isa(mapvolai{sc}.AlphaMap, 'double') && numel(mapvolai{sc}.AlphaMap) == 1
            mapvolai{sc} = mapvolai{sc}.AlphaMap;
        else
            mapvolai{sc} = Inf;
        end
    end
    opts.mapvola = cat(2, mapvolai{:});
end
if numel(opts.mapvola) == 1 && numel(mapvol) > 1
    opts.mapvola = opts.mapvola .* ones(1, numel(mapvol));
elseif numel(opts.mapvola) ~= numel(mapvol)
    opts.mapvola = [];
end
if numel(opts.layers) == 1 && numel(mapvol) > 1
    opts.layers = opts.layers:(opts.layers + numel(mapvol) - 1);
end
if numel(opts.layers) ~= numel(mapvol)
    error('neuroelf:xff:badArgument', 'Number of layers and sampled maps must match.');
end
sval = zeros(size(mapvol));
show3 = (opts.dir == 'a');
oval = zeros(size(mapvol));

% for colored images, keep track of bars
if ~any(opts.type == 'dgl') && ~isempty(opts.rgbbars)
    bars = cell(2, numel(sval));
else
    bars = {};
end

% add gradient option
if any(opts.type == 'dl')
    opts.gradient = true;
    opts.gradtype = ['1' opts.type];
end

% scaling values
if dsc && ~iscell(dscf)
    dscf = repmat({dscf}, 1, nvol);
    dsci = repmat({dsci}, 1, nvol);
end

% iterate over maps to sample
l = opts.layers;
for sc = 1:numel(sval)

    % set mapvol
    opts.mapvol = mapvol(sc);

    % three directions
    if show3

        % sample data
        [sag, cor, tra, scrd, odata, oalpha] = aft_SliceData3D(xo, crd, opts);
        odatasz = size(odata);
        if numel(odatasz) < 3
            odatasz(3) = 1;
        end

        % and get value for sampling
        scrdr = round(scrd);
        if all(scrdr > 0 & scrdr <= odatasz(1:3))
            if dsc
                sval(sc, 1:size(odata, 5)) = dsci{mapvol(sc)} + dscf{mapvol(sc)} .* ...
                    squeeze(odata(scrdr(1), scrdr(2), scrdr(3), :, :));
            else
                sval(sc, 1:size(odata, 5)) = squeeze(odata(scrdr(1), scrdr(2), scrdr(3), :, :));
            end
        end

    % one direction only, but follow same logic
    else
        [sag, scrd, odata, oalpha] = aft_SliceData3D(xo, crd, opts);
        odatasz = size(odata);
        scrdr = round(scrd);
        if all(scrdr > 0 & scrdr <= odatasz(1:3))
            if dsc
                sval(sc, 1:size(odata, 5)) = dsci{mapvol(sc)} + dscf{mapvol(sc)} .* ...
                    squeeze(odata(scrdr(1), scrdr(2), scrdr(3), :, :));
            else
                sval(sc, 1:size(odata, 5)) = squeeze(odata(scrdr(1), scrdr(2), scrdr(3), :, :));
            end
        end
    end

    % for alpha-maps set
    if any(strcmp(ft, {'cmp', 'vmp'})) && ~isempty(opts.mapvola)

        % valid for this particular map
        if ~isinf(opts.mapvola(sc))

            % follow same logic
            opts.mapvol = opts.mapvola(sc);
            atype = bc.Map(opts.mapvol).Type;
            alot = abs(bc.Map(opts.mapvol).LowerThreshold);
            ahit = abs(bc.Map(opts.mapvol).UpperThreshold);
            ahilot = 1 / max(sqrt(eps), abs(ahit - alot));
            if show3
                [saga, cora, traa] = aft_SliceData3D(xo, crd, opts);
                if atype == 1
                    saga = abs(saga);
                    cora = abs(cora);
                    traa = abs(traa);
                end
                saga(sag == 0) = 0;
                cora(cor == 0) = 0;
                traa(tra == 0) = 0;
                oalpha = { ...
                    ne_methods.limitrangec(ahilot .* (saga - alot), 0, 1, 0), ...
                    ne_methods.limitrangec(ahilot .* (cora - alot), 0, 1, 0), ...
                    ne_methods.limitrangec(ahilot .* (traa - alot), 0, 1, 0)};
            else
                saga = aft_SliceData3D(xo, crd, opts);
                if atype == 1
                    saga = abs(saga);
                end
                saga(sag == 0) = 0;
                oalpha = {ne_methods.limitrangec(ahilot .* (saga - alot), 0, 1, 0), [], []};
            end

        % invalid map
        else
            oalpha = 1;
        end
    end

    % for gray scale
    if any(opts.type == 'dgl')

        % scaling
        if isfield(bc.RunTimeVars, 'ScalingWindow') && ...
            numel(bc.RunTimeVars.ScalingWindow) == 2

            % window
            if strcmp(ft, 'vmr') && isfield(opts, 'v16') && islogical(opts.v16) && ...
                numel(opts.v16) == 1 && opts.v16 && ...
                isfield(bc.RunTimeVars, 'ScalingWindow16') && ...
                numel(bc.RunTimeVars.ScalingWindow16) == 2
                win = bc.RunTimeVars.ScalingWindow16;
            else
                opts.v16 = false;
                win = bc.RunTimeVars.ScalingWindow;
            end
            wind = win(2) - win(1);
            win = win(1);

            % adapt window for log-gradient
            if opts.type == 'l'
                win = max(0, win);
                wind = log(1 + wind);
            end

            % for VMR's, only values below 226
            if strcmp(ft, 'vmr') && ~opts.v16

                % round values to begin with, and values > 225 indexed
                scol = sag >= 225.5;
                scolv = round(min(30, max(1, round(sag(scol)) - 225)));
                if show3
                    ccol = cor >= 225.5;
                    ccolv = round(min(30, max(1, round(cor(ccol)) - 225)));
                    tcol = tra >= 225.5;
                    tcolv = round(min(30, max(1, round(tra(tcol)) - 225)));
                end
            end

            % re-compute
            sag = uint8(floor(min(255, max(0, ...
                255 .* ((1 / wind) .* (double(sag) - win))))));
            if show3
                cor = uint8(floor(min(255, max(0, ...
                    255 .* ((1 / wind) .* (double(cor) - win))))));
                tra = uint8(floor(min(255, max(0, ...
                    255 .* ((1 / wind) .* (double(tra) - win))))));
            end

            % for VMR's replace indexed colors
            if strcmp(ft, 'vmr') && ~opts.v16

                % get LUT for values >= 226
                vmrlut = stt_vmrlut;

                % and create XxYx3 uint8 image already (truecolor)
                c_g = sag;
                c_b = sag;
                sag(scol) = vmrlut(scolv, 1);
                c_g(scol) = vmrlut(scolv, 2);
                c_b(scol) = vmrlut(scolv, 3);
                sag = cat(3, sag, c_g, c_b);

                % but still color it?
                if ~isempty(gcol)

                    % get indices
                    scol = find(~scol(:));
                    sagi = uint16(1) + uint16(sag(scol));

                    % now fill planes
                    sagsz = size(sag, 1) * size(sag, 2);
                    sag(scol) = gcol(sagi, 1);
                    sag(scol + sagsz) = gcol(sagi, 2);
                    sag(scol + 2 .* sagsz) = gcol(sagi, 3);
                end

                % for 3-slice sampling
                if show3
                    c_g = cor;
                    c_b = cor;
                    cor(ccol) = vmrlut(ccolv, 1);
                    c_g(ccol) = vmrlut(ccolv, 2);
                    c_b(ccol) = vmrlut(ccolv, 3);
                    cor = cat(3, cor, c_g, c_b);
                    c_g = tra;
                    c_b = tra;
                    tra(tcol) = vmrlut(tcolv, 1);
                    c_g(tcol) = vmrlut(tcolv, 2);
                    c_b(tcol) = vmrlut(tcolv, 3);
                    tra = cat(3, tra, c_g, c_b);
                    if ~isempty(gcol)
                        ccol = find(~ccol(:));
                        cori = uint16(1) + uint16(cor(ccol));
                        corsz = size(cor, 1) * size(cor, 2);
                        cor(ccol) = gcol(cori, 1);
                        cor(ccol + corsz) = gcol(cori, 2);
                        cor(ccol + 2 .* corsz) = gcol(cori, 3);
                        tcol = find(~tcol(:));
                        trai = uint16(1) + uint16(tra(tcol));
                        trasz = size(tra, 1) * size(tra, 2);
                        tra(tcol) = gcol(trai, 1);
                        tra(tcol + trasz) = gcol(trai, 2);
                        tra(tcol + 2 .* trasz) = gcol(trai, 3);
                    end
                end

            % for AVG-VTCs
            elseif strcmp(ft, 'atc') && rtv.AvgVTC

                % condition color tables set
                if isfield(rtv, 'ConditionColorTables') && ...
                    isequal(size(rtv.ConditionColorTables), ...
                    [rtv.NrOfConditions, 8, 256])

                    % get table
                    cctable = rtv.ConditionColorTables;
                else

                    % compute table
                    cctable = uint8(zeros([rtv.NrOfConditions, 8, 256]));
                    for ccc = 1:rtv.NrOfConditions

                        % derive colors and mix
                    end

                    % store table
                    rtv.ConditionColorTables = cctable;
                    thrc = true;
                end

            % colorize
            elseif ~isempty(gcol)
                sag = uint16(1) + uint16(sag);
                ssz = size(sag);
                sag = cat(3, reshape(gcol(sag, 1), ssz), ...
                    reshape(gcol(sag, 2), ssz), reshape(gcol(sag, 3), ssz));
                if show3
                    cor = uint16(1) + uint16(cor);
                    csz = size(cor);
                    cor = cat(3, reshape(gcol(cor, 1), csz), ...
                        reshape(gcol(cor, 2), csz), reshape(gcol(cor, 3), csz));
                    tra = uint16(1) + uint16(tra);
                    tsz = size(tra);
                    tra = cat(3, reshape(gcol(tra, 1), tsz), ...
                        reshape(gcol(tra, 2), tsz), reshape(gcol(tra, 3), tsz));
                end
            end
        end

        % make valid for display
        if isinf(opts.grayalpha)
            saga = single(1);
        else
            saga = (opts.grayalpha / 255) .* single(max(sag, [], 3));
        end
        setlayer(ti(1), l(sc), sag, saga);
        if show3
            if isinf(opts.grayalpha)
                cora = saga;
                traa = saga;
            else
                cora = (opts.grayalpha / 255) .* single(max(cor, [], 3));
                traa = (opts.grayalpha / 255) .* single(max(tra, [], 3));
            end
            setlayer(ti(2), l(sc), cor, cora);
            setlayer(ti(3), l(sc), tra, traa);
        end

    % for RGB display
    elseif opts.type == 'r'

        % override settings from VMP?
        if any(strcmp(ft, {'cmp', 'vmp'}))
            if isempty(opts.rgbcol)
                lut = bc.Map(mapvol(sc)).OverlayColors;
                if isempty(lut)
                    lut = stt_vmrlut;
                end
            end
            if isempty(opts.rgbctails)
                tls = bc.Map(mapvol(sc)).ShowPositiveNegativeFlag;
            end
            if isempty(opts.rgblthr)
                lthr = double(abs(bc.Map(mapvol(sc)).LowerThreshold));
            end
            if isempty(opts.rgbuthr)
                uthr = double(abs(bc.Map(mapvol(sc)).UpperThreshold));
            end
            if iscell(oalpha)
                mapalpha = oalpha;
            elseif bc.Map(mapvol(sc)).TransColorFactor ~= 1
                mapalpha = min(1, alpha * bc.Map(mapvol(sc)).TransColorFactor);
            else
                mapalpha = alpha;
            end

        % for "VMP-like" objects
        elseif any(strcmp(ft, {'ava', 'glm', 'hdr', 'head', 'vtc'}))
            if isempty(opts.rgbcol)
                lut = bc.RunTimeVars.Map(mapvol(sc)).OverlayColors;
                if isempty(lut)
                    lut = stt_vmrlut;
                end
            end
            if isempty(opts.rgbctails)
                tls = bc.RunTimeVars.Map(mapvol(sc)).ShowPositiveNegativeFlag;
            end
            if iscell(oalpha)
                mapalpha = oalpha;
            elseif bc.RunTimeVars.Map(mapvol(sc)).TransColorFactor ~= 1
                mapalpha = min(1, alpha * bc.RunTimeVars.Map(mapvol(sc)).TransColorFactor);
            else
                mapalpha = alpha;
            end

            % auto-thresholding
            if bc.RunTimeVars.Map(mapvol(sc)).Type == 15 && ...
               (isempty(bc.RunTimeVars.Map(mapvol(sc)).LowerThreshold) || ...
                isempty(bc.RunTimeVars.Map(mapvol(sc)).UpperThreshold))

                % winsorize data first
                w = ne_methods.winsorize(ne_methods.limitrangec(odata(:), -1e6, 1e6, 0));

                % then remove zeros
                w(w == 0) = [];

                % and compute range
                w = ne_methods.minmaxmean(abs(w), 5);
                seps = sqrt(eps);
                bc.RunTimeVars.Map(mapvol(sc)).LowerThreshold = ...
                    max(seps, w(3) + 0.25 * sqrt(w(6)) - w(2) * seps);
                bc.RunTimeVars.Map(mapvol(sc)).UpperThreshold = ...
                    w(3) + 3 * sqrt(w(6)) + eps * w(2);
                thrc = true;
            end

            % thresholding
            if isempty(opts.rgblthr)
                lthr = double(abs(bc.RunTimeVars.Map(mapvol(sc)).LowerThreshold));
            end
            if isempty(opts.rgbuthr)
                uthr = double(abs(bc.RunTimeVars.Map(mapvol(sc)).UpperThreshold));
            end

        % for non VMP-files
        else

            % get current thresholds in RunTimeVars
            thrr = thrh(mapvol(sc), :);

            % check cluster size entries
            if any(isnan(thrr(4:5)))
                thrr(4:5) = [1, 0];
                thrh(mapvol(sc), :) = thrr;
                thrc = true;
            end

            % empty colors
            if isempty(opts.rgbcol)

                % colors in object?
                if isfield(bc.RunTimeVars, 'OverlayColors') && ...
                    sc <= numel(bc.RunTimeVars.OverlayColors) && ...
                   ~isempty(bc.RunTimeVars.OverlayColors{sc})
                    lut = bc.RunTimeVars.OverlayColors{sc};
                else
                    lut = stt_vmrlut;
                end
            end

            % no tails selected
            if isempty(opts.rgbctails)
                if isnan(thrr(3))
                    thrr(3) = 3;
                    thrh(mapvol(sc), 3) = 3;
                    thrc = true;
                end
                tls = thrr(3);
            end

            % thresholding
            if isempty(opts.rgblthr) || isempty(opts.rgbuthr)

                % no thresholds defined
                if any(isnan(thrr(1:2)))

                    % winsorize data first
                    w = ne_methods.winsorize(ne_methods.limitrangec(odata(:), -1e6, 1e6, 0));

                    % then remove zeros
                    w(w == 0) = [];

                    % and compute range
                    w = ne_methods.minmaxmean(abs(w), 5);
                    seps = sqrt(eps);
                    thrr(1:2) = [max(seps, w(3) + 0.25 * sqrt(w(6)) - w(2) * seps), ...
                        w(3) + 3 * sqrt(w(6)) + eps * w(2)];
                    thrh(mapvol(sc), 1:2) = thrr(1:2);
                    thrc = true;
                end

                % get threshold
                lthr = thrr(1);
                uthr = thrr(2);
            end
        end

        % threshold images
        sag = ne_methods.threshmapc(double(sag), lthr, uthr, tls);
        if show3
            cor = ne_methods.threshmapc(double(cor), lthr, uthr, tls);
            tra = ne_methods.threshmapc(double(tra), lthr, uthr, tls);
        end

        % full-slice alpha information already
        if iscell(mapalpha)

            % mask alpha
            sagc = mapalpha{1};
            sagc(sag == 0) = 0;
            if show3
                corc = mapalpha{2};
                corc(cor == 0) = 0;
                trac = mapalpha{3};
                trac(tra == 0) = 0;
            end

        % no transparency
        elseif mapalpha >= 1

            % simply show all values ~= 0
            sagc = double(sag ~= 0);
            if show3
                corc = double(cor ~= 0);
                trac = double(tra ~= 0);
            end

        % for positive alpha, set values ~= 0 to that value
        elseif mapalpha >= 0
            sagc = mapalpha .* double(sag ~= 0);
            if show3
                corc = mapalpha .* double(cor ~= 0);
                trac = mapalpha .* double(tra ~= 0);
            end

        % and for negative values
        else

            % multiply with (absolute) stats value
            sagc = min(1, -mapalpha .* abs(sag));
            if show3
                corc = min(1, -mapalpha .* abs(cor));
                trac = min(1, -mapalpha .* abs(tra));
            end
        end

        % colorbars
        if ~isempty(bars)

            % store bars info for positive tail
            if any([1, 3] == tls)
                bars{1, sc} = lut(1:ceil(0.5 .* size(lut, 1)), :);
            end

            % store bars info for negative tail
            if tls > 1
                bars{2, sc} = lut(1+ceil(0.5 .* size(lut, 1)):end, :);
            end
        end

        % set in layer(s)
        setlayer(ti(1), l(sc), ne_methods.threshlutc(sag, lut), sagc);
        if show3
            setlayer(ti(2), l(sc), ne_methods.threshlutc(cor, lut), corc);
            setlayer(ti(3), l(sc), ne_methods.threshlutc(tra, lut), trac);
        end

    % for "real" RGB display
    elseif any(opts.type == '34')

        % scaling
        if isfield(bc.RunTimeVars, 'ScalingWindow') && ...
            numel(bc.RunTimeVars.ScalingWindow) == 2

            % window
            win = bc.RunTimeVars.ScalingWindow;
            wind = win(2) - win(1);
            win = win(1);

            % apply scaling
            sag = uint8(floor(min(255, max(0, ...
                255 .* ((1 / wind) .* (sag - win))))));
            if show3
                cor = uint8(floor(min(255, max(0, ...
                    255 .* ((1 / wind) .* (cor - win))))));
                tra = uint8(floor(min(255, max(0, ...
                    255 .* ((1 / wind) .* (tra - win))))));
            end
        end

        % with alpha
        if size(sag, 3) > 3
            setlayer(ti(1), l(sc), sag(:, :, 1:3), sag(:, :, 4));
            if show3
                setlayer(ti(2), l(sc), cor(:, :, 1:3), cor(:, :, 4));
                setlayer(ti(3), l(sc), tra(:, :, 1:3), tra(:, :, 4));
            end

        % standard RGB
        elseif size(sag, 3) > 2
            setlayer(ti(1), l(sc), sag(:, :, 1:3), 1);
            if show3
                setlayer(ti(2), l(sc), cor(:, :, 1:3), 1);
                setlayer(ti(3), l(sc), tra(:, :, 1:3), 1);
            end

        % something went wrong -> grayscale
        else
            setlayer(ti(1), l(sc), sag(:, :, [1, 1, 1]), 1);
            if show3
                setlayer(ti(2), l(sc), cor(:, :, [1, 1, 1]), 1);
                setlayer(ti(3), l(sc), tra(:, :, [1, 1, 1]), 1);
            end
        end

    % multi-vol display
    else

        % number of volumes
        ndvol = size(sag, 3);
        switch (ndvol)

            % BW+alpha
            case 2
                sagc = sag(:, :, 2);
                sag = sag(:, :, [1, 1, 1]);
                if show3
                    corc = cor(:, :, 2);
                    cor = cor(:, :, [1, 1, 1]);
                    trac = tra(:, :, 2);
                    tra = tra(:, :, [1, 1, 1]);
                end

            % RGB
            case 3
                sagc = 1;
                if show3
                    corc = 1;
                    trac = 1;
                end

            % RGBA
            case 4
                sagc = sag(:, :, 4);
                sag = sag(:, :, 1:3);
                if show3
                    corc = cor(:, :, 4);
                    cor = cor(:, :, 1:3);
                    trac = tra(:, :, 4);
                    tra = tra(:, :, 1:3);
                end

            % Tensor?
            % case 6

            % otherwise
            otherwise
                sagc = 1;
                sag = sum(sag, 3) ./ ndvol;
                if show3
                    corc = 1;
                    cor = repmat(sum(cor, 3) ./ ndvol, [1, 1, 3]);
                    trac = 1;
                    tra = repmat(sum(tra, 3) ./ ndvol, [1, 1, 3]);
                end
        end

        % re-scale
        if show3
            if isempty(opts.xscale)
                mmresc = ne_methods.minmaxmean([sag(:); cor(:); tra(:)]);
            elseif numel(opts.xscale) == 2
                mmresc = opts.xscale;
            else
                mmresc = opts.xscale .* [-1, 1];
            end                
            mmfac = (255 / ((mmresc(2) - mmresc(1)) + eps));
            sag = mmfac .* (sag - mmresc(1));
            cor = mmfac .* (cor - mmresc(1));
            tra = mmfac .* (tra - mmresc(1));
            setlayer(ti(1), l(sc), sag, sagc);
            setlayer(ti(2), l(sc), cor, corc);
            setlayer(ti(3), l(sc), tra, trac);
        else
            if isempty(opts.xscale)
                mmresc = ne_methods.minmaxmean(sag);
            elseif numel(opts.xscale) == 2
                mmresc = opts.xscale;
            else
                mmresc = opts.xscale .* [-1, 1];
            end                
            sag = (255 / ((mmresc(2) - mmresc(1)) + eps)) .* (sag - mmresc(1));
            setlayer(ti(1), l(sc), sag, sagc);
        end

        % set to output
        if show3
        end

    end
end

% show bars
if ~isempty(bars)

    % which bars
    bars = bars(~cellfun('isempty', bars(:)));

    % barpos
    f1 = opts.rgbbars(1);
    f2 = opts.rgbbars(2);
    t1 = opts.rgbbars(3);
    t2 = opts.rgbbars(4);

    % bar height and width
    bh = 1 + t2 - f2;
    bw = (1 + t1 - f1) / numel(bars);

    % render correct layer
    if show3
        ti = ti(2);
        ldata = ne_methods.threshlutc(cor, lut);
        ldataa = corc;
    else
        ti = ti(1);
        ldata = ne_methods.threshlutc(sag, lut);
        ldataa = sagc;
    end
    l = l(numel(sval));

    % iterate over bars
    for sc = 1:numel(bars)

        % get bardata
        bd = bars{sc};
        bd = reshape(bd, [size(bd, 1), 1, size(bd, 2)]);

        % scaling factor
        bs = size(bd, 1);
        bsf = (1 / (bh - 1)) * (bs - 1);

        % target position
        tf = floor(1 + f1 + (sc - 1) * bw);
        tt = tf + floor(bw) - 2;
        if tt < tf
            tt = tf;
        end
        tw = 1 + tt - tf;

        % just one line of color?
        if bsf == 0

            % copy data
            ldata(f2:t2, tf:tt, :) = repmat(bd, bh, tw);
        else

            % interpolate data
            brgb = ne_methods.flexinterpn(bd, ...
                [Inf, Inf, Inf; bs, 1, 1; -bsf, 1, 1; 1, 1, 3]);

            % store
            if tw > 1
                ldata(f2:t2, tf:tt, :) = repmat(round(brgb), 1, tw);
            else
                ldata(f2:t2, tf:tt, :) = round(brgb);
            end
        end

        % set alpha to 1
        ldataa(f2:t2, tf:tt) = 1;
    end

    % update layer
    setlayer(ti, l, ldata, ldataa);
end

% update thresholds
if thrc
    if ~any(strcmp(ft, {'cmp', 'glm', 'hdr', 'head', 'vmp'}))
        bc.RunTimeVars.Thresholds = thrh;
    end
    xo.C = bc;
end
