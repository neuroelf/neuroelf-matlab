function nvox = aft_IntensityMask(xo, mn, mx, opts)
% AFT::IntensityMask  - mask a data object based on intensity values
%
% FORMAT:       [nvox = ] obj.IntensityMask(minvalue [, maxvalue [, opts]])
%
% Input fields:
%
%       minvalue    minimum intensity value
%       maxvalue    maximum intensity  (default: Inf)
%       opts        optional settings
%        .meantype  how to build average for time-variant data, one of
%                    - 1xN double (use regular mean over those timepoints)
%                    - 'bv' (BrainVoyager QX style: 0.5 * d([nvol/2, nvol]))
%                    - {'mean'} (use regular mean over all timepoints)
%                    - 'robmean' (use robust mean over all timepoints)
%        .rangetype how are minimum and maximum values interpreted, one of
%                    - {'absolute'} (direct, as values appear)
%                    - 'hist' / 'percentile' (histogram-based [0 .. 1])
%                    - 'relative' (product with spatial mean)
%                    - 'std' (relative to mean +/- stds)
%
% No output fields. Alters object.
%
% Note: if only minvalue is given and negative, rangetype is set to
%       'std' with the range of +/-abs(minvalue).
%
% TYPES: AMR, FMR, HDR, HEAD, MGH, MTC, NLF, VMR, VTC
%
% Using: findfirst, robustmean.

% Version:  v1.1
% Build:    16031615
% Date:     Mar-16 2016, 3:38 PM EST
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

% neuroelf library
global ne_methods;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ...
   ~xffisobject(xo, true, {'amr', 'fmr', 'hdr', 'head', 'mgh', 'mtc', 'nlf', 'vmr', 'vtc'}) || ...
   ~isa(mn, 'double') || numel(mn) ~= 1 || isinf(mn) || isnan(mn)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 3 || ~isa(mx, 'double') || numel(mx) ~= 1 || isinf(mx) || isnan(mx)
    if mn < 0
        mx = -mn;
    else
        mx = Inf;
    end
end
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
bc = xo.C;
type = lower(xo.S.Extensions{1});

% before testing further options, get data and number of volumes
switch (type)
    case 'amr'
        data = bc.Slice(1).AMRData(:, :);
        data(1, 1, numel(bc.Slice)) = 0;
        for slc = 1:numel(bc.Slice)
            data(:, :, slc) = bc.Slice(slc).AMRData(:, :);
        end
        tdim = [];
    case 'fmr'
        if numel(bc.Slice) == 1
            data = bc.Slice.STCData;
        else
            data = bc.Slice(1).STCData(:, :, :);
            data(1, 1, 1, numel(bc.Slice)) = 0;
            for slc = 1:numel(bc.Slice)
                data(:, :, :, slc) = bc.Slice(slc).STCData(:, :, :);
            end
        end
        tdim = 3;
    case 'hdr'

        % reject complex datatypes
        if any(bc.ImgDim.DataType == [32, 128, 1792, 2048, 2304])
            error('neuroelf:xff:unsupported', ...
                'Intensity masking not supported for complex datatypes.');
        end

        % get data
        data = bc.VoxelData;
        tdim = 4;

        % scaling
        if any([2, 4, 8, 130, 132, 136, 256, 512, 768] == bc.ImgDim.DataType) && ...
           (bc.ImgDim.ScalingIntercept ~= 0 || all([0, 1] ~= bc.ImgDim.ScalingSlope))
            if bc.ImgDim.ScalingSlope ~= 0
                data = bc.ImgDim.ScalingIntercept + bc.ImgDim.ScalingSlope .* double(data);
            else
                data = bc.ImgDim.ScalingIntercept + double(data);
            end
        end

    case 'head'
        data = bc.Brick(1).Data(:, :, :);
        if numel(bc.Brick) > 1
            data(1, 1, 1, numel(bc.Brick)) = 0;
            for brc = 2:numel(bc.Brick)
                data(:, :, :, brc) = bc.Brick(brc).Data(:, :, :);
            end
        end
        tdim = 4;
    case 'mgh'
        data = bc.MGHData;
        tdim = [];
    case 'mtc'
        data = bc.MTCData;
        tdim = 1;
    case 'nlf'
        data = bc.Data;
        tdim = find(bc.DimMeaning == 't');
    case 'vmr'
        data = bc.VMRData;
        tdim = [];
    case 'vtc'
        data = bc.VTCData;
        tdim = 1;
end
if istransio(data)
    data = resolve(data);
end
if ~isempty(tdim)
    nvol = size(data, tdim);
else
    nvol = 1;
end

% now continue further check
if ~isfield(opts, 'meantype') || (~ischar(opts.meantype) && ...
    (~isa(opts.meantype, 'double') || isempty(opts.meantype) || ...
     any(isinf(opts.meantype(:)) | isnan(opts.meantype(:)) | opts.meantype(:) < 1)))
    opts.meantype = 'mean';
end
if isa(opts.meantype, 'double')
    opts.meantype = unique(min(fix(opts.meantype(:)'), nvol));
else
    if ~any(strcmpi(opts.meantype(:)', {'b', 'bv', 'm', 'mean', 'r', 'robmean'}))
        opts.meantype = 'm';
    else
        opts.meantype = lower(opts.meantype(1));
    end
end
if ~isfield(opts, 'rangetype') || ~ischar(opts.rangetype) || ...
   ~any(strcmpi(opts.rangetype(:)', {'a', 'abs', 'absolute', 'h', 'hist', ...
    'p', 'percentile', 'r', 'relative', 's', 'std'}))
    if nargin > 3 || mn >= 0
        opts.rangetype = 'a';
    else
        opts.rangetype = 's';
    end
else
    opts.rangetype = lower(opts.rangetype(1));
end

% get required content
sr = {':'};
sr = sr(ones(1, ndims(data)));
if ~ischar(opts.meantype)
    sr{tdim} = opts.meantype(:)';
    usedata = data(sr{:});
    opts.meantype = 'm';
    nvol = size(usedata, tdim);
else
    usedata = data;
end

% and average
if nvol > 1

    % with required algorithm
    switch (opts.meantype)

        % BrainVoyager QX style
        case 'b'
            sr{tdim} = [floor(0.5 * (nvol + 1)), nvol];
            usedata = 0.5 .* sum(usedata(sr{:}), tdim);

        % normal mean
        case 'm'
            usedata = (1 / nvol) .* sum(usedata, tdim);

        % robust mean
        case 'r'
            usedata = ne_methods.robustmean(double(usedata), tdim);
    end
end

% initial mask ?
imask = (usedata ~= 0);
nvox = sum(imask(:));

% further computation depending on range type
switch (opts.rangetype)

    % histogram-based
    case {'h', 'p'}
        h = hist(usedata(imask), 2000);
        [hc, hp] = cumsum(h);
        mn = ne_methods.findfirst(hc >= (mn * nvox));
        if isempty(mn)
            mn = numel(hp);
        end
        mn = hp(mn);
        mx = ne_methods.findfirst(hc <= (mx * nvox), -1);
        if isempty(mx)
            mx = 1;
        end
        mx = hp(mx);

    % relative to mean
    case 'r'
        mdata = (1 / nvox) .* sum(usedata(imask));
        mn = mdata * mn;
        mx = mdata * mx;

    % relative to mean +/- std's
    case 's'
        mdata = (1 / nvox) .* sum(usedata(imask));
        sdata = std(data(imask));
        mn = mdata + mn * sdata;
        mx = mdata + mx * sdata;
end

% actual mask computation
mask = (imask & usedata >= mn & usedata <= mx);
nvox = sum(mask(:));
mask = ~mask;

% apply mask
rm = ones(1, ndims(data));
rm(tdim) = size(data, tdim);
data(repmat(mask, rm)) = 0;

% now masking depends on type again
switch (type)
    case 'amr'
        for slc = 1:numel(bc.Slice)
            bc.Slice(slc).AMRData = data(:, :, slc);
        end
    case 'fmr'
        if numel(bc.Slice) == 1
            bc.Slice.STCData = data;
        else
            for slc = 1:numel(bc.Slice)
                bc.Slice(slc).STCData = data(:, :, :, slc);
            end
        end
    case 'hdr'
        bc.VoxelData = data;
    case 'head'
        for brc = 1:numel(bc.Brick)
            bc.Brick(brc).Data = data(:, :, :, brc);
        end
    case 'mgh'
        bc.MGHData = data;
    case 'mtc'
        bc.MTCData = data;
    case 'nlf'
        bc.Data = data;
    case 'vmr'
        bc.VMRData = data;
    case 'vtc'
        bc.VTCData = data;
end

% set content
xo.C = bc;
