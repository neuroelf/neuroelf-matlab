function xo = hdr_SliceTiming(xo, opts)
% HDR::SliceTiming  - correct for slice timing differences
%
% FORMAT:       [hdr = ] hdr.SliceTiming([opts]);
%
% Input fields:
%
%       opts        optional settings
%        .interp    interpolation method, one of {'cubic'}, 'lanczos3'
%        .order     slice order, either 1xS spatial scanning order or one
%                   of {'asc'}, 'aint1', 'aint2', 'des', 'dint1', 'dint2'
%        .qweight   use quality measure for weighted interpolation (false)
%        .refslice  reference slice, default: 1
%        .tr        time of repetition in ms (default: 2000)
%        .tshift    if given, either a 1x1 (in ms) or 1xS (in TRs) double 
%
% Output fields:
%
%       hdr         HDR with corrected data
%
% Using: flexinterpn_method, flexinterpnw_method, limitrangec, lsqueeze.

% Version:  v1.1
% Build:    16071416
% Date:     Jul-14 2016, 4:02 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2012, 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'hdr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
nslices = size(bc.VoxelData, 3);
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'qweight') || ~islogical(opts.qweight) || numel(opts.qweight) ~= 1
    opts.qweight = false;
end
if ~isfield(opts, 'interp') || ~ischar(opts.interp) || ~any(strcmpi(opts.interp(:)', ...
    {'linear', 'cubic', 'lanczos2', 'lanczos3', 'lanczos5', 'lanczos8'}))
    opts.interp = 'cubic';
else
    opts.interp = lower(opts.interp(:)');
end
if ~isfield(opts, 'order') || isempty(opts.order) || ~any(strcmpi(class(opts.order), {'char', 'double'}))
    opts.order = 'asc';
end
cord = {'asc', 'aint1', 'aint2', 'des', 'dint1', 'dint2', 'asqr', 'dsqr'};
if ischar(opts.order) && ~any(strcmpi(opts.order(:)', cord))
    opts.order = 0;
elseif ischar(opts.order)
    opts.order = find(strcmpi(opts.order(:)', cord)) - 1;
    if opts.order > 5
        opts.order = opts.order + 14;
    elseif opts.order > 2
        opts.order = opts.order + 7;
    end
end
opts.order = round(opts.order(:)');
opts.order(isinf(opts.order) | isnan(opts.order) | opts.order < 1 | opts.order > nslices) = [];
if ~any(numel(opts.order) == [1, nslices])
    opts.order = 1:nslices;
end
if numel(opts.order) == 1
    switch (opts.order)
        case 1
            opts.order = [1:2:nslices, 2:2:nslices];
        case 2
            opts.order = [2:2:nslices, 1:2:nslices];
        case 10
            opts.order = nslices:-1:1;
        case 11
            opts.order = [nslices:-2:1, (nslices - 1):-2:1];
        case 12
            opts.order = [(nslices - 1):-2:1, nslices:-2:1];
        case 20
            sns = round(sqrt(nslices));
            opts.order = ne_methods.lsqueeze(reshape(1:(sns*(sns+2)), sns, sns + 2)');
            opts.order(opts.order > nslices) = [];
            opts.order = opts.order(:)';
        case 21
            sns = round(sqrt(nslices));
            opts.order = (nslices + 1) - ne_methods.lsqueeze(reshape(1:(sns*(sns+2)), sns, sns + 2)');
            opts.order(opts.order < 1) = [];
            opts.order = opts.order(:)';
        otherwise
            warning('neuroelf:xff:invalidOption', 'Unknown slice-order flag %d.', opts.order);
            opts.order = 1:nslices;
    end
end
if ~isfield(opts, 'refslice') || ~isa(opts.refslice, 'double') || numel(opts.refslice) ~= 1 || ...
   ~any(opts.refslice == (1:nslices))
    opts.refslice = 1;
end
if ~isfield(opts, 'tr') || ~isa(opts.tr, 'double') || numel(opts.tr) ~= 1 || ...
    isinf(opts.tr) || isnan(opts.tr) || opts.tr <= 0
    opts.tr = 2000;
end
if ~isfield(opts, 'tshift') || ~isa(opts.tshift, 'double') || ~any(numel(opts.tshift) == [1, nslices]) || ...
    numel(opts.tshift) ~= max(size(opts.tshift)) || ...
    any(isinf(opts.tshift) | isnan(opts.tshift) | opts.tshift <= -1 | opts.tshift >= 1) || ...
   (numel(opts.tshift) == 1 && (opts.tshift <= 0 || opts.tshift > (opts.tr / nslices)))
    opts.tshift = opts.tr / nslices;
end
if numel(opts.tshift) == 1
    opts.tshift = (opts.tshift / opts.tr) * ((0:(nslices - 1)) - (opts.refslice - 1));
end

% resolve transio if needed
if istransio(bc.VoxelData)
    bc.VoxelData = resolve(bc.VoxelData);
end
stcd = bc.VoxelData;

% patch datatype?
if ~strcmpi(class(stcd), 'single') && ~strcmpi(class(stcd), 'double')
    stcd = single(stcd);
    bc.ImgDim.DataType = 16;
    bc.ImgDim.BitsPerPixel = 32;
end

% apply scaling
if ~any(bc.ImgDim.ScalingSlope == [0, 1])
    if isa(stcd, 'double')
        stcd = double(bc.ImgDim.ScalingSlope) .* stcd;
    else
        stcd = single(bc.ImgDim.ScalingSlope) .* stcd;
    end
end
bc.ImgDim.ScalingSlope = 1;
if bc.ImgDim.ScalingIntercept ~= 0
    if isa(stcd, 'double')
        stcd = double(bc.ImgDim.ScalingIntercept) + stcd;
    else
        stcd = single(bc.ImgDim.ScalingIntercept) + stcd;
    end
end
bc.ImgDim.ScalingIntercept = 0;

% compute voxel-wise quality measure
if opts.qweight

    % weights not already established?
    if ~isfield(bc.RunTimeVars, 'fMRIWeights') || ...
       ~isequal(size(bc.RunTimeVars.fMRIWeights), size(stcd))

        % compute weights
        try
            aft_fMRIWeights(xo);
        catch xfferror
            rethrow(xfferror);
        end
        bc = xo.C;
    end

    % get weights
    stcw = double(bc.RunTimeVars.fMRIWeights);

% otherwise, set to empty array
else
    stcw = [];
end

% extend time course by one at the beginning and the end
szs = size(stcd);
stcd(1, 1, 1, end + 2) = 0;
stcd(:, :, :, 2:end-1) = stcd(:, :, :, 1:end-2);
if ~isempty(stcw)
    stcw(1, 1, 1, end + 2) = 0;
    stcw(:, :, :, 2:end-1) = stcw(:, :, :, 1:end-2);
end
stcx = stcd(:, :, :, 2:6);
flexinterpn_method = ne_methods.flexinterpn_method;
flexinterpnw_method = ne_methods.flexinterpnw_method;
for ic = 1:4
    stcx = flexinterpn_method(stcx, [Inf * ones(1, 4); ...
        1, 1, 1, 0.75; ones(1, 4); size(stcx)], 0, 'lanczos3');
end
stcd(:, :, :, 1) = stcx(:, :, :, 1);
if ~isempty(stcw)
    stcx = stcw(:, :, :, 2:6);
    for ic = 1:4
        stcx = flexinterpn_method(stcx, [Inf * ones(1, 4); ...
            1, 1, 1, 0.75; ones(1, 4); size(stcx)], 0, 'lanczos3');
    end
    stcw(:, :, :, 1) = ne_methods.limitrangec(stcx(:, :, :, 1), 0, 1, 0);
end
stcx = stcd(:, :, :, end-5:end-1);
for ic = 1:4
    stcx = flexinterpn_method(stcx, [Inf * ones(1, 4); ...
         1, 1, 1, 1.25; ones(1, 4); size(stcx) + [0, 0, 0, 1]], 0, 'lanczos3');
end
stcd(:, :, :, end) = stcx(:, :, :, end);
if ~isempty(stcw)
    stcx = stcw(:, :, :, end-5:end-1);
    for ic = 1:4
        stcx = flexinterpn_method(stcx, [Inf * ones(1, 4); ...
             1, 1, 1, 1.25; ones(1, 4); size(stcx) + [0, 0, 0, 1]], 0, 'lanczos3');
    end
    stcw(:, :, :, end) = ne_methods.limitrangec(stcx(:, :, :, end), 0, 1, 0);
end

% patch size argument
szs(3) = 1;
szs(4) = szs(4) + 1;

% iterate over slices
for sc = 1:numel(opts.order)

    % don't do anything for a shift of 0!
    sci = opts.order(sc);

    % interpolate new data
    if isempty(stcw)
        if opts.tshift(sc) == 0
            continue;
        end
        scslc = flexinterpn_method(stcd(:, :, sci, :), ...
            [Inf * ones(1, 4); 1, 1, 1, 2 - opts.tshift(sc); ones(1, 4); szs], ...
            0, opts.interp);

    % for weight-based interpolation
    else

        % patch data first
        stcd(:, :, sci, 2:end-1) = ...
            (stcd(:, :, sci, 1:end-2) .* stcw(:, :, sci, 1:end-2) + ...
            stcd(:, :, sci, 2:end-1) .* stcw(:, :, sci, 2:end-1) + ...
            stcd(:, :, sci, 3:end) .* stcw(:, :, sci, 3:end)) ./ ...
            (stcw(:, :, sci, 1:end-2) + stcw(:, :, sci, 2:end-1) + stcw(:, :, sci, 3:end));
        [scslc, wslc] = flexinterpnw_method(stcd(:, :, sci, :), stcw(:, :, sci, :), ...
            [Inf * ones(1, 4); 1, 1, 1, 2 - opts.tshift(sc); ones(1, 4); szs], ...
            0, opts.interp);
        stcw(:, :, sci, 2:end-1) = wslc;
    end

    % put back into data
    bc.VoxelData(:, :, sci, :) = scslc;
end

% store weights
if ~isempty(stcw)
    bc.RunTimeVars.fMRIWeights = stcw;
end

% set into content
xo.C = bc;
