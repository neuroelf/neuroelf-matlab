function xo = aft_SmoothData3D(xo, k, crd, opts)
% AFT::SmoothData3D  - smooth (partial) spatial data
%
% FORMAT:       [obj = ] obj.SmoothData3D(k, [crd [, opts]]);
%
% Input fields:
%
%       k           1x1 double, smoothing kernel in mm
%       crd         either a 4x3 range definition or Cx3 coordinates
%       opts        optional settings
%        .kcutoff   kernel value cutoff (default: 0.005)
%        .mapvol    required for CMP/VMP objects, defaults to 1
%        .nanrange  data range to exclude from smoothing (default: [])
%        .range     data range to smooth, default [-Inf, Inf]
%        .space     either of 'bvc', 'bvi', 'bvs', 'bvx', {'tal'}
%        .v16       for VMR objects, if true access available V16 data
%
% Output fields:
%
%       data        either XxYxZ or Cx1 data
%
% TYPES: CMP, HDR, HEAD, MGH, NLF, VMP, VMR
%
% Using: bvcoordconv, flexinterpn, rangegrid, smoothdata3, smoothkern.

% Version:  v1.1
% Build:    16031615
% Date:     Mar-16 2016, 3:50 PM EST
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
   ~isa(k, 'double') || numel(k) ~= 1 || isinf(k) || isnan(k) || k < 0 || k > 32
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 3 || ~isa(crd, 'double') || isempty(crd)
    crd = [Inf; -127; 1; 128] * ones(1, 3);
elseif ndims(crd) ~= 2 || ~any(size(crd, 2) == [1, 3])
    error('neuroelf:xff:badArgument', 'Invalid coordinates in call.');
end

% get filetype
bc = xo.C;
ft = lower(xo.S.Extensions{1});

% options
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'kcutoff') || numel(opts.kcutoff) ~= 1 || ~isa(opts.kcutoff, 'double') || ...
    isinf(opts.kcutoff) || isnan(opts.kcutoff) || opts.kcutoff < eps || opts.kcutoff > 0.25
    opts.kcutoff = 0.005;
end
kc = opts.kcutoff;
if ~isfield(opts, 'mapvol') || numel(opts.mapvol) ~= 1 || ~isa(opts.mapvol, 'double') || ...
    isinf(opts.mapvol) || isnan(opts.mapvol) || opts.mapvol < 1
    opts.mapvol = 1;
end
valrange = false;
if ~isfield(opts, 'nanrange') || ~isa(opts.nanrange, 'double') || numel(opts.nanrange) ~= 2 || ...
    any(isnan(opts.nanrange)) || (any(isinf(opts.nanrange)) && opts.nanrange(1) == opts.nanrange(2))
    opts.nanrange = [];
else
    opts.nanrange = [min(opts.nanrange), max(opts.nanrange)];
    valrange = true;
end
if ~isfield(opts, 'range') || ~isa(opts.range, 'double') || numel(opts.range) ~= 2 || ...
    any(isnan(opts.range)) || (any(isinf(opts.range)) && opts.range(1) == opts.range(2))
    opts.range = [-Inf, Inf];
else
    opts.range = [min(opts.range), max(opts.range)];
end
if ~isequal(opts.range, [-Inf, Inf])
    valrange = true;
end
if ~isfield(opts, 'space') || ~ischar(opts.space) || ...
   ~any(strcmpi(opts.space(:)', {'bvc', 'bvi', 'bvs', 'bvx', 'tal'}))
    opts.space = 'tal';
else
    opts.space = lower(opts.space(:)');
end
if ~strcmp(opts.space, 'bvx') && size(crd, 2) ~= 3
    error('neuroelf:xff:badArgument', '1-column coordinates only for ''bvx'' space.');
end
if ~isfield(opts, 'v16') || ~islogical(opts.v16) || numel(opts.v16) ~= 1
    opts.v16 = false;
end

% for a range, unrange (we need this for assigning the data back!)
isrange = false;
if size(crd, 1) == 4 && all(isinf(crd(1, :)))
    isrange = true;
end

% get coordinate frame (needed now and later...)
switch (ft)
    case 'hdr'
        f = hdr_CoordinateFrame(xo);
    case 'head'
        f = head_CoordinateFrame(xo);
    case 'mgh'
        f = mgh_CoordinateFrame(xo);
    case 'nlf'
        f = nlf_CoordinateFrame(xo);
end

% convert coordinates if required
if ~any(strcmp(opts.space, {'bvc', 'bvx'}))

    % for BV types, use bvcoordconv
    if ~any(strcmp(ft, {'hdr', 'head', 'mgh', 'nlf'}))
        trf = ne_methods.bvcoordconv(zeros(0, 3), [opts.space '2bvc'], aft_BoundingBox(xo));

    % for other types
    else

        % only tal allowed
        if ~strcmp(opts.space, 'tal')
            error('neuroelf:xff:badArgument', 'Only TAL (MNI) space for non-BV types.');
        end

        % and convert coordinates
        trf = inv(f.Trf);
    end

    % test TRF for being parallel to standard space axes (e.g. BV types)
    trf3 = trf(1:3, 1:3);
    trf4 = trf(1:3, 4);
    if isrange && all(sum(trf3 ~= 0, 1) == 1) && all(sum(trf3 ~= 0, 2) == 1)

        % axes order -> first axis *IS* first axis
        if (trf3(1, 1) ~= 0)

            % -> second axis *IS* second axis (and third remains third)
            if (trf(2, 2) ~= 0)

                % adapt range
                crd = [crd(1, :); [trf4(1) .* [1; 0; 1] + trf(1, 1) .* crd(2:4, 1), ...
                     trf4(2) .* [1; 0; 1] + trf(2, 2) .* crd(2:4, 2), ...
                     trf4(3) .* [1; 0; 1] + trf(3, 3) .* crd(2:4, 3)]];

            % -> second and third axis swapped
            else

                % adapt range
                crd = [crd(1, :); [trf4(1) .* [1; 0; 1] + trf(1, 1) .* crd(2:4, 1), ...
                     trf4(2) .* [1; 0; 1] + trf(2, 3) .* crd(2:4, 3), ...
                     trf4(3) .* [1; 0; 1] + trf(3, 2) .* crd(2:4, 2)]];
            end

        % -> first axis becomes *SECOND* axis
        elseif (trf(2, 1) ~= 0)

            % -> second axis becomes *THIRD* axis (and third becomes first)
            if (trf(3, 2) ~= 0)

                % adapt range
                crd = [crd(1, :); [trf4(1) .* [1; 0; 1] + trf(1, 3) .* crd(2:4, 3), ...
                     trf4(2) .* [1; 0; 1] + trf(2, 1) .* crd(2:4, 1), ...
                     trf4(3) .* [1; 0; 1] + trf(3, 2) .* crd(2:4, 2)]];

            % -> second axis becomes *FIRST* axis (leaving third the same)
            else

                % adapt range
                crd = [crd(1, :); [trf4(1) .* [1; 0; 1] + trf(1, 2) .* crd(2:4, 2), ...
                     trf4(2) .* [1; 0; 1] + trf(2, 1) .* crd(2:4, 1), ...
                     trf4(3) .* [1; 0; 1] + trf(3, 3) .* crd(2:4, 3)]];
            end

        % -> first axis becomes *THIRD* axis
        else

            % -> third axis becomes *SECOND* axis (and second becomes first)
            if (trf(2, 3) ~= 0)

                % adapt range
                crd = [crd(1, :); [trf4(1) .* [1; 0; 1] + trf(1, 2) .* crd(2:4, 2), ...
                     trf4(2) .* [1; 0; 1] + trf(2, 3) .* crd(2:4, 3), ...
                     trf4(3) .* [1; 0; 1] + trf(3, 1) .* crd(2:4, 1)]];

            % -> third axis becomes *FIRST* axis (leaving second the same)
            else

                % adapt range
                crd = [crd(1, :); [trf4(1) .* [1; 0; 1] + trf(1, 3) .* crd(2:4, 3), ...
                     trf4(2) .* [1; 0; 1] + trf(2, 2) .* crd(2:4, 2), ...
                     trf4(3) .* [1; 0; 1] + trf(3, 1) .* crd(2:4, 1)]];
            end
        end

    % otherwise (either no range or not range-adaptable)
    else

        % we must expand a range!
        if isrange
            crd = ne_methods.rangegrid(crd);
            isrange = false;
        end

        % and then do the multiplication
        crd(:, 4) = 1;
        crd = crd * trf';
        crd(:, 4) = [];
    end

    % then set space to bvc (voxel coordinates)
    opts.space = 'bvc';

end

% get volume and size
odata = aft_GetVolume(xo, opts.mapvol);
vsz = size(odata);

% not range? -> still need to go to bvx
if ~strcmp(opts.space, 'bvx') && ~isrange

    % pre-round off
    crd = round(crd);

    % remove entries out of bounds
    crd(any(crd < 1, 2) | crd(:, 1) > vsz(1) | crd(:, 2) > vsz(2) | crd(:, 3) > vsz(3), :) = [];

    % then convert
    crdx = crd(:, 1) + vsz(1) .* (crd(:, 2) - 1) + (vsz(1) * vsz(2)) .* (crd(:, 3) - 1);

% otherwise
else

    % and further process if not range (BVC)
    if ~isrange

        % just round off
        crdx = round(crd(:));

        % and ensure things are OK
        crdx(crdx < 1 | crdx > prod(vsz)) = [];

        % and create coordinates for sampling (flexinterpn!)
        [crd, cy, cz] = ind2sub(vsz, crdx);
        crd = [crd, cy, cz];

    % for ranges
    else

        % make sure to cover all coordinates
        crd(3, :) = sign(eps + (crd(4, :) - crd(2, :)));
        crd = round(crd);
    end
end

% still range
if isrange

    % ranges must not exceed limits
    crd(2, :) = max(1, min(vsz(1:3), crd(2, :)));
    crd(4, :) = max(1, min(vsz(1:3), crd(4, :)));

    % entire dataset (X)
    if crd(2, 1) == 1 && crd(4, 1) == vsz(1)
        r1 = ':';

    % partial list
    else
        r1 = crd(2, 1):crd(3, 1):crd(4, 1);
    end

    % also for (Y) and (Z)
    if crd(2, 2) == 1 && crd(4, 2) == vsz(2)
        r2 = ':';
    else
        r2 = crd(2, 2):crd(3, 2):crd(4, 2);
    end
    if crd(2, 3) == 1 && crd(4, 3) == vsz(3)
        r3 = ':';
    else
        r3 = crd(2, 3):crd(3, 3):crd(4, 3);
    end
end

% prepare mask for valrange
if valrange

    % preset to true
    wmask = true;

    % if we're using a range
    if isrange

        % fill up to correct size
        wmask = true([numel(r1), numel(r2), numel(r3)]);
        if ischar(r1)
            wmask = wmask(ones(1, vsz(1)), :, :);
        end
        if ischar(r2)
            wmask = wmask(:, ones(1, vsz(2)), :);
        end
        if ischar(r3)
            wmask = wmask(:, :, ones(1, vsz(3)));
        end

        % nan values
        if ~isempty(opts.nanrange)

            % ensure datatype
            if ~isa(odata, 'double') && ~isa(odata, 'single')
                odata = single(odata);
            end

            % set in odata to nan (will be skipped during smoothing)
            odata(odata >= opts.nanrange(1) & odata <= opts.nanrange(2)) = nan;

            % and remove from mask
            if ndims(odata) < 4
                wmask(isnan(odata(r1, r2, r3))) = false;
            else
                wmask(any(isnan(odata(r1, r2, r3, 1, :)), 5)) = false;
            end
        end

        % remove values outside range from mask
        if ~isinf(opts.range(1))
            if ndims(odata) < 4
                wmask(odata(r1, r2, r3) < opts.range(1)) = false;
            else
                wmask(any(odata(r1, r2, r3, 1, :) < opts.range(1), 5)) = false;
            end
        end
        if ~isinf(opts.range(2))
            if ndims(odata) < 4
                wmask(odata(r1, r2, r3) > opts.range(2)) = false;
            else
                wmask(any(odata(r1, r2, r3, 1, :) > opts.range(2), 5)) = false;
            end
        end

        % no actual masking occurred
        if all(wmask(:))
            valrange = false;

        % nothing left to do?
        elseif ~any(wmask(:))
            return;
        end

    % coordinates
    else

        % extend mask
        wmask(1, 1:numel(crdx)) = true;

        % nan values
        if ~isempty(opts.nanrange)

            % set in odata to nan (will be skipped during smoothing)
            odata(odata >= opts.nanrange(1) & odata <= opts.nanrange(2)) = nan;

            % and remove from mask
            if ndims(odata) < 4
                wmask(isnan(odata(crdx))) = false;
            else
                odatas3 = size(odata, 1) * size(odata, 2) * size(odata, 3);
                for d5c = 1:size(odata, 5)
                    wmask(isnan(odata(crdx + (d5c - 1) * odatas3))) = false;
                end
            end
        end

        % remove values outside range from mask
        if ~isinf(opts.range(1))
            if ndims(odata) < 4
                wmask(odata(crdx) < opts.range(1)) = false;
            else
                odatas3 = size(odata, 1) * size(odata, 2) * size(odata, 3);
                for d5c = 1:size(odata, 5)
                    wmask(odata(crdx + (d5c - 1) * odatas3) < opts.range(1)) = false;
                end
            end
        end
        if ~isinf(opts.range(2))
            if ndims(odata) < 4
                wmask(odata(crdx) > opts.range(2)) = false;
            else
                odatas3 = size(odata, 1) * size(odata, 2) * size(odata, 3);
                for d5c = 1:size(odata, 5)
                    wmask(odata(crdx + (d5c - 1) * odatas3) > opts.range(2)) = false;
                end
            end
        end

        % no masking occurred
        if all(wmask)
            valrange = false;

        % remove masked coordinates
        else
            crdx(~wmask) = [];
            if isempty(crdx)
                return;
            end
            crd(~wmask, :) = [];
            wmask = true(1, numel(crdx));
        end
    end
end

% for DMR/FMR/VTC objects, resolve transio, otherwise accept map
switch (ft)

    % CMP/VMP
    case {'cmp', 'vmp'}

        % data field name
        if strcmp(ft, 'cmp')
            if bc.DocumentType ~= 1
                error('neuroelf:xff:invalidSubtype', 'SmoothData3D is only valid for voxel-based CMPs.');
            end
            df = 'CMPData';
        else
            df = 'VMPData';
        end

        % number of volumes
        nvol = numel(bc.Map);
        opts.mapvol = round(min(nvol, opts.mapvol));

        % resolution correction
        k = k / bc.Resolution;

        % resolve transio
        if istransio(bc.Map(opts.mapvol).(df))
            bc.Map(opts.mapvol).(df) = resolve(bc.Map(opts.mapvol).(df));
        end

        % get to work
        sk = ne_methods.smoothkern(k, kc);

        % range
        if isrange

            % smooth entire dataset
            smdata = ne_methods.smoothdata3(odata, [k, k, k], kc);

            % replace data
            if valrange

                % only part
                smdata = smdata(r1, r2, r3);

                % original data
                odata = bc.Map(opts.mapvol).(df)(r1, r2, r3);

                % replace
                smdata(~wmask) = odata(~wmask);

                % put into output
                bc.Map(opts.mapvol).(df)(r1, r2, r3) = smdata;

            % full data
            else

                % put into output
                bc.Map(opts.mapvol).(df)(r1, r2, r3) = smdata(r1, r2, r3);
            end

        % no range, use flexinterpn directly (slower!)
        else
            bc.Map(opts.mapvol).(df)(crdx) = ne_methods.flexinterpn(odata, crd, sk, 1, 0);
        end

    % HDR
    case 'hdr'

        % resolution correction
        if ~isrange
            if ~all(f.Resolution == f.Resolution(1))
                f.Resolution = prod(f.Resolution) ^ (1 / 3);
            else
                f.Resolution = f.Resolution(1);
            end
        end
        k = k ./ f.Resolution;

        % special datatypes
        if any(bc.ImgDim.DataType == [32, 128, 1536, 1792, 2048, 2304])
            switch (bc.ImgDim.DataType)

                % comments, see code below
                case {128, 2304}
                    nvol = numel(size(bc.VoxelDataRGBA, 4));
                    opts.mapvol = round(min(nvol, opts.mapvol));

                    % get data without scaling applyied
                    odata = double(bc.VoxelDataRGBA(:, :, :, opts.mapvol, :));
                    odatasz = size(odata);

                    % working on a range
                    if isrange

                        % smooth entire dataset
                        smdata = reshape(ne_methods.smoothdata3(squeeze(odata), k, kc), odatasz);

                        % replace data
                        if valrange

                            % only part
                            smdata = smdata(r1, r2, r3, 1, :);
                            odata = odata(r1, r2, r3, 1, :);
                            odatas3 = size(odata, 1) * size(odata, 2) * size(odata, 3);
                            wmask = find(~wmask);

                            % replace
                            for d5c = 1:size(odata, 5)
                                smdata(wmask + (d5c - 1) * odatas3) = ...
                                    odata(wmask + (d5c - 1) * odatas3);
                            end

                            % put into output
                            bc.VoxelDataRGBA(r1, r2, r3, opts.mapvol, :) = smdata;

                        % full data
                        else

                            % put into output
                            bc.VoxelDataRGBA(r1, r2, r3, opts.mapvol, :) = ...
                                smdata(r1, r2, r3, 1, :);
                        end

                    % no range, use flexinterpn directly (slower!)
                    else
                        k = ne_methods.smoothkern(k, kc);
                        vsz3 = prod(vsz(1:3));
                        vsz4 = vsz3 * vsz(4);
                        for d5c = 1:size(odata, 5)
                            bc.VoxelDataRGBA(crdx + (opts.mapvol - 1) * vsz3 + (d5c - 1) * vsz4) = ...
                                ne_methods.flexinterpn(odata(:, :, :, d5c), crd, k, 1, 0);
                        end
                    end

                otherwise
                    error('neuroelf:xff:unsupported', 'Smoothing of complex datatypes not supported.');
            end

            % set back and already return
            xo.C = bc;
            return;
        end

        % number of volumes
        nvol = numel(size(bc.VoxelData, 4));
        opts.mapvol = round(min(nvol, opts.mapvol));

        % resolve transio
        if istransio(bc.VoxelData)
            bc.VoxelData = resolve(bc.VoxelData);
        end

        % get data without scaling applyied
        oodata = double(bc.VoxelData(:, :, :, opts.mapvol));

        % but keep nans
        if ~isempty(opts.nanrange)
            oodata(isnan(odata)) = nan;
        end
        odata = oodata;

        % working on a range
        if isrange

            % smooth entire dataset
            smdata = ne_methods.smoothdata3(odata, k, kc);

            % replace data
            if valrange

                % only part
                smdata = smdata(r1, r2, r3);
                odata = odata(r1, r2, r3);

                % replace
                smdata(~wmask) = odata(~wmask);

                % put into output
                bc.VoxelData(r1, r2, r3, opts.mapvol) = smdata;

            % full data
            else

                % put into output
                bc.VoxelData(r1, r2, r3, opts.mapvol) = smdata(r1, r2, r3);
            end

        % no range, use flexinterpn directly (slower!)
        else
            k = ne_methods.smoothkern(k, kc);
            bc.VoxelData(crdx + (opts.mapvol - 1) * prod(vsz)) = ne_methods.flexinterpn( ...
                odata, crd, k, 1, 0);
        end

    % HEAD
    case 'head'

        % resolution correction
        if ~isrange
            if ~all(f.Resolution == f.Resolution(1))
                f.Resolution = prod(f.Resolution) ^ (1 / 3);
            else
                f.Resolution = f.Resolution(1);
            end
        end
        k = k ./ f.Resolution;

        % number of volumes
        nvol = numel(bc.Brick);
        opts.mapvol = round(min(nvol, opts.mapvol));

        % resolve transio
        if istransio(bc.Brick(opts.mapvol).Data)
            bc.Brick(opts.mapvol).Data = resolve(bc.Brick(opts.mapvol).Data);
        end

        % get data without scaling applied
        oodata = double(bc.Brick(opts.mapvol).Data);
        if ~isempty(opts.nanrange)
            oodata(isnan(odata)) = nan;
        end
        odata = oodata;

        % working on a range
        if isrange

            % smooth entire dataset
            smdata = ne_methods.smoothdata3(odata, k, kc);

            % replace data
            if valrange

                % only part
                smdata = smdata(r1, r2, r3);
                odata = odata(r1, r2, r3);

                % replace
                smdata(~wmask) = odata(~wmask);

                % put into output
                bc.Brick(opts.mapvol).Data(r1, r2, r3) = smdata;

            % full data
            else

                % put into output
                bc.Brick(opts.mapvol).Data(r1, r2, r3) = smdata(r1, r2, r3);
            end

        % no range, use flexinterpn directly (slower!)
        else
            k = ne_methods.smoothkern(k, kc);
            bc.Brick(opts.mapvol).Data(crdx) = ne_methods.flexinterpn( ...
                bc.Brick(opts.mapvol).Data, crd, k, 1, 0);
        end

    % MGH
    case 'mgh'

        % resolution correction
        if ~isrange
            if ~all(f.Resolution == f.Resolution(1))
                f.Resolution = prod(f.Resolution) ^ (1 / 3);
            else
                f.Resolution = f.Resolution(1);
            end
        end
        k = k ./ f.Resolution;

        % number of volumes
        opts.mapvol = 1;

        % resolve transio
        if istransio(bc.MGHData)
            bc.MGH.Data = resolve(bc.MGHData);
        end

        % get data without scaling applied
        oodata = double(bc.MGHData);
        if ~isempty(opts.nanrange)
            oodata(isnan(odata)) = nan;
        end
        odata = oodata;

        % working on a range
        if isrange

            % smooth entire dataset
            smdata = ne_methods.smoothdata3(odata, k, kc);

            % replace data
            if valrange

                % only part
                smdata = smdata(r1, r2, r3);
                odata = odata(r1, r2, r3);

                % replace
                smdata(~wmask) = odata(~wmask);

                % put into output
                bc.MGHData(r1, r2, r3) = smdata;

            % full data
            else

                % put into output
                bc.MGHData(r1, r2, r3) = smdata(r1, r2, r3);
            end

        % no range, use flexinterpn directly (slower!)
        else
            k = ne_methods.smoothkern(k, kc);
            bc.MGHData(crdx) = ne_methods.flexinterpn( ...
                bc.MGHData, crd, k, 1, 0);
        end

    % NLF
    case 'nlf'

        % resolution correction
        if ~isrange
            if ~all(f.Resolution == f.Resolution(1))
                f.Resolution = prod(f.Resolution) ^ (1 / 3);
            else
                f.Resolution = f.Resolution(1);
            end
        end
        k = k ./ f.Resolution;

        % number of volumes
        nvol = size(bc.Data, 4);
        opts.mapvol = round(min(nvol, opts.mapvol));

        % resolve transio
        if istransio(bc.Data)
            bc.Data = resolve(bc.Data);
        end
        oodata = bc.Data(:, :, :, opts.mapvol);
        if ~isempty(opts.nanrange)
            oodata(isnan(odata)) = nan;
        end
        odata = oodata;

        % working on a range
        if isrange

            % smooth entire dataset
            smdata = ne_methods.smoothdata3(odata, k, kc);

            % replace data
            if valrange

                % only part
                smdata = smdata(r1, r2, r3);
                odata = odata(r1, r2, r3);

                % replace
                smdata(~wmask) = odata(~wmask);

                % put into output
                bc.Data(r1, r2, r3, opts.mapvol) = smdata;

            % full data
            else

                % put into output
                bc.Data(r1, r2, r3, opts.mapvol) = smdata(r1, r2, r3);
            end

        % no range, use flexinterpn directly (slower!)
        else
            k = ne_methods.smoothkern(k, kc);
            bc.Data(crdx + (opts.mapvol - 1) * prod(vsz)) = ne_methods.flexinterpn( ...
                bc.Data(:, :, :, opts.mapvol), crd, k, 1, 0);
        end

    % VMR
    case 'vmr'

        % resolution correction
        if isrange
            res = [bc.VoxResX, bc.VoxResY, bc.VoxResZ];
        else
            res = (bc.VoxResX * bc.VoxResY * bc.VoxResZ) ^ (1 / 3);
        end
        k = k ./ res;

        % which field
        if isequal(size(bc.VMRData16), size(bc.VMRData)) && opts.v16
            df = 'VMRData16';
        else
            df = 'VMRData';
        end

        % resolve transio
        if istransio(bc.(df))
            bc.(df) = resolve(bc.(df));
        end

        % working on a range
        if isrange

            % patch ranges
            [sr1, sr2, sr3, tr1, tr2, tr3] = ...
                patchrange(r1, r2, r3, size(odata), k, kc, ne_methods.smoothkern);

            % smooth entire dataset
            smdata = ne_methods.smoothdata3(odata(sr1, sr2, sr3), k, kc);

            % replace data
            if valrange

                % only part
                smdata = smdata(tr1, tr2, tr3);
                odata = odata(r1, r2, r3);

                % replace
                smdata(~wmask) = odata(~wmask);

                % put into output
                bc.(df)(r1, r2, r3) = smdata;

            % full data
            else

                % put into output
                bc.(df)(r1, r2, r3) = smdata(tr1, tr2, tr3);
            end

        % no range, use flexinterpn directly (slower!)
        else
            k = ne_methods.smoothkern(k, kc);
            bc.(df)(crdx) = ne_methods.flexinterpn(bc.(df), crd, k, 1, 0);
        end

    otherwise
        error('neuroelf:xff:objectTypeUnsupported', 'Unsupported object type for ::SmoothData3D method.');
end

% set back
xo.C = bc;


% sub-function to patch ranges
function [sr1, sr2, sr3, tr1, tr2, tr3] = patchrange(r1, r2, r3, sz, k, kt, smoothkern)

% ranges don't need patching
if numel(r1) > (sz(1) - 6) && numel(r2) > (sz(2) - 6) && numel(r3) > (sz(3) - 6)

    % ranges remain
    sr1 = r1;
    sr2 = r2;
    sr3 = r3;
    tr1 = 1:numel(r1);
    tr2 = 1:numel(r2);
    tr3 = 1:numel(r3);
    return;
end

% get the size of the kernels
if numel(k) == 1 || all(k == k(1))
    ks1 = 0.5 * (numel(smoothkern(k(1), kt)) - 1);
    ks2 = ks1;
    ks3 = ks1;
else
    ks1 = 0.5 * (numel(smoothkern(k(1), kt)) - 1);
    ks2 = 0.5 * (numel(smoothkern(k(2), kt)) - 1);
    ks3 = 0.5 * (numel(smoothkern(k(3), kt)) - 1);
end

% get range end-points
rp1 = sort(r1([1, end]));
rp2 = sort(r2([1, end]));
rp3 = sort(r3([1, end]));

% the range to smooth is the original range + the kernel size
sr1 = max(1, rp1(1) - ks1):min(sz(1), rp1(2) + ks1);
sr2 = max(1, rp2(1) - ks2):min(sz(2), rp2(2) + ks2);
sr3 = max(1, rp3(1) - ks3):min(sz(3), rp3(2) + ks3);

% the target range sizes then is the difference in range sizes
if r1(1) <= r1(end)
    tr1 = (1 + rp1(1) - sr1(1)):(numel(sr1) - (sr1(end) - rp1(2)));
else
    tr1 = (numel(sr1) - (sr1(end) - rp1(2))):-1:(1 + rp1(1) - sr1(1));
end
if r2(1) <= r2(end)
    tr2 = (1 + rp1(2) - sr2(1)):(numel(sr2) - (sr2(end) - rp2(2)));
else
    tr2 = (numel(sr2) - (sr2(end) - rp2(2))):-1:(1 + rp2(1) - sr2(1));
end
if r3(1) <= r3(end)
    tr3 = (1 + rp3(1) - sr3(1)):(numel(sr3) - (sr3(end) - rp3(2)));
else
    tr3 = (numel(sr3) - (sr3(end) - rp3(2))):-1:(1 + rp3(1) - sr3(1));
end
