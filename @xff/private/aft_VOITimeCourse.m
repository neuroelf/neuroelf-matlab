function [voitc, voin, wr] = aft_VOITimeCourse(xo, voi, opts)
% AFT::VOITimeCourse  - extract VOI time course data
%
% FORMAT:       [voitc, voin] = obj.VOITimeCourse(voi [, opts])
%
% Input fields:
%
%       voi         VOI object or coordinates (e.g. from VOI::BVCoords)
%       opts        optional settings
%        .bvcomp    BV-compatibility flag (default: false)
%        .fliplr    flip left/right for radiological convention
%        .imeth     one of 'cubic', 'lanczos3', 'linear', {'nearest'}
%        .subsel    subject ID of current object (default: from filename)
%        .subvois   if voi is a VOI object, one of 'sub_', '_sub', {'voi'}
%        .voisel    voi selection (either name or number)
%        .weight    voxel weighting (mostly for nearest neighbor)
%                   - scalar 0: unique coordinates
%                   - scalar 1: resample equal coordinates (default)
%                   - scalar 2: to get a cell array of TxC arrays
%                   - scalar 3: cell array of TxC unique coords arrays
%                   - cell array with weights for each selected VOI
%                     * this also then supports giving voxel-wise weights*
%
% Output fields:
%
%       voitc       TxV time course of voi(s) (or 1xV cell with TxC double)
%       voin        1xV voinames (without subject ID)
%
% TYPES: HDR, HEAD, VTC
%
% Note: this implementation supercedes the original VTC::VOITimeCourse,
%       which for compatibility reasons is still available via
%       VTC::VOITimeCourseOrig
%
%       the .bvcomp flag enforces to floor the coordinates (instead of
%       rounding them), in line with BV's internal space (only for VTCs)
%
% Using: applyspmsnc, findfirst, flexinterpn_method, limitrangec.

% Version:  v1.1
% Build:    16032210
% Date:     Mar-22 2016, 10:25 AM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'hdr', 'head', 'vtc'}) || ...
   (~all(xffisobject(voi(:), true, 'voi')) && ~isa(voi, 'double')) || isempty(voi)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin > 2 && isa(opts, 'double') && numel(opts) == 1 && ~isnan(opts)
    opts = struct( 'imeth',  'nearest', 'voisel', [], 'weight', opts);
    if isinf(opts.weight)
        opts.weight = 3;
    end
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
ftyp = lower(xo.S.Extensions{1});
bc = xo.C;
if isa(voi, 'double')
    voicrd = voi;
    voi = xff('new:voi');
    voic = voi.C;
    voic.VOI(1).Voxels = voicrd;
    voic.VOI.NrOfVoxels = size(voicrd, 1);
    voic.VOI.Name = 'temp';
    voi.C = voic;
    fromdouble = true;
    opts.subsel = '';
    opts.subvois = 'voi';
    opts.voisel = 'temp';
else
    fromdouble = false;
    voic = voi.C;
end
vn = {voic.VOI.Name};
[sfp, sf] = fileparts(xo.F);
if ~isfield(opts, 'bvcomp') || ~islogical(opts.bvcomp) || numel(opts.bvcomp) ~= 1
    opts.bvcomp = false;
end
if ~isfield(opts, 'fliplr') || ~islogical(opts.fliplr) || numel(opts.fliplr) ~= 1
    opts.fliplr = false;
end
if ~isfield(opts, 'imeth') || ~ischar(opts.imeth) || ...
   ~any(strcmpi(opts.imeth(:)', {'cubic', 'lanczos3', 'linear', 'nearest'}))
    opts.imeth = 'nearest';
else
    opts.imeth = lower(opts.imeth(:)');
end
if opts.imeth(1) == 'n'
    nearest = true;
else
    nearest = false;
end
if ~isfield(opts, 'subsel') || ~ischar(opts.subsel) || isempty(opts.subsel)
    opts.subsel = '';
end
if ~isfield(opts, 'subvois') || ~ischar(opts.subvois) || ...
   ~any(strcmpi(opts.subvois(:)', {'sub_', '_sub', 'voi'}))
    opts.subvois = 'v';
else
    opts.subvois = lower(opts.subvois(1));
    if opts.subvois ~= 'v'
        if isempty(opts.subsel)
            if isempty(sf) || ~any(sf == '_') || sf(1) == '_'
                error('neuroelf:xff:badArgument', 'Cannot resolve subject ID for VOI access.');
            end
            opts.subsel = regexprep(sf, '^([^_]+)_.*$', '$1');
        end
    end
end
if ~isfield(opts, 'voisel') || isempty(opts.voisel) || ...
   (~isa(opts.voisel, 'double') && ~ischar(opts.voisel) && ~iscell(opts.voisel))
    if opts.subvois == 'v'
        opts.voisel = 1:numel(voic.VOI);
    end
elseif isa(opts.voisel, 'double')
    opts.voisel = opts.voisel(:)';
    if any(isinf(opts.voisel) | isnan(opts.voisel) | ...
        opts.voisel < 1 | opts.voisel > numel(voic.VOI) | opts.voisel ~= fix(opts.voisel)) || ...
        numel(opts.voisel) ~= numel(unique(opts.voisel))
        error('neuroelf:xff:badArgument', 'Invalid numeric VOI selection.');
    end
else
    if ischar(opts.voisel)
        opts.voisel = {opts.voisel(:)'};
    else
        opts.voisel = opts.voisel(:);
    end
    try
        for vc = numel(opts.voisel):-1:1
            switch (opts.subvois)
                case {'v'}
                    opts.voisel{vc} = ne_methods.findfirst(strcmpi(opts.voisel{vc}(:)', vn));
                case {'_'}
                    opts.voisel{vc} = ...
                        ne_methods.findfirst(strcmpi([opts.voisel{vc}(:)' '_' opts.subsel], vn));
                case {'s'}
                    opts.voisel{vc} = ...
                        ne_methods.findfirst(strcmpi([opts.subsel '_' opts.voisel{vc}(:)'], vn));
            end
        end
        opts.voisel = cat(1, opts.voisel{:});
    catch xfferror
        neuroelf_lasterr(xfferror);
        error('neuroelf:xff:badArgument', 'Invalid VOI name selection: %s.', xfferror.message);
    end
end
numvois = numel(opts.voisel);
if ~isfield(opts, 'weight') || ...
   (~isa(opts.weight, 'double') && ~iscell(opts.weight)) || ...
   (isa(opts.weight, 'double') && (numel(opts.weight) ~= 1 || ...
     ~isa(opts.weight, 'double') || isinf(opts.weight) || isnan(opts.weight) || ...
     ~any(opts.weight == (0:3)))) || ...
   (iscell(opts.weight) && ~any(numel(opts.weight) == [1, numvois]))
    opts.weight = 1;
end
if isa(opts.weight, 'double')
    opts.weight = {opts.weight};
end
if numel(opts.weight) ~= numvois
    opts.weight = opts.weight(1, ones(1, numvois));
end
voig = cell(1, numvois);
vois = cell(1, numvois);
cellout = false;
donames = true;
for vc = 1:numvois
    vois{vc} = voi_TalCoords(voi, opts.voisel(vc));
    if opts.fliplr
        vois{vc}(:, 1) = -vois{vc}(:, 1);
    end
    if ~any(numel(opts.weight{vc}) == [1, size(vois{vc}, 1)])
        opts.weight{vc} = 1;
    elseif numel(opts.weight{vc}) == 1 && any(opts.weight{vc} == [2, 3])
        cellout = true;
        if opts.weight{vc} == 3
            donames = false;
        end
    elseif numel(opts.weight{vc}) == size(vois{vc}, 1)
        opts.weight{vc} = ne_methods.limitrangec(opts.weight{vc}(:), -1e6, 1e6, 0);
    end
end
if nargout > 1
    if donames
        voin = cell(1, numvois);
        for vc = 1:numvois
            switch (opts.subvois)
                case 'v'
                    voin{vc} = vn{opts.voisel(vc)};
                case '_'
                    voin{vc} = regexprep(vn{opts.voisel(vc)}, ['_' opts.subsel '^'], '', 'ignorecase');
                case 's'
                    voin{vc} = regexprep(vn{opts.voisel(vc)}, ['$' opts.subsel '_'], '', 'ignorecase');
            end
        end
        if nargout > 2
            wr = cell(1, numvois);
        end
    else
        voin = cell(1, numvois);
        wr = cell(1, numvois);
    end
end

% clear temporary object
if fromdouble
    delete(voi);
end

% depending on filetype
switch (ftyp)
    case 'hdr'
        cfr = hdr_CoordinateFrame(xo);
        trf = inv(cfr.Trf)';
        data = bc.VoxelData;
        vsz = size(data);
        vsz(4:end) = [];
        numtp = size(data, 4);
        sci = double(bc.ImgDim.ScalingIntercept);
        scs = double(bc.ImgDim.ScalingSlope);
        if scs == 0
            scs = 1;
        end
    case 'head'
        cfr = head_CoordinateFrame(xo);
        trf = inv(cfr.Trf)';
        data = {bc.Brick.Data};
        vsz = size(data{1});
        numtp = numel(data);
    case 'vtc'
        bbx = aft_BoundingBox(xo);
        rtv = bc.RunTimeVars;
        trf = bbx.QuatT2B';
        data = bc.VTCData;
        vsz = size(data);
        numtp = vsz(1);
        vsz(1) = [];

        % allow new subject-specific normalization
        if isfield(rtv, 'SPMsn') && isstruct(rtv.SPMsn) && numel(rtv.SPMsn) == 1 && ...
            isfield(rtv.SPMsn, 'VG') && isstruct(rtv.SPMsn.VG) && ...
            isfield(rtv.SPMsn.VG, 'dim') && isfield(rtv.SPMsn.VG, 'mat') && ...
            isfield(rtv.SPMsn, 'VF') && isstruct(rtv.SPMsn.VF) && ...
            isfield(rtv.SPMsn.VF, 'dim') && isfield(rtv.SPMsn.VF, 'mat') && ...
            isfield(rtv.SPMsn, 'Tr') && isa(rtv.SPMsn.Tr, 'double') && ndims(rtv.SPMsn.Tr) == 4 && ...
            isfield(rtv.SPMsn, 'Affine') && isa(rtv.SPMsn.Affine, 'double') && ...
            isequal(size(rtv.SPMsn.Affine), [4, 4])

            % apply SPMsn
            sn = rtv.SPMsn;
            itrf = inv(double(bbx.QuatB2T)) * inv(rtv.TrfPlus);
            ivgm = inv(double(sn.VG(1).mat));
            for vc = 1:numvois
                vois{vc} = ne_methods.applyspmsnc(vois{vc}, sn.Tr, ...
                    sn.VG(1).dim, ivgm, itrf * sn.VF(1).mat * sn.Affine);
            end
            trf = [];
        end
end

% for nearest-neighbor interpolation
if nearest

    vszp = cumprod([1, vsz]);
    vszp(end) = [];
    for vc = 1:numvois
        ncrd = size(vois{vc}, 1);
        if ~isempty(trf)
            if strcmp(ftyp, 'vtc') && opts.bvcomp
                voicrd = floor([vois{vc}, ones(ncrd, 1)] * trf);
            else
                voicrd = round([vois{vc}, ones(ncrd, 1)] * trf);
            end
        else
            voicrd = round(vois{vc});
        end
        voig{vc} = ~(any(voicrd < 1, 2) | voicrd(:, 1) > vsz(1) | ...
            voicrd(:, 2) > vsz(2) | voicrd(:, 3) > vsz(3));
        voicrd(~voig{vc}, :) = 1;
        vois{vc} = 1 + sum((voicrd(:, 1:3) - 1) .* vszp(ones(ncrd, 1), :), 2);
        if isequal(opts.weight{vc}, 0)
            vois{vc}(~voig{vc}) = [];
            voig{vc} = ':';
            vois{vc} = unique(vois{vc});
        elseif isequal(opts.weight{vc}, 3)
            [uvec, voin{vc}, wr{vc}] = unique(vois{vc});
        elseif nargout > 2
            wr(vc) = vois(vc);
        end
    end

% otherwise
else

    % resolve contents for VTCs
    if strcmp(ftyp, 'vtc') && istransio(data)
        data = resolve(data);
    end

    % transform and cut away forth element
    vsz = vsz + 0.5;
    for vc = 1:numvois
        ncrd = size(vois{vc}, 1);
        if ~isempty(trf)
            voicrd = [vois{vc}, ones(ncrd, 1)] * trf;
        else
            voicrd = vois{vc};
        end
        voig{vc} = ~(any(voicrd <= 0.5, 2) | voicrd(:, 1) >= vsz(1) | ...
            voicrd(:, 2) >= vsz(2) | voicrd(:, 3) >= vsz(3));
        vois{vc} = voicrd(:, 1:3);
        if isequal(opts.weight{vc}, 0)
            vois{vc}(~voig{vc}, :) = [];
            voig{vc} = ':';
            vois{vc} = unique(vois{vc}, 'rows');
        elseif isequal(opts.weight{vc}, 3)
            [uvec, voin{vc}, wr{vc}] = unique(vois{vc}, 'rows');
        elseif nargout > 2
            wr(vc) = vois(vc);
        end
    end
end

% prepare output
voitc = cell(1, numvois);
for vc = 1:numvois
    voitc{vc} = zeros(numtp, size(vois{vc}, 1));
end

% for non VTC objects or interpolation
if ~nearest || ~strcmp(ftyp, 'vtc')

    % interpolation kernel
    if ~nearest
        [nulld, ik] = ne_methods.flexinterpn_method(zeros(3, 1), [Inf; 1; 1; 3], opts.imeth);
    end

    % first iterate over volumes
    for tc = 1:numtp

        % get data
        switch (ftyp)
            case 'hdr'
                voldata = data(:, :, :, tc);
                if scs ~= 1 || sci ~= 0
                    voldata = sci + scs .* double(voldata);
                end
            case 'head'
                voldata = data{tc}(:, :, :);
                if ~any(bc.Brick(tc).ScalingFactor == [0, 1])
                    voldata = double(bc.Brick(tc).ScalingFactor) .* double(voldata);
                end
            case 'vtc'
                voldata = squeeze(data(tc, :, :, :));
        end

        % nearest neighbor (hdr/head)
        if nearest

            % access
            for vc = 1:numvois
                voitc{vc}(tc, voig{vc}) = voldata(vois{vc}(voig{vc}));
                if numel(opts.weight{vc}) == size(voitc{vc}, 2)
                    voitc{vc}(tc, :) = opts.weight{vc}(:)' .* voitc{vc}(tc, :);
                end
            end

        % interpolate
        else
            for vc = 1:numvois
                voitc{vc}(tc, voig{vc}) = ...
                    ne_methods.flexinterpn_method(voldata, vois{vc}(voig{vc}, :), ik{:});
                if numel(opts.weight{vc}) == size(voitc{vc}, 2)
                    voitc{vc}(tc, :) = opts.weight{vc}(:)' .* voitc{vc}(tc, :);
                end
            end
        end
    end

% nearest neighbor on VTCs
else

    % iterate over vois
    for vc = 1:numvois

        % don't re-sample for transio!
        if any(voig{vc})
            if istransio(data)
                [voxsu, voxsi, voxsj] = unique(vois{vc}(voig{vc}));
                vvtc = data(:, voxsu);
                voitc{vc}(:, voig{vc}) = double(vvtc(:, voxsj));
            else
                voitc{vc}(:, voig{vc}) = double(data(:, vois{vc}(voig{vc})));
            end
        end
        if numel(opts.weight{vc}) == size(voitc{vc}, 2)
            voitc{vc} = (ones(numtp, 1) * opts.weight{vc}(:)') .* voitc{vc};
        end
    end
end

% average within VOIs
if ~cellout

    for vc = 1:numvois
        voitc{vc} = mean( ...
            voitc{vc}(:, all(~isinf(voitc{vc}) & ~isnan(voitc{vc}) & voitc{vc} ~= 0, 1)), 2);
    end
    voitc = cat(2, voitc{:});
end
