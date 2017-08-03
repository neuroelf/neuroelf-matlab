% FUNCTION ne_smoothstats: smooth stats border (high-q output)
function stvix = ne_smoothstats(varargin)

% Version:  v0.9d
% Build:    14061710
% Date:     Jun-17 2014, 10:36 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, Jochen Weber
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

% global variable
global ne_gcfg;

% this method requires a valid stats var
if nargin < 4
    stvar = ne_gcfg.fcfg.StatsVar;
else
    stvar = varargin{4};
end
if numel(stvar) ~= 1 || ...
   ~isxff(stvar, {'cmp', 'hdr', 'head', 'vmp'})
    return;
end
stvtyp = lower(stvar.Filetype);
if any(strcmp(stvtyp, {'cmp', 'vmp'}))
    res = stvar.Resolution;
else
    res = prod(stvar.CoordinateFrame.Resolution) .^ (1 / 3);
end
if nargin < 5 || ...
   ~isa(varargin{5}, 'double') || ...
    isempty(varargin{5}) || ...
    any(isinf(varargin{5}(:)) | isnan(varargin{5}(:)) | varargin{5}(:) < 1)
    stvix = ne_gcfg.h.StatsVarMaps.Value;
else
    stvix = unique(round(varargin{5}(:)));
end
if isempty(stvix) || ...
    any(stvix(:) > numel(stvar.Map))
    return;
end

% smoothing kernel
if nargin < 3 || ...
   ~isa(varargin{3}, 'double') || ...
    numel(varargin{3}) ~= 1 || ...
    isinf(varargin{3}) || ...
    isnan(varargin{3}) || ...
    varargin{3} <= 0
    smkv = inputdlg({'Smoothking kernel:'}, 'NeuroElf - user input', ...
        1, {'6'});
    if ~iscell(smkv) || ...
        numel(smkv) ~= 1 || ...
        isempty(smkv{1})
        return;
    end
    try
        smkv = str2double(smkv{1});
        if numel(smkv) ~= 1 || ...
            isinf(smkv) || ...
            isnan(smkv) || ...
            smkv <= 0
            return;
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        return;
    end
else
    smkv = varargin{3};
end
smkv = min(64, smkv);
if smkv <= (0.5 * res)
    return;
end
smk = smoothkern(smkv / res, 0.001);

% get selected indices
stvixo = stvix;
if nargin < 6 || ...
   ~isa(varargin{6}, 'double') || ...
    any(isinf(varargin{6}(:)) | isnan(varargin{6}(:)) | varargin{6}(:) < 1 | varargin{6}(:) ~= round(varargin{6}(:))) || ...
    numel(varargin{6}) ~= numel(stvix) || ...
    numel(unique(varargin{6}(:))) ~= numel(varargin{6})
else
    stvixo = varargin{6}(:);
end

% border only
bonly = false;
if nargin > 6 && ...
    islogical(varargin{7}) && ...
    numel(varargin{7}) == 1
    bonly = varargin{7};
end

% for each map
for mc = 1:numel(stvixo)

    % get map, data, and clustering information
    map = stvar.Map(stvixo(mc));
    switch (stvtyp)
        case {'cmp'}
            mapdata = map.CMPData(:, :, :);
            mapdatact = map.CMPDataCT;
        case {'hdr'}
            if istransio(stvar.VoxelData)
                stvar.VoxelData = resolve(stvar.VoxelData);
            end
            mapdata = stvar.VoxelData(:, :, :, stvixo(mc));
            if numel(stvar.VoxelDataCT) >= stvixo(mc)
                mapdatact = stvar.VoxelDataCT{stvixo(mc)};
            else
                mapdatact = logical([]);
            end
            mapdatane = numel(mapdata);
        case {'head'}
            mapdata = stvar.Brick(stvixo(mc)).Data(:, :, :, 1);
            mapdatact = stvar.Brick(stvixo(mc)).DataCT;
        case {'vmp'}
            mapdata = map.VMPData(:, :, :);
            mapdatact = map.VMPDataCT;
    end

    % get border
    if ~isempty(mapdatact);
        brde = mapdatact;
        brd = erode3d(brde);
    else
        brd = (mapdata ~= 0);
        brde = dilate3d(brd);
    end
    brdd = erode3d(brd);

    % find voxels to smooth
    smvi = find(brd(:) & ~brdd(:));
    [smvx, smvy, smvz] = ind2sub(size(brd), smvi);

    % copy map ?
    if stvix(mc) == stvixo(mc)
        stvix(mc) = numel(stvar.Map) + 1;
        stvar.Map(stvix(mc)) = map;
    end

    % remove cluster check?
    if ~isempty(mapdatact)
        stvar.Map(stvix(mc)).EnableClusterCheck = 0;
    end

    % border only ?
    if bonly

        % replace in target
        switch (stvtyp)
            case {'cmp'}
                stvar.Map(stvix(mc)).CMPData = mapdata .* brde;
                stvar.Map(stvix(mc)).CMPData(smvi) = flexinterpn(mapdata, ...
                    [smvx, smvy, smvz], smk, 1, 0);
                stvar.Map(stvix(mc)).CMPDataCT = [];
            case {'hdr'}
                stvar.VoxelData(:, :, :, stvix(mc)) = mapdata .* brde;
                stvar.VoxelData(smvi + mapdatane * (stvix(mc) - 1)) = ...
                    flexinterpn(mapdata, [smvx, smvy, smvz], smk, 1, 0);
                stvar.VoxelDataCT{stvix(mc)} = [];
            case {'head'}
                stvar.Brick(stvix(mc)).Data = mapdata .* brde;
                stvar.Brick(stvix(mc)).Data(smvi) = flexinterpn(mapdata, ...
                    [smvx, smvy, smvz], smk, 1, 0);
                stvar.Brick(stvix(mc)).DataCT = [];
            case {'vmp'}
                stvar.Map(stvix(mc)).VMPData = mapdata .* brde;
                stvar.Map(stvix(mc)).VMPData(smvi) = flexinterpn(mapdata, ...
                    [smvx, smvy, smvz], smk, 1, 0);
                stvar.Map(stvix(mc)).CMPDataCT = [];
        end
        stvar.Map(stvix(mc)).Name = ...
            sprintf('%s (border-smoothed %gmm)', map.Name, smkv);

        % outer border also!
        smvi = find(brde(:) & ~brd(:));
        [smvx, smvy, smvz] = ind2sub(size(brd), smvi);
        switch (stvtyp)
            case {'cmp'}
                stvar.Map(stvix(mc)).CMPData(smvi) = flexinterpn( ...
                    stvar.Map(stvix(mc)).CMPData, [smvx, smvy, smvz], smk, 1, 0);
            case {'hdr'}
                stvar.VoxelData(smvi + mapdatane * (stvix(mc) - 1)) = flexinterpn( ...
                    stvar.VoxelData(:, :, :, stvix(mc)), [smvx, smvy, smvz], smk, 1, 0);
            case {'head'}
                stvar.Brick(stvix(mc)).Data(smvi) = flexinterpn( ...
                    stvar.Brick(stvix(mc)).Data, [smvx, smvy, smvz], smk, 1, 0);
            case {'vmp'}
                stvar.Map(stvix(mc)).VMPData(smvi) = flexinterpn( ...
                    stvar.Map(stvix(mc)).VMPData, [smvx, smvy, smvz], smk, 1, 0);
        end

    % full dataset
    else

        % replace in target
        smk = smkv(1, [1, 1, 1]) ./ res;
        switch (stvtyp)
            case {'cmp'}
                stvar.Map(stvix(mc)).CMPData = single(smoothdata3(mapdata .* brde, smk));
                stvar.Map(stvix(mc)).CMPDataCT = [];
            case {'hdr'}
                stvar.VoxelData(:, :, :, stvix(mc)) = smoothdata3(mapdata .* brde, smk);
                stvar.VoxelDataCT{stvix(mc)} = [];
            case {'head'}
                stvar.Brick(stvix(mc)).Data = smoothdata3(mapdata .* brde, smk);
                stvar.Brick(stvix(mc)).DataCT = [];
            case {'vmp'}
                stvar.Map(stvix(mc)).VMPData = single(smoothdata3(mapdata .* brde, smk));
                stvar.Map(stvix(mc)).CMPDataCT = [];
        end
        stvar.Map(stvix(mc)).Name = ...
            sprintf('%s (smoothed %gmm)', map.Name, smkv);
    end
end

% set new number of maps
if any(strcmp(stvtyp, {'cmp', 'vmp'}))
    stvar.NrOfMaps = numel(stvar.Map);
end

% update
stvar.Browse;
