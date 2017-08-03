% FUNCTION ne_render_smoothstats: smooth stats border (high-q output)
function ne_render_smoothstats(varargin)

% Version:  v0.9d
% Build:    14061710
% Date:     Jun-17 2014, 10:35 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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

% check global fig
if ~isstruct(ne_gcfg.fcfg.Render) || ...
    numel(ne_gcfg.fcfg.Render) ~= 1
    return;
else
    o = ne_gcfg.fcfg.Render;
end

% check options
if ~isfield(o, 'smstat') || ...
   ~isfield(o, 'stvar') || ...
   ~isfield(o, 'stvix') || ...
   ~isfield(o, 'stvixo')
    return;
end

% this method requires a valid stats var
if numel(o.stvar) ~= 1 || ...
   ~isxff(o.stvar, {'cmp', 'hdr', 'head', 'vmp'})
    return;
end

% off-load work
o.stvix = ne_smoothstats(0, 0, o.smstat, o.stvar, o.stvix, o.stvixo, true);

% % for each map
% for mc = 1:numel(stvixo)
%
%     % get map, data, and clustering information
%     map = stvar.Map(stvixo(mc));
%     switch (stvtyp)
%         case {'cmp'}
%             mapdata = map.CMPData;
%             mapdatact = map.CMPDataCT;
%         case {'hdr'}
%             mapdata = stvar.VoxelData(:, :, :, stvixo(mc));
%             if numel(stvar.VoxelDataCT) >= stvixo(mc)
%                 mapdatact = stvar.VoxelDataCT{stvixo(mc)};
%             else
%                 mapdatact = logical([]);
%             end
%             mapdatane = numel(mapdata);
%         case {'head'}
%             mapdata = stvar.Brick(stvixo(mc)).Data;
%             mapdatact = stvar.Brick(stvixo(mc)).DataCT;
%         case {'vmp'}
%             mapdata = map.VMPData;
%             mapdatact = map.VMPDataCT;
%     end
%
%     % get border
%     if ~isempty(mapdatact);
%         brde = mapdatact;
%         brd = erode3d(brde);
%     else
%         brd = (mapdata ~= 0);
%         brde = dilate3d(brd);
%     end
%     brdd = erode3d(brd);
%
%     % find voxels to smooth
%     smvi = find(brd(:) & ~brdd(:));
%     [smvx, smvy, smvz] = ind2sub(size(brd), smvi);
%
%     % copy map ?
%     if stvix(mc) == stvixo(mc)
%         stvix(mc) = numel(stvar.Map) + 1;
%         stvar.Map(stvix(mc)) = map;
%         stvar.Map(stvix(mc)).EnableClusterCheck = 0;
%     end
%
%     % replace in target
%     switch (stvtyp)
%         case {'cmp'}
%             stvar.Map(stvix(mc)).CMPData = mapdata .* brde;
%             stvar.Map(stvix(mc)).CMPData(smvi) = flexinterpn(mapdata, ...
%                 [smvx, smvy, smvz], smk, 1, 0);
%             stvar.Map(stvix(mc)).CMPDataCT = [];
%         case {'hdr'}
%             stvar.VoxelData(:, :, :, stvix(mc)) = mapdata .* brde;
%             stvar.VoxelData(smvi + mapdatane * (stvix(mc) - 1)) = ...
%                 flexinterpn(mapdata, [smvx, smvy, smvz], smk, 1, 0);
%             stvar.VoxelDataCT{stvix(mc)} = [];
%         case {'head'}
%             stvar.Brick(stvix(mc)).Data = mapdata .* brde;
%             stvar.Brick(stvix(mc)).Data(smvi) = flexinterpn(mapdata, ...
%                 [smvx, smvy, smvz], smk, 1, 0);
%             stvar.Brick(stvix(mc)).DataCT = [];
%         case {'vmp'}
%             stvar.Map(stvix(mc)).VMPData = mapdata .* brde;
%             stvar.Map(stvix(mc)).VMPData(smvi) = flexinterpn(mapdata, ...
%                 [smvx, smvy, smvz], smk, 1, 0);
%             stvar.Map(stvix(mc)).CMPDataCT = [];
%     end
%     stvar.Map(stvix(mc)).Name = ...
%         sprintf('%s (smoothed %gmm)', map.Name, smkv);
%
%     % outer border also!
%     smvi = find(brde(:) & ~brd(:));
%     [smvx, smvy, smvz] = ind2sub(size(brd), smvi);
%     switch (stvtyp)
%         case {'cmp'}
%             stvar.Map(stvix(mc)).CMPData(smvi) = flexinterpn( ...
%                 stvar.Map(stvix(mc)).CMPData, [smvx, smvy, smvz], smk, 1, 0);
%         case {'hdr'}
%             stvar.VoxelData(smvi + mapdatane * (stvix(mc) - 1)) = flexinterpn( ...
%                 stvar.VoxelData(:, :, :, stvix(mc)), [smvx, smvy, smvz], smk, 1, 0);
%         case {'head'}
%             stvar.Brick(stvix(mc)).Data(smvi) = flexinterpn( ...
%                 stvar.Brick(stvix(mc)).Data, [smvx, smvy, smvz], smk, 1, 0);
%         case {'vmp'}
%             stvar.Map(stvix(mc)).VMPData(smvi) = flexinterpn( ...
%                 stvar.Map(stvix(mc)).VMPData, [smvx, smvy, smvz], smk, 1, 0);
%     end
% end
%
% % set new number of maps
% if any(strcmp(stvtyp, {'cmp', 'vmp'}))
%     stvar.NrOfMaps = numel(stvar.Map);
% end

% keep track of new assignment
ne_gcfg.fcfg.Render.stvix = o.stvix;

% re-open file (to update list of map names)
ne_openfile(0, 0, o.stvar);
ne_setcstatmap(0, 0, o.stvix);
