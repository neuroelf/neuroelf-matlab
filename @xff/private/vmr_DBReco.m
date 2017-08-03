function srf = vmr_DBReco(xo, opts)
% VMR::DBReco  - direct border reconstruction
%
% FORMAT:       srf = vmr.DBReco([opts])
%
% Input fields:
%
%       opts        optional settings
%        .autoimp   attempt automatic improvement (default: false)
%        .autoivis  visualize improvement attempts (default: true)
%        .onesurf   only create surface for biggest object (default: true)
%        .scol      segmentation color (default: 240)
%        .scolb     border color (default: 235)
%        .tps       triangles-per-surface, {2} or 4
%        .warn      print warnings (default: true)
%
% Output fields:
%
%       srf         reconstructed surface
%
% Using: clustercoordsc, findfirst, floodfill3c, lsqueeze, maxpos,
%        mesh_reconstruct, smoothdata3.

% Version:  v1.1
% Build:    16021413
% Date:     Feb-14 2016, 1:03 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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
clustercoordsc = ne_methods.clustercoordsc;
emptysrf       = ne_methods.emptysrf;
maxpos         = ne_methods.maxpos;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'autoimp') || ~islogical(opts.autoimp) || numel(opts.autoimp) ~= 1
    opts.autoimp = false;
end
if ~isfield(opts, 'autoivis') || ~islogical(opts.autoivis) || numel(opts.autoivis) ~= 1
    opts.autoivis = true;
end
if ~isfield(opts, 'onesurf') || ~islogical(opts.onesurf) || numel(opts.onesurf) ~= 1
    opts.onesurf = true;
end
if ~isfield(opts, 'scol') || numel(opts.scol) ~= 1 || ~isnumeric(opts.scol) || ...
    isinf(opts.scol) || isnan(opts.scol) || opts.scol < 226 || opts.scol > 255 || ...
    opts.scol ~= fix(opts.scol)
    opts.scol = uint8(240);
else
    opts.scol = uint8(opts.scol);
end
if ~isfield(opts, 'scolb') || numel(opts.scolb) ~= 1 || ~isnumeric(opts.scolb) || ...
    isinf(opts.scolb) || isnan(opts.scolb) || opts.scolb < 226 || opts.scolb > 255 || ...
    opts.scolb ~= fix(opts.scolb)
    opts.scolb = uint8(235);
else
    opts.scolb = uint8(opts.scolb);
end
if ~isfield(opts, 'tps') || ~isa(opts.tps, 'double') || numel(opts.tps) ~= 1 || ...
    isinf(opts.tps) || isnan(opts.tps) || ~any(opts.tps == [2, 4])
    opts.tps = 2;
end
if ~isfield(opts, 'warn') || ~isa(opts.warn, 'logical') || numel(opts.warn) ~= 1
    opts.warn = true;
end

% get VMR content
bc = xo.C;
vd = bc.VMRData(:, :, :);
if ~isa(bc.VMRData, 'uint8')
    error('neuroelf:xff:badObject', 'DBReco method only valid for 8-bit VMR data.');
end

% for auto-improvement, all resolutions must be equal
if opts.autoimp && (bc.VoxResX ~= bc.VoxResY || bc.VoxResX ~= bc.VoxResZ)
    error('neuroelf:xff:invalidObject', ...
        'Auto-improvement requires all dimensions having equal resolution.');
end

% un-prepare (border detection is implemented in MEX-file!)
vd(vd == opts.scolb) = opts.scol;
vd = (vd == opts.scol);

% make sure to get biggest chunk
if opts.onesurf
    [cs, cv] = clustercoordsc(vd);

    % nothing in the set?
    if isempty(cs)

        % return empty surface
        srf = emptysrf();
        return;
    end

    % remove dust and speckles
    if numel(cs) > 1
        vd = (cv == maxpos(cs));
    end
end

% but also do that with the background!
bgv = ne_methods.findfirst(ne_methods.lsqueeze(~vd(:, :, 1)));
if isempty(bgv)
    bgv = ne_methods.findfirst(ne_methods.lsqueeze(~vd(:, :, end)));
    if ~isempty(bgv)
        bgv = bgv + size(vd, 1) * size(vd, 2) * (size(vd, 3) - 1);
    end
end
if ~isempty(bgv)
    [bgx, bgy, bgz] = ind2sub(size(vd), bgv);
    [cv, cs] = ne_methods.floodfill3c(~vd, [bgx, bgy, bgz], 2);
    if cs < (0.2 * numel(vd))
        bgv = [];
    end
end
if ~isempty(bgv)
    vd(~cv) = true;
else
    [cs, cv] = clustercoordsc(~vd);
    if numel(cs) > 1

        % set false "holes" in vd to true!
        vd(cv > 0 & cv ~= maxpos(cs)) = true;
    end
end

% now we can run the reconstruction piece
[p, t, n] = ne_methods.mesh_reconstruct(vd, opts.tps);
cc = size(p, 1);
tc = size(t, 1);

% add offset, if necessary
o = [bc.OffsetX, bc.OffsetY, bc.OffsetZ];
if any(o ~= 0)

    % and subtract one, as we require BV-style positions!
    o = o - 1;
    p = p + ones(cc, 1) * o;

% no offset
else

    % compute correct points position -> BV starts at voxel (0, 0, 0)
    o = [-1, -1, -1];
    p = p - 1;
end

% and resolution
r = [bc.VoxResX, bc.VoxResY, bc.VoxResZ];
if any(r ~= 1)
    p = p .* (ones(cc, 1) *  r);
end

% compute Euler characteristic
eul = cc - 0.5 * tc;
if eul ~= 2 && opts.warn
    warning('neuroelf:xff:badObject', ...
        'Probably %d handles/bridges in surface remaining.', round((2 - eul) / 2));
end

% build surface
srf = xff('new:srf');

% fill with vertices and triangles (also calculates normals)
srf_SetCoordsTriangles(srf, p, t, struct('neigh', {n}));

% auto-improvement ?
if opts.autoimp

    % unset autoimp for further (recursive) calls, disable warnings
    opts.autoimp = false;
    opts.warn = false;

    % resolution and offset
    vr = bc.VoxResX;
    vo = -o;

    % store original data in UndoBuffer
    bc.RunTimeVars.UndoBuffer = bc.VMRData(:, :, :);

    % replace with clustered version
    bc.VMRData = uint8(opts.scol) .* uint8(vd);
    xo.C = bc;

    % get SRF content
    srfc = srf.C;

    % and keep backups!
    vdo = vd;
    srfcb = srfc;

    % bridges/handles ?
    hrc = 10;
    if eul ~= 2

        % max counter set to 10
        while hrc > 0

            % visualize SRF
            if opts.autoivis
                aft_Browse(srf);
                drawnow;
            end

            % get original coordinates
            soc = srfc.VertexCoordinate;

            % apply smoothing (gradually more)
            srf_Smooth(srf, 200 - 10 * hrc, 0.05);

            % get content again
            srfc = srf.C;

            % then get new points and normals
            sc = srfc.VertexCoordinate;
            sn = srfc.VertexNormal;

            % just in case, normalize normals
            sn = ((1 ./ sqrt(sum(sn .* sn, 2))) * ones(1, 3)) .* sn;

            % compute version that has one element = 1 (or -1) + length/2
            sno = ((1 ./ max(abs(sn), [], 2)) * ones(1, 3)) .* sn;
            snol = (0.5 * vr) .* sqrt(sum(sno .* sno, 2));

            % compute movement (along normals)
            scd = sc - soc;
            snd = sum(sn .* scd, 2);

            % consider voxels which at least moved to another voxel
            sndin = (snd > snol);
            sndout = (snd < -snol);

            % get vertex positions for moved voxels (plus/minus sno!)
            cin = sc(sndin, :) - sno(sndin, :);
            cout = sc(sndout, :) + sno(sndout, :);

            % and compute the voxel position
            vin = round(ones(size(cin, 1), 1) * vo + (1 / vr) .* cin);
            vout = round(ones(size(cout, 1), 1) * vo + (1 / vr) .* cout);

            % and their actual indices
            vin = sub2ind(size(vd), vin(:, 1), vin(:, 2), vin(:, 3));
            vout = sub2ind(size(vd), vout(:, 1), vout(:, 2), vout(:, 3));

            % for vertices that are in both
            vinout = intersect(vin, vout);

            % remove from either set
            vin = setdiff(vin, vinout);
            vout = setdiff(vout, vinout);

            % then set to new values
            vd(vin) = false;
            vd(vout) = true;

            % and recluster
            [cs, cv] = clustercoordsc(vd);
            if numel(cs > 1)
                vd = (cv == maxpos(cs));
            end
            [cs, cv] = clustercoordsc(~vd);
            if numel(cs) > 1
                vd(cv > 0 & cv ~= maxpos(cs)) = true;
            end

            % visualize old VMR
            if opts.autoivis
                aft_Browse(xo);
                drawnow;
            end

            % re-set in VMRData
            bc.VMRData = uint8(opts.scol) .* uint8(vd);
            xo.C = bc;

            % visualize updated VMR
            if opts.autoivis
                aft_Browse(xo);
                drawnow;
            end

            % *unvisualize SRF*
            if opts.autoivis
                aft_UnBrowse(srf);
            end

            % then re-do the surface
            nsrf = vmr_DBReco(xo, opts);
            nsrfc = nsrf.C;

            % try to reset in old object
            try
                srf.C = nsrfc;

                % and clear new object
                aft_ClearObject(nsrf);

            % if that fails (e.g. closing destroyed the object!)
            catch xfferror
                neuroelf_lasterr(xfferror);

                % simply keep new object
                srf = nsrf;
            end

            % recompute euler characteristic
            srfc = srf.C;
            eul = size(srfc.VertexCoordinate, 1) - 0.5 * size(srfc.TriangleVertex, 1);

            % end this loop already?
            if eul == 2
                break;
            end

            % count down attempts
            hrc = hrc - 1;
        end
    end

    % now create a more smoothed version
    if opts.autoivis
        aft_Browse(srf);
    end

    % record current (unsmoothed) coordinates
    sc = srfc.VertexCoordinate;

    % smooth a bit more
    srf_Smooth(srf, 250, 0.25);

    % and create density map
    dsmp = srf_DensityMap(srf);
    dsmpc = dsmp.C;
    aft_ClearObject(dsmp);

    % *unvisualize SRF*
    if opts.autoivis
        aft_UnBrowse(srf);
    end

    % put original content back in place
    try
        srf.C = srfcb;
    catch xfferror
        neuroelf_lasterr(xfferror);
        srf = xff('new:srf');
        srf.C = srfcb;
    end

    % evaluate density map
    dsmp = dsmpc.Map(1).SMPData;
    dsmp = (isinf(dsmp) | isnan(dsmp) | dsmp > (2 / vr));

    % get those vertices
    dc = sc(dsmp, :);
    ndc = size(dc, 1);

    % and compute potential voxel coordinates
    vc = zeros(8 * ndc, 3);
    vct = 1;
    for xd = 0.5 .* [-vr, vr]
        for yd = 0.5 .* [-vr, vr]
            for zd = 0.5 .* [-vr, vr]
                vc(vct:vct+ndc-1, :) = round(ones(ndc, 1) * vo + (1 / vr) .* ...
                    [dc(:, 1) + xd, dc(:, 2) + yd, dc(:, 3) + zd]);
            end
        end
    end
    vds = size(vd);
    vc(any(vc < 1, 2) | vc(:, 1) > vds(1) | vc(:, 2) > vds(2) | vc(:, 3) > vds(3), :) = [];
    vc = unique(sub2ind(size(vd), vc(:, 1), vc(:, 2), vc(:, 3)));

    % now prepare final set!
    nvd = uint8(opts.scol) .* uint8(vdo);

    % if the bridge/handle detection was run
    if hrc < 10

        % mark all voxels from the bridge detection in green
        brd = vd & ~vdo;
        nvd(brd) = 245;

        % and find further potential candidates (still 0 but smoothed > 1/3)
        brds = (ne_methods.smoothdata3(single(brd), [3, 3, 3]) >= (1 / 3));
        nvd(vc(vd(vc) == 0 & brds(vc))) = 243;

        % then mark outside voxels from bridge detection in bright orange
        nvd(~vd & vdo) = 233;
    end

    % finally, mark all voxels that are segmented but shouldn't be orange
    nvd(vc(nvd(vc) == 240)) = 231;

    % and mark voxels that might better be off segmented dark blue
    brds = (ne_methods.smoothdata3(single(vdo), [3, 3, 3]) >= 0.25);
    nvd(vc(nvd(vc) == 0 & brds(vc))) = 236;
    bc.VMRData = nvd;
    xo.C = bc;

    % update display again?
    if opts.autoivis
        aft_Browse(xo);
    end
end
