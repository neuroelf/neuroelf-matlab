function xo = srf_AddDipoles(xo, bsa, sfh, opts)
% SRF::AddSurfacePoints  - add dipoles from a BSA file
%
% FORMAT:       [srf] = srf.AddDipoles(bsa, sfh [, opts]);
%
% Input fields:
%
%       bsa         BSA object
%       sfh         SFH object (needed for coordinate transformation)
%       opts        struct with optional fields
%        .dipcolor  Dx3 Colors components for SngDip points
%        .dipsize   Dx1 Sizes for SngDip points
%        .srccolor  Sx3 Colors components for RegSrc points
%        .srcsize   Sx1 Sizes for RegSrc points
%
% Output fields:
%
%       srf         altered object
%
% Using: tfmatrix.

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:31 PM EST
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

% neuroelf library and global storage to load unit sphere
global ne_methods xffsngl;
if ~isfield(xffsngl.CONF.loadedobjs, 'srf_AddDipoles')
    cllist = cell(1, 2);
    try
        ufn = neuroelf_path('srf');
        usp = [ufn '/unitsphere320.srf'];
        ulg = [ufn '/unitcylinder4280.srf'];
        usp = xff(usp);
        cllist{1} = usp;
        ulg = xff(ulg);
        cllist{2} = ulg;
        if ~xffisobject(usp, true, 'srf') || ~xffisobject(ulg, true, 'srf')
            error('BAD_SURFACE');
        end
        uspbc = usp.C;
        ulgbc = ulg.C;
        xffsngl.CONF.loadedobjs.srf_AddDipoles = struct;
        xffsngl.CONF.loadedobjs.srf_AddDipols.unitsphere = uspbc;
        xffsngl.CONF.loadedobjs.srf_AddDipols.unitcylinder = ulgbc;
        delete(usp);
        delete(ulg);
    catch xfferror
        neuroelf_lasterr(xfferror);
        clearxffobjects(cllist);
        error('neuroelf:xff:fileNotFound', 'Required element surface not found: ''%s''.', ufn);
    end
else
    try
        uspbc = xffsngl.CONF.loadedobjs.srf_AddDipoles.unitsphere;
        ulgbc = xffsngl.CONF.loadedobjs.srf_AddDipoles.unitcylinder;
    catch xfferror
        neuroelf_lasterr(xfferror);
        error('neuroelf:xff:internalError', 'Element surface(s) removed from global config.');
    end
end

% check arguments
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf') || ...
    numel(bsa) ~= 1 || ~xffisobject(bsa, true, 'bsa') || ...
    numel(sfh) ~= 1 || ~xffisobject(sfh, true, 'sfh')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
srfbc = xo.C;
bsabc = bsa.C;
sfhbc = sfh.C;
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dipcolor') || ~isa(opts.dipcolor, 'double') || size(opts.dipcolor, 2) ~= 3 || ...
   ~any(size(dipopts.color, 1) == [1, bsabc.NrOfSngDip]) || ...
    any(isinf(opts.dipcolor(:)) | isnan(opts.dipcolor(:)) | opts.dipcolor(:) < 0)
    opts.dipcolor = [];
else
    opts.dipcolor = min(255, opts.dipcolor(:, 1:3, 1));
    if all(opts.dipcolor(:) <= 1)
        opts.dipcolor = opts.dipcolor * 255;
    end
    opts.dipcolor = round(opts.dipcolor);
    if size(opts.dipcolor, 1) == 1
        opts.dipcolor = repmat(opts.dipcolor, [sfhbc.NrOfSngDip, 1]);
    end
end
if ~isfield(opts, 'dipsize') || ~isa(opts.dipsize, 'double') || numel(opts.dipsize) ~= length(opts.dipsize) || ...
   ~any(numel(opts.dipsize) == [1, sfhbc.NrOfSngDip]) || ...
    any(isinf(opts.dipsize) | isnan(opts.dipsize) | opts.dipsize <= 0)
    opts.dipsize = [];
else
    opts.dipsize = opts.dipsize(:);
    opts.dipsize(opts.dipsize > 16) = 16;
    if numel(opts.dipsize) == 1
        opts.dipsize = repmat(opts.dipsize, [sfhbc.NrOfSngDip, 1]);
    end
end
if ~isfield(opts, 'srccolor') || ~isa(opts.srccolor, 'double') || size(opts.srccolor, 2) ~= 3 || ...
   ~any(size(srcopts.color, 1) == [1, bsabc.NrOfRegSrc]) || ...
    any(isinf(opts.srccolor(:)) | isnan(opts.srccolor(:)) | opts.srccolor(:) < 0)
    opts.srccolor = [];
else
    opts.srccolor = min(255, opts.srccolor(:, 1:3, 1));
    if all(opts.srccolor(:) <= 1)
        opts.srccolor = opts.srccolor * 255;
    end
    opts.srccolor = round(opts.srccolor);
    if size(opts.srccolor, 1) == 1
        opts.srccolor = repmat(opts.srccolor, [sfhbc.NrOfRegSrc, 1]);
    end
end
if ~isfield(opts, 'srcsize') || ~isa(opts.srcsize, 'double') || numel(opts.srcsize) ~= length(opts.srcsize) || ...
   ~any(numel(opts.srcsize) == [1, sfhbc.NrOfRegSrc]) || ...
    any(isinf(opts.srcsize) | isnan(opts.srcsize) | opts.srcsize <= 0)
    opts.srcsize = [];
else
    opts.srcsize = opts.srcsize(:);
    opts.srcsize(opts.srcsize > 16) = 16;
    if numel(opts.srcsize) == 1
        opts.srcsize = repmat(opts.srcsize, [sfhbc.NrOfRegSrc, 1]);
    end
end

% set to RGB color mode and make copy of unit sphere coordinates
uspbc.VertexColor(:, 1) = NaN;
ulgbc.VertexColor(:, 1) = NaN;
uspc = uspbc.VertexCoordinate;
ulgc = ulgbc.VertexCoordinate;
ulgn = ulgbc.VertexNormal;
nvspc = size(uspc, 1);
nvlgc = size(ulgc, 1);
uspc = uspc - repmat(uspbc.MeshCenter, [nvspc, 1]);
ulgc = ulgc - repmat(ulgbc.MeshCenter, [nvlgc, 1]);

% make new srf
nsrf = xff('new:srf');
nsrfbc = nsrf.C;
nsrfbc.NrOfVertices = 0;
nsrfbc.NrOfTriangles = 0;
nsrfbc.MeshCenter = srfbc.MeshCenter;
nsrfbc.VertexCoordinate = zeros(0, 3);
nsrfbc.VertexNormal = zeros(0, 3);
nsrfbc.VertexColor = zeros(0, 4);
nsrfbc.Neighbors = cell(0, 2);
nsrfbc.TriangleVertex = zeros(0, 3);
nsrfbc.NrOfTriangleStrips = 0;
nsrfbc.TriangleStripSequence = zeros(0, 1);

% get coordinate table of fiducials
bsa = bsabc;
dipt = zeros(numel(bsa.SngDip), 10);
srct = zeros(numel(bsa.RegSrc), 7);
for sc = 1:size(dipt, 1)
    dipt(sc, :) = [bsa.SngDip(sc).Position, bsa.SngDip(sc).Size, bsa.SngDip(sc).Color, bsa.SngDip(sc).Orientation];
end
for sc = 1:size(srct, 1)
    srct(sc, :) = [bsa.RegSrc(sc).Position, bsa.RegSrc(sc).Size, bsa.RegSrc(sc).Color];
end
if ~isempty(opts.dipsize)
    dipt(:, 4) = opts.dipsize;
end
if ~isempty(opts.dipcolor)
    dipt(:, 5:7) = opts.dipcolor;
end
if ~isempty(opts.srcsize)
    srct(:, 4) = opts.srcsize;
end
if ~isempty(opts.srccolor)
    srct(:, 5:7) = opts.srccolor;
end

% get coordinates (in BV fashion)
dipp = -dipt(:, [2, 3, 1]);
dipo = -dipt(:, [9,10, 8]);
srcp = -srct(:, [2, 3, 1]);
dipp(:, 4) = 1;
srcp(:, 4) = 1;

% get matrix to and back from origin
o44 = [[eye(3), (sfhbc.BVMidPoint)']; [0, 0, 0, 1]];
b44 = [[eye(3), (-sfhbc.BVMidPoint)']; [0, 0, 0, 1]];
c44 = [[eye(3), (srfbc.MeshCenter)']; [0, 0, 0, 1]];

% get transformation matrices from SFH
bvt = sfhbc.BVTrf;
bvt(4:6) = bvt(4:6) .* (pi / 180);
r44 = ne_methods.tfmatrix(struct('type', {'r', 'r', 'r'}, ...
    'xyz',  {[-bvt(4), 0, 0], [0, bvt(5), 0], [0, 0, -bvt(6)]}));
t44 = [[eye(3), (bvt(1:3))']; [0, 0, 0, 1]];
s44 = eye(4) * diag([bvt(7:9), 1]);

% apply translation (rotation around (0,0,0), translation,
% then scaling at original midpoint and translation to mesh center
dipp = dipp * o44' * r44' * t44' * b44' * s44' * c44';
srcp = srcp * o44' * r44' * t44' * b44' * s44' * c44';
dipt(:, 1:3) = dipp(:, 1:3);
srct(:, 1:3) = srcp(:, 1:3);
dipo = dipo * r44(1:3, 1:3)';

% first add all RegSrc's to new srf
for sc = 1:size(srct, 1)
    fid = srct(sc, :);
    fidc = repmat(fid(1:3), [nvspc, 1]);
    uspbc.VertexCoordinate = fid(4) * uspc + fidc;
    uspbc.VertexColor(:, 2:4) = repmat(round(fid(5:7)), [nvspc, 1]);
    cnv = size(nsrfbc.VertexCoordinate, 1);
    nsrfbc.VertexCoordinate = [nsrfbc.VertexCoordinate; uspbc.VertexCoordinate];
    nsrfbc.VertexNormal = [nsrfbc.VertexNormal; uspbc.VertexNormal];
    nsrfbc.VertexColor = [nsrfbc.VertexColor; uspbc.VertexColor];
    nnei = [nsrfbc.Neighbors; uspbc.Neighbors];
    nsrfbc.TriangleVertex = [nsrfbc.TriangleVertex; (uspbc.TriangleVertex + cnv)];
    for vc = (cnv+1):size(nnei, 1)
        nnei{vc, 2} = nnei{vc, 2} + cnv;
    end
    nsrfbc.Neighbors = nnei;
end

% then add all SngDip's to new srf
for sc = 1:size(dipt, 1)
    fid = dipt(sc, :);
    fidl = repmat(fid(1:3), [nvlgc, 1]);
    fidc = repmat(fid(1:3), [nvspc, 1]);

    % get orientation right
    ori = dipo(sc, :);
    ori = ori ./ sqrt(sum(ori .* ori));
    if any(ori == 0)
        nori1 = double(ori == 0);
    else
        nori1 = [-ori(2), ori(1), 0];
    end
    nori1 = nori1 ./ sqrt(sum(nori1 .* nori1));
    nori2 = cross(ori, nori1);
    nori2 = nori2 ./ sqrt(sum(nori2 .* nori2));
    n33 = [ori(:), nori1(:), nori2(:)];

    ulgbc.VertexCoordinate = (2.5 * fid(4) * ulgc) * n33' + fidl;
    ulgbc.VertexNormal = ulgn * n33';
    ulgbc.VertexColor(:, 2:4) = repmat(round(fid(5:7)), [nvlgc, 1]);
    uspbc.VertexCoordinate = fid(4) * uspc + fidc;
    uspbc.VertexColor(:, 2:4) = repmat(round(fid(5:7)), [nvspc, 1]);
    cnv = size(nsrfbc.VertexCoordinate, 1);
    nsrfbc.VertexCoordinate = [nsrfbc.VertexCoordinate; ulgbc.VertexCoordinate; uspbc.VertexCoordinate];
    nsrfbc.VertexNormal = [nsrfbc.VertexNormal; ulgbc.VertexNormal; uspbc.VertexNormal];
    nsrfbc.VertexColor = [nsrfbc.VertexColor; ulgbc.VertexColor; uspbc.VertexColor];
    nnei = [nsrfbc.Neighbors; ulgbc.Neighbors; uspbc.Neighbors];
    nsrfbc.TriangleVertex = [nsrfbc.TriangleVertex; (ulgbc.TriangleVertex + cnv); (uspbc.TriangleVertex + cnv + nvlgc)];
    for vc = (cnv+1):(cnv+nvlgc)
        nnei{vc, 2} = nnei{vc, 2} + cnv;
    end
    for vc = (cnv+nvlgc):size(nnei, 1)
        nnei{vc, 2} = nnei{vc, 2} + cnv + nvlgc;
    end
    nsrfbc.Neighbors = nnei;
end

% make nsrf settings
nsrfbc.NrOfVertices = size(nsrfbc.VertexNormal, 1);
nsrfbc.NrOfTriangles = size(nsrfbc.TriangleVertex, 1);
nsrf.C = nsrfbc;

% add spheres to srf
nfile = srf_Combine(xo, nsrf, struct('type', 'wholebrain'));
delete(nsrf);
xo.C = nfile.C;
delete(nfile);
