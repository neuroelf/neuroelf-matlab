function xo = srf_AddSurfacePoints(xo, sfh, opts)
% SRF::AddSurfacePoints  - add points from an SFH file
%
% FORMAT:       [srf] = srf.AddSurfacePoints(sfh [, opts]);
%
% Input fields:
%
%       sfh         SFH object
%       opts        struct with optional fields
%        .color     1x3 or Px3 RGB color, override SFH settings
%        .size      1x1 or Px1 doulbe, override sizes from SFH
%
% Output fields:
%
%       srf         altered object
%
% Using: tfmatrix.

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:19 PM EST
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
if ~isfield(xffsngl.CONF.loadedobjs, 'srf_AddSurfacePoints')
    try
        ufn = [neuroelf_path('srf') '/unitsphere320.srf'];
        usrf = xff(ufn);
        usrfbc = usrf.C;
        if ~xffisobject(usrf, true, 'srf')
            error('BAD_SURFACE');
        end
        xffsngl.CONF.loadedobjs.srf_AddSurfacePoints = struct;
        xffsngl.CONF.loadedobjs.srf_AddSurfacePoints.unitsphere = usrfbc;
        delete(usrf);
    catch xfferror
        neuroelf_lasterr(xfferror);
        error('neuroelf:xff:fileNotFound', ...
            'Required unit sphere surface not found: ''%s''.', ufn);
    end
else
    try
        usrfbc = xffsngl.CONF.loadedobjs.srf_AddSurfacePoints.unitsphere;
    catch xfferror
        neuroelf_lasterr(xfferror);
        error('neuroelf:xff:internalError', ...
            'Unit sphere surface removed from global config.');
    end
end

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf') || ...
    numel(sfh) ~= 1 || ~xffisobject(sfh, true, 'sfh')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
srfbc = xo.C;
sfhbc = sfh.C;
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'color') || ~isa(opts.color, 'double') || size(opts.color, 2) ~= 3 || ...
   ~any(size(opts.color, 1) == [1, sfhbc.NrOfPoints]) || ...
    any(isinf(opts.color(:)) | isnan(opts.color(:)) | opts.color(:) < 0)
    opts.color = [];
else
    opts.color = min(255, opts.color(:, 1:3, 1));
    if all(opts.color(:) <= 1)
        opts.color = opts.color * 255;
    end
    opts.color = round(opts.color);
    if size(opts.color, 1) == 1
        opts.color = repmat(opts.color, [sfhbc.NrOfPoints, 1]);
    end
end
if ~isfield(opts, 'size') || ~isa(opts.size, 'double') || numel(opts.size) ~= length(opts.size) || ...
   ~any(numel(opts.size) == [1, sfhbc.NrOfPoints]) || ...
    any(isinf(opts.size) | isnan(opts.size) | opts.size <= 0)
    opts.size = [];
else
    opts.size = opts.size(:);
    opts.size(opts.size > 16) = 16;
    if numel(opts.size) == 1
        opts.size = repmat(opts.size, [sfhbc.NrOfPoints, 1]);
    end
end

% set to RGB color mode and make copy of unit sphere coordinates
usrfbc.VertexColor(:, 1) = NaN;
usrfc = usrfbc.VertexCoordinate;
nvert = size(usrfc, 1);
usrfc = usrfc - repmat(usrfbc.MeshCenter, [nvert, 1]);

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
fps = sfhbc.Fiducials;
ffps = fieldnames(fps);
fidt = zeros(numel(ffps), 7);
for sc = 1:numel(ffps)
    fidt(sc, :) = fps.(ffps{sc});
end
if ~isempty(opts.size)
    fidt(:, 4) = opts.size;
end
if ~isempty(opts.color)
    fidt(:, 5:7) = opts.color;
end

% get coordinates (in BV fashion)
fidp = -fidt(:, [2, 3, 1]);
fidp(:, 4) = 1;

% get matrix to and back from origin
o44 = [[eye(3), (  sfhbc.BVMidPoint)']; [0, 0, 0, 1]];
b44 = [[eye(3), ( -sfhbc.BVMidPoint)']; [0, 0, 0, 1]];
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
fidp = fidp * o44' * r44' * t44' * b44' * s44' * c44';
fidt(:, 1:3) = fidp(:, 1:3);

% first add all spheres to new srf
ffps = fieldnames(fps);
for sc = 1:numel(ffps)
    fid = fidt(sc, :);
    fidc = repmat(fid(1:3), [nvert, 1]);
    usrfbc.VertexCoordinate = fid(4) * usrfc + fidc;
    usrfbc.VertexColor(:, 2:4) = repmat(round(fid(5:7)), [nvert, 1]);
    cnv = size(nsrfbc.VertexCoordinate, 1);
    nsrfbc.VertexCoordinate = [nsrfbc.VertexCoordinate; usrfbc.VertexCoordinate];
    nsrfbc.VertexNormal = [nsrfbc.VertexNormal; usrfbc.VertexNormal];
    nsrfbc.VertexColor = [nsrfbc.VertexColor; usrfbc.VertexColor];
    nnei = [nsrfbc.Neighbors; usrfbc.Neighbors];
    nsrfbc.TriangleVertex = [nsrfbc.TriangleVertex; (usrfbc.TriangleVertex + cnv)];
    for vc = (cnv+1):size(nnei, 1)
        nnei{vc, 2} = nnei{vc, 2} + cnv;
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
