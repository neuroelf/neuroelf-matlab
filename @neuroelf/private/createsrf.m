function srf = createsrf(coord, tri, opts)
% createsrf  - create a SRF object from coordinates and triangles
%
% FORMAT:       srf = createsrf(coord, tri [, opts])
%
% Input fields:
%
%       coord       Vx3 coordinates of V vertices
%       tri         Tx3 triangles with indices 1...V
%       opts        optional settings in 1x1 struct
%        .axes      1x3 array, axes order and dirs, default: [-2, -3, -1]
%        .offset    1x3 array, added to the coords, default: [128,128,128]
%        .trf       4x4 transformation matrix to apply, default: none
%
% Output fields:
%
%       srf         surface object
%
% Note: internally calls mesh_trianglestoneighbors and SRF::RecalcNormals.

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:15 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% argument check
if nargin < 2 || ...
   ~any(strcmpi(class(coord), {'single', 'double'})) || ...
    size(coord, 2) ~= 3 || ...
    isempty(coord) || ...
    ndims(coord) > 2 || ...
    any(isinf(coord(:)) | isnan(coord(:))) || ...
   ~isa(tri, 'double') || ...
    size(tri, 2) ~= 3 || ...
    isempty(tri) || ...
    ndims(tri) > 2 || ...
    any(isinf(tri(:)) | isnan(tri(:)) | tri(:) < 1 | tri(:) > size(coord, 1))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument.' ...
    );
end

% check options
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'axes') || ...
   ~isa(opts.axes, 'double') || ...
    numel(opts.axes) ~= 3 || ...
    any(isinf(opts.axes) | isnan(opts.axes) | opts.axes < -3 | opts.axes > 3 | opts.axes == 0) || ...
    numel(unique(round(abs(opts.axes)))) ~= 3
    opts.axes = -[2, 3, 1];
end
opts.axsg = sign(opts.axes(:)');
opts.axes = round(abs(opts.axes(:)'));
if ~isfield(opts, 'offset') || ...
   ~isa(opts.offset, 'double') || ...
    numel(opts.offset) ~= 3 || ...
    any(isinf(opts.offset) | isnan(opts.offset) | abs(opts.offset) > 256)
    opts.offset = [128, 128, 128];
end
if ~isfield(opts, 'trf') || ...
   ~isa(opts.trf, 'double') || ...
   ~isequal(size(opts.trf), [4, 4]) || ...
    any(isinf(opts.trf(:)) | isnan(opts.trf(:))) || ...
    opts.trf(4, 4) ~= 1
    opts.trf = [];
end

% try to create neighbors list
try
    n = mesh_trianglestoneighbors(size(coord, 1), tri);
catch ne_eo;
    error( ...
        'neuroelf:InternalError', ...
        'Error creating neighbors list from triangles: %s.', ...
        ne_eo.message ...
    );
end

% create SRF
srf = xff('new:srf');

% create vertices according to axes
srf.NrOfVertices = size(coord, 1);
srf.NrOfTriangles = size(tri, 1);
srf.VertexCoordinate = [ ...
    opts.offset(1) + opts.axsg(1) * coord(:, opts.axes(1)), ...
    opts.offset(2) + opts.axsg(2) * coord(:, opts.axes(2)), ...
    opts.offset(3) + opts.axsg(3) * coord(:, opts.axes(3))];
srf.VertexNormal = zeros(size(coord, 1), 3);
srf.VertexColor = zeros(size(coord, 1), 4);
srf.Neighbors = n;
srf.TriangleVertex = tri;

% normals computation
srf.RecalcNormals;
