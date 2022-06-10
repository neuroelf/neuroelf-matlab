function srf = tom_ConvertToSRF(xo, opts)
% TOM::ConvertToSRF  - convert TOM to SRF format
%
% FORMAT:       srf = tom.ConvertToSRF([opts]);
%
% Input fields:
%
%       opts        optional settings
%        .center    1x3 center coordinate (source system, default: mean)
%        .flipaxes  1x1 boolean flag, flip axes (default: true)
%        .scale     1x1 scale to fit into regular box (default: 0.01)
%
% Output fields:
%
%       srf         SRF type object
%
% Note: uses emptysrf, mesh_trianglestoneighbors.

% Version:  v1.1
% Build:    22061014
% Date:     Jun-10 2022, 2:55 PM EST
% Author:   Jochen Weber, Memorial Sloan Kettering Cancer Center, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2022, Jochen Weber
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
emptysrf = ne_methods.emptysrf;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'tom'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'center') || ~isa(opts.center, 'double') || numel(opts.center) ~= 3 || ...
    any(isinf(opts.center) | isnan(opts.center))
    opts.center = [];
end
if ~isfield(opts, 'flipaxes') || ~islogical(opts.flipaxes) || numel(opts.flipaxes) ~= 1
    opts.flipaxes = true;
end
if ~isfield(opts, 'scale') || ~isa(opts.scale, 'double') || numel(opts.scale) ~= 1 || ...
    isinf(opts.scale) || isnan(opts.scale) || opts.scale <= 0 || opts.scale > 10
    opts.scale = 0.1;
end

% get content
bc = xo.C;

% get coordinates, normals, triangles
c = bc.VertexCoordinate;
t = bc.TriangleVertex(:, [1, 3, 2]);

% compute neighbors
try
    [nei, bn] = ne_methods.mesh_trianglestoneighbors(size(c, 1), t);
catch xfferror;
    neuroelf_lasterr(xfferror);
    rethrow(xfferror);
end

% createew SRF object
srf = emptysrf();

% apply settings
if isempty(opts.center)
    opts.center = mean(c);
end
c = c - ones(size(c, 1), 1) * opts.center(:)';
if opts.scale ~= 1
    c = opts.scale .* c;
end
if opts.flipaxes
    c = c(:, [3, 2, 1]);
    c(:, 2) = -c(:, 2);
end
c = c + 128;

% calc normals
n = -ne_methods.mesh_normals(c, t);

% set if srf
srf.C.NrOfVertices = size(c, 1);
srf.C.NrOfTriangles = size(t, 1);
srf.C.VertexCoordinate = c;
srf.C.VertexNormal = n;
srf.C.VertexColor = zeros(size(c, 1), 4);
srf.C.Neighbors = nei;
srf.C.TriangleVertex = t;
