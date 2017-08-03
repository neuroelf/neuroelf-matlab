function nsph = srf_ApplyTSM(xo, tsm, sphere, opts)
% SRF::ApplyTSM  - apply TSM file on vertex coordinates
%
% FORMAT:       sph = srf.ApplyTSM(tsm, sphere [, opts])
%
% Input fields:
%
%       srf         folded mesh (e.g. RECOSM mesh)
%       tsm         TSM mapping to sphere mesh
%       sphere      sphere mesh SRF used for mapping
%       opts        optional settings
%       .colors     either {'main'} or 'mix'
%       .smforce    smoothing force (default: 0.01)
%       .smiter     number of smoothing steps (default = 5)
%
% Output fields:
%
%       sph         SPH mesh with NrOfVertices of sphere, folded as srf
%
% Using: mesh_morph.

% Version:  v1.1
% Build:    16031813
% Date:     Mar-18 2016, 1:10 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% check arguments
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'fsbf', 'srf'}) || ...
    numel(tsm) ~= 1 || ~xffisobject(tsm, true, 'tsm') || ...
    numel(sphere) ~= 1 || ~xffisobject(sphere, true, {'fsbf', 'srf'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'colors') || ~ischar(opts.colors) || isempty(opts.colors) || ...
   ~any(strcmpi(opts.colors(:)', {'main', 'mix'}))
    opts.colors = 'main';
else
    opts.colors = lower(opts.colors(:)');
end
if ~isfield(opts, 'smforce') || numel(opts.smforce) ~= 1 || ~isa(opts.smforce, 'double') || ...
    isinf(opts.smforce) || isnan(opts.smforce) || opts.smforce <= 0 || opts.smforce >= 2.5
    opts.smforce = 0.01;
end
if ~isfield(opts, 'smiter') || numel(opts.smiter) ~= 1 || ~isa(opts.smiter, 'double') || ...
    isinf(opts.smiter) || isnan(opts.smiter) || opts.smiter < 0 || opts.smiter > 3000
    opts.smiter = 5;
else
    opts.smiter = round(opts.smiter);
end

% get contents
fldc = xo.C;
tsmc = tsm.C;
sphc = sphere.C;

% check numbers
if numel(tsmc.SourceTriangleOfTarget) ~= size(sphc.VertexCoordinate, 1)
    error('neuroelf:xff:badArgument', 'TSM and Sphere mesh must match.');
end
if tsmc.NrOfSourceVertices ~= size(fldc.VertexCoordinate, 1) || ...
   tsmc.NrOfSourceTriangles ~= size(fldc.TriangleVertex, 1)
    error('neuroelf:xff:badArgument', 'TSM and folded mesh must match.');
end
tsml = tsmc.TriangleEdgeLengths;
tsmc = tsmc.SourceTriangleOfTarget;

% create new SRF from folded mesh
nsph = aft_CopyObject(xo);
nspc = nsph.C;

% neighbors and coords
tri = nspc.TriangleVertex(tsmc, :);
crd = nspc.VertexCoordinate;
crb = crd(tri(:, 1), :);
crd = crb + tsml(:, [1, 1, 1]) .* (crd(tri(:, 2), :) - crb) + ...
    tsml(:, [2, 2, 2]) .* (crd(tri(:, 3), :) - crb);

% for coloring, we must know which is the "main" (weighted) vertex
tsmm = 1 - sum(tsml, 2);
vs = tsmm > (max(tsml, [], 2));
mv = (vs + ((~vs) .* (2 + (tsml(:, 1) < tsml(:, 2))))) - 1;

% find deleted vertices (those will be set to deleted on ALL mapped!)
dl = any(reshape(nspc.VertexColor(tri(:), 1), size(tri)) == (2 ^ 32 - 1000), 2);

% which option
switch (opts.colors)

    % main colors
    case 'main'

        % get the main vertex number
        mvt = tri((1:size(tri, 1)) + (size(tri, 1) * mv'));
        col = nspc.VertexColor(mvt, :);

    % mixing
    case 'mix'

        % get color information of vertices
        col1 = nspc.VertexColor(tri(:, 1), :);
        col2 = nspc.VertexColor(tri(:, 1), :);
        col3 = nspc.VertexColor(tri(:, 1), :);
        c1i = ~isnan(col1(:, 1));
        c2i = ~isnan(col2(:, 1));
        c3i = ~isnan(col3(:, 1));
        cai = c1i & c2i & c3i;
        car = ~cai;

        % if all are indexed colors, use main weight after all
        col = zeros(numel(cai), 4);
        col(cai(mv == 0), :) = col1(cai(mv == 0), :);
        col(cai(mv == 1), :) = col1(cai(mv == 1), :);
        col(cai(mv == 2), :) = col1(cai(mv == 2), :);

        % else, nullify indexed colors where RGB is present as well
        tsmm(c1i) = 0;
        tsml(c2i, 1) = 0;
        tsml(c3i, 2) = 0;
        tsms = sum([tsmm, tsml], 2);
        tsmm(car) = tsmm(car) ./ tsms(car);
        tsml(car, :) = tsml(car, :) ./ tsms(car, [1, 1]);

        % mix remaining items
        col(car, :) = round(tsmm(car, [1, 1, 1, 1]) .* col1(car, :) + ...
            tsml(car, [1, 1, 1, 1]) .* col2(car, :) + tsml(car, [2, 2, 2, 2]) .* col3(car, :));

end

% set deleted to deleted!
col(dl, :) = repmat([2^32 - 1000, 0, 0, 0], [sum(dl), 1]);

% perform smoothing
nei = sphc.Neighbors;
tri = sphc.TriangleVertex;
if opts.smiter > 0
    sopts = struct('type', 'smooth', 'force', opts.smforce, 'niter', opts.smiter, 'areac', 0);
    crd = ne_methods.mesh_morph(crd, nei, tri, sopts);
end

% make settings
nspc.ExtendedNeighbors = sphc.ExtendedNeighbors;
nspc.NrOfVertices = sphc.NrOfVertices;
nspc.NrOfTriangles = size(tri, 1);
nspc.VertexCoordinate = crd;
nspc.VertexColor = col;
nspc.Neighbors = nei;
nspc.TriangleVertex = tri;
nspc.NrOfTriangleStrips = sphc.NrOfTriangleStrips;
nspc.TriangleStripSequence = sphc.TriangleStripSequence;
nspc.AutoLinkedSRF = '<none>';

% set back in global memory array
nsph.C = nspc;

% recalc normals
srf_RecalcNormals(nsph);
