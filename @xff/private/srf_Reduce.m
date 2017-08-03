function xo = srf_Reduce(xo, rfac)
% SRF::Reduce  - reduce regular (icosahedron) SRF by factor 4
%
% FORMAT:       [srf = ] srf.Reduce([rfac]);
%
% Input fields:
%
%       rfac        optional face-reduction factor (default 0.25)
%
% Output fields:
%
%       srf         altered object
%
% Using: mesh_trianglestoneighbors.

% Version:  v1.1
% Build:    16021121
% Date:     Feb-11 2016, 9:24 PM EST
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

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 2 || numel(rfac) ~= 1 || ~isa(rfac, 'double') || ...
    isinf(rfac) || isnan(rfac) || rfac <= 0 || (rfac >= 1 && rfac < 20)
    rfac = 0.25;
end
bc = xo.C;
if rfac > bc.NrOfTriangles
    return;
elseif rfac < 1
    rfac = round(bc.NrOfTriangles * rfac);
end

% irregular reduction
pfac = log2((bc.NrOfVertices - 2) / 10) / 2;
ffac = log2(bc.NrOfTriangles / 20) / 2;
if rfac ~= (0.25 * bc.NrOfTriangles) || pfac ~= fix(pfac) || ...
    ffac ~= fix(ffac) || pfac ~= ffac || pfac < 1

    % try to use Matlab's reducepatch function
    try
        nfv = reducepatch(struct('faces', bc.TriangleVertex, 'vertices', bc.VertexCoordinate), rfac);

        % set into
        bc.NrOfVertices = size(nfv.vertices, 1);
        bc.NrOfTriangles = size(nfv.faces, 1);
        bc.VertexCoordinate = nfv.vertices;
        bc.VertexNormal = zeros(bc.NrOfVertices, 3);
        bc.VertexColor = zeros(bc.NrOfVertices, 4);
        bc.Neighbors = ne_methods.mesh_trianglestoneighbors(bc.NrOfVertices, nfv.faces);
        bc.TriangleVertex = nfv.faces;
        bc.NrOfTriangleStrips = 0;
        bc.TriangleStripSequence = zeros(0, 1);
        bc.AutoLinkedSRF = '';
        xo.C = bc;
        srf_RecalcNormals(xo);
        return;
    catch xfferror
        error('neuroelf:xff:reductionFailed', 'Error reducing number of faces: %s.', xfferror.message);
    end
end

% build new size
pnum = (2 ^ (pfac * 2 - 2)) * 10 + 2;
fnum = (2 ^ (pfac * 2 - 2)) * 20;

% set new options
bc.NrOfVertices = pnum;
bc.NrOfTriangles = fnum;
bc.VertexCoordinate(pnum+1:end, :) = [];
bc.VertexNormal(pnum+1:end, :) = [];
bc.VertexColor(pnum+1:end, :) = [];
bc.AutoLinkedSRF = '';

% building new connections
cnx = cell(pnum, 2);
ocn = bc.Neighbors(:, 2);
rfc = 0;
for c = 1:pnum
    oc = ocn{c};
    nc = [];
    ii = [];
    for vc = oc(:)'
        iv = intersect(oc, ocn{vc});
        for ic = iv(:)'
            ii = union(ii, intersect(ocn{vc}, ocn{ic}));
        end
        ii = setdiff(ocn{vc}, union(ii, iv));
        nc = [nc, ii(:)'];
    end
    if any(nc > pnum)
        error('neuroelf:xff:mathError', 'Invalid point in new connection list entry.');
    end
    cnx{c, 1} = length(nc);
    cnx{c, 2} = nc;
    rfc = rfc + length(nc);
end
bc.Neighbors = cnx;

% building new faces
cnx = cnx(:, 2);
nfs = zeros(rfc, 3);
rfc = 1;
for c = 1:pnum
    fcp = [cnx{c} cnx{c}(1)];
    for cc = 1:length(fcp)-1
        if all(fcp(cc:cc+1) > c)
            nfs(rfc, :) = [c fcp(cc:cc+1)];
        elseif fcp(cc) < fcp(cc+1)
            nfs(rfc, :) = [fcp(cc:cc+1) c];
        else
            nfs(rfc, :) = [fcp(cc+1) c fcp(cc)];
        end
        rfc = rfc + 1;
    end
end
bc.TriangleVertex = unique(nfs, 'rows');
if size(bc.TriangleVertex, 1) ~= fnum
    error('neuroelf:xff:mathError', 'Failure reconnecting vertices to faces.');
end

% set back
xo.C = bc;

% recalc normals
srf_RecalcNormals(xo);
