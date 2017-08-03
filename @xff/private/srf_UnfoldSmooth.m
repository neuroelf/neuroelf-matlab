function xo = srf_UnfoldSmooth(xo)
% SRF::UnfoldSmooth  - smooth until all normals point inside
%
% FORMAT:       [srf] = srf.UnfoldSmooth;
%
% No input fileds.
%
% Output fields:
%
%       srf         altered object
%
% Using: mesh_morph, spherecoords.

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:38 PM EST
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
mesh_morph   = ne_methods.mesh_morph;
spherecoords = ne_methods.spherecoords;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end

% get content
bc = xo.C;
crd = bc.VertexCoordinate;
ncd = size(crd);
crd = crd - repmat(bc.MeshCenter, [ncd(1), 1]);
nei = bc.Neighbors;
tri = bc.TriangleVertex;
ntr = size(tri, 1);

% cro = crd; (todo: try to only alter coordinates of mangles triangles and their
% neighbors)

% set morphing opts
opts = struct('areac', 1, 'distc', 3 - eps, 'force', 1e-10, 'niter', 1, 'type', 'smooth');

% continue in steps of 50 with force 0.07 until at most one bad triangle
steps = 0;
nbts = zeros(1, 101);
while true || steps < 5

    % recompute normals (see srf_RecalcNormals.m)
    csa = crd(tri(:, 2), :);
    csb = crd(tri(:, 3), :) - csa;
    csa = csa - crd(tri(:, 1), :);
    crs = cross(csa, csb, 2);
    crl = sqrt(sum(crs .* crs, 2));
    crs = crs ./ crl(:, [1, 1, 1]);
    crs(isnan(crs) | isinf(crs)) = 0;
    nrm = zeros(ncd);
    for tc = 1:ntr
        for vc = 1:3
            nrm(tri(tc, vc), :) = nrm(tri(tc, vc), :) + crs(tc, :);
        end
    end
    crl = sqrt(sum(nrm .* nrm, 2));
    nrm = nrm ./ crl(:, [1, 1, 1]);

    % compute number of badly oriented triangles
    csc = spherecoords(crd);
    nsc = spherecoords(nrm);
    bti = (abs(csc(:, 3) - nsc(:, 3)) < 1.5);
    nbt = sum(bti);
    nbts(steps+1) = nbt;

    % if at most one bad triangle, break
    if nbt < 2 && steps > 5
        break;
    end

    % if already 100 steps, break anyway
    if steps >= 100
        warning('neuroelf:xff:internal', ...
            'More than one inversely oriented triangle after 500 steps.');
        break;
    end

    % perform morph and increase counter
    steps = steps + 1;
    crd = mesh_morph(crd, nei, tri, opts);
end

% if steps > 0, put back
if steps > 0
    nbts = nbts(1:steps+1);
    crd = crd + repmat(bc.MeshCenter, [ncd(1), 1]);
    [fname{1:3}] = fileparts(xo.F);
    fprintf('# of bad triangles in %s:\n', [fname{2} fname{3}]);
    disp([(0:steps)', nbts(:)]);
    bc.VertexCoordinate = crd;
    bc.VertexNormal = nrm;
    xo.C = bc;
end
