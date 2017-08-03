function xo = srf_SetSmoothingInfo(xo)
% SRF::SetSmoothingInfo  - set/update smoothing/smoothness information
%
% FORMAT:       [srf = ] srf.SetSmoothingInfo([reset]);
%
% Input fields:
%
%       reset       if set to 1x1 true, overwrite existing info
%
% No output fields.

% Version:  v1.1
% Build:    16031021
% Date:     Mar-10 2016, 9:40 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2013, 2014, 2016, Jochen Weber
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

% neuroelf methods
global ne_methods;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid object for SRF::NeighborsNDegree');
end
bc = xo.C;
if nargin < 2 || ~islogical(reset) || numel(reset) ~= 1 || ~reset
    if isfield(xo.H, 'SmoothingInfo') && isstruct(xo.H.SmoothingInfo) && numel(xo.H.SmoothingInfo) == 1 && ...
        isfield(xo.H.SmoothingInfo, 'Matrix') && issparse(xo.H.SmoothingInfo.Matrix) && ...
        isequal(size(xo.H.SmoothingInfo.Matrix), [bc.NrOfVertices, bc.NrOfVertices]) && ...
        isfield(xo.H.SmoothingInfo, 'NeighborsList') && size(xo.H.SmoothingInfo.NeighborsList, 1) == bc.NrOfVertices
        return;
    end
end

% create neighbors list and matrix
oneis = bc.Neighbors;
nvert = size(oneis, 1);
numnei = reshape(cellfun('prodofsize', oneis(:, 2)), nvert, 1);
totnei = sum(numnei);
maxnei = max(numnei);
si = zeros(totnei, 1);
si([1; cumsum(numnei(1:end-1, 1)) + 1]) = 1;
si = cumsum(si);
sj = cat(2, oneis{:, 2});
fillnei = (numnei < maxnei);
if sum(fillnei) < (0.01 * nvert)
    fillnei = find(fillnei);
    for fc = fillnei(:)'
        oneis{fc, 2}(end+1:maxnei) = fc;
    end
end
nlist = cat(1, oneis{:, 2});
slist = size(nlist);

% compute distances per vertex->neighbors
crd = xo.C.VertexCoordinate;
vdist = sqrt(sum((repmat(reshape(crd, [nvert, 1, 3]), 1, maxnei) - reshape(crd(nlist, :), [slist, 3])) .^ 2, 3));

% instantiate matrix
neimtx = sparse(si, sj, ones(totnei, 1), nvert, nvert, totnei);

% set into handles
xo.H.SmoothingInfo = struct('Matrix', neimtx, 'NeighborsDist', vdist, 'NeighborsList', nlist);
