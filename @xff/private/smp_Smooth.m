function xo = smp_Smooth(xo, srf, niter, mapno, k)
% SMP::Smooth  - smooth SMP map with neighbors of an SRF
%
% FORMAT:       [smp = ] smp.Smooth(srf [, niter [, mapno [, k]]]);
%
% Input fields:
%
%       srf         SRF object with required neighbor information
%       niter       number of iterations (default: 1)
%       mapno       which map number (default: 1)
%       k           smoothing kernel in mm (only for direct neighbors, 2)
%
% Output fields:
%
%       smp         object with one added, smoothed map
%
% Using: flexinterpn, lsqueeze, mesh_neighborsarray, smoothkern, varc.

% Version:  v1.1
% Build:    16031617
% Date:     Mar-16 2016, 5:04 PM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'smp') || ...
    numel(srf) ~= 1 || ~xffisobject(srf, true, {'fsbf', 'srf'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
sbc = srf.C;
if bc.NrOfVertices ~= sbc.NrOfVertices
    error('neuroelf:xff:badArgument', ...
        'SMP and SRF objects must match in NrOfVertices property.');
end
if isempty(bc.Map)
    error('neuroelf:xff:badArgument', 'SMP does not contain any maps.');
end
if nargin < 3 || ~isa(niter, 'double') || numel(niter) ~= 1 || ...
    isinf(niter) || isnan(niter) || niter < 1 || niter > 10000
    niter = 1;
else
    niter = round(niter);
end
if nargin < 4 || ~isa(mapno, 'double') || isempty(mapno) || ...
    any(isinf(mapno(:)) | isnan(mapno(:)) | mapno(:) < 1 | mapno(:) > numel(bc.Map))
    mapno = 1;
else
    mapno = unique(round(mapno(:)'));
end
if nargin < 5 || ~isa(k, 'double') || numel(k) ~= 1 || isinf(k) || isnan(k) || k <= 0
    k = 2;
else
    k = min(k, 12);
end

% extend map
newmap = numel(bc.Map) + 1;
bc.Map(newmap:end+numel(mapno)) = bc.Map(mapno);

% compute kernel (in um)
kr = 1000 .* ne_methods.smoothkern(1000 * k);

% get central kernel position
kh = 0.5 * (numel(kr) + 1);

% compute 2D weight for own value (square of 1D!)
ks = kr(kh) * kr(kh);

% and compound weight for all other neighbors
ko = (1 - ks);

% get neighbors
nei = ne_methods.mesh_neighborsarray(sbc.Neighbors);
gnei = (nei ~= 0);

% compute neighbor distance (inter-neighbor weighting)
dnei = (1:size(nei, 1))' * ones(1, size(nei, 2));
nei = nei(gnei);
dnei(gnei) = sqrt(sum(...
    (sbc.VertexCoordinate(dnei(gnei), :) - sbc.VertexCoordinate(nei, :)) .^ 2, 2));

% compute weighting (relative to other neighbors)
wnei = zeros(size(dnei));
wnei(gnei) = ne_methods.flexinterpn(kr, kh + 1000 .* ne_methods.lsqueeze(dnei(gnei)));
wnei = wnei .* ((ko ./ sum(gnei .* wnei, 2)) * ones(1, size(dnei, 2)));

% iterate over maps
for mc = mapno

    % get map data
    map = double(bc.Map(mc).SMPData);

    % get overall sum (for final scaling, must be invariant!)
    maps = sum(map);
    mapv = sqrt(ne_methods.varc(map));

    % iterations
    for sc = 1:niter

        % sample data
        dnei(gnei) = map(nei);

        % update map
        map = ks .* map + sum(wnei .* dnei, 2);
    end

    % re-scale (if mean is clear enough away from 0!)
    if (maps / mapv) >= 1
        map = (maps / sum(map)) .* map;
    end

    % then re-store in target
    bc.Map(newmap).SMPData = single(map);
    bc.Map(newmap).Name = [bc.Map(newmap).Name sprintf(' sm%d_%.2fmm', niter, k)];
    newmap = newmap + 1;
end

% update number of maps
bc.NrOfMaps = numel(bc.Map);

% set back in global storage array
xo.C = bc;
