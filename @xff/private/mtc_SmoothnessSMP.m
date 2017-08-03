function smp = mtc_SmoothnessSMP(xo, srf)
% MTC::SmoothnessSMP  - compute inherent smoothness in (residual) data
%
% FORMAT:       smp = mtc.SmoothnessSMP(srf)
%
% Input fields:
%
%       srf         object containing the vertices information
%
% Output fields:
%
%       smp         SMP with estimated smoothness map
%
% Using: cov_nd, mesh_neighborsarray, minmaxmean.

% Version:  v1.1
% Build:    16020917
% Date:     Feb-09 2016, 5:06 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2012, 2014, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'mtc') || ...
    numel(srf) ~= 1 || ~xffisobject(srf, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get object data
mtcc = xo.C;
mtcd = mtcc.MTCData;
srfc = srf.C;

% get neighbors array from surface
na = ne_methods.mesh_neighborsarray(srfc.Neighbors, 1);

% get numbers
nv = size(na, 1);
nn = size(na, 2);
nt = size(mtcd, 1);

% create array for correlation data
cr = zeros(nv, nn);

% process in 100MB steps
vs = ceil(2.5e7 / (nn * nt));
ntc = zeros(vs, nn, nt);
cov_nd = ne_methods.cov_nd;
for vc = 1:vs:nv

    % get end index and number of indices
    ve = min(vc + vs - 1, nv);
    nve = 1 + ve - vc;
    if nve < vs
        ntc = zeros(nve, nn, nt);
    end

    % fill time course array with neighbors' time courses
    for nc = 1:nn
        ntc(:, nc, :) = reshape(mtcd(:, na(vc:ve, nc))', [nve, 1, nt]);
    end

    % compute correlation of each vertex with all (direct) neighbors
    [cv, crp] = cov_nd(repmat(double(reshape(mtcd(:, vc:ve)', [nve, 1, nt])), 1, nn), ntc);
    cr(vc:ve, :) = crp;
end

% remove illegal values
cr(isnan(cr) | isinf(cr)) = 0;

% create array for coordinates of neighbors (to compute distances)
crd = srfc.VertexCoordinate;
ccrd = zeros(nv, 3, nn);

% get neighbors' coordinates
for nc = 1:nn
    ccrd(:, :, nc) = crd(na(:, nc), :);
end

% compute distances
dcrd = squeeze(sqrt(sum((repmat(crd, [1, 1, nn]) - ccrd) .^ 2, 2)));

% compute smoothness in mm, compare resestsmooth code and link
% http://www.fmrib.ox.ac.uk/analysis/techrep/tr00df1/tr00df1/node6.html
sef = sqrt(8 * log(2));
fwhm = sef .* dcrd .* sqrt(-1 ./ (4 .* log(cr)));
fwhm(isinf(fwhm) | isnan(fwhm)) = 0;

% compute distance-weighted average for each vertex
idist = 1 ./ dcrd;
idist(isinf(idist) | isnan(idist)) = 0;
vfwhm = sum(idist .* fwhm, 2) ./ sum(idist, 2);
mfwhm = ne_methods.minmaxmean(vfwhm(vfwhm > 0), 5);

% create map
smp = xff('new:smp');
smpc = smp.C;
smpc.NameOfOriginalSRF = srf.F;
smpc.Map(:) = [];
smpc.Map(1).Type = 1;
smpc.Map.ClusterSize = 25;
smpc.Map.EnableClusterCheck = 0;
smpc.Map.LowerThreshold = mfwhm(3) - 3 * sqrt(mfwhm(6));
smpc.Map.UpperThreshold = mfwhm(3) + 3 * sqrt(mfwhm(6));
smpc.Map.UseValuesAboveThresh = 1;
smpc.Map.DF1 = nt - 1;
smpc.Map.DF2 = 0;
smpc.Map.ShowPositiveNegativeFlag = 3;
smpc.Map.BonferroniValue = 40962;
smpc.Map.RGBLowerThreshPos = [192 0 0];
smpc.Map.RGBUpperThreshPos = [255 192 0];
smpc.Map.RGBLowerThreshNeg = [0 0 128];
smpc.Map.RGBUpperThreshNeg = [0 128 255];
smpc.Map.UseRGBColor = 1;
smpc.Map.LUTName = '<default>';
smpc.Map.TransColorFactor = 1;
smpc.Map.Name = sprintf('Smoothness of %s (%s)', xo.F, srf.F);
smpc.Map.SMPData = single(vfwhm);
smp.C = smpc;
