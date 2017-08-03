function smp = srf_DensityMap(xo, opts)
% SRF::DensityMap  - generate a vertex density map
%
% FORMAT:       dsmp = srf.DensityMap([opts])
%
% Input fields:
%
%       opts        optional settings
%        .area      generate area map (default: true)
%        .dist      generate average distance map (default: true)
%        .mindist   generate minimum distance map (default: true)
%        .stddist   generate distance inhomogeneity map (default: false)
%
% Output fields:
%
%       dsmp        density SMP (in average mm distance to neighbors)
%
% Using: mesh_neighborsarray.

% Version:  v1.1
% Build:    16031617
% Date:     Mar-16 2016, 5:06 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'fsbf', 'srf'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || numel(opts) ~= 1 || ~isstruct(opts)
    opts = struct;
end
if ~isfield(opts, 'area') || ~islogical(opts.area) || numel(opts.area) ~= 1
    opts.area = true;
end
if ~isfield(opts, 'dist') || ~islogical(opts.dist) || numel(opts.dist) ~= 1
    opts.dist = true;
end
if ~isfield(opts, 'mindist') || ~islogical(opts.mindist) || numel(opts.mindist) ~= 1
    opts.mindist = true;
end
if ~isfield(opts, 'stddist') || ~islogical(opts.stddist) || numel(opts.stddist) ~= 1
    opts.stddist = false;
end

% get content
bc = xo.C;
c = bc.VertexCoordinate;
nc = size(c, 1);
nv = ne_methods.mesh_neighborsarray(bc.Neighbors, 1);
ov = repmat((1:nc)', 1, size(nv, 2));
nn = sum(nv ~= ov, 2);
sn = size(nv);

% get average distance to neighbors (also used for density!)
vd = sqrt(sum((cat(3, reshape(c(ov, 1), sn), reshape(c(ov, 2), sn), reshape(c(ov, 3), sn)) - ...
    cat(3, reshape(c(nv, 1), sn), reshape(c(nv, 2), sn), reshape(c(nv, 3), sn))) .^ 2, 3));
ad = sum(vd, 2) ./ nn;

% minimum distance
if opts.mindist
    md = min(vd + max(vd(:)) .* (vd == 0), [], 2);
end

% irregularity across distances
if opts.stddist
    sd = zeros(nc, 1);
    for c = 1:nc
        sd(c) = (1 / (nn(c) - 1)) * sqrt(sum((vd(c, :) - ad(c)) .^ 2));
    end
end

% also compute area
if opts.area
    [a, ta, va] = srf_Area(xo);
end

% create smps
smp = xff('new:smp');
smpc = smp.C;
smpc.NrOfVertices = nc;
smpc.NameOfOriginalSRF = xo.F;
vd = 1 ./ (ad + eps);
smpc.Map.Type = 41;
smpc.Map.DF1 = round(3 * (size(bc.TriangleVertex, 1) / nc));
smpc.Map.LowerThreshold = double(min(vd)) - 16 .* sqrt(eps);
smpc.Map.UpperThreshold = double(mean(vd) + 2 * std(vd));
smpc.Map.ShowPositiveNegativeFlag = 1;
smpc.Map.Name = 'Vertex density w.r.t. direct neighbors';
smpc.Map.SMPData = single(vd);
if opts.dist
    smpc.Map(end+1) = smpc.Map(1);
    smpc.Map(end).Type = 42;
    smpc.Map(end).LowerThreshold = double(min(ad)) - 16 .* sqrt(eps);
    smpc.Map(end).UpperThreshold = double(max(ad));
    smpc.Map(end).Name = 'Average distance to direct neighbors';
    smpc.Map(end).SMPData = single(ad);
end
if opts.mindist
    smpc.Map(end+1) = smpc.Map(1);
    smpc.Map(end).Type = 43;
    smpc.Map(end).LowerThreshold = double(min(md)) - 16 .* sqrt(eps);
    smpc.Map(end).UpperThreshold = double(mean(md) + 2 * std(md));
    smpc.Map(end).Name = 'Minimum distance to a direct neighbor';
    smpc.Map(end).SMPData = single(md);
end
if opts.stddist
    smpc.Map(end+1) = smpc.Map(1);
    smpc.Map(end).Type = 44;
    smpc.Map(end).LowerThreshold = double(min(sd)) - 16 .* sqrt(eps);
    smpc.Map(end).UpperThreshold = double(mean(sd) + 2 * std(sd));
    smpc.Map(end).Name = 'Local distortion (std) in distances to direct neighbors';
    smpc.Map(end).SMPData = sd;
end
if opts.area
    smpc.Map(end+1) = smpc.Map(1);
    smpc.Map(end).Type = 45;
    smpc.Map(end).LowerThreshold = double(min(va)) - 16 .* sqrt(eps);
    smpc.Map(end).UpperThreshold = double(max(va));
    smpc.Map(end).Name = 'Vertex-associated area (mm^2)';
    smpc.Map(end).SMPData = single(va);
end
smpc.NrOfMaps = numel(smpc.Map);
smp.C = smpc;
