function smp = srf_CurvatureMap(xo, opts)
% SRF::CurvatureMap  - compute a curvature map for a surface
%
% FORMAT:       smp = srf.CurvatureMap([opts]);
%
% Input fields:
%
%       opts        optional settings
%        .medrem    also create a median removed map (default: false)
%
% Output fields:
%
%       smp         SMP with one map containing the curvature info

% Version:  v1.1
% Build:    16031617
% Date:     Mar-16 2016, 5:03 PM EST
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

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'fsbf', 'srf'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
fn = xo.F;
ft = lower(xo.S.Extensions{1});
if isempty(fn)
    fn = 'unsaved';
else
    [fp, fn] = fileparts(fn);
end
fn = sprintf('%s.%s (%d vertices, %d triangles)', fn, ft, ...
    size(bc.VertexCoordinate, 1), size(bc.TriangleVertex, 1));
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'medrem') || ~islogical(opts.medrem) || numel(opts.medrem) ~= 1
    opts.medrem = false;
end

% prepare arrays for computation
n = bc.NrOfVertices;
t = bc.TriangleVertex;
vo = double(bc.VertexCoordinate(:, :));
vn = vo - 0.01 * double(bc.VertexNormal(:, :));

% compute area per triangle
tvo = cat(3, vo(t(:, 1), :), vo(t(:, 2), :), vo(t(:, 3), :));
tvo = cat(3, tvo(:, :, 2) - tvo(:, :, 1), tvo(:, :, 3) - tvo(:, :, 2), tvo(:, :, 1) - tvo(:, :, 3));
tvo = squeeze(sqrt(sum(tvo .* tvo, 2)));
so = sum(tvo, 2) ./ 2;
ao = sqrt(so .* (so - tvo(:, 1)) .* (so - tvo(:, 2)) .* (so - tvo(:, 3)));

% the same at 0.01 along normal
tvn = cat(3, vn(t(:, 1), :), vn(t(:, 2), :), vn(t(:, 3), :));
tvn = cat(3, tvn(:, :, 2) - tvn(:, :, 1), tvn(:, :, 3) - tvn(:, :, 2), tvn(:, :, 1) - tvn(:, :, 3));
tvn = squeeze(sqrt(sum(tvn .* tvn, 2)));
sn = sum(tvn, 2) ./ 2;
an = sqrt(sn .* (sn - tvn(:, 1)) .* (sn - tvn(:, 2)) .* (sn - tvn(:, 3)));

% build curvature index (per triangle)
crv = 1 - an ./ ao;

% average for each point of triangle
vv = zeros(n, 1);
vd = zeros(n, 1);
for vc = 1:3
    for tc = 1:size(t, 1)
        tr = t(tc, vc);
        vv(tr) = vv(tr) + crv(tc);
        vd(tr) = vd(tr) + 1;
    end
end

% average correctly
vv = vv ./ vd;
m = median(abs(vv));

% create smp
smp = xff('new:smp');
smc = smp.C;
smc.NrOfVertices = n;
smc.NameOfOriginalSRF = xo.F;
smc.Map(1).Type = 31;
smc.Map = smc.Map(1);
smc.Map.DF1 = 0;
smc.Map.RGBLowerThreshPos = [0, 85, 170];
smc.Map.RGBUpperThreshPos = [64, 128, 255];
smc.Map.RGBLowerThreshNeg = [0, 170, 85];
smc.Map.RGBUpperThreshNeg = [64, 255, 128];
smc.Map.LowerThreshold = 50 * m;
smc.Map.UpperThreshold = 300 * m;
smc.Map.ShowPositiveNegativeFlag = 3;
smc.Map.Name = ['Curvature: ' fn];
smc.Map.SMPData = 100 * vv;
if opts.medrem
    smc.Map(2) = smc.Map;
    smc.Map(2).Type = 32;
    smc.Map(2).RGBLowerThreshPos = [192, 192, 192];
    smc.Map(2).RGBUpperThreshPos = [128, 128, 128];
    smc.Map(2).RGBLowerThreshNeg = [192, 192, 192];
    smc.Map(2).RGBUpperThreshNeg = [240, 240, 240];
    smc.Map(2).SMPData = smc.Map(1).SMPData - mean(smc.Map(1).SMPData(:));
    smc.Map(2).LowerThreshold = sqrt(eps);
    smc.Map(2).UpperThreshold = max(abs(smc.Map(2).SMPData(:)));
end
smp.C = smc;
