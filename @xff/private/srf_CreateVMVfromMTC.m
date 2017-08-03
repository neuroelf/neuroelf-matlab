function vmv = srf_CreateVMVfromMTC(xo, mtc, opts)
% SRF::CreateVMVfromMTC  - create VMV from MTC
%
% FORMAT:       vmv = srf.CreateVMVfromMTC(mtc, opts)
%
% Input fields:
%
%       mtc         MTC object with matching number of vertices
%       opts        1x1 struct with optional fields
%        .sdm       SDM object, regresses out variance
%        .sdmcon    contrast after regression (variance re-inclusion)
%        .smorph    number of morphing steps, default: 8
%        .threshmax maximum threshold, default: 3sd
%        .threshmin minimum threshold, default: 1sd
%
% Output fields:
%
%       vmv         VMV containing the SRF positions and color coding
%
% Using: ztrans.

% Version:  v1.1
% Build:    16021119
% Date:     Feb-11 2016, 7:55 PM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf') || ...
    numel(mtc) ~= 1 || ~xffisobject(mtc, true, 'mtc')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end

% get contents
bc = xo.C;
mtcc = mtc.C;
nrv = bc.NrOfVertices;

% options check
sdmm = [];
if isfield(opts, 'sdm') && numel(opts.sdm) == 1 && xffisobject(opts.sdm, true, 'sdm')
    sdmc = opts.sdm.C;
    if size(sdmc.SDMMatrix, 1) == size(mtcc.MTCData, 1)
        sdmm = sdmc.SDMMatrix;
        sdmcon = [];
        if isfield(opts, 'sdmcon') && isa(opts.sdmcon, 'double') && numel(opts.sdmcon) <= size(sdmm, 2)
            sdmcon = find(opts.sdmcon(:)' ~= 0);
        end
    end
end
if ~isfield(opts, 'smorph') || ~isa(opts.smorph, 'double') || numel(opts.smorph) ~= 1 || ...
    isinf(opts.smorph) || isnan(opts.smorph) || opts.smorph < 1 || opts.smorph ~= fix(opts.smorph)
    opts.smorph = 8;
end
if ~isfield(opts, 'threshmax') || ~isa(opts.threshmax, 'double') || numel(opts.threshmax) ~= 1 || ...
    isinf(opts.threshmax) || isnan(opts.threshmax) || opts.threshmax < 0 || opts.threshmax > 12
    opts.threshmax = 3;
end
if ~isfield(opts, 'threshmin') || ~isa(opts.threshmin, 'double') || numel(opts.threshmin) ~= 1 || ...
    isinf(opts.threshmin) || isnan(opts.threshmin) || opts.threshmin < 0 || opts.threshmin > 12
    opts.threshmin = 1;
end

% checks
if nrv ~= mtcc.NrOfVertices
    error('neuroelf:xff:objectsMismatch', 'NrOfVertices must match between SRF and MTC.');
end

% generate VMV object
vmv = xff('new:vmv');
vmvc = vmv.C;

% fill some fields
nrm = size(mtcc.MTCData, 1);
vmvc.NrOfPositions = nrm + 1;
vmvc.NrOfVertices = nrv;
vmvc.UseViewPoint = 0;
vmvc.UseVertexColor = 1;
vmvc.NameOfOriginateSRF = xo.F;
vmvc.VertexPosition.Coordinates = bc.VertexCoordinate;
vmvc.VertexPosition.Normals = bc.VertexNormal;
vmvc.VertexPosition.Colors = bc.VertexColor;
vmvc.VertexPosition.MorphingSteps = opts.smorph;
vmvc.VertexPosition.NameOfState = 'Rest state';
vmvc.VertexPosition(2:nrm + 1) = vmvc.VertexPosition(1);

% if SDM is given, perform calculus first
if ~isempty(sdmm)

    % get betas and predicted time course
    [rb, ri, rp] = sdm_CalcBetas(opts.sdm, mtcc.MTCData, 1);
    rm = rb(:, end);

    % mask according to mean of confound
    msk = rb(:, end) < (0.5 * mean(rb(:, end)));

    % get res and RSS
    res = mtcc.MTCData - rp;
    flt = sqrt(var(rp) ./ var(double(mtcc.MTCData)));
    rss = sqrt(sum(res .* res) ./ ((size(sdmm, 1) - 1) ./ rm')) ./ flt;

    % if contrast given, only take res + signal of interest
    if ~isempty(sdmcon)
        rb = rb(:, sdmcon);
        sdmm = sdmm(:, sdmcon);
        zt = 100 .* ((res + sdmm * rb') ./ (ones(size(sdmm, 1), 1) * rss));
    else
        zt = 100 .* (double(mtcc.MTCData) - ones(size(sdmm, 1), 1) * rm') ./ (ones(size(sdmm, 1), 1) * rss);
    end
    zt(:, msk) = 0;

% otherwise, simply build z-transformed time course
else
    zt = ne_methods.ztrans(double(mtcc.MTCData));
end
zt(zt < 0 | isnan(zt)) = 0;

% loop over maps
for mc = 1:nrm
    tmap = zt(mc,:);
    tmap = max(0, tmap - opts.threshmin);
    tmap = min(1, tmap ./ (opts.threshmax - opts.threshmin + eps));
    tmi = find(tmap > 0);
    tmap = tmap(tmi);
    tmc = mc + 1;
    vmvc.VertexPosition(tmc).Colors(tmi, 1) = NaN;
    vmvc.VertexPosition(tmc).Colors(tmi, 2) = 255;
    vmvc.VertexPosition(tmc).Colors(tmi, 3) = round(255 * tmap);
    vmvc.VertexPosition(tmc).Colors(tmi, 4) = 0;
    vmvc.VertexPosition(tmc).NameOfState = sprintf('MTC tp %d', mc);
end

% store information back
vmv.C = vmvc;
