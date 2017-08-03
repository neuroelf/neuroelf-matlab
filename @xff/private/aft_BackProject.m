function vox = aft_BackProject(xo, srf, opts)
% AFT::BackProject  - back-project MTC-based data to VTC-space
%
% FORMAT:       vox = obj.BackProject(srf, [opts])
%
% Input fields:
%
%       srf         1x1 xff SRF object (linked coordinates will be used)
%       opts        optional settings
%        .bbox      BV compatible bounding box ([44, 38, 44; 241, 193, 211])
%        .mapno     for SMP, map(s) to back-project (default: all)
%        .res       VTC-space resolution, either 1, 2, or {3}
%
% Output fields:
%
%       vox         VTC-space object with backprojected data
%
% TYPES: FSMF, GLM, MTC, SMP
%
% Note: for BV-objects, normals point inward!

% Version:  v1.1
% Build:    16031611
% Date:     Mar-16 2016, 11:42 AM EST
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

% neuroelf functions
global ne_methods;

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'fsmf', 'glm', 'mtc', 'smp'}) || ...
   (xffisobject(xo, true, 'glm') && xo.C.ProjectType ~= 2) || ...
    numel(srf) ~= 1 || ~xffisobject(srf, true, {'fsbf', 'srf'}) || xo.C.NrOfVertices ~= srf.C.NrOfVertices
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
ot = lower(xo.S.Extensions{1});
bc = xo.C;
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'bbox') || ~isa(opts.bbox, 'double') || ~isequal(size(opts.bbox), [2, 3]) || ...
    any(isinf(opts.bbox(:)) | isnan(opts.bbox(:)) | opts.bbox(:) < 0 | opts.bbox(:) > 256) || ...
    any(opts.bbox(2, :) <= opts.bbox(1, :))
    opts.bbox = [44, 38, 44; 241, 193, 211];
else
    opts.bbox = round(opts.bbox);
end
if ~isfield(opts, 'mapno') || ~isa(opts.mapno, 'double') || isempty(opts.mapno) || ...
    any(isinf(opts.mapno(:)) | isnan(opts.mapno(:)) | opts.mapno(:) < 1)
    if ot(1) == 's'
        opts.mapno = 1:numel(bc.Map);
    else
        opts.mapno = [];
    end
end
if ~isfield(opts, 'res') || ~isa(opts.res, 'double') || numel(opts.res) ~= 1 || ...
    isinf(opts.res) || isnan(opts.res) || ~any(opts.res == [1, 2, 3])
    opts.res = 3;
end
ores = opts.res;
offs = opts.bbox(1, :);
opts.bbox(2, :) = offs + ores .* ceil((opts.bbox(2, :) - offs) ./ ores);
ends = opts.bbox(2, :);
dims = round((ends - offs) ./ ores);

% get SRF content
c = srf.C.VertexCoordinate;
if lower(srf.S.Extensions{1}(1)) == 'f'
    c = 128 - c(:, [2, 3, 1]);
end
nc = size(c, 1);

% create voxel object
if ot(1) == 'g'
    vox = aft_CopyObject(xo);
elseif ot(1) == 'm'
    vox = xff('mew:vtc');
elseif ores == 1
    vox = xff('new:vmp');
    dims = dims + 1;
else
    vox = ne_methods.newnatresvmp(opts.bbox, ores, 1);
end
voxc = vox.C;
voxc.Resolution = ores;
voxc.XStart = offs(1);
voxc.XEnd = ends(1);
voxc.YStart = offs(2);
voxc.YEnd = ends(2);
voxc.ZStart = offs(3);
voxc.ZEnd = ends(3);

% re-compute coordinates
c = (1 / ores) .* (c - ones(nc, 1) * offs);

% dimension checks
tm = ones(nc, 1) * (dims - 2);
dfy = dims(1);
dfz = dfy * dims(2);

% keep the interface responsive
nupi = 1 / 172800;
nup = now + nupi;
tic;

% get lower bound and diff
tcl = floor(c);
gdt = find(all(tcl >= 0 & tcl <= tm, 2));
tc = c(gdt, :) - tcl(gdt, :);
vc = find(gdt);
nnc = numel(vc);

% get base indices
tcxo = 1 + tcl(:, 1) + dfy .* tcl(:, 2) + dfz .* tcl(:, 3);

% value and index arrays
iarr = zeros(8 * nnc, 1);
jarr = zeros(8 * nnc, 1);
varr = zeros(8 * nnc, 1);

% iterate for x, y, z
ti = 1;
for xd = 0:1
    for yd = 0:1
        for zd = 0:1

            % get indices
            tie = ti + nnc - 1;

            % add to sampling matrix
            varr(ti:tie) = abs(xd - tc(:, 1)) .* abs(yd - tc(:, 2)) .* abs(zd - tc(:, 3));
            iarr(ti:tie) = tcxo + (xd + yd * dfy + zd * dfz);
            jarr(ti:tie) = vc;
            ti = ti + nnc;
        end
    end
end

% sparse matrix
pd = prod(dims);
smat = sparse(iarr, jarr, varr, pd, nc, 8 * nnc);

% normalize
ssum = full(sum(smat, 2));
ssum(ssum < 0.25) = Inf;
ssum = 1 ./ ssum;
spv = 1:pd;
smat = sparse(spv(:), spv(:), ssum, pd, pd, pd) * smat;

% sampling
if ot(1) == 'g'
    voxc.ProjectType = 1;
    if bc.ProjectTypeRFX < 1
        
    else
        voxc.GLMData.RFXGlobalMap = reshape(single(full(smat * double(bc.GLMData.RFXGlobalMap))), dims);
        nb = size(bc.GLMData.Subject(1).BetaMaps, 2);
        for sc = 1:numel(bc.GLMData.Subject)
            voxc.GLMData.Subject(sc).BetaMaps = reshape(single(full(smat * ...
                double(bc.GLMData.Subject(sc).BetaMaps))), [dims, nb]);
        end
    end
elseif ot(1) == 'm'
else
    
end

% store back
voxc.RunTimeVars.SourceObject = xo.L;
voxc.RunTimeVars.SamplingSRF = srf.L;
voxc.RunTimeVars.SamplingMatrix = smat;
vox.C = voxc;

% and also store in SRF for sampling
srf.H.BBoxSampling = struct('BBox', [offs; ends], 'DimXYZ', dims, 'FCube', [256, 256, 256], ...
    'RadCnv', 1, 'ResXYZ', ores .* [1, 1, 1], 'SamplingMatrix', smat);
