function neuroelf_makeshenfiles(res, showpbar)
%NEUROELF_MAKESHENFILES  Create surfaces and VOI file from 268-parcel Atlas
%   NEUROELF_MAKESHENFILES creates a whole-brain, right-hemisphere, and
%   left-hemisphere surface as well as a VOI file with all parcels.
%
%   NEUROELF_MAKESHENFILES(RES) uses the internal resolution of RES mm; the
%   default is 2mm.

% Version:  v1.1
% Build:    16061014
% Date:     Jun-10 2016, 2:23 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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

% argument check
if nargin < 1 || ~isa(res, 'double') || numel(res) ~= 1 || isinf(res) || isnan(res) || res <= 0
    res = 2;
end
if nargin < 2 || ~islogical(showpbar) || numel(showpar) ~= 1
    showpbar = true;
end

% progress bar
try
    if showpbar
        pbar = xprogress;
        xprogress(pbar, 'settitle', 'Creating Parcel-based surface...');
        xprogress(pbar, 0, 'Loading and sampling Atlas...', 'visible', 0, 1);
    else
        pbar = [];
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    pbar = [];
end

% load atlas
shenpath = [neuroelf_path('files') '/shenparcel/'];
atlas = xff([shenpath 'shen_1mm_268_parcellation.nii.gz']);

% get label coordinates (in BVI space)
labelcoords = atlas.LabelCoords(struct('space', 'bvi'));

% that's it for the atlas
atlas.ClearObject;

% create VMR
vmr = xff('new:vmr');

% generate empty surface
bh = neuroelf('emptysrf');
bh.ConvexRGBA = [0.7, 0.7, 0.7, 0.5];
bh.ConcaveRGBA = [0.3, 0.3, 0.3, 0.5];
nl = numel(labelcoords);
nv = zeros(nl, 1);
nt = nv;
tv = 1;
tt = 1;

% iterate over labels
for lc = 1:nl

    % progress
    if ~isempty(pbar)
        pbar.Progress(lc/(nl+2), sprintf('Sampling label %d...', lc));
    end

    % set into VMR
    vmr.VMRData(:) = 0;
    c = unique(round((labelcoords{lc} ./ res) - 1), 'rows');
    vmr.VMRData(1 + c(:, 1) + 256 .* c(:, 2) + 65536 .* c(:, 3)) = 240;

    % reconstruct as surface
    srf = vmr.DBReco(struct('warn', false));
    snv = srf.NrOfVertices;
    snt = srf.NrOfTriangles;

    % surface smooth
    srf.Smooth(100, 0.07, struct('show', false));

    % add to combined mesh
    nv(lc) = srf.NrOfVertices;
    nt(lc) = srf.NrOfTriangles;
    bh.VertexCoordinate(tv:tv+snv-1, :) = srf.VertexCoordinate .* res;
    bh.VertexNormal(tv:tv+snv-1, :) = srf.VertexNormal;
    bh.VertexColor(tv:tv+snv-1, :) = srf.VertexColor;
    nei = srf.Neighbors;
    tnei = tv - 1;
    for nc = 1:size(nei, 1)
        nei{nc, 2} = nei{nc, 2} + tnei;
    end
    bh.Neighbors(tv:tv+snv-1, :) = nei;
    bh.TriangleVertex(tt:tt+snt-1, :) = srf.TriangleVertex + tnei;

    % increase counters
    tv = tv + snv;
    tt = tt + snt;

    % clear SRF
    srf.ClearObject;
end

% clear VMR
vmr.ClearObject;

% set in combined mesh
bh.NrOfVertices = tv - 1;
bh.NrOfTriangles = tt - 1;

% generate VOI
if ~isempty(pbar)
    pbar.Progress(1 - 1/nl, 'Creating VOI...');
end

% generate VOI
voi = xff('new:voi');

% space out VOI
voi.VOI(268).Voxels = zeros(0, 3);
ccrd = zeros(268, 3);
pcrd = zeros(268, 3);
tcrd = zeros(268, 3);

% iterate over labels
for lc = 1:numel(labelcoords)

    % get coordinates
    crd = labelcoords{lc};

    % compute mean (for color)
    rgb = round(mean(crd));
    pcrd(lc, :) = rgb - 1;

    % convert to TAL
    crd = unique(round(129 - crd(:, [3, 1, 2])), 'rows');
    ncrd = size(crd, 1);

    % compute mean for sorting
    mcrd = mean(crd);
    ccrd(lc, :) = round(mcrd);
    [scrd, si] = sort(sum((crd - repmat(mcrd, ncrd, 1)) .^ 2, 2));
    crd = crd(si, :);
    tcrd(lc, :) = crd(1, :);

    % look up coordinate closest to COG
    tdl = regexprep(regexprep( ...
        tdclient(round(icbm2tal(crd(1, :))), 'NGM', [13, 1]), '^.*\:\s+', ''), ...
        '\s+.coords.*$', '');
    tdl(tdl == char(10) | tdl == char(13) | tdl == char(9)) = [];
    tdl = splittocell(tdl, ',');
    if numel(tdl) > 4 && numel(tdl{5}) > 1
        if numel(tdl{3}) > 1
            tdl = sprintf('%s (%s)', tdl{3}, tdl{5});
        else
            tdl = tdl{5};
        end
    elseif numel(tdl) > 2 && numel(tdl{3}) > 1
        tdl = tdl{3};
    elseif numel(tdl) > 1 && numel(tdl{2}) > 1
        tdl = tdl{2};
    else
        tdl = 'Brain';
    end

    % add to VOI
    if lc < 134
        voi.VOI(lc).Name = sprintf('RH label %03d: %s', lc, tdl);
    else
        voi.VOI(lc).Name = sprintf('LH label %03d: %s', lc, tdl);
    end
    voi.VOI(lc).Color = rgb;
    voi.VOI(lc).NrOfVoxels = size(crd, 1);
    voi.VOI(lc).Voxels = crd;
end
vnames = voi.VOINames;
bh.RunTimeVars.AutoSave = true;
bh.RunTimeVars.Labels = (1:268)';
bh.RunTimeVars.LabelCoordinates = pcrd;
bh.RunTimeVars.LabelCoordNames = vnames;
bh.RunTimeVars.LabelMeshVertices = nv(:);
bh.RunTimeVars.LabelMeshTriangles = nt(:);
bh.RunTimeVars.LabelTalCoordinates = tcrd;

% create two hemispheres separately
if ~isempty(pbar)
    pbar.Progress(1, 'Creating hemisphere surfaces...');
end
nrhv = sum(nv(1:133));
nrht = sum(nt(1:133));
rh = bh.CopyObject;
rh.NrOfVertices = nrhv;
rh.NrOfTriangles = nrht;
rh.VertexCoordinate(nrhv+1:end, :) = [];
rh.VertexNormal(nrhv+1:end, :) = [];
rh.VertexColor(nrhv+1:end, :) = [];
rh.Neighbors(nrhv+1:end, :) = [];
rh.TriangleVertex(nrht+1:end, :) = [];
rh.RunTimeVars.Labels = (1:133)';
rh.RunTimeVars.LabelCoordinates = pcrd(1:133, :);
rh.RunTimeVars.LabelCoordNames = vnames(1:133);
rh.RunTimeVars.LabelMeshVertices = nv(1:133);
rh.RunTimeVars.LabelMeshTriangles = nt(1:133);
rh.RunTimeVars.LabelTalCoordinates = tcrd(1:133, :);
nlhv = sum(nv(134:end));
nlht = sum(nt(134:end));
lh = bh.CopyObject;
lh.NrOfVertices = nlhv;
lh.NrOfTriangles = nlht;
lh.VertexCoordinate(1:nrhv, :) = [];
lh.VertexNormal(1:nrhv, :) = [];
lh.VertexColor(1:nrhv, :) = [];
lh.Neighbors(1:nrhv, :) = [];
lh.TriangleVertex(1:nrht, :) = [];
lh.TriangleVertex = lh.TriangleVertex - nrhv;
n = lh.Neighbors;
for nc = 1:size(n, 1)
    n{nc, 2} = n{nc, 2} - nrhv;
end
lh.Neighbors = n;
lh.RunTimeVars.Labels = (134:268)';
lh.RunTimeVars.LabelCoordinates = pcrd(134:end, :);
lh.RunTimeVars.LabelCoordNames = vnames(134:end);
lh.RunTimeVars.LabelMeshVertices = nv(134:end);
lh.RunTimeVars.LabelMeshTriangles = nt(134:end);
lh.RunTimeVars.LabelTalCoordinates = tcrd(134:end, :);

% close bar
if ~isempty(pbar)
    closebar(pbar);
end

% save files
disp('Saving Shen et al. 268-parcel surfaces files...');
bh.SaveAs([shenpath 'shen_parcels_BH.srf']);
lh.SaveAs([shenpath 'shen_parcels_LH.srf']);
rh.SaveAs([shenpath 'shen_parcels_RH.srf']);
disp('Saving Shen et al. 268-parcel VOI...');
voi.SaveAs([shenpath 'shen_parcels.voi']);

% clearing objects
voi.ClearObject;
lh.ClearObject;
rh.ClearObject;
bh.ClearObject;
