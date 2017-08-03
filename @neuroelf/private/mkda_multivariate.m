% mkda_multivariate  - coactivation analysis for MKDA/META maps

% code to extract data from a PLP->VMP map series and perform a
% multivariate dimension reduction, followed by parcellation and
% network detection

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:28 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2013, 2016, Jochen Weber
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

% settings
parcelk = 75;

% set transio access for VMP
dispp('Initializing xff for VMP access...');
xffroot = xff();
vmptio = xffroot.TransIOSize('vmp', 1e5);

% load VMP containing the maps with the series of data
dispp('Loading indicator/component maps (transio access)...');
vmp = xff('*.vmp', 'Please select the VMP containing the indicator/kernel maps...');
xffroot.TransIOSize('vmp', vmptio);
nmaps = numel(vmp.Map) - 3;
dispp(sprintf(' -> %d component maps (contrasts) in analysis.', nmaps));

% load VMP containing the thresholded map
dispp('Loading META/MKDA thresholded map...');
tvmp = xff('*.vmp', 'Please select the VMP containing the thresholded map...');
mapsize = size(tvmp.Map(1).VMPData);
mapvnum = prod(mapsize);
tvmpd = [size(tvmp.Map(1).VMPData), 16];
tvmpf = tvmp.FilenameOnDisk;
tvmpm = tvmp.BoundingBox.QuatB2T;
tvmpr = repmat(tvmp.Resolution, 1, 3);

% update parcel-k threshold
parcelk = round(1 + parcelk / (tvmp.Resolution ^ 3));

% get thresholded map
dispp('Clustering supra-threshold voxels...');
tvmp.ClusterTable(1);
tdata = tvmp.Map(1).VMPData .* double(tvmp.Map(1).VMPDataCT);
[ttext, tlist, tvol] = ...
    clustervol(tdata, tvmp.Map(1).LowerThreshold, tvmp.Map(1).ClusterSize);
tmask = find(tvol(:) > 0);
nvox = numel(tmask);
dispp(sprintf(' -> %d supra-threshold voxels (in %d clusters) found.', ...
    nvox, numel(tlist)));

% get thresholded data
dispp('Extracting component map data (transio access)...');
tmaps = zeros(nmaps, nvox);
for mc = 1:nmaps
    tmaps(mc, :) = vmp.Map(mc).VMPData(tmask);
    if mod(mc, 25) == 0
        dispp(sprintf(' -> %d maps done...', mc));
    end
end
dispp(sprintf(' -> %d maps done...', nmaps));

% remove vmp from memory
vmp.ClearObject;

% for now, follow original procedure (which may or may not need clean-up)
dispp('Dimension reduction (SVD)...');
[tcomp, tcS, tcV] = svd(tmaps, 'econ');

% detect number of components (half of the total sum of eigenvalues-main)
dtcS = diag(tcS);
ncomp = findfirst(cumsum(dtcS(2:end))  >(sum(dtcS(2:end)) / exp(1))) + 1;
dispp(sprintf(' -> %d components selected...', ncomp));

% truncate data
scomp = ncomp + 1;
tcomp = tcomp(:, 1:ncomp);
tcS = tcS(1:ncomp, 1:ncomp);
tcV = tcV(:, 1:ncomp)';

% copy VMP and get ncomp + 1 maps worth to show data
cvmp = tvmp.CopyObject;
cvmp.Map = cvmp.Map(1);
cvmp.Map.Type = 2;
cvmp.Map.LowerThreshold = 0.2;
cvmp.Map.UpperThreshold = 0.5;
cvmp.Map.VMPData(:) = 0;
cvmp.Map.EnableClusterCheck = 0;
cvmp.Map.ClusterSize = 10;
cvmp.Map = repmat(cvmp.Map, 1, scomp);
cvmp.Map(scomp).Type = 15;
cvmp.Map(scomp).LowerThreshold = 1;
cvmp.Map(scomp).UpperThreshold = scomp;

% compute rank-transformed correlation (the maps don't need
% transformation, as they are only binary with two "ranks"!)
dispp('Finding best matching component (class)...');
rcomp = ranktrans(tcomp);

% keep track of best correlation and class for each voxel
bcorr = -ones(nvox, 1);
ccorr = zeros(nvox, 1);

% transpose data (for now)
tmaps = tmaps';

% for each component
for mc = 1:ncomp

    % set map name
    cvmp.Map(mc).Name = sprintf('Correlation tmaps ./. component %03d', mc);

    % compute correlation
    [cv, cr] = cov_nd(tmaps, repmat(rcomp(:, mc)', nvox, 1));
    cr(isinf(cr) | isnan(cr)) = 0;

    % store data
    cvmp.Map(mc).VMPData(tmask) = cr;

    % update class (first class signed!)
    if mc == 1
        ccorr(cr > bcorr) = mc;
        bcorr = max([bcorr, cr], [], 2);
    else
        ccorr(abs(cr) > bcorr) = mc;
        bcorr = max([bcorr, abs(cr)], [], 2);
    end
end

% un-transpose data
tmaps = tmaps';

% store best correlation map
cvmp.Map(scomp).VMPData(tmask) = ccorr;
cvmp.Map(scomp + 1) = cvmp.Map(1);
cvmp.Map(end).VMPData(tmask) = bcorr;
cvmp.Map(end).Name = 'Best-matching correlation tmaps ./. any component';
cvmp.NrOfMaps = numel(cvmp.Map);

% cluster maps to get number of components
dispp('Clustering parcels (removing spurious voxels)...');
clist = cell(ncomp, 1);
for mc = 1:ncomp
    [ctext, clist{mc}] = clustervol( ...
        cvmp.Map(mc).VMPData .* double(cvmp.Map(scomp).VMPData == mc), 0.01, parcelk);
end
clist = cat(1, clist{:});
nparcel = numel(clist);

% convert coordinates
coordlist = zeros(1, mapvnum);
coordlist(tmask) = 1:nvox;

% generate study-by-parcel list
dispp('Generating study-by-parcel list...');
actparcel = zeros(nmaps, nparcel);

% extend cluster structure (for Tor's tools)
for pc = 1:nparcel
    clist(pc).coordi = coordlist(sub2ind1(mapsize, clist(pc).coords));
    clist(pc).title = sprintf('Cluster of %d voxels from %s', clist(pc).size, tvmpf);
    clist(pc).threshold = 0.5;
    clist(pc).M = tvmpm;
    clist(pc).dim = tvmpd;
    clist(pc).voxSize = tvmpr;
    clist(pc).name = '';
    clist(pc).Z = clist(pc).values(:)';
    clist(pc).XYZmm = tvmpm * cat(1, clist(pc).coords', ones(1, clist(pc).size));
    clist(pc).XYZmm(4, :) = [];
    clist(pc).XYZ = clist(pc).coords';
    clist(pc).numVox = clist(pc).size;
    clist(pc).mm_center = tvmpm * [clist(pc).center'; 1];
    clist(pc).mm_center = clist(pc).mm_center(1:3)';
    actparcel(:, pc) = any(tmaps(:, clist(pc).coordi), 2);
end

% pass into Tor's function (for now)
MVAR = struct('clusters', clist, 'studyByCluster', actparcel);
% MVAR = meta_SOMclusters('mca',MVAR, 'cluster');
