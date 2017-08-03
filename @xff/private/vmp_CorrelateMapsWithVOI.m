function xo2 = vmp_CorrelateMapsWithVOI(xo, voi, opts)
% VMP::CorrelateMapsWithVOI  - correlate a series of maps with VOI extracts
%
% FORMAT:       rvmp = vmp.ComputeFormula(voi [, opts])
%
% Input fields:
%
%       voi         either Cx3 coordinates (MNI/TAL notation) or 1x1 VOI
%       opts        optional settings
%       .estfwhm    estimate smoothness and store in Map.RunTimeVars (true)
%       .fisherize  fisherize r statistic, either of
%                   {0} - as direct r maps
%                    1  - fisherize without number of datapoints
%                    2  - fisherize with number of datapoints
%       .globsig    HDR/MSK/VMR object with global signal mask (none)
%       .mapsel     selection of maps, either regexpi or Mx1 double (all)
%       .robust     perform regression robustly (default: false)
%       .voicmb     flag, combine multiple VOIs into one (default: false)
%       .voisel     index/indices into VOI (only for object, default: 1)
%
% Output fields:
%
%       rvmp        VMP with correlation map(s)
%
% Using: bvcoordconv, calcbetas, correlinvtstat, fisherr2z,
%        fitrobustbisquare_img, glmtstat, lsqueeze, meannoinfnan,
%        resestsmooth, robustt, ztrans.

% Version:  v1.1
% Build:    16021315
% Date:     Feb-13 2016, 3:21 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2011, 2014, 2016, Jochen Weber
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

% importing from neuroelf library
using(neuroelf, {'bvcoordconv', 'calcbetas', 'correlinvtstat', 'fisherr2z', ...
    'fitrobustbisquare_img', 'glmtstat', 'lsqueeze', 'meannoinfnan', ...
    'resestsmooth', 'robustt', 'ztrans'});

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp') || ...
    isempty(voi) || ((~isa(voi, 'double') || size(voi, 2) ~= 3 || ...
     any(isinf(voi(:)) | isnan(voi(:)) | voi(:) < -128 | voi(:) > 128)) && ...
    (numel(voi) ~= 1 || ~xffisobject(voi, true, {'hdr', 'msk', 'vmr' 'voi'})))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 3 || numel(opts) ~= 1 || ~isstruct(opts)
    opts = struct;
end
bc = xo.C;
maps = bc.Map;
nmaps = numel(maps);
if nmaps < 3
    error('neuroelf:xff:badObject', 'VMP requires at least 3 maps to work.');
end
vmpfile = xo.F;
vmpid = xo.L;
if isempty(vmpfile)
    vmpfile = vmpid;
end
if ~isempty(xo.F)
    [voip, vmpname] = fileparts(xo.F);
else
    vmpname = 'VMP';
end
bbx = aft_BoundingBox(xo);
szmap = size(bc.Map(1).VMPData);
if numel(voi) == 1
    if xffisobject(voi, true, 'voi')
        voic = voi.C;
    else
        voic = xffnewcont('voi');
        voi = (aft_SampleBVBox(voi, bbx) >= 0.5);
        voic.VOI = struct('Name', 'VOI', 'Color', [128, 128, 128], ...
            'NrOfVoxels', sum(voi(:)), 'Voxels', bvcoordconv(find(voi(:)), 'bvx2tal', bbx));
        if ~isempty(voi.F)
            [voip, voic.VOI.Name, voie] = fileparts(voi.F);
            voic.VOI.Name = [voic.VOI.Name voie];
        end
    end
else
    voic = xffnewcont('voi');
    voic.VOI = struct('Name', 'VOI', 'Color', [128, 128, 128], ...
        'NrOfVoxels', size(voi, 1), 'Voxels', voi);
    opts.voisel = 1;
end

% check/parse options
if ~isfield(opts, 'estfwhm') || ~islogical(opts.estfwhm) || numel(opts.estfwhm) ~= 1
    opts.estfwhm = true;
end
if ~isfield(opts, 'fisherize') || ~isa(opts.fisherize, 'double') || numel(opts.fisherize) ~= 1 || ...
    isinf(opts.fisherize) || isnan(opts.fisherize) || opts.fisherize <= 0
    opts.fisherize = 0;
elseif opts.fisherize > 1
    opts.fisherize = 2;
else
    opts.fisherize = 1;
end
if ~isfield(opts, 'globsig') || isempty(opts.globsig) || ...
   (~islogical(opts.globsig) && (numel(opts.globsig) ~= 1 || ...
     ~xffisobject(opts.globsig, true, {'hdr', 'msk', 'vmr'})))
    opts.globsig = false;
elseif numel(opts.globsig) == 1 && xffisobject(opts.globsig, true)
    opts.globsig = (aft_SampleBVBox(opts.globsig, bbx) >= 0.5);
elseif numel(opts.globsig) > 1 && ~isequal(size(opts.globsig), szmap)
    opts.globsig = false;
end
if numel(opts.globsig) > 1 && ~any(opts.globsig(:))
    opts.globsig = false;
end
if ~isfield(opts, 'mapsel') || isempty(opts.mapsel) || ...
   ((~isa(opts.mapsel, 'double') || numel(opts.mapsel) ~= max(size(opts.mapsel)) || ...
     any(isinf(opts.mapsel) | isnan(opts.mapsel) | opts.mapsel < 1 | opts.mapsel > nmaps)) && ...
    ~ischar(opts.mapsel))
    opts.mapsel = (1:nmaps)';
elseif ischar(opts.mapsel)
    mapnames = {maps.Name};
    opts.mapsel = find(~cellfun('isempty', regexpi(mapnames(:), opts.mapsel(:)')));
else
    opts.mapsel = unique(round(opts.mapsel(:)));
end
if numel(opts.mapsel) < 3 || (~isequal(opts.globsig, false) && numel(opts.mapsel) < 4)
    error('neuroelf:xff:badArgument', ...
        'Map selection requires at least 3 (or 4 with globsig) items to work.');
end
maps = bc.Map(opts.mapsel);
nmaps = numel(maps);
if ~isfield(opts, 'robust') || ~islogical(opts.robust) || numel(opts.robust) ~= 1
    opts.robust = false;
end
if ~isfield(opts, 'voicmb') || ~islogical(opts.voicmb) || numel(opts.voicmb) ~= 1
    opts.voicmb = false;
end
if ~isfield(opts, 'voisel') || ~isa(opts.voisel, 'double') || isempty(opts.voisel) || ...
    any(isinf(opts.voisel(:)) | isnan(opts.voisel(:)) | opts.voisel(:) < 1 | opts.voisel(:) > numel(voic.VOI))
    opts.voisel = 1;
else
    opts.voisel = unique(round(opts.voisel(:)));
end

% access VOIs
voic = voic.VOI(opts.voisel);
if opts.voicmb && numel(voic) > 1
    voic(1).Voxels = cat(1, voic(:).Voxels);
    voic = voic(1);
    voic.Name = [voic.Name ' - combined'];
    voic.NrOfVoxels = size(voic.Voxels, 1);
end
nvois = numel(voic);
for vc = 1:nvois
    voic(vc).Voxels = bvcoordconv(voic(vc).Voxels, 'tal2bvx', bbx);
end

% make sure maps are data
for mc = 1:nmaps
    maps(mc).VMPData = double(maps(mc).VMPData);
end
maps = cat(4, maps(:).VMPData);
maps(isinf(maps) | isnan(maps)) = 0;

% access global signal
if isequal(opts.globsig, false)
    opts.globsig = [];
elseif isequal(opts.globsig, true)
    opts.globsig = zeros(nmaps, 1);
    for mc = 1:nmaps
        opts.globsig(mc) = meannoinfnan(lsqueeze(maps(:, :, :, mc)), 1, true);
    end
else
    gsmask = opts.globsig;
    opts.globsig = zeros(nmaps, 1);
    for mc = 1:nmaps
        vd = maps(:, :, :, mc);
        opts.globsig(mc) = meannoinfnan(lsqueeze(vd(gsmask)), 1, true);
    end
    if all(opts.globsig == 0)
        opts.globsig = [];
    end
end
if isempty(opts.globsig)
    opts.globsig = zeros(nmaps, 0);
end

% create and prepare output
xo2 = aft_CopyObject(xo);
vmpc = xo2.C;
vmpc.NrOfTimePoints = 0;
vmpc.NrOfMapParameters = 0;
vmpc.ShowParamsRangeFrom = 0;
vmpc.ShowParamsRangeTo = 0;
vmpc.FingerprintParamsRangeFrom = 0;
vmpc.FingerprintParamsRangeTo = 0;
if ~isempty(vmpc.MapParameter)
    vmpc.MapParameter(:) = [];
end
map1 = vmpc.Map(1);
if opts.fisherize > 0
    map1.Type = 12;
    if opts.fisherize > 1
        map1.LowerThreshold = 1;
        map1.UpperThreshold = 3;
    else
        map1.LowerThreshold = 0.75;
        map1.UpperThreshold = 2;
    end
else
    map1.Type = 2;
    map1.LowerThreshold = 0.5;
    map1.UpperThreshold = 1;
end
map1.RGBLowerThreshPos = [255, 0, 0];
map1.RGBUpperThreshPos = [255, 255, 0];
map1.RGBLowerThreshNeg = [255, 0, 255];
map1.RGBUpperThreshNeg = [0, 0, 255];
map1.UseRGBColor = 0;
map1.LUTName = '<default>';
map1.TransColorFactor = 1;
map1.NrOfLags = 0;
map1.MinLag = 0;
map1.MaxLag = 0;
map1.CCOverlay = 0;
map1.ClusterSize = 10;
map1.EnableClusterCheck = 0;
map1.UseValuesAboveThresh = 1;
map1.DF1 = nmaps - (2 + double(~isempty(opts.globsig)));
map1.DF2 = 0;
map1.ShowPositiveNegativeFlag = 3;
map1.NrOfFDRThresholds = 0;
map1.FDRThresholds = zeros(0, 3);
map1.UnknownValue = 0;
map1.TimePointData = zeros(0, 1);
map1.VMPDataCT = [];
map1.OverlayColors = [];
map1.RunTimeVars = struct('XType', 'VOICorrelation', 'GlobalSignal', opts.globsig, ...
    'SourceVMP', vmpfile, 'SourceVMPID',  vmpid, 'MapSel', opts.mapsel, ...
    'Robust', opts.robust, 'VOIVoxIdx', []);
if opts.estfwhm
    map1.RunTimeVars.FWHMResEst = [Inf, Inf, Inf];
end
if nvois > 1
    vmpc.Map = repmat(map1, 1, nvois);
else
    vmpc.Map = map1;
end

% remove all else from RunTimeVars
vmpc.RunTimeVars = struct('xffID', vmpc.RunTimeVars.xffID, ...
    'AutoSave', true, 'TrfPlus', vmpc.RunTimeVars.TrfPlus);

% fisherr2z argument
if opts.fisherize > 1
    fr2zn = {nmaps - double(~isempty(opts.globsig))};
else
    fr2zn = {};
end

% iterate over VOIs
for vc = 1:nvois

    % set map name
    vmpc.Map(vc).Name = sprintf('VOI correlation (%s w/ %s)', vmpname, voic(vc).Name);

    % extract VOI data from maps
    voisig = zeros(nmaps, 1);
    for mc = 1:nmaps
        vd = maps(:, :, :, mc);
        voisig(mc) = meannoinfnan(lsqueeze(vd(voic(vc).Voxels)), 1, true);
    end

    % create design matrix
    dm = [ztrans([voisig, opts.globsig]), ones(nmaps, 1)];

    % regression
    if ~opts.robust
        [b, irtc, ptc, se] = calcbetas(dm, maps, 4);
        t = glmtstat([1, zeros(1, size(dm, 2) - 1)], b, irtc, se);
        if opts.estfwhm
            [vmpc.Map(vc).RunTimeVars.FWHMResEst, vmpc.Map(vc).RunTimeVars.FWHMResImg] = ...
                resestsmooth(maps - ptc, bc.Resolution);
        end
    else
        [b, w] = fitrobustbisquare_img(dm, maps);
        t = robustt(dm, maps, b, w, [1, zeros(1, size(dm, 2) - 1)]);
        if opts.estfwhm
            ptc = zeros(size(maps));
            for bmc = 1:size(dm, 2)
                ptc = ptc + repmat(b(:, :, :, bmc), [1, 1, 1, size(dm, 1)]) .* ...
                    repmat(reshape(dm(:, bmc), [1, 1, 1, size(dm, 1)]), szmap);
            end
            ptc = w .* ptc + (1 - w) .* maps;
            [vmpc.Map(vc).RunTimeVars.FWHMResEst, vmpc.Map(vc).RunTimeVars.FWHMResImg] = ...
                resestsmooth(maps - ptc, bc.Resolution);
        end
    end
    t(isinf(t) | isnan(t)) = 0;
    t = correlinvtstat(t, nmaps - double(~isempty(opts.globsig)));

    % fisherize
    if opts.fisherize > 0
        t = fisherr2z(t, false, fr2zn{:});
    end

    % store data
    vmpc.Map(vc).VMPData = single(t);

    % update RunTimeVars
    vmpc.Map(vc).RunTimeVars.VOIVoxIdx = voic(vc).Voxels;
end

% set target
vmpc.NrOfMaps = numel(bc.Map);
xo2.C = vmpc;
