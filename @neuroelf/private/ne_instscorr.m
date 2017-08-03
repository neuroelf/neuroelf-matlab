function varargout = ne_instscorr(varargin)
% ne_instscorr  - instantaneous seed-correlation
%
% FORMAT:       ne_instscorr(SRC, EVT, [tcobj [, crd [, rad]]])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       tcobj       optional time-course object (for now VTC only)
%       crd         coordinate (default, current coordinate)
%       rad         radius (default from configuration)
%
% No output fields.
%
% Example:
%
%     ne_instscorr(0, 0, rvtc, [16, -24, -16]);
%
%     this correlates the VTC (rvtc) with the seed from (TAL/MNI)
%     coordinate [16, -24, -16] and stores the output in
%     ne_gcfg.fcfg.StatsVar.

% Version:  v1.1
% Build:    16061500
% Date:     Jun-15 2016, 12:37 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, Jochen Weber
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

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg;
cvmp = cc.StatsVar;
cvmpi = cc.StatsVarIdx;
cini = ne_gcfg.c.ini.Statistics;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get handles
if nargin < 3 || numel(varargin{3}) ~= 1 || ~isxff(varargin{3}, {'hdr', 'head',  'vtc'}) || ...
    varargin{3}.NrOfVolumes < 20
    svar = cc.instscvar;
    if ~isxff(svar, {'hdr', 'head', 'vtc'}) || svar.NrOfVolumes < 20
        return;
    end
else
    svar = varargin{3};
end
if nargin < 4 || ~isa(varargin{4}, 'double') || size(varargin{4}, 2) ~= 3 || ...
    any(isinf(varargin{4}(:)) | isnan(varargin{4}(:)) | abs(varargin{4}(:)) > 128)
    crd = cc.cpos;
else
    crd = varargin{4};
end
if nargin < 5 || ~isa(varargin{5}, 'double') || numel(varargin{5}) ~= 1 || ...
    isinf(varargin{5}) || isnan(varargin{5}) || varargin{5} < 0
    radius = cini.InstantSeedRadius;
else
    radius = varargin{5};
end

% pre-computations necessary
svart = lower(svar.FileType);
svarr = svar.RunTimeVars;
svarh = handles(svar);
if ~isfield(svarh, 'InstCorrTC')

    % get time course information
    if strcmp(svart, 'vtc')
        tc = double(svar.VTCData);
        tr = svar.TR / 1000;
    elseif strcmp(svart, 'hdr')
        tc = double(permute(svar.VoxelData, [4, 1, 2, 3]));
        tr = svar.ImgDim.PixSpacing(5);
        if tr == 0
            tr = 2;
        end
    elseif strcmp(svart, 'head')
        tc = svar.Brick;
        tc = double(permute(cat(4, tc.Data), [4, 1, 2, 3]));
        tr = svar.TimeStep;
        if tr == 0
            tr = 2;
        end
    else
        return;
    end
    tn = size(tc, 1);
    vn = round(numel(tc) / tn);
    tl = tn * tr;

    % mask (implicit)
    if cini.InstantSeedMask && mean(tc(:)) > mean(std(reshape(tc, tn, vn)))
        tcmean = squeeze(mean(tc));
        tc(:, tcmean(:) < (0.5 * mean(tcmean(:)))) = 0;
    end

    % transform, leaving all 0-voxels at 0
    tc = single(ztrans(tc));

    % what needs to be regressed out
    filtstr = struct('nuisreg', zeros(tn, 0));
    if iscell(cini.InstantSeedOutFilter) && numel(cini.InstantSeedOutFilter) > 1 && ...
        ischar(cini.InstantSeedOutFilter{1}) && ...
        any(strcmpi(cini.InstantSeedOutFilter{1}, {'dct', 'fourier', 'poly'})) && ...
        isa(cini.InstantSeedOutFilter{2}, 'double') && numel(cini.InstantSeedOutFilter{2}) == 1 && ...
       ~isinf(cini.InstantSeedOutFilter{2}) && ~isnan(cini.InstantSeedOutFilter{2}) && ...
        cini.InstantSeedOutFilter{2} > 0

        % set filter information
        ftype = lower(cini.InstantSeedOutFilter{1}(1));
        if ftype == 'd'
            filtstr.tempdct = cini.InstantSeedOutFilter{2} / tr;
        elseif ftype == 'f'
            filtstr.tempsc = floor(tl / cini.InstantSeedOutFilter{2});
        else
            filtstr.temppoly = floor(tl / cini.InstantSeedOutFilter{2});
        end
    end
    if cini.InstantSeedOutGlob
        filtstr.nuisreg(:, end+1) = ztrans(mean(mean(mean(tc, 4), 3), 2));
    end
    if cini.InstantSeedOutMP && ...
        isfield(svarr, 'MotionParameters') && ...
        isequal(size(svarr.MotionParameters), [tn, 6])
        mpars = svarr.MotionParameters;
        if cini.InstantSeedOutMPd
            mparsd = 10 .* ( ...
                flexinterpn_method(mpars, [Inf, Inf; 1.05, 1; 1, 1; tn + 0.1, 6]) - ...
                flexinterpn_method(mpars, [Inf, Inf; 0.95, 1; 1, 1; tn, 6]));
            mpars = [mpars, mparsd];
        end
        filtstr.nuisreg(:, end+1:end+size(mpars, 2)) = ...
            ztrans(orthvecs(ztrans(mpars)));
    end

    % apply filter
    if numel(fieldnames(filtstr)) > 1 || ~isempty(filtstr.nuisreg)
        [tc, filtstr] = tempfilter(tc, filtstr);
        tc = single(ztrans(tc));
    else
        filtstr = filtstr.nuisreg;
    end

    % take out first-order AR(1)
    if cini.InstantSeedDiff

        % regress values against following
        xtc = ztrans(tc(1:end-1, :, :, :));
        xtcc = sqrt(1 / (vn - 1)) .* sum(xtc .* tc(2:end, :, :, :), 1);

        % remove (back to front)
        for vc = tn:-1:2
            tc(vc, :, :, :) = tc(vc, :, :, :) - xtcc .* tc(vc-1, :, :, :);
        end
        tc = single(ztrans(tc));
    end

    % store output
    svar.SetHandle('InstCorrTC', tc);
    svar.SetHandle('InstCorrTCFilt', filtstr);

% get filtered time course
else
    tc = svarh.InstCorrTC;
    filtstr = svarh.InstCorrTCFilt;
end
tcsize = size(tc);

% current StatsVar is not a InstCorrMap (VMP)
if svart(1) == 'v'
    bbox = svar.BoundingBox;
else
    bbox = svar.CoordinateFrame;
    bbox.BBox = [zeros(1, 3); bbox.Dimensions(1:3)];
    bbox.DimXYZ = bbox.Dimensions(1:3);
    bbox.FCube = [256, 256, 256];
    bbox.ResXYZ = [1, 1, 1];
    bbox.TrfPlus = bbox.Trf;
end
if numel(cvmp) ~= 1 || numel(cvmpi) ~= 1 || ~isxff(cvmp, 'vmp') || ...
   ~isfield(cvmp.Map, 'RunTimeVars') || ...
   ~isstruct(cvmp.Map(cvmpi).RunTimeVars) || ...
   ~isfield(cvmp.Map(cvmpi).RunTimeVars, 'InstCorrMap') || ...
   ~cvmp.Map(cvmpi).RunTimeVars.InstCorrMap || ...
   ~isequal(bbox.BBox, cvmp.BoundingBox.BBox)

    % create new vmp
    if isfield(bbox, 'QuatB2T')
        cvmp = newnatresvmp(bbox.BBox, bbox.ResXYZ(1));
    else
        cvmp = newnatresvmp([bbox.BBox(1, :); bbox.BBox(2, :) - 1], bbox.ResXYZ(1));
        cvmp.RunTimeVars.TrfPlus = bbox.Trf * cvmp.BoundingBox.QuatT2B;
    end
    df1 = tcsize(1) - (2 + size(filtstr, 2));
    svarname = svar.FilenameOnDisk;
    if isempty(svarname)
         svarname = 'unsaved';
    else
        [svarpath, svarname] = fileparts(svarname);
    end
    cvmp.Map(1).Name = sprintf('Instant-Seed correlation of %s.vtc', svarname);
    cvmp.Map(1).Type = 2;
    cvmp.Map(1).LowerThreshold = correlinvtstat(-sdist('tinv', 0.0025, df1), df1);
    cvmp.Map(1).UpperThreshold = correlinvtstat(-sdist('tinv', 0.00005, df1), df1);
    cvmp.Map(1).DF1 = df1;
    cvmp.Map(1).RunTimeVars = struct('InstCorrMap', true);
    cvmpi = 1;

    % add spatial normalization
    cvmp.RunTimeVars.TrfPlus = svarr.TrfPlus * cvmp.RunTimeVars.TrfPlus;
    if isfield(svarr, 'SPMsn') && numel(svarr.SPMsn) == 1 && ...
        isstruct(svarr.SPMsn) && isfield(svarr.SPMsn, 'VG') && ...
        isfield(svarr.SPMsn, 'VF') && isfield(svarr.SPMsn, 'Tr') && ...
        isfield(svarr.SPMsn, 'Affine')
        cvmp.RunTimeVars.SPMsn = svarr.SPMsn;
    else
        cvmp.RunTimeVars.SPMsn = [];
    end

    % add to GUI
    incb = ne_gcfg.c.incb;
    nuf = ne_gcfg.fcfg.noupdate;
    ne_gcfg.c.incb = false;
    ne_gcfg.fcfg.noupdate = true;
    ne_openfile(0, 0, cvmp, true);
    ne_gcfg.c.incb = incb;
    ne_gcfg.fcfg.noupdate = nuf;
end

% not single coordinate?
if svart(1) == 'v'
    radius = radius / svar.Resolution;
end
if radius >= 1 && size(crd, 1) == 1
    rx = (-floor(radius):floor(radius));
    [rx, ry, rz] = ndgrid(rx, rx, rx);
    rx = [rx(:), ry(:), rz(:)];
    rx(sqrt(sum(rx .* rx, 2)) > radius, :) = [];
    if svart(1) == 'v'
        crd = crd(ones(size(rx, 1), 1), :) + svar.Resolution .* rx;
    end
end

% convert coordinate
snm = cvmp.RunTimeVars.SPMsn;
if isequal(cvmp.RunTimeVars.TrfPlus, eye(4)) && ~isstruct(snm)
    crdv = bvcoordconv(crd, 'tal2bvx', svar.BoundingBox);
    crdv(isnan(crdv)) = [];
else
    ftrf = bvcoordconv(zeros(0, 3), 'tal2bvc', bbox) * ...
        inv(cvmp.RunTimeVars.TrfPlus);
    if ~isstruct(snm)
        crd(:, 4) = 1;
        crdv = crd * ftrf';
        crd(:, 4) = [];
    else
        ivgm = inv(snm.VG(1).mat);
        crdv = applyspmsnc(crd, snm.Tr, snm.VG(1).dim, ivgm, ...
            ftrf * snm.VF.mat * snm.Affine);
    end
    crdv = round(crdv);
    if svart(1) ~= 'v'
        crdv = unique(crdv, 'rows');
    end
    crdv(any(crdv < 1, 2) | crdv(:, 1) > tcsize(2) | crdv(:, 2) > tcsize(3) | crdv(:, 3) > tcsize(4), :) = [];
    if ~isempty(crdv)
        crdv = sub2ind(tcsize(2:4), crdv(:, 1), crdv(:, 2), crdv(:, 3));
    end
end

% disable clustering
cvmp.Map(cvmpi).EnableClusterCheck = 0;
ne_gcfg.h.Stats.UsekThr.Value = 0;

% nothing else to do
if isempty(crdv)

    % clear VMP
    cvmp.Map(cvmpi).VMPData(:) = 0;
    return;
end

% extract time course
xtc = ztrans(mean(tc(:, crdv), 2));

% run correlation and store in result
cvmp.Map(cvmpi).VMPData(:) = single((1 / (tcsize(1) - 1)) .* ...
    xtc' * reshape(tc, tcsize(1), prod(tcsize(2:4))))';

% keep track of coordinates (for later)
cvmp.Map(cvmpi).RunTimeVars.SeedCoords = crd;
