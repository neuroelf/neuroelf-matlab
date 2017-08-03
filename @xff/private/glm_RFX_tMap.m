function map = glm_RFX_tMap(xo, c, opts)
% GLM::RFX_tMap  - calculate a t contrast map
%
% FORMAT:       map = glm.RFX_tMap([c, opts])
%
% Input fields:
%
%       c           CxN contrast vector (default: full model and main eff)
%       opts        structure with optional fields
%        .alphasim  call alphasim with data to assess CDT-based FPs (false)
%        .autobal   auto-balance contrast if pos and abs(neg) == 1 (true)
%        .bbox      bounding box for VMP if necessary (default: full MNI)
%        .brange    acceptible beta-value range (default: [-Inf, Inf])
%        .covs      SxC covariates
%        .estfwhm   estimate smoothness and store in Map.RunTimeVars (true)
%        .estsmap   store estimated smoothness in separate map (false)
%        .groups    Gx2 cell array with group names and subject IDs
%        .imeth     interpolation method if necessary (default: 'cubic')
%        .interp    mesh-based interpolation (default: true)
%        .meanr     add mean regressor of (p>0.25)-voxels of first pass
%        .meanrmsk  mask additionally applied to mean-regressor source
%        .names     map names for contrasts, if not given: list of weights
%        .resvtc    create and open residual VTC in GUI (default: false)
%        .rfxtype   either {'rfx'}, 'vrfx' (global-variance weighted)
%                   or 'wrfx' (requires a GLM created with persubjectglms)
%        .robust    use robust regression for estimation (default: false)
%        .robwmaps  create robust-weight summary maps (default: false)
%        .smk       smoothing kernel for data (prior to regression, 0)
%        .srf       SRF object (required for multi-study MTC-based stats)
%        .subsel    subject selection (numeric list)
%        .swmaps    weight map per subject and regression (default: false)
%
% Output fields:
%
%       map         MAP/VMP/SMP object with C maps
%
% Note: robust regression is only implemented for rfxtype='rfx'
%
% Using: calcbetas, conjval, dilate3d, findfirst, fitrobustbisquare_img,
%        glmtstat, limitrangec, lsqueeze, mtimesnd, resestsmooth (+srf),
%        robustnsamplet_img, robustt, sdist, smoothdata3, varc, ztrans.

% Version:  v1.1
% Build:    16101016
% Date:     Oct-10 2016, 4:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010 - 2016, Jochen Weber
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
alphasim    = ne_methods.alphasim;
calcbetas   = ne_methods.calcbetas;
dilate3d    = ne_methods.dilate3d;
findfirst   = ne_methods.findfirst;
glmtstat    = ne_methods.glmtstat;
limitrangec = ne_methods.limitrangec;
lsqueeze    = ne_methods.lsqueeze;
mtimesnd    = ne_methods.mtimesnd;
robustt     = ne_methods.robustt;
sdist       = ne_methods.sdist;
varc        = ne_methods.varc;
varweights  = ne_methods.varweights;
ztrans      = ne_methods.ztrans;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
bcrtv = bc.RunTimeVars;
glmfile = xo.F;
glmid = xo.L;
if isempty(glmfile)
    glmfile = glmid;
end
if bc.ProjectTypeRFX == 0 && bc.SeparatePredictors ~= 2
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
isrfx = (bc.ProjectTypeRFX > 0);
ffxspred = glm_SubjectPredictors(xo);
ffxsubs = glm_Subjects(xo);
if isrfx
    numsubs = bc.NrOfSubjects;
    numspred = bc.NrOfSubjectPredictors;
else
    ffxpred = bc.Predictor;
    ffxpred = {ffxpred(:).Name2};
    ffxpred = ffxpred(:);
    numsubs = numel(ffxsubs);
    numspred = numel(ffxspred) + 1;
end
if nargin < 2 || ~isa(c, 'double') || isempty(c) || ...
    size(c, 1) > numspred || any(isinf(c(:)) | isnan(c(:)))
    c = [[eye(numspred - 1);zeros(1, numspred - 1)], [ones(numspred - 1, 1); 0]];
elseif numel(c) == size(c, 2)
    if numel(c) > numspred
        c((numspred + 1):end) = [];
    end
    c = real(c');
else
    c = real(c);
end
nummaps = size(c, 2);
if size(c, 1) < numspred
    c = [c; zeros(numspred - size(c, 1), size(c, 2))];
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
ipo = (bc.ProjectType == 2);
if ~isfield(opts, 'alphasim') || ~islogical(opts.alphasim) || numel(opts.alphasim) ~= 1
    opts.alphasim = false;
else
    athrs = (0.001 .* [50, 20, 10, 5, 2, 1, .5, .1, .05, .01, .005, .001])';
end
if ~isfield(opts, 'autobal') || ~islogical(opts.autobal) || numel(opts.autobal) ~= 1
    opts.autobal = true;
end
if ~isfield(opts, 'bbox') || ~isa(opts.bbox, 'double') || ...
   ~isequal(size(opts.bbox), [2, 3]) || ...
    any(isinf(opts.bbox(:)) | isnan(opts.bbox(:)) | ...
        opts.bbox(:) < 0 | opts.bbox(:) > 256 | opts.bbox(:) ~= fix(opts.bbox(:))) || ...
    any(diff(opts.bbox) < 0)
    opts.bbox = [44, 38, 44; 242, 194, 212];
end
if ~isfield(opts, 'brange') || ~isa(opts.brange, 'double') || numel(opts.brange) ~= 2 || ...
    any(isnan(opts.brange)) || opts.brange(1) >= opts.brange(2)
    opts.brange = [-Inf, Inf];
else
    brangedop = false(3, 3, 3);
    brangedop([5, 11, 13, 14, 15, 17, 23]) = true;
end
brange = opts.brange;
if ~isfield(opts, 'covs') || ~isa(opts.covs, 'double') || ndims(opts.covs) ~= 2 || ...
    size(opts.covs, 1) ~= numsubs
    opts.covs = zeros(numsubs, 0);
end
if ~isfield(opts, 'estfwhm') || ~islogical(opts.estfwhm) || numel(opts.estfwhm) ~= 1
    opts.estfwhm = true;
end
if ~isfield(opts, 'estsmap') || ~islogical(opts.estsmap) || numel(opts.estsmap) ~= 1
    opts.estsmap = false;
elseif opts.estsmap
    opts.estfwhm = true;
end
if ~isfield(opts, 'groups') || ~iscell(opts.groups) || size(opts.groups, 2) ~= 2
    ngrp = 0;
    if ~isfield(opts, 'subsel') || ~isa(opts.subsel, 'double') || ...
        isempty(opts.subsel) || any(isinf(opts.subsel(:)) | isnan(opts.subsel(:)))
        ga = ones(numsubs, 1);
    else
        opts.subsel = opts.subsel(:)';
        opts.subsel(opts.subsel < 1 | opts.subsel > numsubs) = [];
        ga = zeros(numsubs, 1);
        ga(unique(round(opts.subsel))) = 1;
    end
    opts.groups = [];
    gamx = 1;
    nval = sum(ga) - 1;
else
    ga = zeros(numsubs, 1);
    for gc = 1:size(opts.groups, 1)
        if ~ischar(opts.groups{gc, 1}) || isempty(opts.groups{gc, 1}) || ...
            isempty(opts.groups{gc, 2}) || ...
            any(isinf(opts.groups{gc, 2}(:)) | isnan(opts.groups{gc, 2}(:))) || ...
            any(opts.groups{gc, 2}(:) < 1 | opts.groups{gc, 2}(:) > numsubs) || ...
            any(ga(round(opts.groups{gc, 2}(:))) > 0)
            error('neuroelf:xff:badArgument', 'Invalid group assignment.');
        end
        opts.groups{gc, 2} = unique(round(opts.groups{gc, 2}(:)));
        ga(opts.groups{gc, 2}) = gc;
    end
    ngrp = size(opts.groups, 1);
    nval = sum(ga > 0) - ngrp;
    if nval < 3
        error('neuroelf:xff:badArgument', 'Too few subjects for grouping.');
    end
    gas = ga(ga > 0);
    gag = tril(ones(ngrp));
    gag(1:(ngrp + 1):end) = 0;
    gamx = sum(gag(:));
    gag(gag > 0) = 1:gamx;
end
gax = find(ga > 0);
numcovs = size(opts.covs, 2);
if numcovs > 0
    if ngrp < 2
        opts.covs = opts.covs(gax, :);
        opts.covs = ztrans(opts.covs);
    else
        for gc = 1:ngrp
            opts.covs(ga == gc, :) = ztrans(opts.covs(ga == gc, :));
        end
        opts.covs = opts.covs(gax, :);
    end
else
    opts.covs = zeros(numel(gax), 0);
end
numsubs = numel(gax);
if ~isfield(opts, 'imeth') || ~ischar(opts.imeth) || isempty(opts.imeth) || ...
    isempty(regexpi(opts.imeth(:)', '^(cubic|lanczos\d|linear|nearest)$'))
    opts.imeth = 'cubic';
else
    opts.imeth = lower(opts.imeth(:)');
end
if isfield(opts, 'interp') && ~isempty(opts.interp) && ...
   (isnumeric(opts.interp) || islogical(opts.interp))
    ipo = ipo && opts.interp(1);
end
if ~isfield(opts, 'meanr') || ~islogical(opts.meanr) || numel(opts.meanr) ~= 1
    opts.meanr = false;
end
if ~isfield(opts, 'meanrmsk') || ~islogical(opts.meanrmsk)
    opts.meanrmsk = [];
end
if ~isfield(opts, 'names') || ~iscell(opts.names) || numel(opts.names) ~= nummaps
    opts.names = {};
else
    for mc = 1:nummaps
        if ~ischar(opts.names{mc}) || isempty(opts.names{mc})
            opts.names = {};
            break;
        end
        opts.names{mc} = opts.names{mc}(:)';
    end
end
if isempty(opts.names)
    opts.names = cell(nummaps, 1);
    for mc = 1:nummaps
        opts.names{mc} = sprintf('%d ', c(1:end-1, mc)');
        opts.names{mc}(end) = [];
    end
end
if ~isfield(opts, 'resvtc') || ~islogical(opts.resvtc) || numel(opts.resvtc) ~= 1
    opts.resvtc = false;
end
resvtc = [];
if ~isfield(opts, 'rfxtype') || ~ischar(opts.rfxtype) || isempty(opts.rfxtype) || ...
   ~any(strcmpi(opts.rfxtype(:)', {'r', 'rfx', 'v', 'vrfx', 'w', 'wrfx'})) || ...
    numsubs < 3 || bc.SeparatePredictors ~= 2
    opts.rfxtype = 'r';
else
    opts.rfxtype = lower(opts.rfxtype(1));
    if opts.rfxtype == 'w'
        if ~isfield(bc.GLMData, 'RunTimeVars') || ~isstruct(bc.GLMData.RunTimeVars) || ...
            numel(bc.GLMData.RunTimeVars) ~= 1 || ...
           ~isfield(bc.GLMData.RunTimeVars, 'Subject') || ...
           ~isstruct(bc.GLMData.RunTimeVars.Subject) || ...
            numel(bc.GLMData.RunTimeVars.Subject) ~= numsubs || ...
           ~isfield(bc.GLMData.RunTimeVars.Subject, 'iXX') || ...
           ~isequal(size(bc.GLMData.RunTimeVars.Subject(1).iXX), [numspred, numspred]) || ...
           ~isfield(bc.GLMData.RunTimeVars.Subject, 'SEMap') || ...
           ~isequal(size(bc.GLMData.RunTimeVars.Subject(1).SEMap), size(bc.GLMData.RFXGlobalMap))
            opts.rfxtype = 'r';
        else
            gsrtvs = bc.GLMData.RunTimeVars.Subject;
            if ~all(cellfun('prodofsize', {gsrtvs.SEMap}) == numel(bc.GLMData.RFXGlobalMap)) || ...
               ~all(cellfun('prodofsize', {gsrtvs.iXX}) == (numspred * numspred))
                opts.rfxtype = 'r';
            else
                gsixx = {gsrtvs.iXX};
                gssem = {gsrtvs.SEMap};
            end
        end
    end
end
if ~isfield(opts, 'robust') || numel(opts.robust) ~= 1 || ~islogical(opts.robust)
    opts.robust = false;
end
if opts.robust
    tmapr = ' (robust)';
else
    tmapr = '';
end
if ~isfield(opts, 'robwmaps') || ~islogical(opts.robwmaps) || numel(opts.robwmaps) ~= 1
    opts.robwmaps = false;
end
if ~isfield(opts, 'smk') || numel(opts.smk) ~= 1 || ~isa(opts.smk, 'double') || ...
    isinf(opts.smk) || isnan(opts.smk) || opts.smk <= (0.5 * bc.Resolution)
    opts.smk = 0;
else
    opts.smk = min(opts.smk, 8 * bc.Resolution);
end
if ~isfield(opts, 'srf') || numel(opts.srf) ~= 1 || ~xffisobject(opts.srf, true, 'srf')
    opts.srf = [];
end
if ~isfield(opts, 'swmaps') || ~islogical(opts.swmaps) || numel(opts.swmaps) ~= 1
    opts.swmaps = false;
end
nmf = 1 + 3 * double(opts.meanr) + (1 + double(opts.meanr)) .* ...
    (2 * double(opts.robwmaps) * double(opts.robust) + ...
    numsubs * double(opts.swmaps) * double(opts.robust)) + ...
    double(opts.estsmap) * (1 + double(opts.meanr));

% collect 1st level SEmaps for non-robust non-rfx stast
if ~opts.robust && opts.rfxtype == 'w'
    ixxs = cat(3, gsixx{:});
    nixxs = size(ixxs, 3);
    semaps = cat(4, gssem{:});
end

% not yet supported for FMR-based GLMs
if bc.ProjectType == 0
    error('neuroelf:xff:notYetSupported', 'FMR->MAP generation not yet supported.');
end

% VTC-based
copymaps = true;
if bc.ProjectType == 1
    mtype = 'v';
    mdim = 4;
    mones = [1, 1, 1];
    map = xff('new:vmp');
    mapc = map.C;
    mapc.Resolution = bc.Resolution;
    if ~isfield(bcrtv, 'SubjectSPMsn') || ~isstruct(bcrtv.SubjectSPMsn) || ...
        isempty(fieldnames(bcrtv.SubjectSPMsn))
        if isrfx
            szmap = size(bc.GLMData.RFXGlobalMap);
        else
            szmap = size(bc.GLMData.MCorrSS);
        end
        mapc.XStart = bc.XStart;
        mapc.XEnd = bc.XEnd;
        mapc.YStart = bc.YStart;
        mapc.YEnd = bc.YEnd;
        mapc.ZStart = bc.ZStart;
        mapc.ZEnd = bc.ZEnd;
    else
        opts.bbox(2, :) = opts.bbox(1, :) + mapc.Resolution .* ...
            ceil((opts.bbox(2, :) - opts.bbox(1, :)) ./ mapc.Resolution - 0.01);
        mapc.XStart = opts.bbox(1, 1);
        mapc.XEnd = opts.bbox(2, 1);
        mapc.YStart = opts.bbox(1, 2);
        mapc.YEnd = opts.bbox(2, 2);
        mapc.ZStart = opts.bbox(1, 3);
        mapc.ZEnd = opts.bbox(2, 3);
        szmap = round((opts.bbox(2, :) - opts.bbox(1, :)) ./ mapc.Resolution);
        copymaps = false;
        sbbox = struct('BBox', opts.bbox, 'ResXYZ', mapc.Resolution);
    end
    nlmap = prod(szmap);
    mfld = 'VMPData';
    mapc.RunTimeVars.TrfPlus = bcrtv.TrfPlus;

    % residuals VTC container
    if opts.resvtc
        resvtc = xff('new:vtc');
        resvtcc = resvtc.C;
        resvtcc.NameOfSourceFMR = glmfile;
        resvtcc.DataType = 2;
        resvtcc.Resolution = mapc.Resolution;
        resvtcc.XStart = mapc.XStart;
        resvtcc.XEnd = mapc.XEnd;
        resvtcc.YStart = mapc.YStart;
        resvtcc.YEnd = mapc.YEnd;
        resvtcc.ZStart = mapc.ZStart;
        resvtcc.ZEnd = mapc.ZEnd;
        resvtcc.TR = 29554; % hex('sr')
    end
    numvert = nlmap;

% MTC-based
elseif bc.ProjectType == 2
    mtype = 's';
    mdim = 2;
    mones = 1;

    % SRF required
    numvert = bc.NrOfVertices;
    if ipo
        if ~isfield(opts, 'srf') || numel(opts.srf) ~= 1 || ...
           ~xffisobject(opts.srf, true, 'srf')
            error('neuroelf:xff:badArgument', ...
                'Missing or bad SRF reference in map options.');
        end
        srfc = opts.srf.C;
        if size(srfc.Neighbors, 1) ~= numvert
            error('neuroelf:xff:badArgument', 'Number of vertices mismatch.');
        end
        nei = srfc.Neighbors(:, 2);
        nnei = numel(nei);
        neil = zeros(nnei, 12);
        neis = zeros(nnei, 12);
        for nc = 1:nnei
            neis(nc) = numel(nei{nc});
            neil(nc, 1:neis(nc)) = nei{nc};
        end
        neii = (neil > 0);
    elseif opts.estfwhm
        if isempty(opts.srf) || opts.srf.C.NrOfVertices ~= numvert
            opts.estfwhm = false;
            opts.estsmap = false;
        end
    end

    % generate output
    map = xff('new:smp');
    mapc = map.C;
    mapc.FileVersion = 4;
    mapc.NrOfVertices = numvert;
    if ipo
        mapc.NameOfOriginalSRF = opts.srf.F;
    end
    szmap = numvert;
    mfld = 'SMPData';
end
szmapo = [szmap, 1];

% preset Map content
mapc.Map.Type = 1;
mapc.Map.LowerThreshold = -sdist('tinv', 0.005, nval);
mapc.Map.UpperThreshold = -sdist('tinv', 0.0001, nval);
mapc.Map.Name = 'Contrast:';
mapc.Map.NrOfLags = [];
mapc.Map.MinLag = [];
mapc.Map.MaxLag = [];
mapc.Map.CCOverlay = [];
mapc.Map.ClusterSize = 25;
mapc.Map.EnableClusterCheck = 0;
mapc.Map.DF1 = nval;
mapc.Map.DF2 = 0;
mapc.Map.BonferroniValue = bc.NrOfVoxelsForBonfCorrection;
mapc.Map.(mfld) = single(zeros(szmapo));
mapc.Map(1:(nummaps*gamx*nmf)) = mapc.Map(1);
mapc.RunTimeVars.AutoSave = true;

% iterate over contrasts
occ = 1;
for cc = 1:size(c, 2)

    % get contrast
    conc = c(:, cc);
    coni = find(conc ~= 0);
    csumc = sum(conc(coni));

    % if 1st level data is required as well
    if ~opts.robust && opts.rfxtype == 'w'
        ixxc = reshape(mtimesnd(mtimesnd(repmat(c', [1, 1, nixxs]), ixxs), ...
            repmat(c, [1, 1, nixxs])), [mones, nixxs]);
        semapc = semaps;
        semapr = 1 ./ (repmat(ixxc, szmapo) .* semaps);
        semapr(isinf(semapr) | isnan(semapr)) = 0;
    end

    % initialize temp map
    if mdim == 4
        tmpmp = zeros([szmap, numsubs]);
    else
        tmpmp = zeros(numvert, numsubs);
    end

    % auto-balancing
    if abs(csumc) < sqrt(eps) && abs(sum(conc(conc > 0)) - 1) < sqrt(eps) && ...
        abs(1 + sum(conc(conc < 0))) < sqrt(eps) && opts.autobal
        autobal = true;
        possum = zeros(numsubs, 1);
        negsum = zeros(numsubs, 1);
        ntmpmp = tmpmp;
    else
        autobal = false;
    end

    % allow for subjects with missing data (FFX)
    keepsubs = true(numsubs, 1);
    gaxk = gax;

    % fill contrast
    for pc = coni(:)'
        for sc = 1:numsubs
            if ~keepsubs(sc)
                continue;
            end
            if isrfx
                if copymaps
                    if mtype == 'v'
                        tmpmpp = bc.GLMData.Subject(gax(sc)).BetaMaps(:, :, :, pc);
                    else
                        tmpmpp = bc.GLMData.Subject(gax(sc)).BetaMaps(:, pc);
                    end
                else
                    tmpmpp = aft_SampleBVBox(xo, sbbox, ...
                        (gax(sc) - 1) * numspred + pc, opts.imeth);
                end
            else
                keepsubi = findfirst(~cellfun('isempty', regexpi(ffxpred, ...
                    sprintf('^subject\\s+%s:\\s*%s', ffxsubs{gax(sc)}, ffxspred{pc}))));
                if ~isempty(keepsubi)
                    if mtype == 'v'
                        tmpmpp = bc.GLMData.BetaMaps(:, :, :, keepsubi);
                    else
                        tmpmpp = bc.GLMData.BetaMaps(:, keepsubi);
                    end
                else
                    keepsubs(sc) = false;
                    continue;
                end
            end
            if ~isinf(brange(1))
                if mtype == 'v'
                    if ~isinf(brange(2))
                        tmpmpp(dilate3d(tmpmpp < (2 .* brange(1)) | tmpmpp > (2 .* brange(2)), brangedop) | tmpmpp < brange(1) | tmpmpp > brange(2)) = NaN;
                    else
                        tmpmpp(dilate3d(tmpmpp < (2 .* brange(1)), brangedop) | tmpmpp < brange(1)) = NaN;
                    end
                else
                    if ~isinf(brange(2))
                        tmpmpp(tmpmpp < brange(1) | tmpmpp > brange(2)) = NaN;
                    else
                        tmpmpp(tmpmpp < brange(1)) = NaN;
                    end
                end
            elseif ~isinf(brange(2))
                if mtype == 'v'
                    tmpmpp(dilate3d(tmpmpp > (2 .* brange(2)), brangedop) | tmpmpp > brange(2)) = NaN;
                else
                    tmpmpp(tmpmpp > brange(2)) = NaN;
                end
            end
            if all(tmpmpp(:) == 0 | isinf(tmpmpp(:)) | isnan(tmpmpp(:)))
                if ~autobal
                    keepsubs(sc) = false;
                end
                continue;
            end
            if autobal
                if conc(pc) > 0
                    possum(sc) = possum(sc) + conc(pc);
                    if mtype == 'v'
                        tmpmp(:, :, :, sc) = tmpmp(:, :, :, sc) + conc(pc) .* tmpmpp;
                    else
                        tmpmp(:, sc) = tmpmp(:, sc) + conc(pc) .* tmpmpp;
                    end
                else
                    negsum(sc) = negsum(sc) + conc(pc);
                    if mtype == 'v'
                        ntmpmp(:, :, :, sc) = ntmpmp(:, :, :, sc) + conc(pc) .* tmpmpp;
                    else
                        ntmpmp(:, sc) = ntmpmp(:, sc) + conc(pc) .* tmpmpp;
                    end
                end
            else
                if mtype == 'v'
                    tmpmp(:, :, :, sc) = tmpmp(:, :, :, sc) + conc(pc) .* tmpmpp;
                else
                    tmpmp(:, sc) = tmpmp(:, sc) + conc(pc) .* tmpmpp;
                end
            end
        end
    end

    % auto-balancing
    if autobal
        if ~all(abs(possum - 1) < sqrt(eps))
            tmpmp = reshape(reshape(tmpmp, numvert, numsubs) * diag(1 ./ possum), [szmap, numsubs]);
        end
        if ~all(abs(1 + negsum) < sqrt(eps))
            ntmpmp = reshape(reshape(ntmpmp, numvert, numsubs) * diag(-1 ./ negsum), [szmap, numsubs]);
        end
        tmpmp = tmpmp + ntmpmp;
        badsum = (possum == 0 | negsum == 0);
        if any(badsum)
            if mtype == 'v'
                tmpmp(:, :, :, badsum) = NaN;
            else
                tmpmp(:, badsum) = NaN;
            end
        end
    end

    % remove bad subjects from list
    if mtype == 'v'
        keepsubs = keepsubs & ~lsqueeze(all(all(all(isnan(tmpmp), 3), 2), 1));
    else
        keepsubs = keepsubs & ~lsqueeze(all(isnan(tmpmp), 1));
    end
    if ~all(keepsubs)
        if mtype == 'v'
            tmpmp(:, :, :, ~keepsubs) = [];
            gaxk(~keepsubs) = [];
            if ~opts.robust && opts.rfxtype == 'w'
                semapc(:, :, :, ~keepsubs) = [];
                semapr(:, :, :, ~keepsubs) = [];
            end
        else
            tmpmp(:, ~keepsubs) = [];
            gaxk(~keepsubs) = [];
            if ~opts.robust && opts.rfxtype == 'w'
                semapc(:, ~keepsubs) = [];
                semapr(:, ~keepsubs) = [];
            end
        end
    end
    goodsubs = sum(keepsubs);
    if ~opts.robust && opts.rfxtype == 'w'
        semapn = sum(semapr ~= 0, mdim);
    end
    if ngrp < 1
        nval = numel(gaxk) - 1;
    else
        nval = numel(gaxk) - ngrp;
    end
    if nval < 1
        continue;
    end
    mapc.Map(occ).DF1 = nval;
    if ~opts.robust && opts.rfxtype == 'w'
        semapw = semapr ./ repmat(max(semapr, [], mdim), [mones, goodsubs]);
        semapw(isinf(semapw) | isnan(semapw)) = 0;
        semapwn = sum(semapw, mdim);
    end

    % smooth data
    if opts.smk > 0 && mdim == 4
        ptc = (isinf(tmpmp) | isnan(tmpmp) | tmpmp == 0);
        tmpmp = ne_methods.smoothdata3(tmpmp, opts.smk / mapc.Resolution);
        tmpmp(ptc) = 0;
    end

    % set additional data
    artv = struct( ...
        'SourceGLM',   glmfile, ...
        'SourceGLMID', glmid, ...
        'Contrast',    conc, ...
        'Covariates',  opts.covs, ...
        'Groups',      {opts.groups}, ...
        'FWHMResEst',  [], ...
        'FWHMResImg',  [], ...
        'GlobSigMap',  [], ...
        'MeanRem',     opts.meanr, ...
        'MeanRemRes',  [], ...
        'RFXGLM',      isrfx, ...
        'Robust',      opts.robust, ...
        'SubPreds',    {ffxspred}, ...
        'SubSel',      gaxk);

    % generate 2nd-level single beta-t-maps
    if ngrp < 1

        % non-robust
        if ~opts.robust

            % no covariates
            if numcovs == 0

                % non first-level weighted
                if opts.rfxtype ~= 'w'

                    % compute OLS mean (beta)
                    mmap = mean(tmpmp, mdim);

                    % regular RFX
                    if opts.rfxtype == 'r'

                        % t-map
                        tmap = sqrt(goodsubs) .* (mmap ./ std(tmpmp, 0, mdim));

                    % variance-weighted
                    else

                        % compute residual
                        rmaps = tmpmp - repmat(mmap, [mones, goodsubs]);

                        % compute variance weights
                        gvvec = varweights(reshape(rmaps, nlmap, goodsubs));

                        % scale to median = 1
                        gvvec = gvvec ./ median(gvvec);
                        gvvec(gvvec > 1) = 1;
                        svvec = sum(gvvec);
                        gvvec = repmat(reshape(gvvec, [mones, goodsubs]), szmapo);

                        % compute weighted t-map
                        mmap = sum(gvvec .* tmpmp, mdim) ./ svvec;
                        mmap(isinf(mmap) | isnan(mmap)) = 0;
                        rmaps = gvvec .* (tmpmp - repmat(mmap, [mones, goodsubs]));
                        tmap = sqrt(svvec) .* mmap ./ sqrt(sum(rmaps .* rmaps, mdim) ./ (svvec - 1));
                    end
                    if opts.estfwhm || opts.resvtc
                        if opts.rfxtype == 'r'
                            ptc = repmat(mmap, [mones, size(tmpmp, mdim)]);
                        elseif opts.rfxtype == 'v'
                            ptc = gvvec .* repmat(mmap, [mones, size(tmpmp, mdim)]) + ...
                                (1 - gvvec) .* tmpmp;
                        else
                            ptc = semapw .* repmat(mmap, [mones, size(tmpmp, mdim)]) + ...
                                (1 - semapw) .* tmpmp;
                        end
                    end
                else
                    mmap = sum(semapw .* tmpmp, mdim) ./ semapwn;
                end
            else
                x = [ones(goodsubs, 1), opts.covs(keepsubs, :)];
                [bmaps, ixx, ptc, se] = calcbetas(x, tmpmp, mdim);
                tmap = glmtstat([1, zeros(1, numcovs)], bmaps, ixx, se);
            end
            if opts.estfwhm && mdim == 4
                [artv.FWHMResEst, artv.FWHMResImg] = ...
                    ne_methods.resestsmooth(tmpmp - ptc, bc.Resolution);
            elseif opts.estfwhm && mdim == 2
                [artv.FWHMResEst, artv.FWHMResImg] = ...
                    ne_methods.resestsmoothsrf(tmpmp - ptc, opts.srf);
            end
            if opts.alphasim && mdim == 4
                [athr, anull, amap, amax] = alphasim(size(tmap), struct( ...
                    'regmaps', tmpmp, 'regmodel', ones(goodsubs, 1), ...
                    'thr', athrs));
                kthrs = zeros(size(athrs));
                for atc = 1:numel(athrs)
                    kthrs(atc) = findfirst(athr{atc}(:, end) < 0.05);
                end
            end
        else
            x = [ones(goodsubs, 1), opts.covs(keepsubs, :)];
            [bmaps, wmaps] = ne_methods.fitrobustbisquare_img(x, tmpmp);
            tmap = robustt(x, tmpmp, bmaps, wmaps, [1, zeros(1, numcovs)]);
            if opts.estfwhm || opts.resvtc
                ptc = zeros(size(tmpmp));
                for bmc = 1:size(x, 2)
                    if mdim == 4
                        ptc = ptc + repmat(bmaps(:, :, :, bmc), [mones, size(x, 1)]) .* ...
                            repmat(reshape(x(:, bmc), [mones, size(x, 1)]), szmapo);
                    else
                        ptc = ptc + repmat(bmaps(:, bmc), [mones, size(x, 1)]) .* ...
                            repmat(reshape(x(:, bmc), [mones, size(x, 1)]), szmapo);
                    end
                end
                ptc = wmaps .* ptc + (1 - wmaps) .* tmpmp;
                if opts.estfwhm && mdim == 4
                    [artv.FWHMResEst, artv.FWHMResImg] = ...
                        ne_methods.resestsmooth(tmpmp - ptc, bc.Resolution);
                elseif opts.estfwhm && mdim == 4
                    [artv.FWHMResEst, artv.FWHMResImg] = ...
                        ne_methods.resestsmoothsrf(tmpmp - ptc, opts.srf);
                end
            end
        end
        tmap(isinf(tmap) | isnan(tmap)) = 0;

        % interpolate (only for MTC-based data anyway)
        if ipo
            rmap = zeros(numvert, 1);
            for nc = 1:size(neii, 2)
                rmap(neii(:, nc)) = rmap(neii(:, nc)) + tmap(neil(neii(:, nc), nc));
            end
            tmap = rmap ./ sum(neii, 2);
        end

        % store and visualize residuals as VTC
        if opts.resvtc && mtype == 'v'
            resvtcc.VTCData = permute(single(tmpmp - ptc), [4, 1, 2, 3]);
            resvtcc.NrOfVolumes = size(resvtcc.VTCData, 1);
            resvtcc.RunTimeVars.Subjects = ffxsubs(gaxk);
            resvtc.C = resvtcc;
            cresvtc = aft_CopyObject(resvtc);
            cresvtc.H.InstCorrTC = single(ztrans(resvtcc.VTCData));
            cresvtc.H.InstCorrTCFilt = zeros(size(resvtcc.VTCData, 1), 0);
            neuroelf_gui('openfile', cresvtc, true);
        end

        % set name and map data
        mapc.Map(occ).Name = sprintf('%s%s', opts.names{cc}, tmapr);
        mapc.Map(occ).(mfld) = single(tmap);
        mapc.Map(occ).RunTimeVars = artv;
        if ~opts.robust && opts.alphasim && mdim == 4
            mapc.Map(occ).RunTimeVars.AlphaSim12 = {0.05, [athrs, kthrs]};
            mapc.Map(occ).RunTimeVars.AlphaSim12zMaps = amap;
            mapc.Map(occ).RunTimeVars.AlphaSimMaxDist = amax;
        end
        occ = occ + 1;

        % add smoothness estimate map
        if opts.estsmap
            mapc.Map(occ).Name = sprintf('%s (smoothness estimate)', opts.names{cc});
            mapc.Map(occ).Type = 83;
            mapc.Map(occ).(mfld) = single(artv.FWHMResImg);
            if mdim == 4
                mapc.Map(occ).LowerThreshold = 2 * bc.Resolution + 0.5 * opts.smk;
                mapc.Map(occ).UpperThreshold = 6 * bc.Resolution + 2.5 * opts.smk;
            else
                mapc.Map(occ).LowerThreshold = 5;
                mapc.Map(occ).UpperThreshold = 20;
            end
            mapc.Map(occ).ShowPositiveNegativeFlag = 1;
            occ = occ + 1;
        end

        % add robust-summary information
        artv.FWHMResEst = [];
        artv.FWHMResImg = [];
        if opts.robust && opts.robwmaps
            mapc.Map(occ).Name = sprintf('%s%s (mean robust weight)', ...
                opts.names{cc}, tmapr);
            mapc.Map(occ).Type = 145;
            mapc.Map(occ).LowerThreshold = 0.75;
            mapc.Map(occ).UpperThreshold = 1;
            mapc.Map(occ).(mfld) = single((1 / size(wmaps, mdim)) .* sum(wmaps, mdim));
            mapc.Map(occ).RunTimeVars = artv;
            occ = occ + 1;
            mapc.Map(occ).Name = sprintf('%s%s (weighted number of outliers)', ...
                opts.names{cc}, tmapr);
            mapc.Map(occ).Type = 146;
            mapc.Map(occ).LowerThreshold = 1;
            mapc.Map(occ).UpperThreshold = 0.5 * size(wmaps, mdim);
            mapc.Map(occ).(mfld) = single(-3 .* sum( ...
                limitrangec(wmaps - 2/3, -1/3, 0, -1/6), mdim));
            mapc.Map(occ).RunTimeVars = artv;
            occ = occ + 1;
        end

        % add robust weight maps to VMP
        if opts.robust && opts.swmaps
            for bmc = 1:size(wmaps, mdim)
                mapc.Map(occ).Name = sprintf('%s%s (%s outlier)', ...
                    opts.names{cc}, tmapr, ffxsubs{gax(bmc)});
                mapc.Map(occ).Type = 144;
                mapc.Map(occ).LowerThreshold = 0.25;
                mapc.Map(occ).UpperThreshold = 1;
                if mdim == 4
                    mapc.Map(occ).(mfld) = single(1 - wmaps(:, :, :, bmc));
                else
                    mapc.Map(occ).(mfld) = single(1 - wmaps(:, bmc));
                end
                mapc.Map(occ).RunTimeVars = artv;
                occ = occ + 1;
            end
        end

        % create additional map without residual
        if opts.meanr

            % compute std==1-scaled data
            resmp = tmpmp;
            gsmap = 1 ./ sqrt(varc(tmpmp, mdim));
            gsmap(isinf(gsmap) | isnan(gsmap) | gsmap > 1) = 1;
            for sc = 1:goodsubs
                if mdim == 4
                    resmp(:, :, :, sc) = resmp(:, :, :, sc) .* gsmap;
                else
                    resmp(:, sc) = resmp(:, sc) .* gsmap;
                end
            end
            resmp(isinf(resmp) | isnan(resmp)) = 0;
            artv.GlobSigMap = single(gsmap);

            % average according to tmap
            tmin = -sdist('tinv', 0.25, nval);
            tmax = -sdist('tinv', 0.05, nval);
            wmp = limitrangec(tmax - (1 / (tmax - tmin)) .* tmap, 0.25, 1, 0);
            wmp = wmp .* wmp;
            if isequal(size(wmp), size(opts.meanrmsk))
                wmp = wmp .* opts.meanrmsk;
            end
            wmp(tmap == 0) = 0;
            for sc = 1:goodsubs
                if mdim == 4
                    resmp(:, :, :, sc) = resmp(:, :, :, sc) .* wmp;
                else
                    resmp(:, sc) = resmp(:, sc) .* wmp;
                end
            end
            if mdim == 4
                resmp = lsqueeze(sum(sum(sum(resmp, 1), 2), 3)) ./ sum(wmp(:));
            else
                resmp = lsqueeze(sum(resmp, 1)) ./ sum(wmp(:));
            end

            % non-robust remodeling
            x = [ones(goodsubs, 1), ztrans(resmp)];
            if ~opts.robust

                % use calcbetas and glmtstat
                [bmaps, ixx, ptc, se] = calcbetas(x, tmpmp, mdim);
                tmap = glmtstat([1, 0], bmaps, ixx, se);
                if opts.estfwhm && mdim == 4
                    [artv.FWHMResEst, artv.FWHMResImg] = ...
                        ne_methods.resestsmooth(tmpmp - ptc, bc.Resolution);
                elseif opts.estfwhm && mdim == 2
                    [artv.FWHMResEst, artv.FWHMResImg] = ...
                        ne_methods.resestsmoothsrf(tmpmp - ptc, opts.srf);
                end

            % robust remodeling
            else
                [bmaps, wmaps] = ne_methods.fitrobustbisquare_img(x, tmpmp);
                tmap = robustt(x, tmpmp, bmaps, wmaps, [1, 0]);
                if opts.estfwhm
                    ptc = zeros(size(tmpmp));
                    for bmc = 1:size(x, 2)
                        if mdim == 4
                            ptc = ptc + repmat(bmaps(:, :, :, bmc), [mones, size(x, 1)]) .* ...
                                repmat(reshape(x(:, bmc), [mones, size(x, 1)]), szmapo);
                        else
                            ptc = ptc + repmat(bmaps(:, bmc), [mones, size(x, 1)]) .* ...
                                repmat(reshape(x(:, bmc), [mones, size(x, 1)]), szmapo);
                        end
                    end
                    ptc = wmaps .* ptc + (1 - wmaps) .* tmpmp;
                    if mdim == 4
                        [artv.FWHMResEst, artv.FWHMResImg] = ...
                            ne_methods.resestsmooth(tmpmp - ptc, bc.Resolution);
                    else
                        [artv.FWHMResEst, artv.FWHMResImg] = ...
                            ne_methods.resestsmoothsrf(tmpmp - ptc, opts.srf);
                    end
                end
            end
            tmap(isinf(tmap) | isnan(tmap)) = 0;

            % set name and map data
            mapc.Map(occ).Name = sprintf('%s%s (without average res.)', ...
                opts.names{cc}, tmapr);
            mapc.Map(occ).(mfld) = single(tmap);
            mapc.Map(occ).DF1 = nval - 1;
            mapc.Map(occ).LowerThreshold = -sdist('tinv', 0.005, nval - 1);
            mapc.Map(occ).UpperThreshold = -sdist('tinv', 0.0001, nval - 1);
            mapc.Map(occ).RunTimeVars = artv;
            occ = occ + 1;

            % add smoothness estimate map
            if opts.estsmap
                mapc.Map(occ).Name = sprintf('%s (smoothness estimate w/o avg res.)', opts.names{cc});
                mapc.Map(occ).Type = 83;
                mapc.Map(occ).(mfld) = single(artv.FWHMResImg);
                if mdim == 4
                    mapc.Map(occ).LowerThreshold = 2 * bc.Resolution + 0.5 * opts.smk;
                    mapc.Map(occ).UpperThreshold = 6 * bc.Resolution + 2.5 * opts.smk;
                else
                    mapc.Map(occ).LowerThreshold = 5;
                    mapc.Map(occ).UpperThreshold = 20;
                end
                mapc.Map(occ).ShowPositiveNegativeFlag = 1;
                occ = occ + 1;
            end

            % compute contrast for actual average residual
            if ~opts.robust
                tmap = glmtstat([0, 1], bmaps, ixx, se);
            else
                tmap = robustt(x, tmpmp, bmaps, wmaps, [0, 1]);
            end
            tmap(isinf(tmap) | isnan(tmap)) = 0;

            % set name and map data
            mapc.Map(occ).Name = sprintf('%s%s (corr. average res.)', ...
                opts.names{cc}, tmapr);
            mapc.Map(occ).(mfld) = single(tmap);
            mapc.Map(occ).DF1 = nval - 1;
            mapc.Map(occ).LowerThreshold = -sdist('tinv', 0.005, nval - 1);
            mapc.Map(occ).UpperThreshold = -sdist('tinv', 0.0001, nval - 1);
            mapc.Map(occ).RunTimeVars = artv;
            occ = occ + 1;

            % compute contrast vs. correlation with residual map
            if ~opts.robust
                tmap = ne_methods.conjval(glmtstat([1, 0.5], bmaps, ixx, se), ...
                    glmtstat([1, -0.5], bmaps, ixx, se));
            else
                tmap = ne_methods.conjval(robustt(x, tmpmp, bmaps, wmaps, [1, 0.5]), ...
                    robustt(x, tmpmp, bmaps, wmaps, [1, -0.5]));
            end
            tmap(isinf(tmap) | isnan(tmap)) = 0;

            % set name and map data
            mapc.Map(occ).Name = sprintf('%s%s (w/o+xmsk average res.)', ...
                opts.names{cc}, tmapr);
            mapc.Map(occ).(mfld) = single(tmap);
            mapc.Map(occ).DF1 = nval - 1;
            mapc.Map(occ).LowerThreshold = -sdist('tinv', 0.005, nval - 1);
            mapc.Map(occ).UpperThreshold = -sdist('tinv', 0.0001, nval - 1);
            mapc.Map(occ).RunTimeVars = artv;
            occ = occ + 1;

            % add robust-summary information
            artv.FWHMResEst = [];
            artv.FWHMResImg = [];
            artv.GlobSigMap = [];
            if opts.robust && opts.robwmaps
                mapc.Map(occ).Name = sprintf('%s%s (mean robust weight, w/o global mean)', ...
                    opts.names{cc}, tmapr);
                mapc.Map(occ).Type = 145;
                mapc.Map(occ).LowerThreshold = 0.75;
                mapc.Map(occ).UpperThreshold = 1;
                mapc.Map(occ).(mfld) = single((1 / size(wmaps, mdim)) .* sum(wmaps, mdim));
                mapc.Map(occ).RunTimeVars = artv;
                occ = occ + 1;
                mapc.Map(occ).Name = sprintf('%s%s (weighted number of outliers, w/o global mean)', ...
                    opts.names{cc}, tmapr);
                mapc.Map(occ).Type = 146;
                mapc.Map(occ).LowerThreshold = 1;
                mapc.Map(occ).UpperThreshold = 0.5 * size(wmaps, mdim);
                mapc.Map(occ).(mfld) = single(-3 .* sum( ...
                    limitrangec(wmaps - 2/3, -1/3, 0, -1/6), mdim));
                mapc.Map(occ).RunTimeVars = artv;
                occ = occ + 1;
            end

            % add robust weight maps to VMP
            if opts.robust && opts.swmaps
                for bmc = 1:size(wmaps, mdim)
                    mapc.Map(occ).Name = sprintf('%s%s (%s outlier, w/o global mean)', ...
                        opts.names{cc}, tmapr, ffxsubs{gax(bmc)});
                    mapc.Map(occ).Type = 144;
                    mapc.Map(occ).LowerThreshold = 0.25;
                    mapc.Map(occ).UpperThreshold = 1;
                    if mdim == 4
                        mapc.Map(occ).(mfld) = single(1 - wmaps(:, :, :, bmc));
                    else
                        mapc.Map(occ).(mfld) = single(1 - wmaps(:, bmc));
                    end
                    mapc.Map(occ).RunTimeVars = artv;
                    occ = occ + 1;
                end
            end
        end

    % for multiple groups
    else

        % create additional map without residual
        if opts.meanr

            % compute std==1-scaled data
            resmp = tmpmp;
            gsmap = 1 ./ sqrt(varc(tmpmp, mdim));
            gsmap(isinf(gsmap) | isnan(gsmap) | gsmap > 1) = 1;
            for sc = 1:goodsubs
                if mdim == 4
                    resmp(:, :, :, sc) = resmp(:, :, :, sc) .* gsmap;
                else
                    resmp(:, sc) = resmp(:, sc) .* gsmap;
                end
            end
            resmp(isinf(resmp) | isnan(resmp)) = 0;
            artv.GlobSigMap = single(gsmap);
        end

        if ~opts.robust
            for gc1 = 1:ngrp
                ga1 = find(ga(gaxk) == gc1);
                if mdim == 4
                    m1 = (1 / numel(ga1)) .* sum(tmpmp(:, :, :, ga1), mdim);
                    s1 = (1 / numel(ga1)) .* var(tmpmp(:, :, :, ga1), [], mdim);
                else
                    m1 = (1 / numel(ga1)) .* sum(tmpmp(:, ga1), mdim);
                    s1 = (1 / numel(ga1)) .* var(tmpmp(:, ga1), [], mdim);
                end
                for gc2 = (gc1 + 1):ngrp
                    ga2 = find(ga(gaxk) == gc2);
                    if mdim == 4
                        m2 = (1 / numel(ga2)) .* sum(tmpmp(:, :, :, ga2), mdim);
                        s2 = (1 / numel(ga2)) .* var(tmpmp(:, :, :, ga2), [], mdim);
                    else
                        m2 = (1 / numel(ga2)) .* sum(tmpmp(:, ga2), mdim);
                        s2 = (1 / numel(ga2)) .* var(tmpmp(:, ga2), [], mdim);
                    end
                    t12 = (m1 - m2) ./ sqrt(s1 + s2);
                    df12 = ((s1 + s2) .^ 2) ./ ...
                        ((1 / (numel(ga1) - 1)) .* s1 .* s1 + ...
                         (1 / (numel(ga2) - 1)) .* s2 .* s2);
                    badv = find(isinf(t12) | isnan(t12) | isnan(df12) | df12 < 1);
                    t12(badv) = 0;
                    df12(badv) = 1;
                    rdf12 = sum(ga == gc1 | ga == gc2) - 2;
                    mapc.Map(occ).Name = sprintf('%s (%s > %s)', ...
                        opts.names{cc}, opts.groups{gc1, 1}, ...
                        opts.groups{gc2, 1});
                    mapc.Map(occ).(mfld) = single(...
                        sdist('tinv', sdist('tcdf', t12, df12), rdf12));
                    mapc.Map(occ).DF1 = rdf12;
                    mapc.Map(occ).RunTimeVars = artv;
                    mapc.Map(occ).RunTimeVars.SubSel(ga(gaxk) ~= gc1 & ga(gaxk) ~= gc2) = [];
                    if opts.estfwhm
                        gresmp = tmpmp;
                        if mdim == 4
                            gresmp(:, :, :, ga1) = gresmp(:, :, :, ga1) - ...
                                m1(:, :, :, ones(1, numel(ga1)));
                            gresmp(:, :, :, ga2) = gresmp(:, :, :, ga2) - ...
                                m2(:, :, :, ones(1, numel(ga2)));
                        else
                            gresmp(:, ga1) = gresmp(:, ga1) - m1(:, ones(1, numel(ga1)));
                            gresmp(:, ga2) = gresmp(:, ga2) - m2(:, ones(1, numel(ga2)));
                        end
                        if mdim == 4
                            [mapc.Map(occ).RunTimeVars.FWHMResEst, ...
                             mapc.Map(occ).RunTimeVars.FWHMResImg] = ...
                                ne_methods.resestsmooth(gresmp(:, :, :, [ga1(:)', ga2(:)']), bc.Resolution);
                        else
                            [mapc.Map(occ).RunTimeVars.FWHMResEst, ...
                             mapc.Map(occ).RunTimeVars.FWHMResImg] = ...
                                ne_methods.resestsmoothsrf(gresmp(:, [ga1(:)', ga2(:)']), opts.srf);
                        end
                    end
                    occ = occ + 1;

                    % remove residual and re-create maps as well
                    if opts.meanr
                        tmin = -sdist('tinv', 0.25, rdf12);
                        tmax = -sdist('tinv', 0.05, rdf12);
                        wmp = limitrangec(tmax - (1 / (tmax - tmin)) .* ...
                            double(mapc.Map(occ-1).(mfld)), 0.25, 1, 0);
                        wmp = wmp .* wmp;
                        if isequal(size(wmp), size(opts.meanrmsk))
                            wmp = wmp .* opts.meanrmsk;
                        end
                        wmp(mapc.Map(occ-1).(mfld) == 0) = 0;
                        gresmp = resmp;
                        for sc = [ga1(:)', ga2(:)']
                            if mdim == 4
                                gresmp(:, :, :, sc) = gresmp(:, :, :, sc) .* wmp;
                            else
                                gresmp(:, sc) = gresmp(:, sc) .* wmp;
                            end
                        end
                        if mdim == 4
                            gresmp = lsqueeze(sum(sum(sum(gresmp, 1), 2), 3)) ./ sum(wmp(:));
                        else
                            gresmp = lsqueeze(sum(gresmp, 1)) ./ sum(wmp(:));
                        end
                        x = zeros(numel(ga1) + numel(ga2), 4);
                        x(1:numel(ga1), 1:2) = [ones(numel(ga1), 1), ztrans(gresmp(ga1))];
                        x(numel(ga1)+1:end, 3:4) = [ones(numel(ga2), 1), ztrans(gresmp(ga2))];
                        if mdim == 4
                            bmaps = calcbetas(x, tmpmp(:, :, :, [ga1(:)', ga2(:)']), mdim);
                            bmaps(:, :, :, [1, 3]) = 0;
                            gresmp = tmpmp(:, :, :, [ga1(:)', ga2(:)']) - ...
                                reshape((reshape(bmaps, numel(df12), 4) * x'), ...
                                [size(df12), size(x, 1)]);
                            m1rm = (1 / numel(ga1)) .* sum(gresmp(:, :, :, 1:numel(ga1)), mdim);
                            s1rm = (1 / numel(ga1)) .* var(gresmp(:, :, :, 1:numel(ga1)), [], mdim);
                            m2rm = (1 / numel(ga2)) .* sum(gresmp(:, :, :, numel(ga1)+1:end), mdim);
                            s2rm = (1 / numel(ga2)) .* var(gresmp(:, :, :, numel(ga1)+1:end), [], mdim);
                        else
                            bmaps = calcbetas(x, tmpmp(:, [ga1(:)', ga2(:)']), mdim);
                            bmaps(:, [1, 3]) = 0;
                            gresmp = tmpmp(:, [ga1(:)', ga2(:)']) - ...
                                reshape((reshape(bmaps, numel(df12), 4) * x'), ...
                                [numel(df12), size(x, 1)]);
                            m1rm = (1 / numel(ga1)) .* sum(gresmp(:, 1:numel(ga1)), mdim);
                            s1rm = (1 / numel(ga1)) .* var(gresmp(:, 1:numel(ga1)), [], mdim);
                            m2rm = (1 / numel(ga2)) .* sum(gresmp(:, numel(ga1)+1:end), mdim);
                            s2rm = (1 / numel(ga2)) .* var(gresmp(:, numel(ga1)+1:end), [], mdim);
                        end
                        t12 = (m1rm - m2rm) ./ sqrt(s1rm + s2rm);
                        df12 = ((s1rm + s2rm) .^ 2) ./ ...
                            ((1 / (numel(ga1) - 1)) .* s1rm .* s1rm + ...
                             (1 / (numel(ga2) - 1)) .* s2rm .* s2rm);
                        badv = find(isinf(t12) | isnan(t12) | isnan(df12) | df12 < 1);
                        t12(badv) = 0;
                        df12(badv) = 1;
                        rdf12 = sum(ga == gc1 | ga == gc2) - 4;
                        mapc.Map(occ).Name = sprintf('%s (%s > %s, without average res.)', ...
                            opts.names{cc}, opts.groups{gc1, 1}, ...
                            opts.groups{gc2, 1});
                        mapc.Map(occ).(mfld) = single(...
                            sdist('tinv', sdist('tcdf', t12, df12), rdf12));
                        mapc.Map(occ).DF1 = rdf12;
                        mapc.Map(occ).RunTimeVars = artv;
                        mapc.Map(occ).RunTimeVars.SubSel(ga(gaxk) ~= gc1 & ga(gaxk) ~= gc2) = [];
                        if opts.estfwhm && mdim == 4
                            gresmp(:, :, :, 1:numel(ga1)) = gresmp(:, :, :, 1:numel(ga1)) - ...
                                m1rm(:, :, :, ones(1, numel(ga1)));
                            gresmp(:, :, :, numel(ga1)+1:end) = gresmp(:, :, :, numel(ga1)+1:end) - ...
                                m2rm(:, :, :, ones(1, numel(ga2)));
                            [mapc.Map(occ).RunTimeVars.FWHMResEst, ...
                             mapc.Map(occ).RunTimeVars.FWHMResImg] = ...
                                ne_methods.resestsmooth(gresmp, bc.Resolution);
                        elseif opts.estfwhm && mdim == 2
                            gresmp(:, 1:numel(ga1)) = gresmp(:, 1:numel(ga1)) - ...
                                m1rm(:, ones(1, numel(ga1)));
                            gresmp(:, numel(ga1)+1:end) = gresmp(:, numel(ga1)+1:end) - ...
                                m2rm(:, ones(1, numel(ga2)));
                            [mapc.Map(occ).RunTimeVars.FWHMResEst, ...
                             mapc.Map(occ).RunTimeVars.FWHMResImg] = ...
                                ne_methods.resestsmoothsrf(gresmp, opts.srf);
                        end
                        occ = occ + 1;
                    end
                end
            end

        % robust multi-group regression
        else
            [tmap, wmaps, noot, bmaps] = ...
                ne_methods.robustnsamplet_img(tmpmp, ga(gaxk), struct('wvols', true));
            if opts.estfwhm
                ptc = zeros(size(tmpmp));
                for bmc = 1:size(bmaps, mdim)
                    if mdim == 4
                        ptc(:, :, :, ga(gaxk) == bmc) = ...
                            repmat(bmaps(:, :, :, bmc), [mones, sum(ga(gaxk) == bmc)]);
                    else
                        ptc(:, ga(gaxk) == bmc) = ...
                            repmat(bmaps(:, bmc), [mones, sum(ga(gaxk) == bmc)]);
                    end
                end
                ptc = wmaps .* ptc + (1 - wmaps) .* tmpmp;
            end
            for gmc = 1:gamx
                [g2, g1] = find(gag == gmc);
                mapc.Map(occ).Name = sprintf('%s%s (%s > %s)', ...
                    opts.names{cc}, tmapr, opts.groups{g1, 1}, ...
                    opts.groups{g2, 1});
                if mdim == 4
                    mapc.Map(occ).(mfld) = single(tmap(:, :, :, gmc));
                else
                    mapc.Map(occ).(mfld) = single(tmap(:, gmc));
                end
                mapc.Map(occ).DF1 = sum(ga == g1 | ga == g2) - 2;
                mapc.Map(occ).RunTimeVars = artv;
                mapc.Map(occ).RunTimeVars.SubSel(ga(gaxk) ~= g1 & ga(gaxk) ~= g2) = [];
                if opts.estfwhm && mdim == 4
                    [mapc.Map(occ).RunTimeVars.FWHMResEst, ...
                     mapc.Map(occ).RunTimeVars.FWHMResImg] = ...
                     ne_methods.resestsmooth(tmpmp - ptc, bc.Resolution);
                elseif opts.estfwhm && mdim == 2
                    [mapc.Map(occ).RunTimeVars.FWHMResEst, ...
                     mapc.Map(occ).RunTimeVars.FWHMResImg] = ...
                     ne_methods.resestsmoothsrf(tmpmp - ptc, opts.srf);
                end
                occ = occ + 1;

                % leave one map empty for the mean-removed maps
                if opts.meanr
                    occ = occ + 1;
                end
            end

            % save counter
            socc = occ;

            % add robust-summary information
            if opts.robust && opts.robwmaps
                artv.FWHMResEst = [];
                artv.FWHMResImg = [];
                mapc.Map(occ).Name = sprintf('%s%s (mean robust weight)', ...
                    opts.names{cc}, tmapr);
                mapc.Map(occ).Type = 145;
                mapc.Map(occ).LowerThreshold = 0.75;
                mapc.Map(occ).UpperThreshold = 1;
                mapc.Map(occ).(mfld) = single((1 / size(wmaps, mdim)) .* sum(wmaps, mdim));
                mapc.Map(occ).RunTimeVars = artv;
                occ = occ + 1;
                mapc.Map(occ).Name = sprintf('%s%s (weighted number of outliers)', ...
                    opts.names{cc}, tmapr);
                mapc.Map(occ).Type = 146;
                mapc.Map(occ).LowerThreshold = 1;
                mapc.Map(occ).UpperThreshold = 0.5 * size(wmaps, mdim);
                mapc.Map(occ).(mfld) = single(-3 .* sum( ...
                    limitrangec(wmaps - 2/3, -1/3, 0, -1/6), mdim));
                mapc.Map(occ).RunTimeVars = artv;
                occ = occ + 1;
            end

            % add robust weight maps to VMP
            if opts.robust && opts.swmaps
                for bmc = 1:size(wmaps, mdim)
                    mapc.Map(occ).Name = sprintf('%s%s (%s outlier)', ...
                        opts.names{cc}, tmapr, ffxsubs{gax(bmc)});
                    mapc.Map(occ).Type = 144;
                    mapc.Map(occ).LowerThreshold = 0.25;
                    mapc.Map(occ).UpperThreshold = 1;
                    if mdim == 4
                        mapc.Map(occ).(mfld) = single(1 - wmaps(:, :, :, bmc));
                    else
                        mapc.Map(occ).(mfld) = single(1 - wmaps(:, bmc));
                    end
                    mapc.Map(occ).RunTimeVars = artv;
                    occ = occ + 1;
                end
            end
            socc2 = occ;

            % for mean-removed maps
            if opts.meanr
                occ = socc - 2 * gamx(end);
                wmp = ones(size(mapc.Map(occ).(mfld)));
                for gmc = 1:gamx
                    [g2, g1] = find(gag == gmc);
                    rdf12 = sum(ga == g1 | ga == g2) - 4;
                    tmin = -sdist('tinv', 0.25, rdf12);
                    tmax = -sdist('tinv', 0.05, rdf12);
                    wmp = min(wmp, limitrangec(tmax - (1 / (tmax - tmin)) .* ...
                        double(mapc.Map(occ + 2 * gmc - 2).(mfld)), 0.25, 1, 0));
                    wmp(mapc.Map(occ + 2 * gmc - 2).(mfld) == 0) = 0;
                end
                wmp = wmp .* wmp;
                if isequal(size(wmp), size(opts.meanrmsk))
                    wmp = wmp .* opts.meanrmsk;
                end
                gresmp = resmp;
                for sc = 1:goodsubs
                    if mdim == 4
                        gresmp(:, :, :, sc) = gresmp(:, :, :, sc) .* wmp;
                    else
                        gresmp(:, sc) = gresmp(:, sc) .* wmp;
                    end
                end
                artv.GlobSigMap = single(wmp);
                if mdim == 4
                    gresmp = lsqueeze(sum(sum(sum(gresmp, 1), 2), 3)) ./ sum(wmp(:));
                else
                    gresmp = lsqueeze(sum(gresmp, 1)) ./ sum(wmp(:));
                end
                x = zeros(size(gresmp, 4), 2 * ngrp);
                gs = zeros(1, ngrp);
                for gmc = 1:ngrp
                    gs(gmc) = sum(ga(gaxk) == gmc);
                    x(ga(gaxk) == gmc, 2*gmc-1:2*gmc) = ...
                        [ones(gs(gmc), 1), ztrans(gresmp(ga == gmc))];
                end
                [bmaps, wmaps] = ne_methods.fitrobustbisquare_img(x, tmpmp);
                if opts.estfwhm
                    ptc = zeros(size(tmpmp));
                    for bmc = 1:size(x, 2)
                        if mdim == 4
                            ptc = ptc + repmat(bmaps(:, :, :, bmc), [mones, size(x, 1)]) .* ...
                                repmat(reshape(x(:, bmc), [mones, size(x, 1)]), szmapo);
                        else
                            ptc = ptc + repmat(bmaps(:, bmc), [mones, size(x, 1)]) .* ...
                                repmat(reshape(x(:, bmc), [mones, size(x, 1)]), szmapo);
                        end
                    end
                    ptc = wmaps .* ptc + (1 - wmaps) .* tmpmp;
                end
                wmaps = sqrt(wmaps);
                rmaps = tmpmp - reshape((reshape(bmaps, numel(wmp), size(bmaps, mdim)) * x'), ...
                    [size(wmp), size(tmpmp, mdim)]);
                for g1c = 1:ngrp
                    if mdim == 4
                        g1w = wmaps(:, :, :, ga(gaxk) == g1c);
                    else
                        g1w = wmaps(:, ga(gaxk) == g1c);
                    end
                    g1ws = sum(g1w .* g1w, mdim);
                    g1w = g1w .* repmat((gs(g1c) ./ g1ws), [mones, size(g1w, mdim)]);
                    if mdim == 4
                        g1r = rmaps(:, :, :, ga(gaxk) == g1c) .* g1w;
                    else
                        g1r = rmaps(:, ga(gaxk) == g1c) .* g1w;
                    end
                    g1rv = varc(g1r, mdim, true) .* (gs(g1c) ./ (g1ws .* g1ws));
                    for g2c = (g1c+1):ngrp
                        if mdim == 4
                            g2w = wmaps(:, :, :, ga(gaxk) == g2c);
                        else
                            g2w = wmaps(:, ga(gaxk) == g2c);
                        end
                        g2ws = sum(g2w .* g2w, mdim);
                        g2w = g2w .* repmat((gs(g2c) ./ g2ws), [mones, size(g2w, mdim)]);
                        if mdim == 4
                            g2r = rmaps(:, :, :, ga(gaxk) == g2c) .* g2w;
                        else
                            g2r = rmaps(:, ga(gaxk) == g2c) .* g2w;
                        end
                        g2rv = varc(g2r, mdim, true) .* (gs(g2c) ./ (g2ws .* g2ws));
                        if mdim == 4
                            tstat = (bmaps(:, :, :, 2*g1c-1) - ...
                                bmaps(:, :, :, 2*g2c-1)) ./ sqrt(g1rv + g2rv);
                        else
                            tstat = (bmaps(:, 2*g1c-1) - ...
                                bmaps(:, 2*g2c-1)) ./ sqrt(g1rv + g2rv);
                        end
                        dfstat = ((g1rv + g2rv) .^ 2) ./ ...
                            (g1rv .* g1rv ./ (g1ws - 1) + g2rv .* g2rv ./ (g2ws - 1));
                        rdf12 = sum(ga == g1c | ga == g2c) - 4;
                        badv = (isinf(tstat) | isnan(tstat) | isnan(dfstat) | dfstat < 1);
                        if any(badv(:))
                            tstat(badv) = 0;
                            dfstat(badv) = 1;
                        end

                        % don't override non-mean removed maps
                        occ = occ + 1;

                        % then store correct map
                        mapc.Map(occ).Name = ...
                            sprintf('%s (%s > %s, robust, without average res.)', ...
                            opts.names{cc}, opts.groups{g1c, 1}, opts.groups{g2c, 1});
                        mapc.Map(occ).(mfld) = single(...
                            sdist('tinv', sdist('tcdf', tstat, dfstat), rdf12));
                        mapc.Map(occ).DF1 = rdf12;
                        mapc.Map(occ).RunTimeVars = artv;
                        mapc.Map(occ).RunTimeVars.SubSel(ga(gaxk) ~= g1c & ga(gaxk) ~= g2c) = [];
                        if opts.estfwhm && mdim == 4
                            [mapc.Map(occ).RunTimeVars.FWHMResEst, ...
                                mapc.Map(occ).RunTimeVars.FWHMResImg] = ...
                                ne_methods.resestsmooth(tmpmp(:, :, :, ga(gaxk) == g1c | ga(gaxk) == g2c) - ...
                                ptc(:, :, :, ga(gaxk) == g1c | ga(gaxk) == g2c), bc.Resolution);
                        elseif opts.estfwhm && mdim == 4
                            [mapc.Map(occ).RunTimeVars.FWHMResEst, ...
                                mapc.Map(occ).RunTimeVars.FWHMResImg] = ...
                                ne_methods.resestsmoothsrf(tmpmp(:, ga(gaxk) == g1c | ga(gaxk) == g2c) - ...
                                ptc(:, ga(gaxk) == g1c | ga(gaxk) == g2c), opts.srf);
                        end
                        occ = occ + 1;
                    end
                end
                wmaps = wmaps .* wmaps;
            end

            % add robust-summary information
            occ = socc2;
            artv.FWHMResEst = [];
            artv.FWHMResImg = [];
            artv.GlobSigMap = [];
            if opts.robust && opts.robwmaps
                mapc.Map(occ).Name = sprintf('%s%s (mean robust weight, w/o global mean)', ...
                    opts.names{cc}, tmapr);
                mapc.Map(occ).Type = 145;
                mapc.Map(occ).LowerThreshold = 0.75;
                mapc.Map(occ).UpperThreshold = 1;
                mapc.Map(occ).(mfld) = single((1 / size(wmaps, mdim)) .* sum(wmaps, mdim));
                mapc.Map(occ).RunTimeVars = artv;
                occ = occ + 1;
                mapc.Map(occ).Name = sprintf('%s%s (weighted number of outliers, w/o global mean)', ...
                    opts.names{cc}, tmapr);
                mapc.Map(occ).Type = 146;
                mapc.Map(occ).LowerThreshold = 1;
                mapc.Map(occ).UpperThreshold = 0.5 * size(wmaps, mdim);
                mapc.Map(occ).(mfld) = single(-3 .* sum( ...
                    limitrangec(wmaps - 2/3, -1/3, 0, -1/6), mdim));
                mapc.Map(occ).RunTimeVars = artv;
                occ = occ + 1;
            end

            % add robust weight maps to VMP
            if opts.robust && opts.swmaps
                for bmc = 1:size(wmaps, mdim)
                    mapc.Map(occ).Name = sprintf('%s%s (%s outlier, w/o global mean)', ...
                        opts.names{cc}, tmapr, ffxsubs{gax(bmc)});
                    mapc.Map(occ).Type = 144;
                    mapc.Map(occ).LowerThreshold = 0.25;
                    mapc.Map(occ).UpperThreshold = 1;
                    if mdim == 4
                        mapc.Map(occ).(mfld) = single(1 - wmaps(:, :, :, bmc));
                    else
                        mapc.Map(occ).(mfld) = single(1 - wmaps(:, bmc));
                    end
                    mapc.Map(occ).RunTimeVars = artv;
                    occ = occ + 1;
                end
            end
        end
    end
end
if occ <= numel(mapc.Map)
    mapc.Map(occ:end) = [];
end

% put back
mapc.NrOfMaps = numel(mapc.Map);
map.C = mapc;

% clear residual VTC
delete(resvtc);
