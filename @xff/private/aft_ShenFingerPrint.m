function [p, ptc] = aft_ShenFingerPrint(xo, opts)
% AFT::ShenFingerPrint  - compute Shen et al. based finger print for data
%
% FORMAT:       [p, ptc] = vtc.FingerPrint([opts])
%
% Input fields:
%
%       opts        settings
%        .conds     for GLM objects: select conditions over which to run
%        .glmcdim   GLM-based correlation dim, either {'cond'} or 'subj'
%        .networks  instead of the 268x268 matrix for full parcellation,
%                   return a 8x8 matrix (networks, default: false)
%        .robcorr   use robcorrcoef (instead of corrcoef, false)
%        .srf       for MTC-based files, surface with coordinates
%        .subsel    for GLM objects: subjects to extract data from (all)
%
% Output fields:
%
%       p           268x268 or 8x8 connectivity finger print(s)
%       ptc         NrOfDataPointsx268 or NrOfDataPointsx8 matrices
%
% TYPES: GLM, HEAD, HDR, MTC, VTC

% Version:  v1.1
% Build:    16061611
% Date:     Jun-16 2016, 11:33 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% NeuroElf library
global ne_methods;
multimatch = ne_methods.multimatch;
robcorrcoef = ne_methods.robcorrcoef;
splittocellc = ne_methods.splittocellc;

% persistent fingerprint NII file
persistent fpnii;
if numel(fpnii) ~= 1 || ~xffisobject(fpnii, true, 'hdr')
    try

        % use temp object, so the un-gziped file won't survive in /tmp
        fpnii_o = xff([neuroelf_path('files') filesep 'shenparcel' filesep 'shen_1mm_268_parcellation.nii.gz']);
        aft_LoadTransIOData(fpnii_o);
        fpnii = aft_CopyObject(fpnii_o);
        aft_ClearObject(fpnii_o);

        % also load network labels
        fpnlo = xff([neuroelf_path('files') filesep 'shenparcel' filesep 'shen_268_parcellation_networklabels.csv'], 'ntt');
        fpnii.H.NetworkLabels = fpnlo.C.Data;
        aft_ClearObject(fpnlo);
    catch xfferror
        error('neuroelf:xff:fileLoadError', 'Error loading finger-print dataset: %s', xfferror.message);
    end

    % init cache for sparse spatial averaging matrices for repeated access
    fpnii.H.FPSampling = {'NONE', [], [], [], [], []};
end

% check input
if nargin < 1 || ~xffisobject(xo, true, {'glm', 'hdr', 'head', 'mtc', 'vtc'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
xotype = lower(xo.S.Extensions{1});
bc = xo.C;
subids = {'FFX'};
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if strcmp(xotype, 'glm')
    spreds = glm_SubjectPredictors(xo);
    if strcmpi(spreds{end}, 'constant')
        spreds(end) = [];
    end
    subids = glm_Subjects(xo);
else
    spreds = {};
    subids = {};
end
if ~isfield(opts, 'conds') || ~iscell(opts.conds) || isempty(opts.conds) || ...
   ~all(cellfun(@ischar, opts.conds(:)))
    opts.conds = spreds;
else
    opts.conds = opts.conds(:);
    opts.conds(cellfun('isempty', opts.conds)) = [];
end
if ~isfield(opts, 'glmcdim') || ~ischar(opts.glmcdim) || isempty(opts.glmcdim) || ...
   ~any(lower(opts.glmcdim(1)) == 'cs')
    opts.glmcdim = 'c';
else
    opts.glmcdim = lower(opts.glmcdim(1));
end
if strcmp(xotype, 'glm')
    if isequal(opts.conds, spreds)
        condi = (1:numel(spreds))';
    else
        condi = multimatch(opts.conds, spreds);
        if all(condi < 1)
            condi = find(multimatch(spreds, opts.conds, true) > 0);
        end
    end
    if opts.glmcdim == 'c'
        if any(condi < 1)
            error('neuroelf:xff:badArgument', 'Invalid condition selection.');
        end
    elseif any(condi < 1) || isempty(condi)
        condi = zeros(numel(spreds), numel(opts.conds));
        for cc = 1:numel(opts.conds)
            condt = opts.conds{cc};
            if any(condt == '>')
                condt = splittocellc(condt, ' > ');
                condn = condt{2};
                condt = condt{1};
            else
                condn = '';
            end
            condt = splittocellc(condt, ' + ');
            condti = multimatch(condt, spreds);
            if all(condti < 0)
                condti = find(multimatch(spreds, condt, true) > 0);
            elseif any(condti < 0)
                error('neuroelf:xff:badArgument', 'Invalid condition selection.');
            end
            condi(condti, cc) = 1 ./ numel(condti);
            if ~isempty(condn)
                condt = splittocellc(condn, ' + ');
                condti = multimatch(condt, spreds);
                if all(condti < 0)
                    condti = find(multimatch(spreds, condt, true) > 0);
                elseif any(condti < 0)
                    error('neuroelf:xff:badArgument', 'Invalid condition selection.');
                end
                condi(condti, cc) = -1 ./ numel(condti);
            end
        end
    end
end
if ~isfield(opts, 'networks') || ~islogical(opts.networks) || numel(opts.networks) ~= 1
    opts.networks = false;
end
if ~isfield(opts, 'robcorr') || ~islogical(opts.robcorr) || numel(opts.robcorr) ~= 1
    opts.robcorr = false;
end
if ~isfield(opts, 'srf') || numel(opts.srf) ~= 1 || ~xffisobject(opts.srf, true ,'srf')
    opts.srf = [];
end
if strcmp(xotype, 'glm')
    if ~isfield(opts, 'subsel') || isempty(opts.subsel) || ...
       (~isa(opts.subsel, 'double') && ~iscell(opts.subsel))
        opts.subsel = (1:numel(subids))';
    end
    opts.subsel = opts.subsel(:);
    if iscell(opts.subsel)
        if ~all(cellfun(@ischar, opts.subsel)) || any(cellfun('isempty', opts.subsel))
            error('neuroelf:xff:badArgument', 'Invalid subject selection.');
        end
        opts.subsel = multimatch(opts.subsel, subids);
        if any(opts.subsel == 0)
            error('neuroelf:xff:badArgument', 'Invalid subject selection.');
        end
    elseif any(isinf(opts.subsel) | isnan(opts.subsel) | opts.subsel < 1 | opts.subsel > numel(subids))
        error('neuroelf:xff:badArgument', 'Invalid subject selection.');
    elseif ~isempty(opts.subsel)
        opts.subsel = unique(round(opts.subsel));
    end
else
    opts.subsel = [];
end

% bounding boxes mismatch
if strcmp(xotype, 'vtc') || (strcmp(xotype, 'glm') && bc.ProjectType == 1)
    bb = aft_BoundingBox(xo);
else
    bb = aft_Layout(xo);
end
if ~isequal(fpnii.H.FPSampling{1}, bb)

    % sample fpnii for VTC/GLM
    if isstruct(bb)
        fpvtc = round(aft_SampleBVBox(fpnii, bb, 1, 'nearest'));

    % for HDR
    elseif strcmp(xotype, 'hdr')
        cfr = hdr_CoordinateFrame(xo);
        vsz = [size(bc.VoxelData), 1, 1, 1];
        [cx, cy, cz] = ndgrid(1:vsz(1), 1:vsz(2), 1:vsz(3));
        cx = [cx(:), cy(:), cz(:), ones(numel(cx), 1)];
        cx = cx * cfr.Trf';
        fpvtc = round(aft_SampleData3D(fpnii, cx(:, 1:3), struct('mapvol', 1, 'method', 'nearest')));

    % for HEAD
    elseif strcmp(xotype, 'head')
        cfr = head_CoordinateFrame(xo);
        vsz = [size(bc.Brick(1).Data), 1, 1, 1];
        [cx, cy, cz] = ndgrid(1:vsz(1), 1:vsz(2), 1:vsz(3));
        cx = [cx(:), cy(:), cz(:), ones(numel(cx), 1)];
        cx = cx * cfr.Trf';
        fpvtc = round(aft_SampleData3D(fpnii, cx(:, 1:3), struct('mapvol', 1, 'method', 'nearest')));

    % not for MTC-based objects
    elseif isempty(opts.srf) || bc.NrOfVertices ~= size(opts.srf.C.VertexCoordinate, 1)
        error('neuroelf:xff:missingOrBadSurface', 'Missing or bad surface supplied.');

    % for MTC-based objects
    else
        fpvtc = round(aft_SampleData3D(128 - opts.srf.C.VertexCoordinate(:, [3, 1, 2]), ...
            struct('mapvol', 1, 'method', 'nearest')));
    end
    fpvtc = double(fpvtc(:));
    fpi = find(fpvtc > 0);
    hpvtc = histc(fpvtc, 1:268);
    nlab = fpnii.H.NetworkLabels(:, 2);
    hpntc = histc(nlab(fpvtc(fpi)), 1:8);
    nsvtc = numel(fpi);

    % create sparse matrices
    spfp = sparse(fpi, fpvtc(fpi), true(nsvtc, 1), numel(fpvtc), 268, nsvtc);
    snfp = sparse(fpi, nlab(fpvtc(fpi)), true(nsvtc, 1), numel(fpvtc), 8, nsvtc);
    wpfp = sparse((1:268)', (1:268)', 1 ./ hpvtc, 268, 268, 268);
    wnfp = sparse((1:8)', (1:8)', 1 ./ hpntc, 8, 8, 8);

    % store data
    fpnii.H.FPSampling = {bb, fpvtc, spfp, snfp, wpfp, wnfp};
end

% transpose TC data
tpflag = false;

% get data (code for expression, etc. below) for VTC
if strcmp(xotype, 'vtc')
    vd = reshape(double(bc.VTCData), size(bc.VTCData, 1), prod(bb.DimXYZ));

% for MTC
elseif strcmp(xotype, 'mtc')
    vd = double(bc.MTCData);

% for HDR (Analyze/NIftI)
elseif strcmp(xotype, 'hdr')
    tpflag = true;
    if prod(vsz(1:4)) == prod(vsz)
        vd = double(bc.VoxelData);
    else
        vd = squeeze(double(bc.VoxelData(:, :, :, 1, :)));
    end
    vd = reshape(vd, prod(vsz(1:3)), vsz(4));
    if ~any(bc.ImgDim.ScalingSlope == [0, 1]) || bc.ImgDim.ScalingIntercept ~= 0
        if ~any(bc.ImgDim.ScalingSlope == [0, 1])
            if bc.ImgDim.ScalingIntercept ~= 0
                vd = bc.ImgDim.ScalingSlope .* vd + bc.ImgDim.ScalingIntercept;
            else
                vd = bc.ImgDim.ScalingSlope .* vd;
            end
        else
            vd = vd + single(bc.ImgDim.ScalingIntercept);
        end
    end

% for HEAD (AFNI)
elseif strcmp(xotype, 'head')
    nvox = prod(vsz(1:3));
    vd = zeros(numel(bc.Brick), nvox);
    for brc = 1:numel(bc.Brick)
        scf = bc.Brick(brc).ScalingFactor;
        if numel(scf) ~= 1 || isinf(scf) || isnan(scf) || scf == 0 || scf == 1
            vd(brc, :) = reshape(double(bc.Brick(brc).Data), 1, nvox);
        else
            vd(brc, :) = scf .* reshape(double(bc.Brick(brc).Data), 1, nvox);
        end
    end

% for GLMs
else

    % dimension is conditions
    nsubs = numel(opts.subsel);
    if any(condi(:) < 1)
        ucvec = true;
        nconds = size(condi, 2);
    else
        ucvec = false;
        nconds = 1;
    end

    % dimension is subjects
    if opts.glmcdim == 's'

        % create output
        np = cell(1, nconds);
        ntc = cell(1, nconds);
        p = cell(1, nconds);
        ptc = cell(1, nconds);

        % get data
        for cc = 1:nconds
            if ucvec
                cvec = condi(:, cc);
            else
                cvec = zeros(numel(spreds), 1);
                cvec(condi) = 1 ./ numel(cvec);
            end
            cdata = double(glm_RFX_conmaps(xo, cvec));
            if ndims(cdata) > 3
                cdsz = size(cdata);
                vd = reshape(cdata, [prod(cdsz(1:3)), cdsz(4)])';
            end

            % compute values
            ntc{cc} = (vd * fpnii.H.FPSampling{4}) * fpnii.H.FPSampling{6};
            if opts.robcorr
                np{cc} = robcorrcoef(ntc{cc});
            else
                np{cc} = corrcoef(ntc{cc});
            end
            np{cc}(1:9:end) = 0;
            ptc{cc} = (vd * fpnii.H.FPSampling{3}) * fpnii.H.FPSampling{5};
            if opts.robcorr
                p{cc} = robcorrcoef(ptc{cc});
            else
                p{cc} = corrcoef(ptc{cc});
            end
            p{cc}(1:269:end) = 0;
        end
        if numel(ntc) == 1
            ntc = ntc{1};
            np = np{1};
            ptc = ptc{1};
            p = p{1};
        end

    % dimension is conditions
    else

        % create output
        np = cell(nsubs, nconds);
        ntc = cell(nsubs, nconds);
        p = cell(nsubs, nconds);
        ptc = cell(nsubs, nconds);

        % iterate over subject
        ptrfx = (bc.ProjectTypeRFX > 0);
        for subc = 1:nsubs

            % for RFX
            if ptrfx
            else
            end
        end
    end

    % return
    xo.C.RunTimeVars.FingerPrintOpts = opts;
    xo.C.RunTimeVars.FingerPrintNetworksTC = ntc;
    xo.C.RunTimeVars.FingerPrintNetworks = np;
    xo.C.RunTimeVars.FingerPrintParcelsTC = ptc;
    xo.C.RunTimeVars.FingerPrintParcels = p;
    if opts.networks
        ptc = ntc;
        p = np;
    end
    return;
end

% compute (both, so cheap!)
if tpflag
    ntc = (fpnii.H.FPSampling{4}' * vd)' * fpnii.H.FPSampling{6};
else
    ntc = (vd * fpnii.H.FPSampling{4}) * fpnii.H.FPSampling{6};
end
xo.C.RunTimeVars.FingerPrintNetworksTC = ntc;
if opts.robcorr
    np = robcorrcoef(ntc);
else
    np = corrcoef(ntc);
end
np(1:9:end) = 0;
xo.C.RunTimeVars.FingerPrintNetworks = np;
if tpflag
    ptc = (fpnii.H.FPSampling{3}' * vd)' * fpnii.H.FPSampling{5};
else
    ptc = (vd * fpnii.H.FPSampling{3}) * fpnii.H.FPSampling{5};
end
xo.C.RunTimeVars.FingerPrintParcelsTC = ptc;
if opts.robcorr
    p = robcorrcoef(ptc);
else
    p = corrcoef(ptc);
end
p(1:269:end) = 0;
xo.C.RunTimeVars.FingerPrintParcels = p;

% what to return
if opts.networks
    ptc = ntc;
    p = np;
end
