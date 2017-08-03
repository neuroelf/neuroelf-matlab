function [vb, vbv, vbi] = glm_VOIBetas(xo, vfile, opts)
% GLM::VOIBetas  - returns a table of VOI betas (per subjects)
%
% FORMAT:       [vb, vbv, vbi] = glm.VOIBetas(voi [, opts]);
%
% Input fields:
%
%       voi         VOI object
%       opts        optional settings
%        .c         CxP contrasts (defaults: each beta map on its own)
%        .interp    either of {'nearest'}, 'linear', 'cubic'
%        .refspace  ReferenceSpace for VOI if not well defined ('bvi')
%        .rmean     remove mean (of map) first (default: false)
%        .robust    perform robust mean estimation of average (for VOIs
%                   with at least 4 voxels, default: false)
%        .rfx       return RFX table if possible (default: true)
%        .vl        indices of VOIs in object to sample (default: all)
%
% Output fields:
%
%       vb          SxCxV double table with data
%       vbv         SxCxV cell array with source data
%       vbi         1xV or SxV cell array with source voxel indices
%
% Using: applyspmsnc, bvcoordconv, findfirst, fitrobustbisquare,
%        flexinterpn_method, lsqueeze, makelabel, meannoinfnan.

% Version:  v1.1
% Build:    16020315
% Date:     Feb-03 2016, 3:59 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2012, 2014, 2016, Jochen Weber
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
applyspmsnc        = ne_methods.applyspmsnc;
bvcoordconv        = ne_methods.bvcoordconv;
findfirst          = ne_methods.findfirst;
fitrobustbisquare  = ne_methods.fitrobustbisquare;
flexinterpn_method = ne_methods.flexinterpn_method;
meannoinfnan       = ne_methods.meannoinfnan;

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm') || ...
   ((numel(vfile) ~= 1 || ~xffisobject(vfile, true, 'voi')) && ...
    (~isa(vfile, 'double') || size(vfile, 2) ~= 3 || isempty(vfile) || ...
     any(isinf(vfile(:)) | isnan(vfile(:)) | vfile(:) < -128 | vfile(:) > 128)))
    error('neuroelf:xff:badArgument', 'Invalid object handle in call.');
end
bc = xo.C;
if bc.ProjectType ~= 1
    error('neuroelf:xff:badArgument', 'Only valid for VTC based GLM files.');
end
switch (bc.ProjectTypeRFX)
    case 0
        ns = 1;
        np = bc.NrOfPredictors;
    case 1
        ns = numel(bc.GLMData.Subject);
        np = bc.NrOfSubjectPredictors;
    otherwise
        error('neuroelf:xff:badArgument', 'Invalid/unsupported ProjectTypeRFX flag in file.');
end
rtv = bc.RunTimeVars;
if ~isa(vfile, 'double')
    vc = vfile.C;
else
    vc = struct('ReferenceSpace', 'tal', 'VOI', struct('Voxels', vfile));
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'c') || ~isa(opts.c, 'double') || isempty(opts.c) || ...
   (size(opts.c, 2) ~= np && size(opts.c, 2) ~= (np - 1)) || ...
    any(isinf(opts.c(:)) | isnan(opts.c(:))) || any(any(opts.c < 0, 2) & sum(opts.c, 2) ~= 0)
    opts.c = eye(np);
end
nc = size(opts.c, 1);
if ~isfield(opts, 'interp') || ~ischar(opts.interp) || isempty(opts.interp) || ...
   ~any(lower(opts.interp(1)) == 'clns')
    opts.interp = 'n';
else
    opts.interp = lower(opts.interp(1));
end
if ~isfield(opts, 'refspace') || ~ischar(opts.refspace) || ~any(strcmpi(opts.refspace(:)', {'bvi', 'bvs'}))
    opts.refspace = 'bvi';
else
    opts.refspace = lower(opts.refspace(:)');
end
if ~isfield(opts, 'rmean') || ~islogical(opts.rmean) || numel(opts.rmean) ~= 1
    opts.rmean = false;
end
if ~isfield(opts, 'rfx') || ~islogical(opts.rfx) || numel(opts.rfx) ~= 1
    opts.rfx = true;
end
if ~isfield(opts, 'robust') || ~islogical(opts.robust) || numel(opts.robust) ~= 1
    opts.robust = false;
end
if ~isfield(opts, 'vl') || ~isa(opts.vl, 'double') || isempty(opts.vl) || ...
    any(isinf(opts.vl(:)) | isnan(opts.vl(:)))
    opts.vl = 1:numel(vc.VOI);
else
    opts.vl = intersect(opts.vl(:)', 1:numel(vc.VOI));
end

% prepare output
nv = numel(opts.vl);
vb = zeros(ns, nc, nv);
if nargout > 1
    vbv = cell(size(vb));
end

% interpolation option
switch (opts.interp)
    case 'c'
        ipo = 'cubic';
    case 'l'
        ipo = 'linear';
    case 'n'
        ipo = 'nearest';
    case 's'
        ipo = 'lanczos3';
end

% no subject SPMsn information
if bc.ProjectTypeRFX == 0 || ~isfield(rtv, 'SubjectSPMsn') || ...
   ~isstruct(rtv.SubjectSPMsn) || isempty(fieldnames(rtv.SubjectSPMsn))

    % get bounding box
    bb = aft_BoundingBox(xo);

    % convert coordinates
    voi = cell(1, numel(opts.vl));
    voiui = voi;
    if strcmpi(vc.ReferenceSpace, 'tal')
        convtype = 'tal2bvc';
    else
        convtype = [opts.refspace '2bvc'];
    end
    if opts.interp == 'n'
        if bc.ProjectTypeRFX > 0
            bms = size(bc.GLMData.RFXGlobalMap);
        else
            bms = size(bc.GLMData.MultipleRegressionR);
        end
        convtype(end) = 'x';
    end
    for c = 1:numel(voi)
        voic = vc.VOI(opts.vl(c)).Voxels;
        if isempty(voic)
            if opts.interp ~= 'n'
                voi{c} = [-1, -1, -1];
            end
            continue;
        end
        voi{c} = bvcoordconv(voic, convtype, bb);
        if opts.interp == 'n'
            voi{c}(isnan(voi{c})) = [];
        end
        if convtype(end) == 'x'
            [voi{c}, voiu, voiui{c}] = unique(voi{c}(:));
        end
    end

    % voxel indices
    if nargout > 2
        vbi = voi(:)';
        for vc = 1:nv
            vbi{vc} = ne_methods.lsqueeze(vbi{vc}(voiui{vc}));
        end
    end

    % for RFX GLMs
    if bc.ProjectTypeRFX == 1

        % iterate over subjects
        for sc = 1:ns

            % for non-interpolated data, it's very easy
            if opts.interp == 'n'

                % get size offset
                so = prod(bms(1:3));

                % iterate over VOIs first
                for vc = 1:nv

                    % reject empty VOIs
                    if isempty(voi{vc})
                        continue;
                    end

                    % get sampled values per contrast
                    for cc = 1:nc
                        conval = zeros(numel(voiui{vc}), 1);
                        for pc = find(opts.c(cc, :) ~= 0)
                            xval = bc.GLMData.Subject(sc).BetaMaps((pc - 1) * so + voi{vc});
                            conval = conval + opts.c(cc, pc) .* xval(voiui{vc});
                        end
                        if nargout > 1
                            vbv{sc, cc, vc} = conval;
                        end
                        conval(conval == 0) = [];
                        if ~opts.robust || numel(conval) < 4
                            vb(sc, cc, vc) = sum(conval) / numel(conval);
                        else
                            vb(sc, cc, vc) = fitrobustbisquare( ...
                                ones(numel(conval), 1), double(conval));
                        end
                    end
                end

            % for interpolated data it's more complicated
            else

                % iterate over contrasts, predictors and VOIs
                for cc = 1:nc
                    for pc = find(opts.c(cc, :) ~= 0)
                        cval = opts.c(cc, pc);
                        bm = bc.GLMData.Subject(sc).BetaMaps(:, :, :, pc);
                        bn = (bm ~= 0);
                        bm(~bn) = Inf;
                        if opts.rmean && pc < np
                            bm(bn) = bm(bn) - mean(bm(bn));
                        end
                        for vc = 1:nv
                            ipv = flexinterpn_method(bm, voi{vc}, ipo);
                            if nargout > 1
                                if isempty(vbv{sc, cc, vc})
                                    vbv{sc, cc, vc} = zeros(size(ipv));
                                end
                                vbv{sc, cc, vc} = vbv{sc, cc, vc} + cval .* ipv;
                            end
                            ipv(isinf(ipv) | isnan(ipv) | ipv == 0) = [];
                            if ~opts.robust || numel(ipv) < 4
                                vb(sc, cc, vc) = vb(sc, cc, vc) + cval * mean(ipv);
                            else
                                vb(sc, cc, vc) = vb(sc, cc, vc) + cval * ...
                                    fitrobustbisquare(ones(numel(ipv), 1), double(ipv));
                            end
                        end
                    end
                end
            end
        end

    % for non-RFX GLMs
    else

        % for non-interpolated data, it's very easy
        if opts.interp == 'n'

            % get size offset
            so = prod(bms(1:3));

            % iterate over VOIs first
            for vc = 1:nv

                % reject empty VOIs
                if isempty(voi{vc})
                    continue;
                end

                % get sampled values per contrast
                for cc = 1:nc
                    conval = zeros(numel(voiui{vc}), 1);
                    for pc = find(opts.c(cc, :) ~= 0)
                        xval = bc.GLMData.BetaMaps((pc - 1) * so + voi{vc});
                        conval = conval + opts.c(cc, pc) .* xval(voiui{vc});
                    end
                    if nargout > 1
                        vbv{1, cc, vc} = conval;
                    end
                    conval(conval == 0) = [];
                    if ~opts.robust || numel(conval) < 4
                        vb(1, cc, vc) = sum(conval) / numel(conval);
                    else
                        vb(1, cc, vc) = fitrobustbisquare( ...
                            ones(numel(conval), 1), double(conval));
                    end
                end
            end

        % for interpolated data it's more complicated
        else

            % iterate over contrasts, predictors and VOIs
            for cc = 1:nc
                for pc = find(opts.c(cc, :) ~= 0)
                    cval = opts.c(cc, pc);
                    bm = bc.GLMData.BetaMaps(:, :, :, pc);
                    bn = (bm ~= 0);
                    bm(~bn) = Inf;
                    if opts.rmean && pc < np
                        bm(bn) = bm(bn) - mean(bm(bn));
                    end
                    for vc = 1:nv
                        ipv = flexinterpn_method(bm, voi{vc}, ipo);
                        if nargout > 1
                            if isempty(vbv{1, cc, vc})
                                vbv{1, cc, vc} = zeros(size(ipv));
                            end
                            vbv{1, cc, vc} = vbv{1, cc, vc} + cval .* ipv;
                        end
                        ipv(isinf(ipv) | isnan(ipv) | ipv == 0) = [];
                        if ~opts.robust || numel(ipv) < 4
                            vb(1, cc, vc) = vb(1, cc, vc) + cval * mean(ipv);
                        else
                            vb(1, cc, vc) = vb(1, cc, vc) + cval * ...
                                fitrobustbisquare(ones(numel(ipv), 1), double(ipv));
                        end
                    end
                end
            end
        end

        % try to re-shape into RFX betas
        if opts.rfx && bc.SeparatePredictors == 2

            % create copy of vb
            vbc = vb;
            if nargout > 1
                vbvc = vbv;
            end

            % get subjects and subject predictors and full predictor names
            subjects = glm_Subjects(xo);
            subpreds = glm_SubjectPredictors(xo);
            fpnames = {bc.Predictor.Name2};

            % generate new vb
            vb = NaN .* zeros(numel(subjects), numel(subpreds), size(vbc, 3));
            if nargout > 1
                vbv = cell(size(vb));
            end

            % look up indices
            for sc = 1:numel(subjects)
                for pc = 1:numel(subpreds)

                    % copy data
                    vbi = findfirst(~cellfun('isempty', regexpi(fpnames, ...
                        sprintf('^subject\\s+%s\\:\\s+%s$', subjects{sc}, subpreds{pc}))));
                    if ~isempty(vbi)
                        vb(sc, pc, :) = vbc(1, vbi, :);
                        if nargout > 1
                            vbv(sc, pc, :) = vbvc(1, vbi, :);
                        end
                    end
                end
            end
        end
    end

% with Subject-speficic normalization (RFX-only!)
else

    % subject-based indices!
    if nargout > 2 || ipo(1) == 'n'
        vbi = cell(ns, nv);
    end

    % ensure VOI is TAL space
    if lower(vc.ReferenceSpace(1)) ~= 't'
        for vvc = 1:numel(vc.VOI)
            vc.VOI(vvc).Voxels = 128 - vc.VOI(vvc).Voxels;
        end
    end

    % get number of subject predictors
    mbx = aft_BoundingBox(xo);
    cfr = struct('Trf', mbx.QuatB2T);
    tmatc = inv(double(cfr.Trf));
    if isfield(rtv, 'TrfPlus') && isequal([4, 4], size(rtv.TrfPlus)) && ...
        any(any(rtv.TrfPlus ~= eye(4)))
        rtvtrfplus = inv(rtv.TrfPlus);
    else
        rtvtrfplus = eye(4);
    end

    % iterate over subjects
    sids = ne_methods.makelabel(glm_Subjects(xo));
    sidsn = rtv.SubjectSPMsn;
    if isfield(rtv, 'SubjectTrfPlus') && isstruct(rtv.SubjectTrfPlus) && ...
        numel(rtv.SubjectTrfPlus) == 1 && ...
        isequal(fieldnames(rtv.SubjectTrfPlus), fieldnames(sidsn))
        sidtrfpl = rtv.SubjectTrfPlus;
    else
        sidtrfpl = struct;
        sidsnf = fieldnames(sidsn);
        for sc = 1:numel(sidsnf)
            sidtrfpl.(sidsnf{sc}) = eye(4);
        end
    end
    indexvol = reshape(1:numel(bc.GLMData.RFXGlobalMap), size(bc.GLMData.RFXGlobalMap));
    indexnum = numel(indexvol);
    for sc = 1:ns

        % test for and get subject-specific information
        if ~isfield(sidsn, sids{sc}) || ~isfield(sidtrfpl, sids{sc})
            error('neuroelf:xff:badObject', 'SPMsn information not available for all subjects.');
        end
        sn = sidsn.(sids{sc});
        trfplus = sidtrfpl.(sids{sc});

        % computations
        tmat = tmatc * inv(trfplus) * rtvtrfplus;
        ivgm = inv(sn.VG(1).mat);

        % iterate over contrasts, predictors and VOIs
        for cc = 1:nc
            for pc = find(opts.c(cc, :) ~= 0)
                cval = opts.c(cc, pc);
                if ipo(1) ~= 'n' || opts.rmean
                    bm = bc.GLMData.Subject(sc).BetaMaps(:, :, :, pc);
                    bn = (bm ~= 0);
                    bm(~bn) = Inf;
                    if opts.rmean && pc < np
                        bm(bn) = bm(bn) - meannoinfnan(bm(bn));
                    end
                end
                for vvc = 1:nv
                    if (nargout > 2 || ipo(1) == 'n') && cc == 1 && pc == 1
                        voxcrd = applyspmsnc(vc.VOI(opts.vl(vvc)).Voxels, sn.Tr, ...
                            sn.VG(1).dim, ivgm, tmat * sn.VF.mat * sn.Affine);
                        vbi{sc, vvc} = flexinterpn_method(indexvol, voxcrd, 0, 'nearest');
                    end
                    if ipo(1) ~= 'n' || opts.rmean
                        ipv = flexinterpn_method(bm, ...
                            applyspmsnc(vc.VOI(opts.vl(vvc)).Voxels, sn.Tr, ...
                            sn.VG(1).dim, ivgm, tmat * sn.VF.mat * sn.Affine), 0, ipo);
                    else
                        [upv, upi1, upi2] = unique(vbi{sc, vvc});
                        ipv = zeros(size(upv));
                        ipv(upv > 0) = bc.GLMData.Subject(sc).BetaMaps((pc - 1) * indexnum + upv(upv > 0));
                        ipv = ipv(upi2);
                    end
                    if nargout > 1
                        if isempty(vbv{sc, cc, vvc})
                            vbv{sc, cc, vvc} = zeros(size(ipv));
                        end
                        vbv{sc, cc, vvc} = vbv{sc, cc, vvc} + cval .* ipv;
                    end
                    ipv(isinf(ipv) | isnan(ipv) | ipv == 0) = [];
                    if ~opts.robust || numel(ipv) < 4
                        vb(sc, cc, vvc) = vb(sc, cc, vvc) + cval * mean(ipv);
                    else
                        vb(sc, cc, vvc) = vb(sc, cc, vvc) + cval * ...
                            fitrobustbisquare(ones(numel(ipv), 1), double(ipv));
                    end
                end
            end
        end
    end
end
