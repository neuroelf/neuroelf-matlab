function [cxo, condsel] = aft_PrepForConnectivity(xo, varargin)
% AFT::PrepForConnectivity  - prepare 4D data for connectivity
%
% FORMAT:       cobj = obj.PrepForConnectivity([opts])
%
% Input fields:
%
%       opts        settings
%        .condsel   condition selection using regexp (char or cell, GLM only)
%        .globsigd  also add diff of global signals as nuisance regressors
%        .globsigs  add global signals as confound, one of (default: 1)
%                   0 - none
%                   1 - entire dataset (voxels with > .1 median variance)
%                   2 or more, perform PCA of time courses and first N
%                   xff object(s), extract average time course from masks
%        .motpars   regress out motion (from RunTimeVars, default: true)
%        .motparsd  also add diff of motion parameters (default: true)
%        .motparsq  also add squared motion parameters (default: false)
%        .subsel    subject selection (GLM only)
%        .tfiltbp   band-pass settings for Butterworth filter ([0.01 0.1])
%        .tfiltbwo  Butterworth filter order (default: 2, see
%                   http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3849338/)
%        .tfilter   filtering cut-off (in seconds, default: 120)
%        .tfiltlp   low-pass (smoothing kernel) in seconds (default: 4)
%                   only applied if tfilttyp not 'bw' (or 'none')
%        .tfilttyp  filter type, one of 'bw', {'dct'}, 'fourier', or 'poly'
%        .trans     perform either 'psc' or 'z' transform (default: 'none')
%        .xconfound additional confounds, if PRT object, regress out task
%
% Output fields:
%
%       cvtc        VTC prepared for connectivity analyses
%
% TYPES: GLM, HDR, HEAD, VTC
%
% Notes: a VTC will be forced to FileVersion 3 / DataType 2!
%        most options will be ignored for GLM objects (filtering, etc.)
%
% Using: calcbetas, multimatch, psctrans, tempfilter, ztrans.

% Version:  v1.1
% Build:    16081911
% Date:     Aug-19 2016, 11:01 AM EST
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

% neuroelf library
global ne_methods;
multimatch = ne_methods.multimatch;
ztrans = ne_methods.ztrans;

% argument check
if nargin < 1 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'glm', 'hdr', 'head', 'vtc'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get object, type, and RunTimeVars
bc = xo.C;
xotype = lower(xo.S.Extensions{1});
rtv = bc.RunTimeVars;

% create opts if necessary
if nargin < 2 || numel(varargin{end}) ~= 1 || ~isstruct(varargin{end})
    opts = struct;
else
    opts = varargin{end};
end

% get data (VTC)
if strcmp(xotype, 'vtc')

    % get VTCData
    vtcd = double(bc.VTCData);
    vtcs = size(vtcd);
    voxs = vtcs(2:4);
    vtcd = reshape(vtcd, vtcs(1), prod(vtcs(2:end)));

% HDR
elseif strcmp(xotype, 'hdr')

    % requires multiple volumes
    if size(bc.VoxelData, 4) < 20
        error('neuroelf:xff:badObject', 'Connectivity requires at least 20 volumes.');
    end

    % get VoxelData
    vtcd = double(bc.VoxelData);
    if size(vtcd, 5) > 1
        error('neuroelf:xff:badObject', 'Connectivity requires a 4D-only dataset.');
    end
    vtcs = size(vtcd);
    voxs = vtcs(1:3);
    vtcd = reshape(vtcd, prod(vtcs(1:3)), vtcs(4))';
    vtcs = size(vtcd);

    % requires scaling
    if ~any(bc.ImgDim.ScalingSlope == [0, 1])

        % and intercept
        if bc.ImgDim.ScalingIntercept ~= 0
            vtcd = bc.ImgDim.ScalingIntercept + bc.ImgDim.ScalingSlope .* vtcd;

        % scaling only
        else
            vtcd = bc.ImgDim.ScalingSlope .* vtcd;
        end

    % only add intercept
    elseif bc.ImgDim.ScalingIntercept ~= 0
        vtcd = vtcd + bc.ImgDim.ScalingIntercept;
    end

% HEAD
elseif strcmp(xotype, 'head')

    % requires multiple volumes
    if numel(bc.Brick) < 20
        error('neuroelf:xff:badObject', 'Connectivity requires at least 20 volumes.');
    end

    % generate data first
    vsz = size(bc.Brick(1).Data);
    if numel(vsz) < 3
        vsz(3) = 1;
    end
    voxs = vsz;
    vtcd = zeros([numel(bc.Brick), vsz]);
    for vc = 1:numel(bc.Brick)
        vtcd(vc, :) = bc.Brick(vc).ScalingFactor .* ...
            reshape(double(bc.Brick(vc).Data), 1, prod(vsz));
    end
    vtcs = size(vtcd);

% GLM
else

    % requires VTC-based RFX GLM with at least 12 regular regressors
    if bc.ProjectTypeRFX ~= 1 || bc.ProjectType ~= 1 || bc.NrOfSubjectPredictors < 13
        error('neuroelf:xff:badArgument', 'Only available for VTC-based RFX GLMs.');
    end

    % disable filter settings
    opts.motpars = false;
    opts.tfilter = Inf;
    opts.tfiltbp = [0, 1];
    opts.tfiltlp = 0;
    opts.tfilttyp = 'none';

    % conditions
    if ~isfield(opts, 'condsel') || isempty(opts.condsel) || ...
       (~ischar(opts.condsel) && ~iscell(opts.condsel))
        opts.condsel = {'.*'};
    elseif ischar(opts.condsel)
        opts.condsel = {opts.condsel(:)'};
    else
        opts.condsel = opts.condsel(:);
    end
    condsel = zeros(1, numel(opts.condsel));
    spred = glm_SubjectPredictors(xo);
    if strcmpi(spred{end}, 'constant')
        spred(end) = [];
    end
    for vc = 1:numel(condsel)
        opts.condsel{vc} = find(multimatch(spred, opts.condsel{vc}(:)', true) > 0);
        condsel(vc) = numel(opts.condsel{vc});
    end
    if any(condsel == 0) || sum(condsel) ~= numel(cat(1, opts.condsel{:}))
        error('neuroelf:xff:badArgument', 'Invalid condition selection.');
    end

    % subject selection (for this function, only double supported)
    subids = glm_Subjects(xo);
    numsubs = numel(subids);
    if ~isfield(opts, 'subsel') || isempty(opts.subsel) || ...
       ((~isa(opts.subsel, 'double') || ...
        any(isinf(opts.subsel(:)) | isnan(opts.subsel(:))) || ...
        numel(unique(round(opts.subsel(:)))) ~= numel(opts.subsel) || ...
        any(opts.subsel(:) < 1 | opts.subsel(:) > numsubs)) && ...
        (~iscell(opts.subsel) || ~all(cellfun(@ischar, opts.subsel(:))) || ...
         any(cellfun('isempty', opts.subsel(:)))))
        opts.subsel = 1:numsubs;
    elseif isa(opts.subsel, 'double')
        opts.subsel = round(opts.subsel(:)');
    else
        opts.subsel = find(multimatch(subids, opts.subsel(:)) > 0);
    end
    subsel = opts.subsel;
    numsubs = numel(subsel);

    % vtcd
    vtcs = [sum(condsel), size(bc.GLMData.RFXGlobalMap), numsubs];
    vtcd = zeros(vtcs);
    voxs = vtcs(2:4);

    % for each selection, get data
    vci = 1;
    for sc = 1:numsubs
        for vc = 1:numel(condsel)
            vtcd(vci:vci+condsel(vc)-1, :, :, :, sc) = permute( ...
                bc.GLMData.Subject(subsel(sc)).BetaMaps(:, :, :, opts.condsel{vc}), ...
                [4, 1, 2, 3]);
        end
    end
    vtcd = reshape(vtcd, [sum(condsel), prod(vtcs(2:4)), numsubs]);
    vtcs = size(vtcd);
end

% number of volumes (per subject)
nvol = size(vtcd, 1);
if ~strcmp(xotype, 'glm')
    condsel = nvol;
end
if ~isfield(opts, 'globsigd') || ~islogical(opts.globsigd) || numel(opts.globsigd) ~= 1
    opts.globsigd = false;
end
cleargs = false;
if ~isfield(opts, 'globsigs')
    opts.globsigs = 1;
end
if isfield(opts, 'globsigs') && ischar(opts.globsigs) && ...
    any(strcmpi(opts.globsigs(:)', {'cb2', 'cb23', 'cb3', 'cb33'}))
    mp = neuroelf_path('masks');
    switch (lower(opts.globsigs(:)'))
        case {'cb2'}
            if exist([mp '/colin_brain_ICBMnorm_brain2mm.msk'], 'file') == 2
                try
                    opts.globsigs = {xff([mp '/colin_brain_ICBMnorm_brain2mm.msk'])};
                    cleargs = true;
                catch xfferror
                    fprintf('Error loading 2mm Colin brain mask: %s.\n', xfferror.message);
                end
            end
        case {'cb23'}
            if exist([mp '/colin_brain_ICBMnorm_gray2mm.msk'], 'file') == 2 && ...
                exist([mp '/colin_brain_ICBMnorm_white2mm.msk'], 'file') == 2 && ...
                exist([mp '/colin_brain_ICBMnorm_csfx2mm.msk'], 'file') == 2
                try
                    opts.globsigs = {xff([mp '/colin_brain_ICBMnorm_gray2mm.msk']); ...
                        xff([mp '/colin_brain_ICBMnorm_white2mm.msk']); ...
                        xff([mp '/colin_brain_ICBMnorm_csfx2mm.msk'])};
                    cleargs = true;
                catch xfferror
                    fprintf('Error loading 2mm Colin brain masks: %s.\n', xfferror.message);
                end
            end
        case {'cb3'}
            if exist([mp '/colin_brain_ICBMnorm_brain3mm.msk'], 'file') == 2
                try
                    opts.globsigs = {xff([mp '/colin_brain_ICBMnorm_brain3mm.msk'])};
                    cleargs = true;
                catch xfferror
                    fprintf('Error loading 3mm Colin brain mask: %s.\n', xfferror.message);
                end
            end
        case {'cb33'}
            if exist([mp '/colin_brain_ICBMnorm_gray3mm.msk'], 'file') == 2 && ...
                exist([mp '/colin_brain_ICBMnorm_white3mm.msk'], 'file') == 2 && ...
                exist([mp '/colin_brain_ICBMnorm_csfx3mm.msk'], 'file') == 2
                try
                    opts.globsigs = {xff([mp '/colin_brain_ICBMnorm_gray3mm.msk']); ...
                        xff([mp '/colin_brain_ICBMnorm_white3mm.msk']); ...
                        xff([mp '/colin_brain_ICBMnorm_csfx3mm.msk'])};
                    cleargs = true;
                catch xfferror
                    fprintf('Error loading 3mm Colin brain masks: %s.\n', xfferror.message);
                end
            end
    end
    if ischar(opts.globsigs)
        warning('neuroelf:xff:fileNorFound', 'Required mask file(s) not found.');
        opts.globsigs = [];
    end
end
if ~isfield(opts, 'globsigs') || ((~isa(opts.globsigs, 'double') || ...
     numel(opts.globsigs) ~= 1 || isinf(opts.globsigs) || isnan(opts.globsigs) || opts.globsigs < 0) && ...
    (~iscell(opts.globsigs) || isempty(opts.globsigs) || (~ischar(opts.globsigs{1}) && ...
     ~xffisobject(opts.globsigs{1}, true, {'hdr', 'msk', 'vmr', 'voi'}))) && ...
    (~xffisobject(opts.globsigs, true, {'hdr', 'msk', 'vmr', 'voi'})))
    opts.globsigs = 0;
elseif isa(opts.globsigs, 'double')
    opts.globsigs = floor(opts.globsigs);
elseif xffisobject(opts.globsigs, true)
    opts.globsigs = {opts.globsigs};
else
    opts.globsigs = opts.globsigs(:)';
    gcc = numel(opts.globsigs);
    if isa(opts.globsigs{end}, 'double')
        gcc = gcc - 1;
    end
    for gc = gcc:-1:1
        if (~ischar(opts.globsigs{gc}) || isempty(opts.globsigs{gc}) || ...
            ~isfield(rtv, 'GlobSigs') || ~iscell(rtv.GlobSigs) || ...
             all(cellfun('isempty', regexpi(rtv.GlobSigs(:, 3), opts.globsigs{gc}(:)')))) && ...
            (numel(opts.globsigs{gc}) ~=1 || ...
             ~xffisobject(opts.globsigs{gc}, true, {'hdr', 'msk', 'vmr', 'voi'}))
            fprintf('Invalid globsigs option (entry %d).\n', gc);
            opts.globsigs(gc) = [];
        end
    end
    if numel(opts.globsigs) == 1 && isa(opts.globsigs{1}, 'double')
        opts.globsigs = 0;
    end
end
if ~isfield(opts, 'motpars')
    opts.motpars = true;
end
if isfield(opts, 'motpars') && islogical(opts.motpars) && numel(opts.motpars) == 1 && ...
    opts.motpars
    if isfield(rtv, 'MotionParameters') && isnumeric(rtv.MotionParameters) && ...
        isequal(size(rtv.MotionParameters), [nvol, 6])
        opts.motpars = rtv.MotionParameters;
    else
        fprintf('Warning: RunTimeVars do not contain valid MotionParameters.\n');
        opts.motpars = [];
    end
end
if ~isfield(opts, 'motparsd') || ~islogical(opts.motparsd) || numel(opts.motparsd) ~= 1
    opts.motparsd = true;
end
if ~isfield(opts, 'motparsq') || ~islogical(opts.motparsq) || numel(opts.motparsq) ~= 1
    opts.motparsq = false;
end
if ~isfield(opts, 'tfiltbp') || ~isa(opts.tfiltbp, 'double') || numel(opts.tfiltbp) ~= 2 || ...
    any(isinf(opts.tfiltbp) | isnan(opts.tfiltbp) | opts.tfiltbp <= 0 | opts.tfiltbp > 1) || ...
    opts.tfiltbp(2) <= opts.tfiltbp(1)
    opts.tfiltbp = [0.01, 0.1];
end
if ~isfield(opts, 'tfiltbwo') || ~isa(opts.tfiltbwo, 'double') || numel(opts.tfiltbwo) ~= 1 || ...
    isinf(opts.tfiltbwo) || isnan(opts.tfiltbwo) || opts.tfiltbwo <= 0 || opts.tfiltbwo > 9
    opts.tfiltbwo = 2;
end
if ~isfield(opts, 'tfilter') || ~isa(opts.tfilter, 'double') || numel(opts.tfilter) ~= 1 || ...
    isnan(opts.tfilter) || opts.tfilter < 30
    opts.tfilter = 120;
end
if ~isfield(opts, 'tfiltlp') || ~isa(opts.tfiltlp, 'double') || numel(opts.tfiltlp) ~= 1 || ...
    isinf(opts.tfiltlp) || isnan(opts.tfiltlp) || opts.tfiltlp < 0
    opts.tfiltlp = 4;
end
if ~isfield(opts, 'tfilttype') || ~ischar(opts.tfilttype) || isempty(opts.tfilttype) || ...
   ~any(strcmpi(opts.tfilttype(:)', {'bw', 'dct', 'fourier', 'none', 'poly'}))
    opts.tfilttype = 'dct';
else
    opts.tfilttype = lower(opts.tfilttype(:)');
    if opts.tfilttype(1) == 'n'
        opts.tfilter = Inf;
    end
end
if ~isfield(opts, 'trans') || ~ischar(opts.trans) || isempty(opts.trans) || ...
   ~any(lower(opts.trans(1)) == 'npz')
    opts.trans = 'n';
else
    opts.trans = lower(opts.trans(1));
end
if ~isfield(opts, 'xconfound') || isempty(opts.xconfound)
    opts.xconfound = [];
elseif ischar(opts.xconfound)
elseif numel(opts.xconfound) == 1 && xffisobject(opts.xconfound, true, 'prt')
else
    error('neuroelf:xff:invalidOption', 'Invalid xconfound option.');
end

% create design matrix for regression
X = ones(nvol, 1);
xfs = 1;
bwb = [];

% filter content
if opts.tfilter > 0 && ~isinf(opts.tfilter) && opts.tfilttype(1) ~= 'b'

    % prepare tempfilter options
    topts = opts;
    topts.spat = false;
    topts.tdim = 1;
    topts.temp = true;
    topts.tempdt = false;
    if opts.tfilttype(1) == 'd'
        topts.tempdct = opts.tfilter / (0.001 * bc.TR);
        topts.temppoly = 0;
        topts.tempsc = 0;
    elseif opts.tfilttype(1) == 'f'
        topts.tempdct = Inf;
        topts.temppoly = 0;
        topts.tempsc = floor(nvol * 0.001 * bc.TR / opts.tfilter);
    else
        topts.tempdct = Inf;
        topts.temppoly = floor(2 * nvol * 0.001 * bc.TR / opts.tfilter - 1);
        topts.tempsc = 0;
    end

    % temp filter data of first object
    [null, Xf] = ne_methods.tempfilter(zeros(nvol, 1), topts);

    % and temp filter regressors
    X = [X, Xf];
    xfs = size(X, 2);

% butterworth filter
elseif opts.tfilter > 0 && ~isinf(opts.tfilter)

    % correct for TR and butterworth
    bpset = min(1, opts.tfiltbp .* (0.002 .* bc.TR));

    % get filter coefficient (for filter function)
    [bwb, bwa] = ne_methods.butterworth(bpset, ceil(opts.tfiltbwo / 2));

    % compute and remove mean
    mvtc = mean(vtcd, 1);
    vtcd = vtcd - repmat(mvtc, nvol, 1);

    % apply filter
    vtcd = filter(bwb, bwa, vtcd, [], 1);

    % add back mean
    vtcd = vtcd + repmat(mvtc, nvol, 1);
end

% add regressors
for ac = 1:nargin-1

    % numeric matrix
    if isnumeric(varargin{ac}) && size(varargin{ac}, 1) == nvol

        % add to X
        X = [X, varargin{ac}(:, :)];

    % RTC/SDM
    elseif numel(varargin{ac}) == 1 && xffisobject(varargin{ac}, true, 'sdm')

        % get content
        sdmc = varargin{ac}.C;

        % size check
        if size(sdmc.SDMMatrix, 1) ~= nvol
            error('neuroelf:xff:badArgument', 'SDM size mismatch with VTCData.');
        end

        % add to design
        X = [X, sdmc.SDMMatrix];

    % text/mat file
    elseif ischar(varargin{ac}) && ~isempty(varargin{ac})

        % try to load
        try
            xa = load(varargin{ac}(:)');
        catch xfferror
            neuroelf_lasterr(xfferror);
            try
                xa = xff(varargin{ac}(:)');
                if xffisobject(xa, true)
                    xac = xa.C;
                    delete(xa);
                else
                    xac =[];
                end
                if isstruct(xac) && isfield(xac, 'SDMMatrix')
                    xa = xac.SDMMatrix;
                    if size(xa, 1) ~= nvol
                        error('neuroelf:xff:badArgument', 'SDM size mismatch with VTCData.');
                    end
                else
                    xa = [];
                end
            catch xfferror
                neuroelf_lasterr(xfferror);
                xa = [];
            end
        end
        if isnumeric(xa) && size(xa, 1) == nvol
            X = [X, xa(:, :)];
        elseif isstruct(xa) && numel(fieldnames(xa)) == 1
            xf = fieldnames(xa);
            xa = xa.(xf{1});
            if isnumeric(xa) && size(xa, 1) == nvol
                X = [X, xa(:, :)];
            end
        end
    end
    if size(X, 2) > xfs && ~isempty(bwb)
        X(:, xfs+1:end) = filter(bwb, bwa, X(:, xfs+1:end), [], 1);
    end
end
X(:, any(isnan(X) | isinf(X))) = [];
xs = (sum(abs(diff(X))) < sqrt(eps));
xs(1) = false;
if any(xs)
    X(:, xs) = [];
end
if size(X, 2) > xfs
    X(:, xfs+1:end) = X(:, xfs+1:end) - X(:, 1:xfs) * (((X(:, 1:xfs)' * X(:, 1:xfs)) \ X(:, 1:xfs)') * X(:, xfs+1:end));
    X(:, any(isnan(X) | isinf(X))) = [];
    xs = (sum(abs(diff(X))) < sqrt(eps));
    xs(1) = false;
    if any(xs)
        X(:, xs) = [];
    end
    xfs = size(X, 2);
end

% add motion parameters
if ~isempty(opts.motpars) && size(opts.motpars, 1) == nvol
    mp = ztrans(opts.motpars);
    if opts.motparsd
        mpd = ztrans([zeros(1, size(mp, 2)); diff(mp)]);
    else
        mpd = zeros(nvol, 0);
    end
    if opts.motparsq
        mpq = ztrans(mp .* mp);
    else
        mpq = zeros(nvol, 0);
    end
    X = [X, mp, mpd, mpq];
    X(:, any(isnan(X) | isinf(X))) = [];
    xs = (sum(abs(diff(X))) < sqrt(eps));
    xs(1) = false;
    if any(xs)
        X(:, xs) = [];
    end
    if ~isempty(bwb)
        X(:, xfs+1:end) = filter(bwb, bwa, X(:, xfs+1:end), [], 1);
    end
    X(:, xfs+1:end) = X(:, xfs+1:end) - X(:, 1:xfs) * (((X(:, 1:xfs)' * X(:, 1:xfs)) \ X(:, 1:xfs)') * X(:, xfs+1:end));
    X(:, any(isnan(X) | isinf(X))) = [];
    xs = (sum(abs(diff(X))) < sqrt(eps));
    xs(1) = false;
    if any(xs)
        X(:, xs) = [];
    end
end

% global signals
if isa(opts.globsigs, 'double') && numel(opts.globsigs) == 1 && opts.globsigs > 0

    % for all but GLMs
    if ~strcmp(xotype, 'glm')

        % load mask
        try
            gsmask = xff([neuroelf_path('masks') filesep 'colin_brain_ICBMnorm_brain2mm.msk']);
        catch xfferror
            error('neuroelf:xff:maskNotReadable', 'Global signal mask not readable: %s', xfferror.message);
        end

        % create grid for sampling
        [msx, msy, msz] = ndgrid(1:voxs(1), 1:voxs(2), 1:voxs(3));
        if strcmp(xotype, 'vtc')
            msx = [msx(:), msy(:), msz(:)];
            msx = ne_methods.bvcoordconv(msx, 'bvc2tal', aft_BoundingBox(xo));
        elseif strcmp(xotype, 'hdr')
            msx = [msx(:), msy(:), msz(:), ones(numel(msx))];
            msy = hdr_CoordinateFrame(xo);
            msx = msx * msy.Trf';
        else
            msx = [msx(:), msy(:), msz(:), ones(numel(msx))];
            msy = head_CoordinateFrame(xo);
            msx = msx * msy.Trf';
        end
        msz = (aft_SampleData3D(gsmask, msx(:, 1:3)) < 0.5);
        msz = msz(:)';
        aft_ClearObject(gsmask);

        % perform PCA of voxels most contributing to variance
        vvtc = var(vtcd);
        vvtc(isinf(vvtc) | isnan(vvtc) | msz) = 0;
        vtcim = (vvtc > (0.1 * median(vvtc(vvtc > 0))));

        % get data
        gsd = ztrans(vtcd(:, vtcim));

        % remove known influences already
        gsd = gsd - X * (((X' * X) \ X') * gsd);
        pcvecs = ne_methods.ne_fastica(gsd, struct('step', 'pca'));

        % keep requested components
        X = [X, pcvecs(:, end:-1:end+1-opts.globsigs)];

    % GLMs
    else
    end

% global signal objects
elseif iscell(opts.globsigs)

    % get VTC bounding box
    bb = aft_BoundingBox(xo);

    % number of components
    if numel(opts.globsigs) > 1 && isa(opts.globsigs{end}, 'double') && ...
        numel(opts.globsigs{end}) == 1 && ~isinf(opts.globsigs{end}) && ...
       ~isnan(opts.globsigs{end}) && opts.globsigs{end} >= 1
        ngs = round(opts.globsigs{end});
        opts.globsigs(end) = [];
    else
        ngs = 1;
    end

    % for each object, find coordinates
    for oc = 1:numel(opts.globsigs)
        
        % already extracted
        if ischar(opts.globsigs{oc})
            gsidx = find(~cellfun('isempty', regexpi(rtv.GlobSigs(:, 3), opts.globsigs{gc}(:)')));
            for gic = 1:numel(gsidx)
                X = [X, ztrans(rtv.GlobSigs{gsidx(gic), 2})];
            end
            continue;
        end

        % get object's content
        gso = opts.globsigs{oc}.C;

        % for VOIs
        if xffisobject(opts.globsigs{oc}, true, 'voi')

            % get all coordinates as one
            vtcvox = unique(cat(1, gso.VOI.Voxels), 'rows');
            
            % recompute into VTC voxels
            if upper(gso.ReferenceSpace(1)) == 'T'
                vtcvox = ne_methods.bvcoordconv(vtcvox, 'tal2bvx', bb);
            else
                vtcvox = ne_methods.bvcoordconv(vtcvox, 'bvi2bvx', bb);
            end
            vtcvox(isnan(vtcvox) | vtcvox < 1) = [];
            vtcvox = unique(vtcvox);
            
        % for all other objects
        else

            % sample box
            vtcvox = find(ne_methods.lsqueeze(aft_SampleBVBox(opts.globsigs{oc}, bb, 1)) >= 0.5);
        end
        
        % get data
        gsd = ztrans(vtcd(:, vtcvox));
        gsd(:, any(isinf(gsd) | isnan(gsd)) | all(gsd == 0)) = [];

        % remove known influences already
        gsd = gsd - X * (((X' * X) \ X') * gsd);
        if ngs > 1
            pcvecs = ne_methods.ne_fastica(gsd, struct('step', 'pca'));
        else
            pcvecs = mean(gsd, 2);
        end

        % keep requested components
        X = [X, pcvecs(:, end:-1:end+1-ngs)];
    end
end

% cleanup
X(:, any(isnan(X) | isinf(X))) = [];
X(:, sum(abs(diff(X))) < sqrt(eps)) = [];

% constant
X(:, end + 1) = 1;

% clear global signal masks
if cleargs
    clearxffobjects(opts.globsigs);
end

% calcbetas
b = ne_methods.calcbetas(X, vtcd);

% set constant to 0 (unless z-trans)
if opts.trans ~= 'z'
    b(:, end) = 0;
end

% regress out
vtcd = vtcd - X * b';

% LP filter (smoothing)
if opts.tfiltlp > 0 && ~any('bn' == opts.tfilttype(1))
    vtcd = ne_methods.flexinterpn(vtcd, [Inf, Inf; ones(2, 2); size(vtcd)], ...
        {ne_methods.smoothkern(opts.tfiltlp / (0.001 * bc.TR), 1e-6), [0; 1 ;0]}, {1, 1});
end

% transformation
if opts.trans == 'p'
    vtcd = ne_methods.psctrans(vtcd);
elseif opts.trans == 'z'
    vtcd = ztrans(vtcd);
end

% remove infs/nans
vtcd(repmat(any(isinf(vtcd) | isnan(vtcd)), nvol, 1)) = 0;

% copy object
cxo = aft_CopyObject(xo);

% set content
if strcmp(xotype, 'vtc')
    cxo.C.FileVersion = max(3, bc.FileVersion);
    cxo.C.DataType = 2;
    cxo.C.VTCData = single(reshape(vtcd, [nvol, vtcs(2:end)]));
elseif strcmp(xotype, 'hdr')
    cxo.C.ImgDim.DataType = 16;
    cxo.C.ImgDim.BitsPerPixel = 32;
    cxo.C.ImgDim.ScalingSlope = 1;
    cxo.C.ImgDim.ScalingIntercept = 0;
    cxo.C.VoxelData = reshape(single(vtcd'), [vtcs(2:end), nvol]);
elseif strcmp(xotype, 'head')
    vtcd = reshape(vtcd, [nvol, vtcs(2:end)]);
    for vc = 1:nvol
        cxo.C.Brick(vc).DataType = 3;
        cxo.C.Brick(vc).Data = single(squeeze(vtcd(vc, :, :, :)));
        cxo.C.Brick(vc).ScalingFactor = 1;
    end
else
    vtcd = reshape(vtcd, [nvol, vtcs(2:4), numsubs]);
    for sc = 1:numsubs
        if istransio(cxo.C.GLMData.Subject(subsel(sc)).BetaMaps)
            cxo.C.GLMData.Subject(subsel(sc)).BetaMaps = ...
                resolve(cxo.C.GLMData.Subject(subsel(sc)).BetaMaps);
        end
        for vc = 1:numel(condsel)
             bc.GLMData.Subject(subsel(sc)).BetaMaps(:, :, :, opts.condsel{vc}) = ...
                 permute(vtcd(vci:vci+condsel(vc)-1, :, :, :, sc), [2, 3, 4, 1]);
        end
    end
end
