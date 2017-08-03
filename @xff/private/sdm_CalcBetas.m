function [betas, irtc, ptc, se, weights, rsm, tcd] = sdm_CalcBetas(xo, tcd, dim)
% SDM::CalcBetas  - perform GLM calculcation
%
% FORMAT:       [betas, irtc, ptc, se, w, rsm, res] = sdm.CalcBetas(tc [, tdim])
%
% Input fields:
%
%       tc          time course data (numeric data or FMR/VTC/MTC)
%       tdim        required for numeric tc data, temporal dimension
%                   alternatively, a struct can be given with options
%        .drop      time point indices to drop prior to regression
%        .mask      if tc is a xff object, apply corresponding mask
%        .maxiter   maximum number of iterations for robust regression
%        .pbar      progress bar object (otherwise created)
%        .prange    progress range (default: [0, 1])
%        .regdiff   regress discreet derivative instead
%        .reshpf    residual high-pass filter (default: 0.02Hz)
%        .reslpf    residual low-pass filter (default: 0.1Hz)
%        .resrmgsig remove global signal from residual (true; automask)
%        .robust    perform robust regression instead of OLS
%        .rwreuse   attempt to re-use robust weights (only for objects)
%        .singlew   single vector of weights (must match tc in size)
%        .sngtpool  pool remaining single-trials (lower model DF, false)
%        .tdim      temporal dimension (otherwise auto-detect)
%        .thresh    value, which mean(tc) must be greater as (default: 0)
%        .tmaps     compute t-stats from betas (e.g. for RSA, false)
%        .trans     either of {'none'}, 'psc', or 'z'
%        .tssm      either SSM or TSM object with mappings
%        .writeres  write residual TC, e.g. '%s_res.vtc' (only for objects)
%
% Output fields:
%
%       betas       betas, in time course data dimension
%       irtc        inverse design matrix
%       ptc         predicted tc
%       se          standard error
%       w           robust-regression weights
%       rsm         1x2 cell array with residual smoothness estimate
%       res         residual time course
%
% Note: if no constant term is in the matrix, it will be added as the
%       last column
%       also, for robust regression, the standard error, if requested
%       is adapted to the nominal d.f.
%
% Using: calcbetas, findfirst, fitrobustbisquare_img, flexinterpn, glmtstat,
%        lsqueeze, meannoinfnan, poolnonsingletrial, psctrans, resestsmooth,
%        robustt, smoothkern, tempfilter, varc, ztrans.

% Version:  v1.1
% Build:    16093014
% Date:     Sep-30 2016, 2:38 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2015, 2016, Jochen Weber
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

% importing from neuroelf library into workspace
using(neuroelf, {'calcbetas', 'findfirst', 'fitrobustbisquare_img', 'flexinterpn', ...
    'glmtstat', 'lsqueeze', 'meannoinfnan', 'poolnonsingletrial', 'psctrans', ...
    'resestsmooth', 'robustt', 'smoothkern', 'tempfilter', 'varc', 'ztrans'});

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'sdm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
nout = nargout;

% get design matrix right
bc = xo.C;
sdm = bc.SDMMatrix;

% further arguments
dopbar = true;
drop = [];
mask = [];
maxiter = 30;
mypbar = true;
pbar = [];
prange = [0, 1];
regdiff = false;
reshpf = 0.02;
reslpf = 0.1;
resrmgsig = true;
robust = false;
rwreuse = false;
singlew = [];
sngtpool = false;
thresh = 0;
tmaps = false;
trans = 'n';
tssm = [];
writeres = '';
if nargin > 2 && isstruct(dim) && numel(dim) == 1
    if isfield(dim, 'drop') && isa(dim.drop, 'double') && ~isempty(dim.drop) && ...
       ~any(isinf(dim.drop(:)) | isnan(dim.drop(:)) | dim.drop(:) < 1)
        drop = dim.drop(:)';
    end
    if isfield(dim, 'mask') && islogical(dim.mask) && ~isempty(dim.mask)
        mask = dim.mask;
    end
    if isfield(dim, 'maxiter') && isa(dim.maxiter, 'double') && numel(dim.maxiter) == 1 && ...
       ~isinf(dim.maxiter) && ~isnan(dim.maxiter) && dim.maxiter > 1
        maxiter = min(100, round(dim.maxiter));
    end
    if isfield(dim, 'pbar') && isempty(dim.pbar)
        dopbar = false;
    end
    if isfield(dim, 'pbar') && numel(dim.pbar) == 1 && ...
        any(strcmpi(class(dim.pbar), {'xfigure', 'xprogress'}))
        pbar = dim.pbar;
        mypbar = false;
    end
    if isfield(dim, 'prange') && numel(dim.prange) == 2 && isa(dim.prange, 'double') && ...
       ~any(isinf(dim.prange) | isnan(dim.prange) | dim.prange < 0 | dim.prange > 1) && ...
        dim.prange(2) > dim.prange(1)
        prange = dim.prange;
    end
    if isfield(dim, 'regdiff') && islogical(dim.regdiff) && numel(dim.regdiff) == 1
        regdiff = dim.regdiff;
    end
    if isfield(dim, 'reshpf') && isa(dim.reshpf, 'double') && numel(dim.reshpf) == 1 && ...
       ~isnan(dim.reshpf) && dim.reshpf < 0.1
        reshpf = dim.reshpf;
    end
    if isfield(dim, 'reslpf') && isa(dim.reslpf, 'double') && numel(dim.reslpf) == 1 && ...
       ~isnan(dim.reslpf) && dim.reslpf > 0.05
        reslpf = dim.reslpf;
    end
    if isfield(dim, 'resrmgsig') && islogical(dim.resrmgsig) && numel(dim.resrmgsig) == 1
        resrmgsig = dim.resrmgsig;
    end
    if isfield(dim, 'robust') && islogical(dim.robust) && numel(dim.robust) == 1
        robust = dim.robust;
    end
    if isfield(dim, 'rwreuse') && islogical(dim.rwreuse) && numel(dim.rwreuse) == 1
        rwreuse = dim.rwreuse;
    end
    if isfield(dim, 'singlew') && isa(dim.singlew, 'double') && numel(dim.singlew) == size(sdm, 1)
        singlew = dim.singlew(:);
        singlew(isinf(singlew) | isnan(singlew) | singlew < 0) = 0;
        singlew(singlew > 1) = 1;
    end
    if isfield(dim, 'sngtpool') && islogical(dim.sngtpool) && numel(dim.sngtpool) == 1
        sngtpool = dim.sngtpool;
    end
    if isfield(dim, 'thresh') && isa(dim.thresh, 'double') && numel(dim.thresh) == 1 && ...
       ~isinf(dim.thresh) && ~isnan(dim.thresh) && dim.thresh > 0
        thresh = dim.thresh;
    end
    if isfield(dim, 'tmaps') && islogical(dim.tmaps) && numel(dim.tmaps) == 1
        tmaps = dim.tmaps;
    end
    if isfield(dim, 'trans') && ischar(dim.trans) && ~isempty(dim.trans) && ...
        any('pz' == lower(dim.trans(1)))
        trans = lower(dim.trans(1));
    end
    if isfield(dim, 'tssm') && numel(dim.tssm) == 1 && ...
       (xffisobject(dim.tssm, true, 'ssm') || xffisobject(dim.tssm, true, 'tsm')) && ...
        numel(tcd) == 1 && xffisobject(tcd, true, 'mtc')
        tssm = dim.tssm;
    end
    if isfield(dim, 'tdim') && isa(dim.tdim, 'double') && numel(dim.tdim) == 1 && ...
        any((1:5) == dim.tdim)
        dim = dim.tdim;
    end
    if isfield(dim, 'writeres') && ischar(dim.writeres) && numel(dim.writeres) > 5 && ...
        sum(dim.writeres(:)' == '%') == 1 && numel(tcd) == 1 && ...
        xffisobject(tcd, true, {'fmr', 'vtc', 'mtc'}) && ...
        strcmpi(dim.writeres(end-3:end), ['.' aft_Filetype(tcd)])
        writeres = dim.writeres(:)';
    else
        writeres = '';
    end
end
if nargin < 3 || ~isa(dim, 'double') || numel(dim) ~= 1 || isinf(dim) || isnan(dim) || ...
    fix(dim) ~= dim || dim < 1 || dim > 5
    rdim = [];
else
    rdim = dim;
end
if ~dopbar
    pbar = [];
end

% add all-1 column for confound
mrtc = mean(sdm);
vrtc = var(sdm);
if ~any(mrtc ~= 0 & vrtc == 0) && size(sdm, 2) < size(sdm, 1)
    sdm(:, end + 1) = 1;
end

% number of rows and betas
numrows = size(sdm, 1);
numbets = size(sdm, 2);

% use derivatives?
if regdiff
    sdmco = (sum(abs(sdm)) > 0 & var(sdm) == 0);
    sdm = diff(sdm);
    sdm(:, sdmco) = 1;
end

% get betas that can be ~= 0
usebets = (sum(abs(sdm)) ~= 0);

% create single trial arguments
stpu = find(usebets);
nstpu = numel(stpu);
if sngtpool

    % if less than 3 _00X condition
    pn = bc.PredictorNames(:);
    stpn = pn(stpu, 1);
    stpm = (~cellfun('isempty', regexp(stpn, '_T\d+$')));
    if sum(stpm) < 3

        % disable, as those two need to be separated anyway
        sngtpool = false;
    else

        % find conditions that need to be "re-pooled"
        stpc = find(stpm);
        nstpc = numel(stpc);

        % and then, for each unique condition
        stpx = cell(numel(stpc), 3);
        for pc = 1:numel(stpc)

            % combine SDM accordingly
            stpx{pc, 1} = poolnonsingletrial(sdm(:, stpu), stpn, stpc(pc));
        end
    end
end

% check second argument -> xff object
stcwo = {[]};
tcc = [];
if numel(tcd) == 1 && xffisobject(tcd, true)

    % get struct version
    stcd = tcd.C;
    stcf = tcd.F;

    % write residuals?
    if ~isempty(writeres)

        % create copy of object
        tcc = aft_CopyObject(tcd);

        % filename
        resorig = stcf;
        writeres = strrep(writeres, '%s', regexprep(resorig, '\.([a-zA-Z])+$', ''));
        tcc.F = '';

        % also patch Prefix for FMRs
        if strcmpi(writeres(end-2:end), 'fmr')
            [tccp, tcc.C.Prefix] = fileparts(writeres);
        end

        % make sure ptc is computed
        if nout < 3
            nout = 3;
        end
    end

    % robust weights
    if ~isempty(stcf) && stcf(end-3) == '.' && robust && rwreuse
        stcwf = [stcf(1:end-4) '_robweights' stcf(end-3:end)];
        try
            stcwo{1} = xff(stcwf);
            if ~xffisobject(stcwo{1}, true, tcd.S.Extensions{1})
                error('neuroelf:xff:internalError', 'Invalid xff object for weights.');
            end
            stcwc = stcwo{1}.C;
        catch xfferror
            neuroelf_lasterr(xfferror);
            clearxffobjects(stcwo);
            stcwo{1} = [];
        end
    end

    % get filetype
    ftype = lower(tcd.S.Extensions{1});

    % switch over filetype
    switch lower(ftype)

        % FMR
        case {'fmr'}

            % rdim MUST be 3
            rdim = 3;

            % STCs loaded?
            if isempty(stcd.Slice) || ~isstruct(stcd.Slice) || ...
               ~isfield(stcd.Slice, 'STCData')
                try
                    fmr_LoadSTCs(tcd);
                    stcd = tcd.C;
                catch xfferror
                    rethrow(xfferror);
                end
            end

            % check number of volumes
            nvol = numel(stcd.Slice(1).STCData(1, 1, :, 1));
            if nvol ~= numrows
                error('neuroelf:xff:internalError', 'Invalid number of time points in FMR/STC.');
            end

            % get STC data
            switch (stcd.DataStorageFormat)

                % old FMR format (with multiple STCs)
                case {1}
                    nslc = numel(stcd.Slice);
                    tcd = stcd.Slice(1).STCData(:, :, :);
                    tcd(end, end, end, nslc) = 0;
                    for sc = 2:nslc
                        tcd(:, :, :, sc) = stcd.Slice(sc).STCData(:, :, :);
                    end

                % new FMR format (one STC file)
                case {2}
                    tcd = stcd.Slice.STCData(:, :, :, :);

                % unknown FMR format (?)
                otherwise
                    error('neuroelf:xff:invalidObject', 'Unsupported FMR DataStorageFormat.');
            end

            % test robust weights
            if ~isempty(stcwo{1})
                if ~isequal(size(stcwc.Slice.STCData), size(tcd)) || ...
                   ~isa(stcwc.Slice.STCData(1), 'single')
                    clearxffobjects(stcwo);
                    stcwo{1} = [];
                end
            end

        % MTC
        case {'mtc'}

            % rdim MUST be 1
            rdim = 1;

            % check MTC
            if size(stcd.MTCData, 1) ~= numrows
                error('neuroelf:xff:internalError', 'Invalid number of time points in MTC.');
            end

            % get MTC data
            tcd = stcd.MTCData(:, :);

            % apply SSM/TSM ?
            if ~isempty(tssm)

                % get content
                tssm = tssm.C;

                % SSM
                if ~isfield(tssm, 'NrOfSourceTriangles')

                    % get unique indices
                    [ssmui, ssmua, ssmub] = unique(tssm.SourceOfTarget(:));

                    % access data (uniquely first!)
                    tcd = tcd(:, ssmui);
                    tcd = tcd(:, ssmub);

                % TSM
                else

                    % NO SRF OBJECT DEFINED!! MDM NEEDS ADDITIONAL WORK!!

                    % get elements
                    srv = srfc.TriangleVertex(:, :);
                    srv = srv(tssm.SourceTriangleOfTarget(:), :);
                    srw2 = tssm.TriangleEdgeLengths(:, :)';
                    srw1 = 1 - sum(srw2, 1);
                    srw3 = srw2(2, :);
                    srw2 = srw2(1, :);

                    % re-sample data
                    otcd = tcd(:, :);
                    tcd = single(0);
                    tcd(size(otcd, 1), numel(srw1)) = 0;
                    for vfrom = 1:size(otcd, 1)
                        tcd(vfrom, :) = srw1 .* otcd(vfrom, srv(:, 1)) + ...
                            srw2 .* otcd(vfrom, srv(:, 2)) + ...
                            srw3 .* otcd(vfrom, srv(:, 3));
                    end
                end
            end

        % VTC
        case {'vtc'}

            % rdim MUST be 1
            rdim = 1;

            % check VTC
            if size(stcd.VTCData, 1) ~= numrows
                error('neuroelf:xff:internalError', 'Invalid number of volumes in VTC.');
            end

            % get VTC data
            tcd = stcd.VTCData(:, :, :, :);

        % otherwise
        otherwise
            error('neuroelf:xff:typeNotSupported', ...
                'xff object of type %s not supported.', ftype);
    end

    % update drop
    if isfield(stcd.RunTimeVars, 'Discard')
        drop = stcd.RunTimeVars.Discard(:)';
    end

% for other numerics
elseif isnumeric(tcd) && isa(tcd, 'transio')

    % create double
    tcd = resolve(tcd);
end

% check dim
if isempty(rdim)
    rdim = findfirst(size(tcd) == size(sdm, 1));
    if isempty(rdim)
        error('neuroelf:xff:internalError', ...
            'Invalid data matrix size vs. given data.');
    end
elseif rdim < 1 || rdim > length(size(tcd)) || size(tcd, rdim) ~= numrows
    error('neuroelf:xff:internalError', ...
        'Invalid data matrix vs. selected dimension.');
end
if rdim > 1
    neword = [rdim, 1:length(size(tcd))];
    newodr = find(neword == rdim);
    neword(newodr(2)) = [];
    tcd = permute(tcd, neword);
    [newsrt, oldord] = sort(neword);
else
    oldord = [];
end
drop = unique(round(min(numrows, max(1, drop))));
keep = 1:numrows;
if ~isempty(drop)
    keep(drop) = [];
    if numel(keep) <= numbets
        error('neuroelf:xff:BadArgument', 'Too many dropped time points!');
    end
end

% other things required for tmaps!
if ~robust && tmaps
    nout = max(4, nout);
end

% reshaping data to comply
tcds = size(tcd);
numvox = prod(tcds(2:end));
tcd = reshape(tcd, [tcds(1), numvox]);
tcds(1) = [];
tcdrs = tcds;
if length(tcdrs) < 2
    tcdrs(2) = 1;
end

% pre-allocate output
betas = zeros(numvox, numbets);
if nout > 3
    se = zeros(tcdrs);
end
if ~robust
    if nout > 2
        ptc = zeros(numrows, numvox);
        if ~isempty(drop)
            ptc(drop, :) = tcd(drop, :);
        end
    end
    if nout > 4
        weights = ones(numrows, 1);
        if ~isempty(drop)
            weights(drop) = 0;
        end
    end
else
    if nout > 2
        ptc = zeros(numrows, numvox);
        if ~isempty(drop)
            ptc(drop, :) = tcd(drop, :);
        end
        weights = zeros(numrows, numvox);
    end

    % for progress
    copts = struct('maxiter', maxiter, 'pbar', pbar, 'prange', prange);
end

% adapt numrows if required
if ~isempty(drop)
    numrows = numrows - numel(drop);
end

% reject incorrect mask
if ~isempty(mask) && numel(mask) ~= numvox
    mask = [];
else
    mask = mask(:)';
end

% implicit masking
if trans ~= 'n'
    if ~isempty(mask)
        mask(any(isnan(tcd) | isinf(tcd), 1)) = false;
    else
        tcd(:, any(isnan(tcd) | isinf(tcd), 1)) = 0;
    end
end

% thresholding
if thresh > 0

    % relative thresholding for values up to (and including) 2
    if thresh <= 2

        % compute threshold as relative mean of all voxels
        maski = mean(tcd);
        maski(isinf(maski) | isnan(maski)) = 0;
        thresh = thresh * mean(maski(:));
    end

    % then depending on existing mask
    if ~isempty(mask)
        mask = mask & (sum(tcd) >= (numrows * thresh));
    else
        mask = (sum(tcd) >= (numrows * thresh));
    end
end

% mask indices
if ~isempty(mask)
    maski = find(mask(:));
    numvox = numel(maski);
else
    maski = 1:numvox;
end

% for robust regression
if robust

    % precompute h and perm and iXX
    [robiXX, robh, robperm] = qr(sdm(keep, usebets), 0);
    robh = {robh, robperm};
    robiXX = inv(sdm(keep, usebets)' * sdm(keep, usebets));
end

% progress bar
if numvox >= 128 && dopbar
    pmin = prange(1);
    pmax = prange(2);
    pdif = eps + pmax - pmin;
    if mypbar
        try
            pbar = xprogress;
            xprogress(pbar, 'setposition', [80, 200, 640, 36]);
            xprogress(pbar, 'settitle', 'Calculating betas ...');
            xprogress(pbar, pmin, sprintf('Estimating %d samples...', ...
                numvox), 'visible', pmin, pmax);
            if robust
                copts.pbar = pbar;
            end
        catch xfferror
            neuroelf_lasterr(xfferror);
            pbar = [];
        end
    else
        pbar.Progress(pmin);
    end
else
    mypbar = false;
    pbar = [];
end

% compute necessary additions
if sngtpool
    for pc = 1:size(stpx, 1)
        stpx{pc, 1} = stpx{pc, 1}(keep, :);
        if ~robust
            stpx{pc, 2} = stpx{pc, 1}';
            stpx{pc, 3} = inv(stpx{pc, 2} * stpx{pc, 1});
            stpx{pc, 2} = stpx{pc, 3} * stpx{pc, 2};
        end
    end
end
if tmaps
    tmapc = eye(size(sdm, 2));
end

% looping over voxels
vstep = max(128, floor(1e7 / (numrows * numbets))) - 1;
vfrom = 1;
while vfrom <= numvox

    % next block
    vnext = min(vfrom + vstep, numvox);
    maskn = maski(vfrom:vnext);

    % progress
    if ~isempty(pbar)
        pbar.Progress(pmin + pdif * vfrom / numvox);
    end

    % get time course snippet
    if isempty(drop)
        tcs = double(tcd(:, maskn));
    else
        tcs = double(tcd(keep, maskn));
    end

    % transformation
    if trans == 'p'
        tcs = psctrans(tcs);
    elseif trans == 'z'
        tcsv = (varc(tcs) > sqrt(eps));
        if any(tcsv)
            tcs(:, tcsv) = ztrans(tcs(:, tcsv));
        end
        tcs(:, ~tcsv) = 0;
    end

    % compute discreet derivative first?
    if regdiff
        tcs = diff(tcs);
    end

    % robust
    if robust
        if ~isempty(pbar)
            copts.prange = pmin + (pdif / numvox) .* [vfrom, vnext];
        end
        if nout < 3 && ~tmaps
            betas(maskn, usebets) = fitrobustbisquare_img( ...
                sdm(keep, usebets), tcs, [], [], copts, robh, robiXX);
            if sngtpool
                for pc = 1:nstpc
                    sbetas = fitrobustbisquare_img(stpx{pc, 1}, tcs, [], [], copts);
                    betas(maskn, stpc(pc)) = sbetas(:, 1);
                end
            end
        else
            if isempty(drop)
                [betas(maskn, usebets), weights(:, maskn)] = fitrobustbisquare_img( ...
                    sdm(:, usebets), tcs, [], [], copts, robh, robiXX);
                if sngtpool
                    for pc = 1:nstpc
                        sbetas = fitrobustbisquare_img(stpx{pc, 1}, tcs, [], [], copts);
                        betas(maskn, stpc(pc)) = sbetas(:, 1);
                    end
                end
                ptc(:, maskn) = (sdm(:, usebets) * betas(maskn, usebets)') .* ...
                    weights(:, maskn) + (tcs .* (1 - weights(:, maskn)));
            else
                [betas(maskn, usebets), weights(keep, maskn)] = fitrobustbisquare_img( ...
                    sdm(keep, usebets), tcs, [], [], copts, robh, robiXX);
                if sngtpool
                    for pc = 1:nstpc
                        sbetas = fitrobustbisquare_img(stpx{pc, 1}, tcs, [], [], copts);
                        betas(maskn, stpc(pc)) = sbetas(:, 1);
                    end
                end
                ptc(keep, maskn) = (sdm(keep, usebets) * betas(maskn, usebets)') .* ...
                    weights(keep, maskn) + (tcs .* (1 - weights(keep, maskn)));
            end
            if tmaps
                sbetas = betas(maskn, usebets);
                for pc = 1:nstpu
                    if ~all(sdm(keep, stpu(pc)) == 1)
                        rtcv = zeros(1, nstpu);
                        rtcv(pc) = 1;
                        sbetas(:, pc) = robustt(sdm(keep, usebets), tcs, ...
                            betas(maskn, usebets), weights(keep, maskn), rtcv);
                    end
                end
                betas(maskn, usebets) = sbetas;
            end
            if nout > 3
                if isempty(drop)
                    se(maskn) = squeeze(sqrt((1 / (numrows - numbets)) .* ...
                       (sum(weights(:, maskn)) - 1)) .* std(tcs - ptc(:, maskn), 0));
                else
                    se(maskn) = squeeze(sqrt((1 / (numrows - numbets)) .* ...
                       (sum(weights(keep, maskn)) - 1)) .* std(tcs - ptc(keep, maskn), 0));
                end
            end
        end

    % or OLS
    else
        switch (nout)
            case {1, 2}
                betas(maskn, usebets) = calcbetas(sdm(keep, usebets), tcs, 1, [], singlew);
                if sngtpool
                    for pc = 1:nstpc
                        sbetas = stpx{pc, 2} * tcs;
                        betas(maskn, stpc(pc)) = sbetas(1, :)';
                    end
                end
            case {3}
                if isempty(drop)
                    [betas(maskn, usebets), irtc, ptc(:, maskn)] = ...
                        calcbetas(sdm(:, usebets), tcs, 1, [], singlew);
                else
                    [betas(maskn, usebets), irtc, ptc(keep, maskn)] = ...
                        calcbetas(sdm(keep, usebets), tcs, 1, [], singlew);
                end
                if sngtpool
                    for pc = 1:nstpc
                        sbetas = stpx{pc, 2} * tcs;
                        betas(maskn, stpc(pc)) = sbetas(1, :)';
                    end
                    ptc(:, maskn) = sdm(:, usebets) * betas(maskn, usebets)';
                end
            otherwise
                [betas(maskn, usebets), irtc, ptc(keep, maskn), se(maskn)] = ...
                    calcbetas(sdm(keep, usebets), tcs, 1, [], singlew);
                if sngtpool
                    for pc = 1:nstpc
                        sbetas = stpx{pc, 2} * tcs;
                        betas(maskn, stpc(pc)) = sbetas(1, :)';
                    end
                    ptc(:, maskn) = sdm(:, usebets) * betas(maskn, usebets)';
                    se(maskn) = sqrt((numrows - 1) / (numrows - numbets)) .* ...
                        lsqueeze(std(tcs - ptc(:, maskn), 0));
                end
        end
        if tmaps
            for pc = 1:nstpu
                betas(maskn, stpu(pc)) = glmtstat(tmapc(stpu(pc), usebets), ...
                    betas(maskn, usebets), irtc, se(maskn));
            end
        end
    end
    vfrom = vfrom + vstep + 1;
end
if nout > 3
    se(isinf(se) | isnan(se)) = Inf;
end

% also build standard irtc ?
if nout > 1
    if sum(usebets) ~= numbets
        sdm(:, ~usebets) = 0;
        irtc = pinv(sdm' * sdm);
        irtc(~usebets, :) = 0;
        irtc(:, ~usebets) = 0;
    elseif exist('irtc', 'var') == 0
        irtc = pinv(sdm' * sdm);
    end
end

% reshape betas
betas = reshape(betas, [tcds, numbets]);

% compute residual
if ~isempty(tcc) || nargout > 5
    if trans == 'p'
        tcd = psctrans(double(tcd));
    elseif trans == 'z'
        tcd = ztrans(double(tcd));
    end
    tcd = double(tcd) - ptc;

    % compute residual image
    rsm = cell(1, 2);
    if nargout > 5 && numel(tcdrs) == 3
        [rsm{1}, rsm{2}] = resestsmooth( ...
            reshape(tcd, [size(tcd, 1), tcdrs]), [1, 1, 1], struct('tdim', 1));
    end
end

% write out residual
if ~isempty(tcc)

    % progress
    if ~isempty(pbar)
        pbar.Progress(pmin + pdif * 0.99, 'Writing residuals...');
    end

    % set additional fields
    tcc.C.RunTimeVars.AutoSave = true;
    tcc.C.RunTimeVars.IsResidual = true;
    tcc.C.RunTimeVars.MaskIndex = maski;
    tcc.C.RunTimeVars.OriginalFile = resorig;
    tcc.C.RunTimeVars.PredictorNames = bc.PredictorNames;
    tcc.C.RunTimeVars.ResHPF = reshpf;
    tcc.C.RunTimeVars.ResLPF = reslpf;
    tcc.C.RunTimeVars.ResRemoveGlobSig = true;
    tcc.C.RunTimeVars.Robust = robust;
    tcc.C.RunTimeVars.SDMMatrix = sdm;
    tcc.C.RunTimeVars.SDMKeptPoints = keep;

    % temporally filter residual
    reshpf = reshpf * (tcc.C.TR(1) / 1000);
    reslpf = reslpf * (tcc.C.TR(1) / 1000);
    if reslpf < 1

        % compute smoothing kernel
        smk = smoothkern(0.5 / reslpf, 0, false, 'lanczos8');

        % reduce size
        smk(smk < max(smk)*sqrt(eps)) = [];

        % smooth data
        tcd = single(flexinterpn(tcd, [Inf, Inf; 1, 1; 1, 1; size(tcd)], ...
            {smk, [0; 1; 0]}, {1, 1}, 0));
    end
    if reshpf > (1.5 / numrows)
        tcd = single(tempfilter(tcd, struct('tempdct', 1 / reshpf)));
    end
    if resrmgsig
        gsigmean = mean(tcd(:));
        if gsigmean >= 50
            gsigmask = (mean(tcd, 1) > gsigmean);
        else
            gsigmask = ':';
        end
        gsig = double(meannoinfnan(tcd(:, gsigmask), 2, true));
        tcc.C.RunTimeVars.ResGlobalSignal = gsig;
        gsig = gsig - mean(gsig);
        gsigb = calcbetas([gsig, ones(size(gsig))], double(tcd));
        tcd = tcd - gsig * gsigb(:, 1)';
    end

    % depending on type
    switch lower(tcc.S.Extensions{1})
        case {'fmr'}
            tcdosz = size(tcd);
            tcd = reshape(tcd, [size(tcd, 1), tcc.C.ResolutionX, tcc.C.ResolutionY, tcc.C.NrOfSlices]);
            tcd = permute(tcd, [2, 3, 1, 4]);
            tcc.C.FileVersion = 6;
            tcc.C.DataType = 2;
            if tcc.C.DataStorageFormat == 1
                for tcsc = 1:numel(tcc.Slice)
                    tcc.C.Slice(tcsc).STCData = single(tcd(:, :, :, tcsc));
                end
            else
                tcc.C.Slice.STCData = single(tcd);
            end
            tcd = permute(tcd, [3, 1, 2, 4]);
            tcd = reshape(tcd, tcdosz);
        case {'mtc'}
        case {'vtc'}

            % set data
            tcc.C.DataType = 2;
            tcc.C.VTCData = reshape(tcd, size(tcc.C.VTCData));
            [rsmfwm, rsmfwi] = resestsmooth( ...
                reshape(tcd, [size(tcd, 1), tcdrs]), [1, 1, 1], struct('tdim', 1));
            tcc.C.RunTimeVars.FWHMResEst = rsmfwm;
            tcc.C.RunTimeVars.FWHMResImg = single(rsmfwi);
    end

    % try saving, then clear
    try
        aft_SaveAs(tcc, writeres);
        aft_SaveRunTimeVars(tcc);
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
    try
        delete(tcc);
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end

% reshape ptc if necessary
if nout > 2
    ptc = reshape(ptc, [numrows + numel(drop), tcds]);
    if ~isempty(oldord)
        if nout > 4 && size(weights, 2) > 1
            weights = reshape(weights, size(ptc));
            weights = permute(weights, oldord);
        end
        ptc = permute(ptc, oldord);
    end
    if nout > 3
        se = reshape(se, tcdrs);
    end
end

% progress bar
if mypbar && ~isempty(pbar)
    closebar(pbar);
end
