function [varargout] = alphasim(ddim, opts)
% alphasim  - simulate noise data to estimate cluster threshold
%
% FORMAT:       [at = ] alphasim(ddim [, opts])
%
% Input fields:
%
%       ddim        data dimension (1x3 integer values)
%       opts        optional settings
%        .clconn    connectivity of clusters, ('face', {'edge'}, 'vertex')
%        .conj      conjunction simulation (1x1 double, number of maps)
%        .fftconv   boolean flag, use FFT convolution (default: true)
%        .fwhm      FWHM kernel sizes (default: [4, 4, 4])
%        .mask      boolean mask (size must be == ddim!, default: none)
%        .niter     number of iterations, default: 1000 (5000 with model)
%        .pbar      either xprogress or xfigure:XProgress object
%        .reglink   link regressions (if number of subjects match, true)
%        .regmaps   regression maps (e.g. betas, contrasts)
%        .regmodel  regression model (if empty, all-1s/mean-test model)
%        .regmodsc  test contrast (multiple rows, conjunction or F-test)
%        .regrank   rank-transform data before useing regression
%        .srf       optional surface (perform surface-based simulation)
%        .srfsmp    surface sampling (from, step, to, along normals,
%                   default: [-3, 1, 1])
%        .srftrf    transformation required to sample surface coordinates
%                   derived from bvcoordconv
%        .stype     1x1 or 1x2 statistics type, default: [1, 2], meaning
%                   that one tail of a two-tailed statistic is taken
%                   a single 1 is one tail of a one-tailed statistic (F)
%                   a single 2 is both tails of a two-tailed statistic (t)
%        .tdf       simulate actual t-stats (for 2-tailed stats only)
%        .thr       applied (raw) threshold(s), default: p<0.001
%        .zshift    shift normal distribution by this Z value (default: 0)
%
% Output fields:
%
%       at          optional output table
%
% Note: other than AFNI's AlphaSim, the data is considered to be
%       iso-voxel for the default kernel, but that can be altered
%       accordingly by changing the kernel!
%
%       to simulate specific regression results, both options, .regmaps
%       .regmodel must be set; if only .regmaps is given, random numbers
%       (using randn) will be generated instead of permuting the predictor

% Version:  v1.1
% Build:    16101016
% Date:     Oct-10 2016, 4:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% argument check
if nargin < 1 || ~isa(ddim, 'double') || numel(ddim) > 3 || ...
    any(isinf(ddim) | isnan(ddim) | ddim < 1 | ddim > 256)
    error('neuroelf:general:badArgument', 'Missing or invalid ddim argument.');
else
    ddim = round(ddim);
end
ndim = numel(ddim) + 1;
if nargin < 2 || isempty(opts)
    opts = struct;
elseif ~isstruct(opts) || numel(opts) ~= 1
    error('neuroelf:general:badArgument', 'Invalid opts argument.');
end
if ~isfield(opts, 'clconn')
    opts.clconn = 'edge';
elseif ~ischar(opts.clconn) || ...
   ~any(strcmpi(opts.clconn(:)', {'edge', 'face', 'vertex'}))
    error('neuroelf:general:badArgument', 'Invalid opts.clconn field.');
else
    opts.clconn = lower(opts.clconn(:)');
end
switch (opts.clconn(1))
    case {'e'}
        clconn = 2;
    case {'f'}
        clconn = 1;
    otherwise
        clconn = 3;
end
givenmaps = false;
multimaps = false;
if ~isfield(opts, 'regmaps') || isempty(opts.regmaps) || ...
   ((~isnumeric(opts.regmaps) || ~any(ndims(opts.regmaps) ~= [2, 4]) || ...
     size(opts.regmaps, ndims(opts.regmaps)) < 3) && ...
    (~iscell(opts.regmaps) || ~any(ndims(opts.regmaps{1}) ~= [2, 4]) || ...
     size(opts.regmaps{1}, ndims(opts.regmaps{1})) < 3))
    opts.regmaps = {[]};
    regnsub = 0;
else
    givenmaps = true;
    if isnumeric(opts.regmaps)
        if ~isa(opts.regmaps, 'double')
            opts.regmaps = {double(opts.regmaps)};
        else
            opts.regmaps = {opts.regmaps};
        end
    elseif numel(opts.regmaps) > 1
        multimaps = true;
        opts.regmaps = opts.regmaps(:)';
    end
    ddim = size(opts.regmaps{1});
    ndim = numel(ddim);
    regnsub = ddim(ndim);
    ddim(ndim) = [];
    pddim = [1, prod(ddim)];
    for rc = 1:numel(opts.regmaps)
        tdim = size(opts.regmaps{rc});
        if ~isequal(tdim(1:end-1), ddim)
            error('neuroelf:general:badArgument', 'Data must match across regressions.');
        end
        regnsub(rc) = tdim(ndim);
        if ~isa(opts.regmaps{rc}, 'double')
            opts.regmaps{rc} = double(opts.regmaps{rc});
        end
        opts.regmaps{rc} = reshape(opts.regmaps{rc}, pddim(2), regnsub(rc))';
    end
end
nregmaps = numel(opts.regmaps);
ndim = ndim - 1;
if any(regnsub < 3) || ~isfield(opts, 'regmodel') || isempty(opts.regmodel) || ...
   ((~isa(opts.regmodel, 'double') || size(opts.regmodel, 1) ~= regnsub(1) || ...
     any(isinf(opts.regmodel(:)) | isnan(opts.regmodel(:)))) && ...
    (~iscell(opts.regmodel) || ~isa(opts.regmodel{1}, 'double') || ...
     size(opts.regmodel{1}, 1) ~= regnsub(1) || ...
     any(isinf(opts.regmodel{1}(:)) | isnan(opts.regmodel{1}(:)))))
    opts.regmodel = cell(1, nregmaps);
    for rc = 1:nregmaps
        opts.regmodel{rc} = ones(regnsub(rc), 1);
    end
    regnreg = regnsub - 1;
else
    if isa(opts.regmodel, 'double')
        opts.regmodel = repmat({opts.regmodel}, 1, nregmaps);
    end
    regnreg = zeros(1, nregmaps);
    for rc = 1:nregmaps
        if size(opts.regmodel{rc}, 1) ~= regnsub(rc)
            error('neuroelf:general:badArgument', 'Invalid data/model combination.');
        end
        if ~any(all(opts.regmodel{rc} == 1))
            opts.regmodel{rc}(:, end+1) = 1;
        end
        regnreg(rc) = regnsub(rc) - size(opts.regmodel{rc}, 2);
    end
end
regnsfc = sqrt((regnsub - 1) ./ regnreg);
modttest = true(size(opts.regmodel));
if ~isfield(opts, 'regmodsc') || isempty(opts.regmodsc) || ...
   ((~isa(opts.regmodsc, 'double') || ...
     ~any(size(opts.regmodsc, 2) == (size(opts.regmodel{1}, 2) - [0, 1]))) && ...
    (~iscell(opts.regmodsc) || ~isa(opts.regmodsc{1}, 'double') || ...
     ~any(size(opts.regmodsc{1}, 2) == (size(opts.regmodel{1}, 2) - [0, 1]))))
    opts.regmodsc = repmat({[]}, 1, nregmaps);
elseif isa(opts.regmodsc, 'double')
    opts.regmodsc = repmat({opts.regmodsc}, 1, nregmaps);
end
redmod = opts.regmodel;
redixx = opts.regmodel;
for rc = 1:nregmaps
    if isempty(opts.regmodsc{rc})
        opts.regmodsc{rc} = [1, zeros(1, size(opts.regmodel{rc}, 2) - 1)];
    elseif size(opts.regmodsc{rc}, 1) > 1
        if all(opts.regmodsc{rc}(:) == 0 | opts.regmodsc{rc}(:) == 1)
            modttest(rc) = false;
            opts.regmodsc{rc} = double(any(opts.regmodsc > 0, 1));
            redmod{rc} = opts.regmodel{rc}(:, opts.regmodsc{rc} == 0);
            if ~isempty(redmod)
                redixx{rc} = inv(redmod{rc}' * redmod{rc});
            end
        else
            opts.conj = size(opts.regmodsc{rc}, 1);
        end
    end
end
if ~isempty(opts.regmodel{1})
    modixx = opts.regmodel;
    newmod = opts.regmodel;
    for rc = 1:nregmaps
        modixx{rc} = inv(opts.regmodel{rc}' * opts.regmodel{rc});
    end
end
if ~isfield(opts, 'reglink') || numel(opts.reglink) ~= 1 || ~islogical(opts.reglink)
    opts.reglink = all(regnsub == regnsub(1));
end
if ~isfield(opts, 'conj')
    nconj = 1;
elseif ~isa(opts.conj, 'double') || numel(opts.conj) ~= 1 || ...
    isinf(opts.conj) || isnan(opts.conj) || opts.conj < 1 || opts.conj > 5
    error('neuroelf:general:badArgument', 'Invalid opts.conj field.');
else
    nconj = floor(opts.conj);
end
if ~isfield(opts, 'fftconv')
    opts.fftconv = true;
elseif ~islogical(opts.fftconv) || numel(opts.fftconv) ~= 1
    error('neuroelf:general:badArgument', 'Invalid opts.fftconv field.');
end
fftconv = opts.fftconv;
if ~isfield(opts, 'fwhm')
    opts.fwhm = (12 / ndim) .* ones(1, ndim);
elseif ~isa(opts.fwhm, 'double') || ...
    numel(opts.fwhm) ~= ndim || ...
    any(isinf(opts.fwhm) | isnan(opts.fwhm) | opts.fwhm <= 0 | opts.fwhm(:)' > ddim)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid opts.fwhm field.' ...
    );
else
    opts.fwhm = opts.fwhm(:)';
end
if ndim == 3
    kcell = { ...
        smoothkern(opts.fwhm(1), 0), ...
        smoothkern(opts.fwhm(2), 0), ...
        smoothkern(opts.fwhm(3), 0)};
    kern = {zeros(numel(kcell{1}), numel(kcell{2}), numel(kcell{3}))};
    kern = kern(1, [1, 1, 1]);
    kern{1}(:, (numel(kcell{2}) + 1) / 2, (numel(kcell{3}) + 1) / 2) = kcell{1};
    kern{2}((numel(kcell{1}) + 1) / 2, :, (numel(kcell{3}) + 1) / 2) = kcell{2};
    kern{3}((numel(kcell{2}) + 1) / 2, (numel(kcell{2}) + 1) / 2, :) = kcell{3};
end
if ~isfield(opts, 'mask')
    opts.mask = ([] > 0);
elseif ~islogical(opts.mask) || ...
   ~isequal(size(opts.mask), ddim)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid opts.mask field.' ...
    );
end
mask = opts.mask;
if ~isempty(mask)
    summask = sum(mask(:));
    if summask == 0
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid opts.mask field.' ...
        );
    end
    msktxt = sprintf(' in %d-voxel mask', summask);
    if any(regnsub > 0)
        rma = ones(1, ndim + 1);
        for rc = 1:nregmaps
            rma(end) = regnsub(rc);
            opts.regmaps{rc}(~repmat(mask, rma)) = 0;
        end
    end
else
    msktxt = '';
end
if ~isfield(opts, 'niter')
    if isempty(opts.regmodel)
        niter = 1000;
    else
        niter = 5000;
    end
elseif ~isa(opts.niter, 'double') || ...
    numel(opts.niter) ~= 1 || ...
    isinf(opts.niter) || ...
    isnan(opts.niter) || ...
    opts.niter < 1 || ...
    opts.niter > 1e6
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid opts.niter field.' ...
    );
else
    niter = round(opts.niter);
end
if ~isfield(opts, 'pbar') || ...
    numel(opts.pbar) ~= 1 || ...
   ~any(strcmpi(class(opts.pbar), {'xfigure', 'xprogress'}))
    opts.pbar = [];
end
if ~isfield(opts, 'regrank') || ...
   ~islogical(opts.regrank) || ...
    numel(opts.regrank) ~= 1
    opts.regrank = false;
end
if ~isempty(opts.regmodel) && ...
    opts.regrank
    for rc = 1:nregmaps
        opts.regmodel{rc} = ztrans(ranktrans(opts.regmodel{rc}, 1));
        opts.regmodel{rc}(:, all(opts.regmodel{rc}) == 0, 1) = 1;
    end
end
if ~isfield(opts, 'srf') || ...
    numel(opts.srf) ~= 1 || ...
   ~isxff(opts.srf, 'srf')
    opts.srf = [];
end
if ~isfield(opts, 'srfsmp') || ...
   ~isa(opts.srfsmp, 'double') || ...
    numel(opts.srfsmp) ~= 3 || ...
    any(isinf(opts.srfsmp) | isnan(opts.srfsmp) | abs(opts.srfsmp) > 12) || ...
    opts.srfsmp(1) > opts.srfsmp(3) || ...
    isempty(opts.srfsmp(1):opts.srfsmp(2):opts.srfsmp(3))
    opts.srfsmp = -3:1;
else
    opts.srfsmp = opts.srfsmp(1):opts.srfsmp(2):opts.srfsmp(3);
    while numel(opts.srfsmp) > 12
        opts.srfsmp = opts.srfsmp(1:2:end);
    end
end
if ~isempty(opts.srf) && ...
   (~isfield(opts, 'srftrf') || ...
    ~isa(opts.srftrf, 'double') || ...
    ~isequal(size(opts.srftrf), [4, 4]) || ...
     any(isinf(opts.srftrf(:)) | isnan(opts.srftrf(:))) || ...
     any(opts.srftrf(4, 1:3) ~= 0)) && ...
    (ndim == 3 || ...
     ddim(1) ~= opts.srf.NrOfVertices)
    opts.srf = [];
end
if ~isempty(opts.srf)
    tri = opts.srf.TriangleVertex;
    crd = opts.srf.VertexCoordinate;
    try
        [nei, bn, trb] = mesh_trianglestoneighbors(size(crd, 1), tri);
        if ~isempty(bn)
            warning( ...
                'neuroelf:BadSurface', ...
                'Cluster sizes potentially flawed. %d bad neighborhoods!', ...
                numel(bn) ...
            );
        end
        if isempty(nei{end}) || ...
            isempty(trb{end})
            error('BAD_SURFACE');
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error( ...
            'neuroelf:BadSurface', ...
            'Invalid surface, neighborhood references invalid.' ...
        );
    end
    nei = nei(:, 2);
    nrm = opts.srf.VertexNormal;
    tsa = sqrt(sum((crd(tri(:, 1), :) - crd(tri(:, 2), :)) .^ 2, 2));
    tsb = sqrt(sum((crd(tri(:, 1), :) - crd(tri(:, 3), :)) .^ 2, 2));
    tsc = sqrt(sum((crd(tri(:, 2), :) - crd(tri(:, 3), :)) .^ 2, 2));
    tss = 0.5 * (tsa + tsb + tsc);
    tra = sqrt(tss .* (tss - tsa) .* (tss - tsb) .* (tss - tsc));
    smp = opts.srfsmp;
    if numel(smp) == 1
        opts.srf = crd + smp .* nrm;
    else
        nrm = [lsqueeze(nrm(:, 1) * smp), ...
               lsqueeze(nrm(:, 2) * smp), ...
               lsqueeze(nrm(:, 3) * smp)];
        opts.srf = repmat(crd, numel(smp), 1) + nrm;
    end
    opts.srfsmp = [size(crd, 1), numel(smp)];
    opts.srf(:, 4) = 1;
    opts.srf = opts.srf * opts.srftrf';
    opts.srf(:, 4) = [];
end
srf = opts.srf;
if ~isfield(opts, 'stype') || ...
   ~isa(opts.stype, 'double') || ...
   ~any(numel(opts.stype) == [1, 2]) || ...
    any(isinf(opts.stype) | isnan(opts.stype)) || ...
    any(opts.stype ~= 1 & opts.stype ~= 2)
    opts.stype = [1, 2];
elseif numel(opts.stype) == 1
    opts.stype = opts.stype .* ones(1, 2);
else
    opts.stype = opts.stype(:)';
end
opts.stype(1) = min(opts.stype);
if all(opts.stype == 1)
    stypes = ' (1-tailed statistic)';
elseif opts.stype(1) == 1
    stypes = ' (1 tail of a 2-tailed statistic)';
else
    stypes = ' (2 tails of a 2-tailed statistic)';
end
if ~isfield(opts, 'tdf') || ...
   ~isa(opts.tdf, 'double') || ...
    numel(opts.tdf) ~= 1 || ...
    isinf(opts.tdf) || ...
    isnan(opts.tdf) || ...
    opts.tdf < 1
    opts.tdf = 1;
elseif regnsub ~= 0
    opts.tdf = regnsub;
else
    opts.tdf = min(240, round(opts.tdf));
end
if ~isfield(opts, 'thr')
    thr = 0.001;
elseif ~isa(opts.thr, 'double') || ...
    isempty(opts.thr) || ...
    numel(opts.thr) > 100 || ...
    any(isinf(opts.thr(:)) | isnan(opts.thr(:)) | opts.thr(:) <= 0 | opts.thr(:) > 0.5)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid opts.thr field.' ...
    );
else
    thr = opts.thr(:)';
end
nthr = numel(thr);
if ~isfield(opts, 'zshift')
    zshift = 0;
elseif ~isa(opts.zshift, 'double') || ...
    numel(opts.zshift) ~= 1 || ...
    isinf(opts.zshift) || ...
    isnan(opts.zshift)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid opts.zshift field.' ...
    );
else
    zshift = opts.zshift;
end

% tails
onet = (opts.stype(2) == 1);
botht = all(opts.stype == 2);

% for simulated datamaps
if regnsub == 0

    % one-tailed statistic
    if onet
        if opts.tdf == 1
            zthr = sdist('norminv', 0.5 .* thr, 0, 1);
        else
            zthr = sdist('tinv', 0.5 .* thr, opts.tdf - 1);
        end
        zthr = zthr .* zthr;

    % two-tailed with both tails
    elseif botht
        if opts.tdf == 1
            zthr = -sdist('norminv', 0.5 .* thr, 0, 1);
        else
            zthr = -sdist('tinv', 0.5 .* thr, opts.tdf - 1);
        end

    % one of a two-tailed
    else
        if opts.tdf == 1
            zthr = -sdist('norminv', thr, 0, 1);
        else
            zthr = -sdist('tinv', thr, opts.tdf - 1);
        end
    end

% for actual data
else
    
    % compute thresholds
    zthr = zeros(nregmaps, numel(thr));
    for rc = 1:nregmaps

        % t-test
        if modttest(rc)
            % both tails
            if botht
                zthr(rc, :) = -sdist('tinv', 0.5 .* thr, regnreg(rc));

            % one of a two-tailed
            else
                zthr(rc, :) = -sdist('tinv', thr, regnreg(rc));
            end

        % F-test
        else
            zthr(rc, :) = sdist('finv', thr, sum(opts.regmodsc), regnreg(rc), true);
        end
    end
    
    % smallest value
    if rc > 1
        zthr = min(zthr, [], 1);
    end
    
    % compute point estimate
    if nargout > 2
        newtc = [];
        for rc = 1:nregmaps
            if modttest(rc)
                newb = modixx{rc} * newmod{rc}' * opts.regmaps{rc};
                newe = regnsfc(rc) .* sqrt(varc(opts.regmaps{rc} - newmod{rc} * newb));
                for tc = 1:size(opts.regmodsc{rc}, 1)
                    cv = opts.regmodsc{rc}(tc, :);
                    newt = reshape((cv * newb) ./ (sqrt(cv * modixx{rc} * cv') .* newe), ddim);
                    newt(isinf(newt(:)) | isnan(newt(:))) = 0;
                    if isempty(newtc)
                        newtc = newt;
                    else
                        newtc = conjval(newtc, newt);
                    end
                end

            % F-test (todo)
            else
            end
        end
        pe = newt;
        zcnt = zeros([size(pe), 2]);
    end
end
scc = 0;

% create counting arrays
cc = zeros(nthr, 1000);
fc = zeros(nthr, niter);
maxstd = zeros(niter, 1 + double(botht));

% compute scaling factor
kern = smoothkern(opts.fwhm, 0, false, 'linear');
scf = sum(abs(kern(:))) / sqrt(sum(kern(:) .* kern(:)));

% prepare convolution FFT kernel if required
if regnsub(1) == 0 && ...
    fftconv
    kdim = size(kern);
    rsd = 1 + round(0.5 .* (kdim - ddim));
    if kdim(1) >= ddim(1)
        kern = kern(rsd(1):rsd(1)+2*(floor(0.5*(ddim(1)-1))), :, :);
    end
    if kdim(2) >= ddim(2)
        kern = kern(:, rsd(2):rsd(2)+2*(floor(0.5*(ddim(2)-1))), :);
    end
    if kdim(3) >= ddim(3)
        kern = kern(:, :, rsd(3):rsd(3)+2*(floor(0.5*(ddim(3)-1))));
    end
    kdim = size(kern);
    kdimh = floor(kdim ./ 2);
    fftkern = zeros(ddim);
    ddh = round((ddim + 1) / 2);
    fftkern(ddh(1)-kdimh(1):ddh(1)+kdimh(1), ...
            ddh(2)-kdimh(2):ddh(2)+kdimh(2), ...
            ddh(3)-kdimh(3):ddh(3)+kdimh(3)) = kern;
    fftkern = fftn(fftkern);
end

% extend mask
if ~isempty(mask) && ...
    nconj > 1
    mask = mask(:, :, :, ones(1, nconj));
end

% test xprogress
if niter >= 50
    if isempty(opts.pbar)
        try
            pbar = xprogress;
            xprogress(pbar, 'settitle', 'Running alphasim...');
            xprogress(pbar, 0, sprintf('0/%d iterations, %d thresholds%s...', ...
                niter, nthr, msktxt), 'visible', 0, 1);
            pbarn = '';
            pst = niter / 100;
            psn = pst;
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            pbar = [];
            psn = Inf;
        end
    else
        pst = ceil(niter / 200);
        psn = pst;
        pbar = opts.pbar;
        pbar.Progress(0, sprintf('alphasim: 0/%d iterations, %d thresholds%s...', niter, nthr, msktxt));
        pbarn = 'alphasim: ';
    end
else
    pbar = [];
    psn = Inf;
end

% extend ddim
ddim(4) = nconj;
ddim(5) = opts.tdf;
nconjdf = nconj * opts.tdf;

% run loop
for n = 1:niter

    % simulated data
    if regnsub(1) == 0

        % create data
        r = randn(ddim);

        % voxel-space convolution
        if ~fftconv
            for nc = 1:nconjdf
                r(:, :, :, nc) = conv3d(conv3d(conv3d(r(:, :, :, nc), ...
                    kern{1}), kern{2}), kern{3});
            end

        % frequency-space convolution
        else
            for nc = 1:nconjdf
                rf = fftn(r(:, :, :, nc));
                rf = rf .* fftkern;
                r(:, :, :, nc) = fftshift(ifftn(rf));
            end
        end

        % re-scale to unity variance again
        if opts.tdf == 1
            for nc = 1:nconj
                scc = scc + 1;

                % re-scale map
                r(:, :, :, nc) = scf .* r(:, :, :,nc);
            end

        % or compute test statistic
        else
            r = sqrt(opts.tdf) .* mean(r, 5) ./ std(r, [], 5);
            if nconj > 1
                r = conjval(r, 4);
            end
        end

        % one-tailed statistic
        if onet

            % square maps
            r = r .* r;
        end

        % shift towards requested end
        if zshift ~= 0
            r = r + zshift;
        end

        % mask
        if ~isempty(mask)
            r = r .* mask;
        end

        % conjunction of different signed tails?
        if nconj > 1 && ...
           ~onet

            % get sign of main map
            rs = sign(r(:, :, :, 1));
        end

        % iterate over other maps
        for nc = 2:nconj

            % for one-tailed
            if onet

                % simple minimum
                r(:, :, :, 1) = min(r(:, :, :, 1), r(:, :, :, nc));

            % otherwise
            else

                % absolute minimum where direction is the same
                r(:, :, :, 1) = rs .* (rs == sign(r(:, :, :, nc))) .* ...
                     min(abs(r(:, :, :, 1)), abs(r(:, :, :, nc)));
            end
        end

        % make sure we end up with one map
        if nconj > 1
            r = r(:, :, :, 1);
        end

        % with surface
        if ~isempty(srf)

            % sample volume at coordinates
            smp = (1 / opts.srfsmp(2)) .* sum(reshape(limitrangec( ...
                flexinterpn_method(r, srf, 'linear'), -1e10, 1e10, 0), ...
                opts.srfsmp), 2);
        end

        % global estimate
        if nargout > 3
            if isempty(srf)
                if botht
                    maxstd(n, 1) = max(r(:));
                    maxstd(n, 2) = -min(r(:));
                else
                    maxstd(n) = max(abs(r(:)));
                end
            else
                if botht
                    maxstd(n, 1) = max(smp);
                    maxstd(n, 2) = -min(smp);
                else
                    maxstd(n) = max(abs(smp));
                end
            end
        end

        % do for each threshold
        for tc = 1:nthr

            % for volume-based output
            if isempty(srf)

                % both tails
                if botht

                    % compute cluster frequency for both tails!
                    cf = [lsqueeze(clustercoordsc(r >= zthr(tc), clconn)); ...
                        lsqueeze(clustercoordsc(r <= -zthr(tc), clconn))];
                else

                    % just for positive tail
                    cf = clustercoordsc(r >= zthr(tc), clconn);
                end

            % for surface-based output
            else

                % then cluster surface maps
                if botht
                    cf = [lsqueeze(ceil(clustermeshmapbin(smp >= zthr(tc), ...
                            nei, crd, tra, trb, 0, 1))); ...
                          lsqueeze(ceil(clustermeshmapbin(smp <= -zthr(tc), ...
                            nei, crd, tra, trb, 0, 1)))];
                else
                    cf = ceil(clustermeshmapbin(smp >= zthr(tc), ...
                        nei, crd, tra, trb, 0, 1));
                end
            end

            % largest cluster
            if ~isempty(cf)
                mc = max(cf);
            else
                mc = 0;
            end

            % extend array if necessary
            if mc > size(cc, 2)
                cc(1, mc + ceil(size(cc, 2) / 12)) = 0;
            end

            % put into frequency arrays
            fc(tc, n) = mc;
            for nc = 1:numel(cf)
                cc(tc, cf(nc)) = cc(tc, cf(nc)) + 1;
            end
        end

    % actual data supplied
    else

        % iterate over models
        newtc = [];
        if opts.reglink
            rsgn = sign(randn(regnsub(1), 1));
            [rdt, rord] = sort(rand(regnsub(1), 1));
        end
        for rc = 1:nregmaps

            % generate new model
            if size(opts.regmodel{rc}, 2) == 1
                if opts.reglink
                    opts.regmaps{rc} = ...
                        repmat(rsgn, pddim) .* opts.regmaps{rc};
                else
                    opts.regmaps{rc} = ...
                        repmat(sign(randn(regnsub(rc), 1)), pddim) .* opts.regmaps{rc};
                end
            else
                if opts.reglink
                    neword = rord;
                else
                    [rdt, neword] = sort(rand(regnsub(rc), 1));
                end
                newmod{rc} = opts.regmodel{rc}(neword, :);
            end

            % perform regression and compute t-stats (the fast way)
            if modttest(rc)
                newb = modixx{rc} * newmod{rc}' * opts.regmaps{rc};
                newe = regnsfc(rc) .* sqrt(varc(opts.regmaps{rc} - newmod{rc} * newb));
                for tc = 1:size(opts.regmodsc{rc}, 1)
                    cv = opts.regmodsc{rc}(tc, :);
                    newt = reshape((cv * newb) ./ (sqrt(cv * modixx{rc} * cv') .* newe), ddim);
                    newt(isinf(newt(:)) | isnan(newt(:))) = 0;
                    if isempty(newtc)
                        newtc = newt;
                    else
                        newtc = conjval(newtc, newt);
                    end
                end

            % F-test
            else % STILL NEEDS IMPLEMENTATION
                newmod(:, 2) = 1;
                newi = invnd(newmod' * newmod);
                newb = newi * newmod' * opts.regmaps;
                newe = regnsfc .* sqrt(varc(opts.regmaps - newmod * newb));
                newt = reshape(newb(1, :) ./ (sqrt(newi(1)) .* newe), ddim);
                newt(isinf(newt(:)) | isnan(newt(:))) = 0;
            end
        end
        newt = newtc;

        % with surface
        if ~isempty(srf)

            % sample volume at coordinates
            smp = (1 / opts.srfsmp(2)) .* sum(reshape(limitrangec( ...
                flexinterpn_method(newt, srf, 'linear'), -1e10, 1e10, 0), ...
                opts.srfsmp), 2);
        end

        % z-score?
        if nargout > 2
            zcnt(:, :, :, 1) = zcnt(:, :, :, 1) + (newt < pe) + 0.5 .* (newt == pe);
            if all(modttest)
                zcnt(:, :, :, 2) = zcnt(:, :, :, 2) + (newt > pe) + 0.5 .* (newt == pe);
            end
            if nargout > 3
                if isempty(srf)
                    if botht
                        maxstd(n, 1) = max(newt(:));
                        maxstd(n, 2) = -min(newt(:));
                    else
                        maxstd(n) = max(abs(newt(:)));
                    end
                else
                    if botht
                        maxstd(n, 1) = max(smp);
                        maxstd(n, 2) = -min(smp);
                    else
                        maxstd(n) = max(abs(smp));
                    end
                end
            end
        end

        % do for each threshold
        for tc = 1:nthr

            % for volume-based output
            if isempty(srf)

                % both tails
                if botht

                    % compute and combine both tails' clusters
                    cf = [lsqueeze(clustercoordsc(newt >= zthr(tc), clconn)); ...
                        lsqueeze(clustercoordsc(newt <= -zthr(tc), clconn))];

                % only positive tail
                else

                    % cluster frequency
                    cf = clustercoordsc(newt >= zthr(tc), clconn);
                end

            % for surface-based output
            else

                % both tails
                if botht
                    cf = [lsqueeze(ceil(clustermeshmapbin(smp >= zthr(tc), ...
                            nei, crd, tra, trb, 0, 1))); ...
                          lsqueeze(ceil(clustermeshmapbin(smp <= -zthr(tc), ...
                            nei, crd, tra, trb, 0, 1)))];

                % just positive tail
                else
                    cf = ceil(clustermeshmapbin(smp >= zthr(tc), ...
                        nei, crd, tra, trb, 0, 1));
                end
            end

            % largest cluster
            if ~isempty(cf)
                mc = max(cf);
            else
                mc = 0;
            end

            % extend array if necessary
            if mc > size(cc, 2)
                cc(1, mc + ceil(size(cc, 2) / 12)) = 0;
            end

            % put into frequency arrays
            fc(tc, n) = mc;
            for nc = 1:numel(cf)
                cc(tc, cf(nc)) = cc(tc, cf(nc)) + 1;
            end
        end
    end

    % update progress bar
    if n >= psn && ~isempty(pbar)
        pbar.Progress(n / niter, sprintf('%s%d/%d iterations, %d thresholds%s...', ...
            pbarn, n, niter, nthr, msktxt));
        pbar.Visible = 'on';
        psn = psn + pst;
    end
end

% close progress bar
if ~isempty(pbar) && isempty(opts.pbar)
    closebar(pbar);
end

% get size and data
mf = max(1, max(fc, [], 2));
cc = cc(:, 1:max(mf));

% prepare output
tout = cell(nthr, 1);
for tc = 1:nthr
    hf = hist(fc(tc, :), 1:mf(tc));
    hfs = cumsum(hf(:));
    hx = cc(tc, 1:mf(tc)) .* (1:mf(tc));
    hxs = [0; hx(:)];
    hxs(end) = [];
    ht = sum(hx);
    sc = sum(cc(tc, 1:mf(tc)));
    ccx = cc(tc, :)';
    tout{tc} = [(1:mf(tc))', ccx(1:mf(tc)), cumsum(ccx(1:mf(tc))) ./ sc, ...
        thr(tc) .* (sum(hx) - cumsum(hxs(:))) ./ ht, hf(1:mf(tc))', ...
        1 - ([0; hfs(1:mf(tc)-1)]) ./ niter];
end

% output variable or table
if nargout < 1
    for tc = 1:nthr
        disp(' ');
        disp(sprintf('Uncorrected threshold: p < %f%s', thr(tc), stypes));
        disp('------------------------------------------------------------');
        disp(' Cl Size  Frequency  CumProbCl  p / Voxel  MaxFreq   Alpha ');
        stout = tout{tc};
        stout(stout(:, 2) == 0, :) = [];
        if isempty(srf)
            disp(sprintf(' %7d  %9d  %9.7f  %9.7f  %8d  %7.5f\n', lsqueeze(stout')));
        else
            disp(sprintf('%5dmm2  %9d  %9.7f  %9.7f  %8d  %7.5f\n', lsqueeze(stout')));
        end
    end
else
    if nthr == 1
        varargout{1} = tout{1};
    else
        varargout{1} = tout;
    end
    if nargout > 1
        varargout{2} = scf;
        if nargout > 2
            zcnt(zcnt == 0) = 0.5;
            zcnt(zcnt == niter) = (niter - 0.5);
            varargout{3} = sdist('norminv', zcnt ./ niter, 0, 1);
            if nargout > 3
                varargout{4} = maxstd;
            end
        end
    end
end
