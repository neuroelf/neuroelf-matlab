function [p, se, t, df, fr] = mediationpset(x, m, y, varargin)
% mediationpset  - compute paths, std. errors, t-scores of a mediation model
%
% FORMAT:       [p, se, t] = mediationpset(x, m, y [, c] [, opts])
%
% Input fields:
%
%       x           independent variable (X)
%       m           mediation variable(s) (M)
%       y           outcome variable (Y)
%       c           potential covariate(s) (C)
%       opts        optional settings
%        .bootmeth  either of 'bca', {'perc'}, 'var'
%        .bootn     number of bootstrapping samples (default: 5000)
%        .bootre    re-use sampling across multiple cases (default: false)
%        .bootsmax  bootstrap tuning parameter, if given must be set to a
%                   number up to which any given sample can (re-) occur in
%                   the resampling to be considered (default: N)
%        .bootsmin  bootstrap tuning parameter, if given must be set to a
%                   percentile of values which must be present in each
%                   resampling to be considered (default: 0)
%        .mcmama    MCMAM alpha level (normative value, default: 0.05)
%        .progress  either 1x1 xprogress of xfigure:Progress object
%        .robust    robust regression of coefficients (and s(a), s(b))
%        .sabmeth   s(ab) method, one of {'boot'}, 'mcmam', 'sobel'
%        .ztrans    apply z-transform to model and data (default: false)
%
% Output fields:
%
%       p           path coefficients
%       se          standard error of path coefficients
%       t           t statistic (p / se(p)); order of all outputs is
%                   - V(c) (X->Y)
%                   - V(a) (X->M1, ..., X->Mn)
%                   - V(b) (M1->Y, ..., Mn->Y)
%                   - V(ab) (M1, ..., Mn)
%                   - V(c') (X->Y w/o M*->Y)
%                   for now, no statistic is computed for the covariates!
%                   the path coefficients for a, b, and ab paths are stored
%                   interspersedly (a1, b1, ab1, ..., aN, bN, abN)
%
% Note: the bootsmax/min options can be used to force a minimal coverage
%       (bootsmin) of the sample (e.g. 0.6 for 60% of values must at
%       least be sampled once) as well as to restrict the maximal number
%       of draws per sample (e.g. 3 so that each sample can be drawn up to
%       3 times at most) which can be beneficial to diminish the effects of
%       a highly unusual distribution as well as outliers on the resampling

% Version:  v0.9b
% Build:    13022213
% Date:     Apr-09 2011, 11:08 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the
%       distribution.
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

% preliminary argument check
if nargin < 3 || ...
   ~isnumeric(x) || ...
   ~isnumeric(m) || ...
   ~isnumeric(y) || ...
    numel(x) < 2 || ...
    numel(m) < 2 || ...
    numel(y) < 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
c = [];
opts = struct;
if nargin > 3
    if isnumeric(varargin{1}) && ...
        numel(varargin{1}) > 1
        c = double(varargin{1});
    end
    if isstruct(varargin{end}) && ...
        numel(varargin{end}) == 1
        opts = varargin{end};
    end
end
if isfield(opts, 'robust') && ...
    islogical(opts.robust) && ...
    numel(opts.robust) == 1
    dorobust = opts.robust;
else
    dorobust = false;
end

% try to establish required dimensions
xs = size(x);
xs(xs < 2) = [];
ms = size(m);
ms(ms < 2) = [];
ys = size(y);
ys(ys < 2) = [];
cs = size(c);
cs(cs < 2) = [];

% get number of cases (dim along which regression is performed)
ncases = xs(end);

% check number of cases
if (ncases ~= ms(end) && ...
    (numel(ms) < 2 || ...
     ncases ~= ms(end-1))) || ...
    ncases ~= ys(end)
    error( ...
        'neuroelf:BadArgument', ...
        'X, M, and Y must share the size of their last non-singleton dimension.' ...
    );
end
if ~isempty(cs) && ...
   (ncases ~= cs(end) && ...
    (numel(cs) < 2 || ...
     ncases ~= cs(end-1)))
    error( ...
        'neuroelf:BadArgument', ...
        'Any given covariate must match the number of cases.' ...
    );
end

% number of mediators
if numel(ms) > 1 && ...
    ncases ~= ms(end)
    nmedi = ms(end);
else
    nmedi = 1;
end

% number of covariates
if ~isempty(cs) && ...
    ncases ~= cs(end)
    ncovs = cs(end);
elseif ~isempty(cs)
    ncovs = 1;
else
    ncovs = 0;
end

% get output size
if numel(xs) > 1
    os = xs(1:end-1);
elseif numel(ys) > 1
    os = ys(1:end-1);
elseif numel(ms) > 1 && ...
    ms(end) == ncases
    os = ms(1:end-1);
elseif numel(ms) > 2 && ...
    ms(end-1) == ncases
    os = ms(1:end-2);
elseif numel(cs) > 1 && ...
    cs(end) == ncases
    os = cs(1:end-1);
elseif numel(cs) > 2 && ...
    cs(end-1) == ncases
    os = cs(1:end-2);
else
    os = [];
end
on = prod(os);

% check number of elements for all arrays
if ~any(numel(x) == (ncases .* [1, on])) || ...
   ~any(numel(m) == ((ncases * nmedi) .* [1, on])) || ...
   ~any(numel(y) == (ncases .* [1, on])) || ...
   ~any(numel(c) == ((ncases * ncovs) .* [1, on]))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid number of elements in one of the inputs.' ...
    );
end

% parse options
if ~isfield(opts, 'bootmeth') || ...
   ~ischar(opts.bootmeth) || ...
   ~any(strcmpi(opts.bootmeth(:)', {'bca', 'perc', 'var'}))
    opts.bootmeth = 'perc';
else
    opts.bootmeth = lower(opts.bootmeth(:)');
end
if ~isfield(opts, 'bootn') || ...
   ~isa(opts.bootn, 'double') || ...
    numel(opts.bootn) ~= 1 || ...
    isinf(opts.bootn) || ...
    isnan(opts.bootn) || ...
    opts.bootn < 1
    opts.bootn = 5000;
else
    opts.bootn = round(min(10000000, max(100, opts.bootn)));
end
if ~isfield(opts, 'bootre') || ...
   ~islogical(opts.bootre) || ...
    numel(opts.bootre) ~= 1
    opts.bootre = false;
end
if ~isfield(opts, 'bootsmax') || ...
   ~isa(opts.bootsmax, 'double') || ...
    numel(opts.bootsmax) ~= 1 || ...
    isinf(opts.bootsmax) || ...
    isnan(opts.bootsmax) || ...
    opts.bootsmax < 1
    opts.bootsmax = ncases;
else
    opts.bootsmax = round(min(ncases, max(1, opts.bootsmax)));
end
if ~isfield(opts, 'bootsmin') || ...
   ~isa(opts.bootsmin, 'double') || ...
    numel(opts.bootsmin) ~= 1 || ...
    isinf(opts.bootsmin) || ...
    isnan(opts.bootsmin) || ...
    opts.bootsmin < 0 || ...
    opts.bootsmin > 1
    opts.bootsmin = 0;
end
if ~isfield(opts, 'mcmama') || ...
   ~isa(opts.mcmama, 'double') || ...
    numel(opts.mcmama) ~= 1 || ...
    isinf(opts.mcmama) || ...
    isnan(opts.mcmama) || ...
    opts.mcmama < (2 / opts.bootn) || ...
    opts.mcmama > 0.1
    opts.mcmama = max(0.05, 2 / opts.bootn);
end
if ~isfield(opts, 'robust') || ...
   ~islogical(opts.robust) || ...
    numel(opts.robust) ~= 1
    opts.robust = false;
end
if ~isfield(opts, 'progress') || ...
    numel(opts.progress) ~= 1 || ...
   (~isa(opts.progress, 'xprogress') && ...
    ~isa(opts.progress, 'xfigure'))
    opts.progress = [];
end
if ~isfield(opts, 'sabmeth') || ...
   ~ischar(opts.sabmeth) || ...
    isempty(opts.sabmeth) || ...
   ~any(strcmpi(opts.sabmeth(:)', {'boot', 'mcmam', 'sobel'}))
    opts.sabmeth = 'boot';
else
    opts.sabmeth = lower(opts.sabmeth(:)');
end
if strcmp(opts.sabmeth, 'boot')
    opts.robust = dorobust;
end
if ~isfield(opts, 'ztrans') || ...
   ~islogical(opts.ztrans) || ...
    numel(opts.ztrans) ~= 1
    opts.ztrans = false;
end

% total number of statistics
nstats = 2 + 3 * nmedi;

% degrees of freedom
if opts.robust
    df = zeros(on, nstats);
else
    df = zeros(nstats, 1);
    df(1) = ncases - (2 + ncovs);
    df(2:3:end) = ncases - (2 + ncovs);
    df(3:3:end) = ncases - (2 + nmedi + ncovs);
    df(4:3:end) = ncases - (2 + nmedi + ncovs);
    df(end) = ncases - (2 + nmedi + ncovs);
end

% permute and reshape X
if numel(xs) > 1
    pd = findfirst(size(x) == ncases, -1);
    pa = 1:ndims(x);
    pa(pd) = [];
    pa = [pd, pa];
    x = permute(double(x), pa);
else
    x = double(x(:));
end
x = reshape(x, ncases, numel(x) / ncases);
xs = (numel(x) == ncases);

% permute and reshape M
if numel(ms) > 1
    pd = findfirst(size(m) == ncases, -1);
    pa = 1:ndims(m);
    if ms(end) == ncases
        pa(pd) = [];
        pa = [pd, pa];
    else
        pa(pd:end) = [];
        pa = [pd:ndims(m), pa];
    end
    m = permute(double(m), pa);
else
    m = double(m(:));
end
m = reshape(m, [ncases, nmedi, numel(m) / (ncases * nmedi)]);
ms = (numel(m) == (nmedi * ncases));

% permute and reshape Y
if numel(ys) > 1
    pd = findfirst(size(y) == ncases, -1);
    pa = 1:ndims(y);
    pa(pd) = [];
    pa = [pd, pa];
    y = permute(double(y), pa);
else
    y = double(y(:));
end
y = reshape(y, [ncases, 1, numel(y) / ncases]);
ys = (numel(y) == ncases);

% permute and reshape C (if necessary)
if numel(cs) > 1
    pd = findfirst(size(c) == ncases, -1);
    pa = 1:ndims(c);
    if cs(end) == ncases
        pa(pd) = [];
        pa = [pd, pa];
    else
        pa(pd:end) = [];
        pa = [pd:ndims(c), pa];
    end
    c = permute(c, pa);
    c = reshape(c, [ncases, ncovs, numel(c) / (ncases * ncovs)]);
else
    c = c(:);
end
if ~isempty(c)
    usecov = true;
else
    usecov = false;
end
cs = (numel(c) == (ncovs * ncases));

% transform data ?
if opts.ztrans

    % transform each component
    x = ztrans(x, 1);
    if any(isinf(x(:)) | isnan(x(:)))
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid X data for z transformation.' ...
        );
    end
    m = ztrans(m, 1);
    if any(isinf(m(:)) | isnan(m(:)))
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid M data for z transformation.' ...
        );
    end
    y = ztrans(y, 2);
    if any(isinf(y(:)) | isnan(y(:)))
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid Y data for z transformation.' ...
        );
    end
    if usecov
        c = ztrans(c, 1);
        if any(isinf(c(:)) | isnan(c(:)))
            error( ...
                'neuroelf:BadArgument', ...
                'Invalid C data for z transformation.' ...
            );
        end
    end
end

% first model predictor and mediator
if opts.ztrans
    fmp = 1;
else
    fmp = 2;
end
fmed = fmp + ncovs + 1;

% model is the same for all datapoints
if xs && ...
    cs && ...
    ms

    % ones required?
    if fmp > 1

        % have ones first
        models = [ones(ncases, fmp - 1), x, c, m];

    % no ones required
    else

        % catenate model
        models = [x, c, m];
    end
    mpred = fmed:size(models, 2);

% model differs for data
else

    % prepare model space
    models = ones(ncases, fmp + ncovs + nmedi, on);
    mpred = fmed:size(models, 2);

    % fill models with X
    if xs
        models(:, fmp, :) = x(:, :, ones(1, on));
    else
        models(:, fmp, :) = reshape(x, [ncases, 1, on]);
    end

    % covariates in model?
    if usecov

        % covariate predictors
        covpred = fmp+1:fmp+ncovs;

        % same covariates for all models or not
        if cs
            models(:, covpred, :) = c(:, :, ones(1, on));
        else
            models(:, covpred, :) = c;
        end
    end

    % same mediators for all models or not
    if ms
        models(:, mpred, :) = m(:, :, ones(1, on));
    else
        models(:, mpred, :) = m;
    end
end

% get model predictor lists and targets in stats
apred = 1:(size(models, 2) - nmedi);
tga = 2:3:1+3*nmedi;
tgb = 3:3:1+3*nmedi;
tgab = 4:3:1+3*nmedi;

% number of outputs and reserve space (for paths and standard errors)
p = zeros(on, nstats);

% robust regression ?
if opts.robust

    % robustly regress X (+C) -> Y/M (total effect: c, and a path coefficients)
    if xs && ...
        cs
        if ys
            [obsc, robw, ixxwom, sec] = ...
                fitrobustbisquare_multi(models(:, apred, 1), y);
        else
            [obsc, robw, ixxwom, sec] = ...
                fitrobustbisquare_img(models(:, apred, ones(1, on)), y);
            obsc = squeeze(obsc);
            sec = squeeze(sec);
        end
        if ms
            if numel(m) == ncases
                [obsa, robw, robixx, sea] = ...
                   fitrobustbisquare_multi(models(:, apred, 1), m);
            else
                [obsa, robw, robixx, sea] = ...
                   fitrobustbisquare_img(models(:, apred, 1), reshape(m, [ncases, 1, nmedi]));
                obsa = squeeze(obsa);
                sea = squeeze(sea);
            end
        else
            [obsa, robw, robixx, sea] = ...
               fitrobustbisquare_multi(models(:, apred, ones(1, on)), m);
        end
    else
        [obsc, robw, ixxwom, sec] = ...
            fitrobustbisquare_multi(models(:, apred, :), y);
        [obsa, robw, robixx, sea] = ...
            fitrobustbisquare_multi(models(:, apred, :), m);
    end

    % robustly regress X+M (+C) -> Y (c' and b path coefficients)
    if size(models, 3) == 1
        [obsb, robw, ixxwm, seb] = ...
            fitrobustbisquare_multi(models(:, :, ones(1, on)), y);
    else
        [obsb, robw, ixxwm, seb] = ...
            fitrobustbisquare_multi(models, y);
    end

% OLS regression
else

    % OLS regress X (+C) -> Y/M (total effect: c path coefficient)
    [ixxwom, obsc, obsa] = mmregress(models(:, apred, :), y, m);

    % regress X+M (+C) -> Y (c' and b path coefficients)
    [ixxwm, obsb] = mmregress(models, y);
end

% expand if necessary
if on > 1
    if size(obsc, 3) == 1
        obsc = obsc(:, :, ones(1, on));
    end
    if size(obsa, 3) == 1
        obsa = obsa(:, :, ones(1, on));
    end
end

% store into paths
p(:, 1) = squeeze(obsc(fmp, :, :));
p(:, tga) = reshape(obsa(fmp, :, :), nmedi, on)';
p(:, tgb) = reshape(obsb(mpred, :, :), nmedi, on)';
p(:, end) = squeeze(obsb(fmp, :, :));

% compute ab product
p(:, tgab) = p(:, tga) .* p(:, tgb);

% output size
os = [os, nstats, 1];

% only one output, end here
if nargout < 2
    p = reshape(p, os);
    return;
end

% allocate space for SE computation
se = zeros(on, nstats);

% for robust path
if opts.robust

    % expand if necessary
    if on > 1
        if size(sec, 3) == 1
            sec = sec(:, :, ones(1, on));
        end
        if size(sea, 3) == 1
            sea = sea(:, :, ones(1, on));
        end
    end

    % store standard errors
    se(:, 1) = squeeze(sec(fmp, :, :));
    se(:, tga) = reshape(sea(fmp, :, :), nmedi, on)';
    se(:, tgb) = reshape(seb(mpred, :, :), nmedi, on)';
    se(:, end) = squeeze(seb(fmp, :, :));

% for OLS path
else

    % compute residuals of X->Y
    if ys
        resxy = varc(y(:, :, ones(1, on)) - ...
            transmul(models(:, apred, :), obsc), 1, 1);
    elseif size(ixxwom, 3) == 1
        resxy = varc(y - transmul(repmat(models(:, apred, 1), [1, 1, on]), obsc), 1, 1);
    else
        resxy = varc(y - transmul(models(:, apred, :), obsc), 1, 1);
    end

    % compute residuals of X->M (unless a is bootstrapped)
    if size(models, 3) == 1
        resxm = varc(models(:, mpred) - transmul(models(:, apred), obsa(:, :, 1)), 1, 1);
        if on > 1
            resxm = resxm(1, 1, ones(1, on));
        end
    else
        resxm = varc(models(:, mpred, :) - transmul(models(:, apred, :), obsa), 1, 1);
    end

    % compute residuals of (X+M)->Y
    if ys
        fr = y(:, :, ones(1, on)) - transmul(models, obsb);
    elseif size(ixxwm, 3) == 1
        fr = y - transmul(repmat(models(:, :, 1), [1, 1, on]), obsb);
    else
        fr = y - transmul(models, obsb);
    end
    resxmy = varc(fr, 1, 1);

    % compute canonical SEs
    sem = 2 + nmedi;
    lcov = fmp + ncovs;

    % -> c path
    se(:, 1) = squeeze(sqrt(((ncases - 1) / (ncases - 2)) .* ixxwom(fmp, fmp, :) .* resxy));

    % depending on what is bootstrapped
    for mc = 1:nmedi
        se(:, tga(mc)) = sqrt(((ncases - 1) / (ncases - 2)) .* ...
            ixxwom(fmp, fmp, :) .* resxm(:, mc, :));
    end
    for mc = 1:nmedi
        se(:, tgb(mc)) = sqrt(((ncases - 1) / (ncases - sem)) .* ...
            ixxwm(mc+lcov, mc+lcov, :) .* resxmy);
    end

    % -> c' path
    se(:, end) = sqrt(((ncases - 1) / (ncases - sem)) .* ixxwm(fmp, fmp, :) .* resxmy);
end

% for bootstrapping / MCMAM -> initialize progress bar
pbar = [];
if (strcmp(opts.sabmeth, 'boot') && ...
    (on * opts.bootn) > 1e6) || ...
   (strcmp(opts.sabmeth, 'mcmam') && ...
    (on * opts.bootn) > 5e6)
    try
        if isempty(opts.progress)
            pbar = xprogress;
            xprogress(pbar, 'setposition', [80, 200, 640, 36]);
            xprogress(pbar, 'settitle', ...
                sprintf('Running %dx%d mediation on %d samples...', ...
                size(models, 1), size(models, 2), on));
            xprogress(pbar, 0, 'Mediation analysis...', 'visible', 0, 1);
        else
            pbar = opts.progress;
            pbar.Progress(0, sprintf('Mediation analysis: sample %d/%d', ...
                1, on));
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        pbar = [];
    end
end

% algorithm depends on semetho d (first)
switch (opts.sabmeth)

    % bootstrapping
    case {'boot'}

        % bootstrapping options
        if on < 2
            opts.bootre = false;
        end
        btmeth = opts.bootmeth;
        btuniq = ~opts.bootre;
        numsmp = opts.bootn;
        yrs = [ncases, 1, numsmp];
        bso = struct( ...
            'maxsmp',  opts.bootsmax, ...
            'minsmp',  opts.bootsmin * ncases, ...
            'numparm', 1, ...
            'numsmp',  numsmp, ...
            'perm',    false);

        % create jack-knife sampling
        if strcmp(btmeth, 'bca')
            jki = reshape(jksample(ncases), [ncases - 1, 1, ncases]);
            jkn = size(jki, 3);
            jks = true;
            ncasej = ncases - 1;
            yrj = [ncasej, 1, jkn];

            % the X-M path does not to be jack-knifed every time
            if size(ixxwm, 3) == 1
                obsaj = cell(1, nmedi);
                for mc = 1:nmedi
                    obsaj{mc} = reshape(indexarray(m, jki, mc), [ncasej, 1, jkn]);
                end
                [ixxt, obsaj{:}] = mmregress(reshape( ...
                    indexarray(models(:, apred, 1), jki, 1:numel(apred)), ...
                    [ncasej, numel(apred), jkn]), obsaj{:});
                obsaj = cat(2, obsaj{:});
                obsaj = reshape(obsaj(fmp, :, :), nmedi, jkn)';
            end
        else
            jks = false;
        end

        % re-use the same bootstrapping sampling?
        if opts.bootre || ...
            size(ixxwm, 3) == 1

            % simply create random indices
            bss = bssample(ncases, bso);
        end

        % the X-M path does not need to be re-done every time
        if size(ixxwm, 3) == 1

            % bootstrap a path coefficients once per mediator
            obsab = cell(1, nmedi);
            for mc = 1:nmedi
                obsab{mc} = reshape( ...
                indexarray(m, bss, mc), [ncases, 1, numsmp]);
            end
            [ixxt, obsab{:}] = mmregress(reshape( ...
                indexarray(models(:, apred, 1), bss, 1:numel(apred)), ...
                [ncases, numel(apred), numsmp]), obsab{:});
            obsab = cat(2, obsab{:});
            obsab = reshape(obsab(fmp, :, :), nmedi, numsmp)';
        end

        % compute progress update interval
        if ~isempty(pbar)
            pint = ceil(20000 / numsmp);
        end

        % iterate over number of outputs
        for oc = 1:on

            % draw new bootstrap sampling
            if btuniq
                bss = bssample(ncases, bso);
            end

            % do the a path coefficients need to be bootstrapped?
            if size(ixxwm, 3) > 1

                % bootstrap a
                [ixxt, obsab] = mmregress(reshape( ...
                    indexarray(models(:, apred, oc), bss, 1:numel(apred)), ...
                    [ncases, numel(apred), numsmp]), reshape( ...
                    indexarray(models(:, mpred, oc), bss, 1:nmedi), ...
                    [ncases, nmedi, numsmp]));
                obsab = reshape(obsab(fmp, :, :), nmedi, numsmp)';

                % also create jack-knife samples
                if jks
                    [ixxt, obsaj] = mmregress(reshape( ...
                        indexarray(models(:, apred, oc), jki, 1:numel(apred)), ...
                        [ncasej, numel(apred), jkn]), reshape( ...
                        indexarray(models(:, mpred, oc), jki, 1:nmedi), ...
                        [ncasej, nmedi, jkn]));
                    obsaj = reshape(obsaj(fmp, :, :), nmedi, jkn)';
                end
            end

            % bootstrap b
            if size(models, 3) == 1
                [ixxt, obsbb] = mmregress(reshape( ...
                    indexarray(models, bss, 1:size(models, 2)), ...
                    [ncases, size(models, 2), numsmp]), reshape( ...
                    indexarray(y(:, :, oc), bss, 1), yrs));
                if jks
                    [ixxt, obsbj] = mmregress(reshape( ...
                        indexarray(models, jki, 1:size(models, 2)), ...
                        [ncasej, size(models, 2), jkn]), reshape( ...
                        indexarray(y(:, :, oc), jki, 1), yrj));
                end
            elseif size(y, 3) == 1
                [ixxt, obsbb] = mmregress(reshape( ...
                    indexarray(models(:, :, oc), bss, 1:size(models, 2)), ...
                    [ncases, size(models, 2), numsmp]), reshape( ...
                    indexarray(y, bss, 1), yrs));
                if jks
                    [ixxt, obsbj] = mmregress(reshape( ...
                        indexarray(models(:, :, oc), jki, 1:size(models, 2)), ...
                        [ncasej, size(models, 2), jkn]), reshape( ...
                        indexarray(y, jki, 1), yrj));
                end
            else
                [ixxt, obsbb] = mmregress(reshape( ...
                    indexarray(models(:, :, oc), bss, 1:size(models, 2)), ...
                    [ncases, size(models, 2), numsmp]), reshape( ...
                    indexarray(y(:, :, oc), bss, 1), yrs));
                if jks
                    [ixxt, obsbj] = mmregress(reshape( ...
                        indexarray(models(:, :, oc), jki, 1:size(models, 2)), ...
                        [ncasej, size(models, 2), jkn]), reshape( ...
                        indexarray(y(:, :, oc), jki, 1), yrj));
                end
            end
            obsbb = reshape(obsbb(mpred, :, :), nmedi, numsmp)';

            % build products
            obsabb = obsab .* obsbb;

            % remove bad values
            obsabb(isinf(obsabb(:)) | isnan(obsabb(:))) = 0;

            % depending on the se strategy
            switch (btmeth)

                % BCa method
                case {'bca'}

                    % get jack-knifed ab estimate
                    obsabj = obsaj .* reshape(obsbj(mpred, :, :), nmedi, jkn)';
                    obsabji = any(isinf(obsabj) | isnan(obsabj));
                    if any(obsabji)
                        obsabj(:, obsabji) = 0;
                    end

                    % get z-score with jack-knife estimate
                    se(oc, tgab) = ...
                        p(oc, tgab) ./ bstrapbca(p(oc, tgab), obsabj, obsabb, 0);

                % straight forward percentile method
                case {'perc'}

                    % look up percentile z-score
                    se(oc, tgab) = ...
                        p(oc, tgab) ./ bstrappct(p(oc, tgab), obsabb);

                % compute simple stdev over ab products
                case {'var'}

                    % compute standard deviation of sampled statistics
                    se(oc, tgab) = sqrt(varc(obsabb, 1, 1));
            end

            % progress
            if ~isempty(pbar) && ...
                mod(oc, pint) == 0
                pbar.Progress(oc / on, ...
                    sprintf('Mediation analysis: sample %d/%d', oc, on));
            end
        end

    % MCMAM
    case {'mcmam'}

        % create random data (with known normal distribution)
        mcm = ztrans(randn(opts.bootn, 2));

        % uncorrelate data (totally)
        mcm(:, 2) = ztrans(mcm(:, 2) - ((1 ./ sum(mcm(:, 1) .* mcm(:, 1))) * ...
            sum(mcm(:, 1) .* mcm(:, 2))) .* mcm(:, 1));

        % determine size/replication argument
        if nmedi == 1
            mco = 1;
        else
            mco = ones(opts.bootn, 1);
        end

        % LCL and UCL indices and factors, as well as median index
        cif = 1 / (-sdist('tinv', 0.5 * opts.mcmama, 1e6));
        lci = max(1, floor(0.5 * opts.mcmama * opts.bootn));
        uci = min(opts.bootn, ceil((1 - 0.5 * opts.mcmama) * opts.bootn) + 1);
        if mod(opts.bootn, 2) > 0
            mdi = ceil(0.5 * opts.bootn);
            mdn = 1;
        else
            mdi = 0.5 * opts.bootn + [0, 1];
            mdn = 0.5;
        end

        % compute progress update interval
        if ~isempty(pbar)
            pint = ceil(100000 / opts.bootn);
        end

        % iterate over number of outputs
        for oc = 1:on

            % get positive/negative path (for test != 0)
            posab = p(oc, tgab) >= 0;

            % sort standard error of product
            mcms = sort( ...
                (mco * p(oc, tga) + mcm(:, 1) * se(oc, tga)) .* ...
                (mco * p(oc, tgb) + mcm(:, 2) * se(oc, tgb)));

            % for positive ab path coefficient
            if any(posab)

                % determine the re-worked SE
                se(oc, tgab(posab)) = cif .* ...
                   (mdn .* sum(mcms(mdi, posab), 1) - ...
                    mcms(lci, posab));

            % for negative ab path coefficient
            elseif ~all(posab)

                % determine the re-worked SE
                se(oc, tgab(~posab)) = cif .* ...
                   (mcms(uci, ~posab) - ...
                    mdn .* sum(mcms(mdi, ~posab), 1));
            end

            % progress
            if ~isempty(pbar) && ...
                mod(oc, pint) == 0
                pbar.Progress(oc / on, ...
                    sprintf('Mediation analysis: sample %d/%d', oc, on));
            end
        end

    % Sobel-test
    case {'sobel'}

        % use Sobel formula
        se(:, tgab) = sqrt( ...
            p(:, tga) .* p(:, tga) .* se(:, tgb) .* se(:, tgb) + ...
            p(:, tgb) .* p(:, tgb) .* se(:, tga) .* se(:, tga));
end

% clear progress bar
if ~isempty(pbar) && ...
    isempty(opts.progress)
    closebar(pbar);
end

% reshape outputs
p = reshape(p, os);
se = reshape(se, os);

% compute t-scores
if nargout > 2

    % free up some memory (just in case)
    models(:) = [];

    % computation
    t = p ./ se;
end
