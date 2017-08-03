function [acc, cacc, tacc] = slnbc(data, cond, run, opts)
%SLNBC  Search-light implementation for Naive Bayesian Classifier
%   ACC = SLNBC(DATA, COND) computes the leave-one-out classification
%   accuracy (overall) for the data in DATA using the condition information
%   in COND with, by default, at most a total of 5000 pair drawings.
%
%   ACC = SLNBC(DATA, COND, RUN) performs the computation with different
%   runs configured. For a single run (empty run argument, or with one
%   unique value), the function uses a leave-one-out scheme. When run
%   information is given, the function computes the classification using
%   out-of-run data only.
%
%   ACC = SLNBC(DATA, COND, RUN, OPTS) allows to set optional settings as
%   a 1x1 struct with fields:
%
%       .npairs     number of pairs (default: 5000)
%       .outtype    either 'hits', 'logodds', or 'scaledz', {'z'}
%
%   Notes: outtype 'hits' averages binary classification accuracy (to a
%   fraction of 1, with a 2-class, balanced design typically reaching 0.5).
%   Outtype 'logodds' computes the log-odds for each sample and averages
%   this value across pairings from the other class(es).

% argument check
if nargin < 2 || (~isa(data, 'double') && ~isa(data, 'single')) || size(data, 1) < 4 || ...
    isempty(data) || numel(cond) ~= size(data, 1) || any(isinf(cond(:)) | isnan(cond(:)) < cond(:) < 1)
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
data = double(data);
cond = round(cond(:));
nt = numel(cond);
nf = size(data, 2);

% runs?
if nargin < 3 || isempty(run)
    run = 1;
elseif ~isa(run, 'double') || numel(run) ~= nt || any(isinf(run(:)) | isnan(run(:)) | run(:) < 1)
    error('neuroelf:general:badArgument', 'Bad RUN argument.');
else
    run = fix(real(run(:)));
end
ur = unique(run);
nur = numel(ur);

% recode conditions (to 1:N)
uc = unique(cond);
nuc = numel(uc);
if nuc < 2
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
uci = zeros(max(uc), 1);
uci(uc) = 1:nuc;
cond = uci(cond);
cond = cond(:);

% options
npairs = 5000;
outtype = 'z';
if nargin > 3 && isstruct(opts) && numel(opts) == 1
    if isfield(opts, 'npairs') && isa(opts.npairs, 'double') && ...
        numel(opts.npairs) == 1 && ~isinf(opts.npairs) && ~isnan(opts.npairs)
        npairs = min(100000, ceil(opts.npairs));
    end
    if isfield(opts, 'outtype') && ischar(opts.outtype) && ~isempty(opts.outtype) && ...
        any(lower(opts.outtype(1)) == 'hlsz')
        outtype = lower(opts.outtype(1));
    end
end

% init number-per-class and weighting (averaging) matrix
nc = zeros(1, nuc);
w = zeros(nt, nuc);

% single run and leave-one-out
if nur == 1

    % allocate data for means and variance
    m = data;
    v = data;

    % compute means and vars (using leave-one-out scheme)
    ci = cell(1, nuc);
    um = cell(1, nuc);
    uv = cell(1, nuc);
    for cc = 1:nuc
        ci{cc} = find(cond == cc);
        nc(cc) = numel(ci{cc});
        w(ci{cc}, cc) = 1 / nc(cc);
        [m(ci{cc}, :), v(ci{cc}, :)] = loomeanvar(data(ci{cc}, :));
        um{cc} = (1 / nc(cc)) .* sum(data(ci{cc}, :), 1);
        uv{cc} = var(data(ci{cc}, :), [], 1);
    end

% for multiple runs
else

    % allocate data for means and variance
    re = [1, 1, nuc];
    m = repmat(data, re);
    v = repmat(data, re);

    % recode runs
    uri = zeros(max(ur), 1);
    uri(ur) = 1:nur;
    run = uri(run);
    run = run(:);
    
    % we need runs according to compute leave-run-out data
    um = zeros(nur, nf);
    uv = um;
    for cc = 1:nuc
        ci = (cond == cc);
        nc(cc) = sum(ci);
        w(ci, cc) = (1 / sum(ci));
        for rc = 1:nur
            ri = find(ci & run ~= rc);
            um(rc, :) = (1 / numel(ri)) .* sum(data(ri, :), 1);
            uv(rc, :) = var(data(ri, :), [], 1);
        end
        m(:, :, cc) = um(run, :);
        v(:, :, cc) = uv(run, :);
    end
end

% total number of pairs
if nur == 1
    cpairs = squeeze(prod(pairs(nc(:)), 3));
    tpairs = sum(cpairs);
else
    tpairs = nt;
end

% don't do more than we need to
npairs = max(nt, min(npairs, tpairs));

% for leave-one-out (within-run) use class pairs
if nur == 1

    % numbers of pairs to draw
    if npairs ~= tpairs
        dpairs = cpairs ./ (tpairs / npairs);
        while sum(round(dpairs)) < npairs
            dx = mod(dpairs, 1) + 0.5;
            [dx, dy] = max(dx);
            dpairs(dy) = dpairs(dy) + 0.5;
        end
        while sum(round(dpairs)) > npairs
            dx = mod(dpairs, 1) - 0.5;
            [dx, dy] = min(dx);
            dpairs(dy) = dpairs(dy) - 0.5;
        end
        di = 1;
        dpairs = round(dpairs);
    end

    % pre-allocate trial-accuracy data (and trial indices)
    tacc = zeros(npairs, 2);
    ttst = tacc;
    ti = 1;

    % iterate over possible class pairs
    for c1 = 1:(nuc-1)
        for c2 = (c1+1):nuc

            % all indices to test
            if npairs == tpairs

                % use ndgrid
                [ci1, ci2] = ndgrid(1:nc(c1), 1:nc(c2));

            % only some index pairs
            else

                % use random draws
                ci1 = ceil(nc(c1) .* rand(dpairs(di), 1));
                ci2 = ceil(nc(c2) .* rand(dpairs(di), 1));
                di = di + 1;
            end

            % get mean-removed data (and matching variances) for pairs
            ci1 = ci{c1}(ci1(:));
            ci2 = ci{c2}(ci2(:));
            d11 = data(ci1, :) - m(ci1, :);
            d12 = data(ci1, :) - m(ci2, :);
            d21 = data(ci2, :) - m(ci1, :);
            d22 = data(ci2, :) - m(ci2, :);
            v1 = v(ci1, :);
            v2 = v(ci2, :);
            lv1 = log(v1);
            lv2 = log(v2);

            % compute Naive-Bayesian Classification sum(log(p))
            d11 = squeeze(sum(lv2 + d12 .* d12 ./ v2, 2)) - ...
                squeeze(sum(lv1 + d11 .* d11 ./ v1, 2));
            d22 = squeeze(sum(lv1 + d21 .* d21 ./ v1, 2)) - ...
                squeeze(sum(lv2 + d22 .* d22 ./ v2, 2));

            % store accuracy data (log-odds-ratio)
            tacc(ti:ti+numel(d11)-1, :) = [d11, d22];
            ttst(ti:ti+numel(d11)-1, :) = [ci1(:), ci2(:)];
            ti = ti + numel(d11);
        end
    end

% winner
else

    % compute all data - means combinations (class: 3rd dim)
    data = repmat(data, re) - m;

    % compute log odds (without 2*pi factor)
    data = squeeze(sum(log(v) + data .* data ./ v, 2));

    % get value for correct class
    ttst = (1:nt)';
    cdi = ttst + nt .* (cond - 1);
    cdata = data(cdi);

    % get next lowest value
    data(cdi) = -Inf;
    data = sort(data, 2);
    tacc = data(:, 2) - cdata;
end

% convert to hits?
if outtype == 'h'
    tacc = double(tacc > 0);

% convert to Z?
elseif any('sz' == outtype)
    tacc = sign(tacc) .* sqrt(abs(2 .* tacc));
    if outtype == 's'
        tacc = tacc ./ (1 + 1 ./ (tacc .* tacc));
    end
end

% return averages
tst = histcountn(ttst(:), 1, nt + 1, 1);
tacc = histcountn(ttst(:), 1, nt + 1, 1, tacc(:)) ./ tst;
cacc = w' * tacc;
acc = mean(cacc);
