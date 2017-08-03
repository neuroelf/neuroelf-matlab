function [msq, msqn, cm, df, f] = mrmanova(y, w, b, opts)
% mrmanova  - mixed, repeated measures ANOVA
%
% FORMAT:       [msq, msqn, cm, df, f] = mrmanova(y, w, b [, opts])
%
% Input fields:
%
%       y           dependent variable(s) (last two dims BxW)
%       w           LxW factor level assignments (can be empty)
%       b           SxL subject-to-group assignments (can be empty)
%       opts        optional settings (1x1 struct)
%        .cov       covariates (must be BxC in size)
%        .nonull    reject null-values
%
% Output fields:
%
%       msq         mean sum-of-squares
%       msqn        name of MSS value/map
%       f           f-stats
%       df          Sx2 d.f. for f-stats
%       cm          cell means
%
% Notes: currently implemented models are 1-within and 1-within-1-between
%        planned implementation for up to 3-within and 2-between

% TODO:
%        .robust    perform estimation robustly (experimental!!)

% Version:  v0.9c
% Build:    13032013
% Date:     May-31 2011, 5:53 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
if nargin < 3 || ...
   ~isnumeric(y) || ...
    isempty(y) || ...
   (isempty(w) && ...
    isempty(b)) || ...
   ~isa(w, 'double') || ...
    ndims(w) > 2 || ...
    size(w, 1) > 3 || ...
   (~isempty(w) && ...
    (size(w, 2) ~= size(y, ndims(y)) || ...
     any(isinf(w(:)) | isnan(w(:))))) || ...
   ~isa(b, 'double') || ...
    ndims(b) > 2 || ...
    size(b, 2) > 2 || ...
   (~isempty(b) && ...
    (size(b, 1) ~= size(y, ndims(y) - 1) || ...
     any(isinf(b(:)) | isnan(b(:)))))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument (only up to 3-within/2-between allowed).' ...
    );
end

% get size of data
nd = ndims(y);
sz = size(y);

% subject-dimension is second to last dim
sd = nd - 1;
ns = sz(sd);

% number of within/between factors
nwf = size(w, 1);
ngf = size(b, 2);

% error if too few subjects
if ns < 3
    error( ...
        'neuroelf:BadArgument', ...
        'At least 3 subject observations required.' ...
    );
end

% parse options
if nargin < 4 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end

% covariates (not yet fully implemented!)
if ~isfield(opts, 'cov') || ...
    isempty(opts.cov) || ...
   ~isa(opts.cov, 'double') || ...
    ndims(opts.cov) > 2 || ...
    size(opts.cov, 1) ~= sd || ...
    any(isinf(opts.cov(:)) | isnan(opts.cov(:))) || ...
    any(sum(abs(diff(opts.cov)), 1) == 0)
    opts.cov = [];
end

% flag to interact covariates with group (allow per-group differences)
if ~isfield(opts, 'covgrpi') || ...
   ~islogical(opts.covgrpi) || ...
    numel(opts.covgrpi) ~= 1
    opts.covgrpi = false;
end

% regress out covariates prior to estimation
if ~isfield(opts, 'covregi') || ...
   ~islogical(opts.covregi) || ...
    numel(opts.covregi) ~= 1
    opts.covregi = true;
end

% reject null (0) as missing values
if ~isfield(opts, 'nonull') || ...
   ~islogical(opts.nonull) || ...
    numel(opts.nonull) ~= 1
    opts.nonull = false;
end

% perform robust MSQ (not tested/partially implemented)
if ~isfield(opts, 'robust') || ...
   ~islogical(opts.robust) || ...
    numel(opts.robust) ~= 1
    opts.robust = false;
end

% indexing arguments (to access data = y, use ':' to access all data)
yi = repmat({':'}, 1, nd);
yt = yi;

% generate "repmat arguments"
rmt = sz;
rmt(1:(sd-1)) = 1;
rmtm = sz;
rmtm(sd:nd) = 1;

% and source index (as many as dimensions)
rsi = ones(1, nd);

% no groups defined
if isempty(b)

    % take as one group
    b = ones(sz(sd), 1);
    gib = b;
    gwi = 1;

% group-factors
else

    % find unique rows (groups with lowest interactions)
    [g, gi, gib] = unique(b, 'rows');

    % and determine number of factor levels (of un-interacted groups)
    gwi = ones(1, size(b, 2));
    for gfc = 1:size(b, 2)
        gwi(gfc) = numel(unique(b(:, gfc)));
    end
end

% group degrees of freedom
gdf = gwi - 1;
gdft = prod(gdf);

% check number of total groups (fully balanced/crossed?)
ng = max(gib);
if prod(gwi) ~= ng
    error( ...
        'neuroelf:BadArgument', ...
        'Grouping not fully crossed.' ...
    );
end

% one group only
if ng == 1

    % disable group-covariate interactions
    opts.covgrpi = false;
end

% for each group (-interaction)
for gc = 1:ng

    % at least three subjects?
    if sum(gib == gc) < 3
        error( ...
            'neuroelf:BadArgument', ...
            'Grouping requires at least 3 entries per group.' ...
        );
    end
end

% determine the number of subjects in each group
gsz = zeros([gwi, 1]);
ug1 = unique(b(:, 1));

% only one grouping factor (or none)
if ngf <= 1

    % create cell array with lists of indices (to access data)
    gfi = cell(gwi(1), 1);

    % for each level (group)
    for gc1 = 1:gwi(1)

        % find indices
        gfi{gc1} = find(b(:, 1) == ug1(gc1));

        % and store size
        gsz(gc1) = numel(gfi{gc1});
    end

    % normalizing factor is the harmonic mean of group sizes
    gmt = repmat(harmmean(gsz, 1), gwi(1), 1);

% two factors
else

    % create G1-by-G2 cell array for indices
    gfi = cell(gwi(1), gwi(2));

    % find indices
    ug2 = unique(b(:, 2));
    for gc1 = 1:gwi(1)
        for gc2 = 1:gwi(2)

            % store indices and get size
            gfi{gc1, gc2} = find(b(:, 1) == ug1(gc1) & b(:, 2) == ug2(gc2));
            gsz(gc1, gc2) = numel(gfi{gc1, gc2});
        end
    end

    % store normalizing factors
    gm1 = repmat(harmmean(gsz, 1), gwi(1), 1);
    gm2 = repmat(harmmean(gsz, 2), 1, gwi(2));
    gmt = repmat(harmmean(gsz(:), 1), gwi(1), gwi(2));
end
gci = zeros(size(gfi));

% no condition assignments
if isempty(w)

    % assume one-factor with DIM levels
    w = 1:sz(end);
    wib = w;
    wwi = numel(w);

% otherwise
else

    % get unique columns (rows of transpose)
    [g, gi, wib] = unique(w', 'rows');

    % for each factor
    wwi = ones(1, size(w, 1));
    for wfc = 1:size(w, 1)

        % get number of levels
        wwi(wfc) = numel(unique(w(wfc, :)'));
    end
end

% within-degrees of freedom
wdf = wwi - 1;
if any(wdf < 1)
    error( ...
        'neuroelf:BadArgument', ...
        'Within-factors must have at least 2 levels.' ...
    );
end

% total degrees of freedom
wdft = prod(wdf);

% check that all levels are covered
nw = size(w, 2);
if prod(wwi) ~= nw || ...
    max(wib) ~= nw
    error( ...
        'neuroelf:BadArgument', ...
        'Conditions not fully or correctly crossed.' ...
    );
end

% get factor crossings (interactions)
uw1 = unique(w(1, :)');

% for one factor
if nwf == 1

    % find the column for levels 1:L
    wci = zeros(wwi(1), 1);
    for wc1 = 1:wwi(1)
        wci(wc1) = find(w(1, :) == uw1(wc1));
    end

% for two factors
elseif nwf == 2

    % find the columns for level combinations
    wci = zeros(wwi(1), wwi(2));
    uw2 = unique(w(2, :)');
    for wc1 = 1:wwi(1)
        for wc2 = 1:wwi(2)
            wci(wc1, wc2) = find(w(1, :) == uw1(wc1) & w(2, :) == uw2(wc2));
        end
    end

% for three factors
else

    % find the columns for level combinations
    wci = zeros(wwi(1), wwi(2), wwi(3));
    uw2 = unique(w(2, :)');
    uw3 = unique(w(3, :)');
    for wc1 = 1:wwi(1)
        for wc2 = 1:wwi(2)
            for wc3 = 1:wwi(3)
                wci(wc1, wc2, wc3) = find( ...
                    w(1, :) == uw1(wc1) & w(2, :) == uw2(wc2) & w(3, :) == uw3(wc3));
            end
        end
    end
end

% keep weights
if any(strcmpi(class(y), {'double', 'single'}))
    cnan = true;
else
    cnan = false;
end

% no weights
if ~opts.robust && ...
   ~cnan
    yw = 1;

% otherwise, set weights to 1 (individually!)
else
    yw = ones(size(y));
end

% make sure y uses a good datatype
if ~isa(y, 'double')
    y = double(y);
end

% care-about NaN's ?
if cnan

    % set weights for invalid to 0
    yw(isinf(y)) = 0;
    yw(isnan(y)) = 0;

    % also for 0-values in Y
    if opts.nonull
        yw(y == 0) = 0;
    end

    % no invalid values at all and not robust
    if all(yw(:)) == 1 && ...
       ~opts.robust

        % set global weight
        yw = 1;

    % otherwise
    else

        % replace invalid values (Inf/NaN's) with 0
        y(yw == 0) = 0;
    end
end

% generate cell means and msq/msqn arrays (for up to 100 combinations)
cm = zeros([sz(1:nd-2), ng, nw]);
msq = cell(100, 1);
msqn = cell(100, 1);
f = cell(100, 1);
df = zeros(100, 2);

% initialize output counter
msc = 1;

% overall df (number of subjects - number of groups) * (number-within - 1)
odf = (ns - ng) * max(1, (nw - 1));

% remove mean and recompute
tm = mean(mean(y, sd), nd);
y = y - repmat(tm, rmt);

% total sum of squares (without anything computed yet)
tss = sum(sum(y .* y, sd), nd);

% covariates matrix
if ~isempty(opts.cov)

    % per-group covariates
    if opts.covgrpi

        % number of covariates (+ intercept)
        cvis = size(opts.cov, 2) + 1;

        % covariates and intercepts
        Xc = zeros(ns, ng * cvis);

        % fill regression matrix (with transformed covariates)
        for cc = 1:cvis:size(Xc, 2)
            Xci = find(gib == cc);
            Xc(Xci, 1+(cc-1)*cvis:cc*cvis) = ...
                [opts.cov(Xci, :) - ones(numel(Xci), 1) * mean(opts.cov(Xci, :), 1), ...
                 ones(numel(Xci), 1)];
        end

    % fixed-effects covariates
    else

        % transform covariate
        Xc = opts.cov - ones(size(opts.cov, 1), 1) * mean(opts.cov, 1);
    end
end

% remove covariate
if ~isempty(opts.cov)

    % remove variability due to this
    if opts.covregi

        % reshape (fold within and between into one dim)
        y = reshape(y, [prod(sz(1:sd-1)), ns * nw]);

        % calculate betas
        b = calcbetas(Xc, y, 2);

        % individual regression per group
        if opts.covgrpi

            % keep intercepts intact
            b(:, cvis:cvis:end) = 0;

            % adapt DF
            df(msc, :) = [size(Xc, 2) - ng, ng * nw - (size(Xc, 2) + 1)];

        % fixed-effects covariate
        else

            % adapt DF
            df(msc, :) = [size(Xc, 2), ng * nw - (size(Xc, 2) + 1)];
        end

        % remove variance associated with covariate(s) and reshape back
        y = reshape(y - (Xc * b')', sy);

    % explain variance (and store)
    else

        % reshape data (not folding within and between)
        y = reshape(y, [prod(sz(1:sd-1)), ns, nw]);

        % no weights
        if numel(yw) == 1

            % simple mean (across within factors, no interaction!)
            yr = (1 / nw) .* sum(y, 2);

        % with weights
        else

            % weighted mean (across within factors, no interaction!)
            yr = (1 / sum(yw, 2)) .* sum(y .* yw, 2);
        end

        % compute (common variance across within factor levels!)
        b = calcbetas(Xc, yr, 2);

        % for individual covariates
        if opts.covgrpi

            % keep intercept(s) intact
            b(:, cvis:cvis:end) = 0;
        end

        % remove variance associated with covariate(s) and reshape back
        y = reshape(y - repmat((Xc * b')', [1, 1, nw]), sy);
    end

    % compute difference in total SS
    rss = sum(sum(y .* y, sd), nd);

    % store attributable to Covariate
    msq{msc} = (1 / df(msc, 1)) .* (tss - rss);
    msqn{msc} = 'Covariate(s)';

    % replace total (to-be-explained) variance
    tss = rss;

    % adjust DF
    odf = odf - df(msc, 1);

    % next output
    msc = msc + 1;
end

% compute cell means (and initialize counter)
cmgc = 1;

% no/one grouping factor
if ngf <= 1

    % for each group
    for gc1 = 1:gwi(1)

        % store output counter
        gci(gc1) = cmgc;

        % get indices for group (which subjects)
        yi{sd} = gfi{gc1};

        % and access output with target indices (counter)
        yt{sd} = cmgc;

        % for each level-combination (within-interaction)
        for wc1 = 1:nw

            % get and set index
            yi{nd} = wc1;
            yt{nd} = wc1;

            % no weighting
            if numel(yw) == 1

                % compute plain mean (across subjects)
                cm(yt{:}) = (1 / numel(gfi{gc1})) .* sum(y(yi{:}), sd);

            % with weights
            else

                % compute weighted mean
                cm(yt{:}) = (1 / sum(yw(yi{:}))) .* sum(y(yi{:}) .* yw(yi{:}), sd);
            end
        end

        % increment counter
        cmgc = cmgc + 1;
    end

% for two group factors (interaction)
else
    for gc2 = 1:gwi(2)
        for gc1 = 1:gwi(1)
            gci(gc1, gc2) = cmgc;
            yi{sd} = gfi{gc1, gc2};
            yt{sd} = cmgc;
            cmgc = cmgc + 1;
            for wc1 = 1:nw
                yi{nd} = wc1;
                yt{nd} = wc1;
                if numel(yw) == 1
                    cm(yt{:}) = (1 / numel(gfi{gc1, gc2})) .* sum(y(yi{:}), sd);
                else
                    cm(yt{:}) = (1 / sum(yw(yi{:}))) .* sum(y(yi{:}) .* yw(yi{:}), sd);
                end
            end
        end
    end
end

% compute subject SS
if numel(yw) == 1
    sss = (1 / nw) .* sum(sum(y, nd) .^ 2, sd);
else
    sss = (1 ./ sum(sum(yw, nd) .^ 2, sd)) .* ...
        sum(sum(y .* yw, nd) .^ 2, sd);
end

% between-subjects factor/s
yi{nd} = ':';
if ngf > 0
    bss = zeros(size(tm));
    for gc1 = 1:ng
        yi{sd} = gci(gc1);
        bss = bss + (1 / (nw * gsz(gc1))) .* ((gsz(gc1) .* sum(cm(yi{:}), nd)) .^ 2);
    end
    if ngf > 1
        bss1 = zeros(size(tm));
        bss2 = bss1;
        for gc1 = 1:gwi(1)
            yi{nd} = gci(gc1, :);
            rsi(sd) = numel(yi{nd});
            bss1 = bss1 + (1 / (nw * sum(gsz(gc1, :)))) .* (sum((cm(yi{:}) .* ...
                repmat(reshape(gsz(gc1, :), rsi), rmtm)), nd) .^ 2);
        end
        for gc2 = 1:gwi(2)
            yi{nd} = gci(:, gc2);
            rsi(sd) = numel(yi{nd});
            bss2 = bss2 + (1 / nw) .* (sum((1 / sum(gsz(:, gc2))) .* (cm(yi{:}) .* ...
                repmat(reshape(gsz(:, gc2), rsi), rmtm)), nd) .^ 2);
        end
        df(msc, :) = [gdft(1), ns - ng];
        msq{msc} = (1 / gdft(1)) .* bss1;
        msqn{msc} = 'Between-subjects factor 1';
        msc = msc + 1;
        df(msc, :) = [gdft(2), ns - ng];
        msq{msc} = (1 / gdft(2)) .* bss2;
        msqn{msc} = 'Between-subjects factor 2';
        msc = msc + 1;
        df(msc, :) = [gdft, ns - ng];
        msq{msc} = (1 / gdft) .* (bss - (bss1 + bss2));
        msqn{msc} = 'Between-subjects factors 1x2';
        msc = msc + 1;
    else
        df(msc, :) = [gdft, ns - ng];
        msq{msc} = (1 / gdft) .* bss;
        msqn{msc} = 'Between-subjects factor';
        msc = msc + 1;
    end
    df(msc, :) = [ns - ng, odf];
    msq{msc} = (1 / (ns - ng)) .* (sss - bss);
    msqn{msc} = 'Between-subjects error';
    msc = msc + 1;
else
    yi{sd} = 1;
    bss = (ns / nw) .* (sum(cm(yi{:}), nd) .^ 2);
    df(msc, :) = [ns - 1, odf];
    msq{msc} = (1 / (ns - 1)) .* (sss - bss);
    msqn{msc} = 'Between-subjects error';
    msc = msc + 1;
end

% within-subjects factor/s
if nwf > 0
    wss = zeros(size(tm));
    wsu = zeros(size(tm));
    yi{sd} = ':';
    rsi(sd) = numel(gsz);
    for wc1 = 1:nw
        yi{nd} = wc1;
        wss = wss + (1 / sum(gmt(:))) .* ((gmt(1) .* sum(cm(yi{:}), sd)) .^ 2);
        wsu = wsu + (1 / sum(gsz(:))) .* (sum(repmat(reshape(gsz, rsi), rmtm) .* cm(yi{:}), sd) .^ 2);
    end
    if nwf > 1
        error( ...
            'neuroelf:NotYetImplemented', ...
            'Not yet implemented.' ...
        );
    else
        df(msc, :) = [wdft, odf];
        msq{msc} = (1 / wdft) .* wss;
        msqn{msc} = 'Within-subjects factor';
        msc = msc + 1;
    end
else
    wss = 0;
    wsu = 0;
end

% uncorrected version of within-subjects factor for calculating interaction
if ngf > 0 && ...
    nwf > 0
    rmtm(nd) = nw;
    switch (ngf)
        case {1}
            switch (nwf)
                case {1}
                    iss = sum(sum(repmat(reshape(gsz, rsi), rmtm) .* ...
                        (cm .^ 2), sd), nd) - (bss + wsu);
                    df(msc, :) = [(ng - 1) * (nw - 1), odf];
                    msq{msc} = (1 / df(msc, 1)) .* iss;
                    msqn{msc} = 'Interaction between-subjects factor X within-subjects factor';
                    msc = msc + 1;
                case {2}
                case {3}
            end
        case {2}
            switch (nwf)
                case {1}
                case {2}
                case {3}
            end
    end
else
    iss = 0;
end

% within subjects error map
df(msc, :) = [odf, odf];
msq{msc} = (1 / odf) .* (tss - (sss + wsu + iss));
msqn{msc} = 'Within-subjects error';
df(msc+1:end, :) = [];
msq(msc+1:end) = [];
msqn(msc+1:end) = [];

% compute F maps
if nargout > 4
    bsem = findfirst(strcmp(msqn, 'Between-subjects error'));
    for c = 1:msc
        if ~isempty(regexpi(msqn{c}, 'error'))
            f{c} = zeros(size(tm));
            continue;
        end
        if ~isempty(regexpi(msqn{c}, 'within'))
            f{c} = msq{c} ./ msq{msc};
        else
            f{c} = msq{c} ./ msq{bsem};
        end
    end
    f = cat(sd, f{:});
    f(isinf(f) | isnan(f)) = 0;
end

% concatenate
msq = cat(sd, msq{:});
msq(isinf(msq) | isnan(msq)) = 0;

% add global mean to cell means
if nargout > 2
    rmt(sd:nd) = [ng, nw];
    cm = cm + repmat(tm, rmt);
end
