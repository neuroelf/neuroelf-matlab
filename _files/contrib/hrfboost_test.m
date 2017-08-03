% test program for HRF-boosting algorithms

% set design type (e = event-related, dur=0; d = duration, b = block)
% designtype = 'b'; % block design: every 45 seconds for 30 seconds
% designtype = 'd'; % duration-based design: every ~12 seconds with duration
designtype = 'e'; % event-related design: every ~6 seconds with dur = 0

% for duration-based designs
duration = 4; % in seconds

% number of simulated time courses
numberoftcs = 10000;

% length of time courses (in seconds)
tclength = 800;

% amount of noise (as opposed to signal being unit = 1)
noise = 1.5;

% rescale HRF (not in SPM!)
scalehrf = true;

% initialize run-time
tic;

% get canonical HRF basis functions (in 1/16s resolution and 2 derivatives)
if scalehrf
    nt = 'a';
else
    nt = 's';
end
h = hrf('twogamma', [], [], [], [], [], [], [], [], 2, nt);

% generate resampled version
hs = h(1:32:end, :);

% initialize stimulus function (300s in 16 bins/s resolution)
sflength = tclength * 16;
sf = zeros(sflength, 1);

% depending on design type
switch (designtype)

    % block design
    case {'b'}

        % generate indices for block design
        sfi = 0:720:(sflength-1);

        % generate time-bin resolution stimulus function
        for bc = 1:480
            sf(sfi + bc) = 1;
        end

    % duration-based design
    case {'d'}

        % generate indices for duration stimulus function (around every 12s)
        sfi = ceil(cumsum(144 + 96 * rand(ceil(tclength / 8), 1))) - 144;

        % remove onsets within and after the last 4 seconds
        sfi(sfi > (sflength - 80)) = [];

        % generate time-bin resolution stimulus function
        for bc = 1:(duration * 16)
            sf(sfi + bc) = 1;
        end

    % event-related design
    case {'e'}

        % generate indices for ER (dur=0) stimulus function (around every 6s)
        sfi = ceil(cumsum(48 + 96 * rand(ceil(tclength / 4), 1))) - 48;

        % remove onsets within and after the last 4 seconds
        sfi(sfi > (sflength - 64)) = [];

        % generate time-bin resolution stimulus function
        sf(sfi) = 1;

    otherwise
        error( ...
            'neuroelf:BadSetting', ...
            'Unknown design type: %s.', ...
            designtype ...
        );
end

% convolve stimulus function (with canonical HRF)
csf = conv(sf, h(:, 1));
csf((sflength+1):end) = [];

% generate regressors (convolving stimulus function with HRF + derivatives)
rbf1 = conv(h(:, 1), sf);
rbf2 = conv(h(:, 2), sf);
rbf3 = conv(h(:, 3), sf);

% scale noise to be relative to amplitude of HRF basis function (max. conv)
noise = noise * max(rbf1(:, 1));

% generate 10000 sample time courses (with 300s in 16 bins/s resolution)
if noise <= 0
    tc = zeros(sflength, numberoftcs);
else
    tc = noise .* randn(sflength, numberoftcs);
end

% store first 35 seconds of all basis functions
hrs = zeros(560, numberoftcs);

% generate parameters
timepospeak =  4 + rand(numberoftcs, 1) + rand(numberoftcs, 1); % time to positive peak between 4 and 6 secs
timenegpeak = 13 + 2 .* rand(numberoftcs, 1) + 2 .* rand(numberoftcs, 1); % time to negative peak between 13 and 17 secs
ratioofpeak =  6 + 0.5  .* randn(numberoftcs, 1); % ratio between positive and negative peak around 6 +/- 0.5
hrfonsshift =      0.25 .* randn(numberoftcs, 1); % onset of HRF 0 +/- 0.25
disppospeak =  1 + 0.1  .* randn(numberoftcs, 1); % dispersion of positive peak 1 +/- 0.1
dispnegpeak =  1 + 0.05 .* randn(numberoftcs, 1); % dispersion of negative peak 1 +/- 0.15

% fill time courses (with differently shaped HRF functions)
for c = 1:numberoftcs

    % generate different HRF shape
    hr = hrf( ...
        'twogamma',     ... % two-gamma
        [],             ... % default sampling frequency!
        timepospeak(c), ...
        timenegpeak(c), ...
        ratioofpeak(c), ...
        hrfonsshift(c), ...
        disppospeak(c), ...
        dispnegpeak(c), ...
        [], 0, nt);

     % store for later inspection
     hrs(1:min(numel(hr), 560), c) = hr(1:min(numel(hr), 560));

     % convolve with stimulus function
     hc = conv(sf, hr);

     % store in time courses array
     tc(:, c) = tc(:, c) + hc(1:sflength);
end

% create TR-sampled version (150 TRs with 2s each)
tcs = tc(1:32:end, :);

% generate design matrix (just one condition, truncated to 300s * 16 bins)
X = [rbf1, rbf2, rbf3, ones(numel(rbf1), 1)];
X((sflength+1):end, :) = [];

% create resampled version as well
Xs = X(1:32:end, :);

% also create version with canonical-only HRF
Xh = X(:, [1, 4]);
Xhs = Xs(:, [1, 4]);

% compute regression (full and resampled data, full and HRF-only models)
b   = pinv(X'   * X)   * X'   * tc;
bs  = pinv(Xs'  * Xs)  * Xs'  * tcs;
bh  = pinv(Xh'  * Xh)  * Xh'  * tc;
bhs = pinv(Xhs' * Xhs) * Xhs' * tcs;

% compute how much variance is explained in each combination
tcv    = var(tc);              % full variance
tcsv   = var(tcs);             % resampled variance
tcrv   = var(tc  - X   * b);   % residual of full data after full model
tcsrv  = var(tcs - Xs  * bs);  % residual of resampled data after full model
tchrv  = var(tc  - Xh  * bh);  % residual of full data after HRF-only model
tcshrv = var(tcs - Xhs * bhs); % residual of resampled data after HRF-only model

% histograms over residuals (fractional)
figure;
ax = subplot(2, 2, 1);
hist(tcrv ./ tcv, 0:0.005:1);
set(ax, 'XLim', [0, 1], 'YLim', [0, 1000]);
title(ax, sprintf('FRV full model/full data (median = %.3g)', median(tcrv ./ tcv)));

ax = subplot(2, 2, 2);
hist(tcsrv ./ tcsv, 0:0.005:1);
set(ax, 'XLim', [0, 1], 'YLim', [0, 1000]);
title(ax, sprintf('FRV full model/resampled data (median = %.3g)', median(tcsrv ./ tcsv)));

ax = subplot(2, 2, 3);
hist(tchrv ./ tcv, 0:0.005:1);
set(ax, 'XLim', [0, 1], 'YLim', [0, 1000]);
title(ax, sprintf('FRV HRF-only model/full data (median = %.3g)', median(tchrv ./ tcv)));

ax = subplot(2, 2, 4);
hist(tcshrv ./ tcsv, 0:0.005:1);
set(ax, 'XLim', [0, 1], 'YLim', [0, 1000]);
title(ax, sprintf('FRV HRF-only model/resampled data (median = %.3g)', median(tcshrv ./ tcsv)));

% compute HRF-boost with different settings
bcal   = hrfboost(b(1:3,  :)', struct('bf', h,  'comp', 'boost')); % Calhoun
bscal  = hrfboost(bs(1:3, :)', struct('bf', h,  'comp', 'boost')); % Calhoun
bscals = hrfboost(bs(1:3, :)', struct('bf', hs, 'comp', 'boost')); % Calhoun
bmax   = hrfboost(b(1:3,  :)', struct('bf', h,  'comp', 'max'));   % MAX ampl.
bsmax  = hrfboost(bs(1:3, :)', struct('bf', h,  'comp', 'max'));   % MAX ampl.
bsmaxs = hrfboost(bs(1:3, :)', struct('bf', hs, 'comp', 'max'));   % MAX ampl.
bauc   = hrfboost(b(1:3,  :)', struct('bf', h,  'comp', 'auc'));   % AUC
bsauc  = hrfboost(bs(1:3, :)', struct('bf', h,  'comp', 'auc'));   % AUC
bsaucs = hrfboost(bs(1:3, :)', struct('bf', hs, 'comp', 'auc'));   % AUC
baup   = hrfboost(b(1:3,  :)', struct('bf', h,  'comp', 'pos'));   % pos. AUC
bsaup  = hrfboost(bs(1:3, :)', struct('bf', h,  'comp', 'pos'));   % pos. AUC
bsaups = hrfboost(bs(1:3, :)', struct('bf', hs, 'comp', 'pos'));   % pos. AUC

% histograms over HRf-only betas
fig = figure;
set(fig, 'Position', [200, 80, 800, 1000]);
ax = subplot(5, 3, 1);
hist(bh(1, :), 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('HRF-only (median = %.3g, std = %.3g)', ...
    median(bh(1, :), 2), std(bh(1, :), [], 2)));
ax = subplot(5, 3, 3);
hist(bhs(1, :), 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('resampled data+hrf (median = %.3g, std = %.3g)', ...
    median(bhs(1, :), 2), std(bhs(1, :), [], 2)));

% histograms for Calhoun-based boosted values
ax = subplot(5, 3, 4);
hist(bcal, 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('Calhoun-boost (median = %.3g, std = %.3g)', median(bcal), std(bcal)));
ax = subplot(5, 3, 5);
hist(bscal, 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('resampled data (median = %.3g, std = %.3g)', median(bscal), std(bscal)));
ax = subplot(5, 3, 6);
hist(bscals, 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('resampled data+hrf (median = %.3g, std = %.3g)', median(bscals), std(bscals)));

% histograms for MAX-amplitude boosted values
ax = subplot(5, 3, 7);
hist(bmax, 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('MAX-ampl.-boost (median = %.3g, std = %.3g)', median(bmax), std(bmax)));
ax = subplot(5, 3, 8);
hist(bsmax, 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('resampled data (median = %.3g, std = %.3g)', median(bsmax), std(bsmax)));
ax = subplot(5, 3, 9);
hist(bsmaxs, 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('resampled data+hrf (median = %.3g, std = %.3g)', median(bsmaxs), std(bsmaxs)));

% histograms for full AUC-boosted values
ax = subplot(5, 3, 10);
hist(bauc, 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('AUC-boost (median = %.3g, std = %.3g)', median(bauc), std(bauc)));
ax = subplot(5, 3, 11);
hist(bsauc, 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('resampled data (median = %.3g, std = %.3g)', median(bsauc), std(bsauc)));
ax = subplot(5, 3, 12);
hist(bsaucs, 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('resampled data+hrf (median = %.3g, std = %.3g)', median(bsaucs), std(bsaucs)));

% histograms for positive-only AUC-boosted values
ax = subplot(5, 3, 13);
hist(baup, 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('pos-AUC-boost (median = %.3g, std = %.3g)', median(baup), std(baup)));
ax = subplot(5, 3, 14);
hist(bsaup, 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('resampled data (median = %.3g, std = %.3g)', median(bsaup), std(bsaup)));
ax = subplot(5, 3, 15);
hist(bsaups, 0.5:0.01:2);
set(ax, 'XLim', [0.5, 2], 'YLim', [0, 500]);
title(ax, sprintf('resampled data+hrf (median = %.3g, std = %.3g)', median(bsaups), std(bsaups)));

% print run-time
toc
