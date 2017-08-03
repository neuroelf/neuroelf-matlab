% neuroelf object
ne = neuroelf;

% data size
szd = [18, 24, 50000]; % number of trials, subjects, and voxels (samples)

% fraction of missing values (random location)
msvals = 2 / szd(1);

% effect size (of interest)
esize = 0.000; % in SDs

% within subject noise levels
wsize = 0.500; % randn
wsiz3 = 0.000; % rand .* randn .^ 3

% between subject noise level
bsize = 0.500; % randn
bsiz3 = 0.000; % rand(S) .* randn .^ 3

% generate data (keeping the original)
sdata = esize + wsize .* randn(szd) + ...
    repmat(bsize .* rand(1, szd(2), szd(3)), szd(1), 1) .* randn(szd);
sdata = sdata + wsiz3 .* randn(szd) .^ 3 + ...
    repmat(bsiz3 .* rand(1, szd(2), szd(3)), szd(1), 1) .* (randn(szd) .^ 3);
sodata = sdata;

% and with missing values
sdata(ceil(numel(sdata) .* rand(ceil(msvals * prod(szd)), 1))) = NaN;

% initialize outputs (temp variables for computation)
ns = zeros(szd(3), 1);
t = ns;

% total number of data points
sn = ns;

% total (ffx) sum and sum-of-squares
st = ns;
sst = ns;

% total (w-rfx) sum and sum-of-squares
ts = ns;
tss = ns;

% total (w-rfx) variance
tv = ns;
tvm = ns;
tvs = ns;

% iterate over samples
for vc = 1:szd(3)

    % subjects
    for c2 = 1:szd(2)

        % and trials
        n = 0;
        s = 0;
        ss = 0;
        for c1 = 1:szd(1)

            % get data
            sc = sdata(c1, c2, vc);

            % don't process invalid data
            if isnan(sc)
                continue;
            end

            % add to sum and sum-of-squares (and counter)
            s = s + sc;
            ss = ss + sc * sc;
            n = n + 1;
        end

        % compute weight: (within-subject) 1/variance (temporary value)
        sv = (n - 1) / (ss - (s * s / n));

        % add to total N, SUM, and SUM-OF-SQUARES (for FFX)
        sn(vc) = sn(vc) + n;
        st(vc) = st(vc) + s;
        sst(vc) = sst(vc) + ss;

        % add to weighted estimate of total N
        ns(vc) = ns(vc) + n * sv;

        % add to (weighted) sum
        t(vc) = t(vc) + sv * s / n;

        % and weighted sum of squares
        ts(vc) = ts(vc) + sv * s;
        tss(vc) = tss(vc) + sv * s * s;

        % and also keep track of total variance
        tv(vc) = tv(vc) + sv;
        tvm(vc) = max(tvm(vc), sv);
        tvs(vc) = tvs(vc) + sv * sv;
    end
end

% overall mean
ns = ns ./ tv;
tm = t ./ tv;
tf = tv ./ szd(2);
tsv = ((tf .* tss - ts .* ts ./ szd(2)) ./ (tf .* tf .* (szd(2) - 1))) ./ (ns .* ns);
ttm = squeeze(mean(mean(sodata), 2));
ttsv = squeeze(var(mean(sodata), [], 2));
ttt = ttm ./ sqrt(ttsv ./ szd(2));

% ffx variance
ttv = (sst - (st .* st ./ sn)) ./ (sn - 1);
tttv = var(reshape(sodata, szd(1) * szd(2), szd(3)))';

% rfx variance-variance
tvv = (tvs - (tv .* tv ./ szd(2))) ./ (szd(2) - 1);

% estimate w-rfx DF and error
df = tv ./ (tv ./ szd(2) + 1.96 .* sqrt(tvv ./ szd(2)));
% df = tv ./ tvm;
se = sqrt(tsv ./ df);

% recompute as nominal z-values
wrfxz = -sign(tm) .* ne.sdist('norminv', ne.sdist('tcdf', -abs(tm ./ se), df - 1), 0, 1);
grfxz = -sign(ttm) .* ne.sdist('norminv', ne.sdist('tcdf', -abs(ttt), szd(2) - 1), 0, 1);

% visualize results, t-scores
fh = figure;

% means
ah = subplot(2, 2, 1);
r = corrcoef(ttm, tm);
scatter(ttm, tm);
title(sprintf('mean comparison (r = %5.3f)', r(2)));
xlabel('mean (OLS, all data)');
ylabel('mean (wRFX)');
xl = get(ah, 'XLim');
yl = get(ah, 'YLim');
xlyl = [min(xl(1), yl(1)), max(xl(2), yl(2))];
hold(ah, 'on');
line(xlyl, xlyl);
set(ah, 'XLim', xlyl, 'YLim', xlyl);

% variances
ah = subplot(2, 2, 2);
r = corrcoef(log(1+ttsv), log(1+tsv));
scatter(log(1+ttsv), log(1+tsv));
title(sprintf('variance comparison (r = %5.3f)', r(2)));
xlabel('log(1+variance) (OLS, all data)');
ylabel('log(1+variance) (wRFX)');
xl = get(ah, 'XLim');
yl = get(ah, 'YLim');
xlyl = [min(xl(1), yl(1)), max(xl(2), yl(2))];
hold(ah, 'on');
line(xlyl, xlyl);
set(ah, 'XLim', xlyl, 'YLim', xlyl);

% z-scores
ah = subplot(2, 2, 3);
scatter(grfxz, wrfxz, 24, 4 .* (ceil(df) - min(ceil(df))));
title('Z-score comparison');
xlabel('z-score (OLS, all data)');
ylabel('z-score (wRFX)');
xl = get(ah, 'XLim');
yl = get(ah, 'YLim');
xlyl = [min(xl(1), yl(1)), max(xl(2), yl(2))];
hold(ah, 'on');
line(xlyl, xlyl);
set(ah, 'XLim', xlyl, 'YLim', xlyl);

% wRFX z-scores
ah = subplot(2, 2, 4);
hist(wrfxz, 100);
z025 = sum(wrfxz > -ne.sdist('norminv', 0.025, 0, 1)) ./ numel(wrfxz);
z001 = sum(wrfxz > -ne.sdist('norminv', 0.001, 0, 1)) ./ numel(wrfxz);
title(sprintf('wRFX z-scores (%6.3f%% p<.025, %6.3f%% p<.001)', 100 * z025, 100 * z001));
xlabel('Z-value');
ylabel('count');
yl = get(ah, 'YLim');
hold(ah, 'on');
line([1.96, 1.96], [0, yl(2)]);
line([3.09, 3.09], [0, yl(2)]);
