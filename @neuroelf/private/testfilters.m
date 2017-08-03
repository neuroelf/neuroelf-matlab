function [fh, t] = testfilters(opts)

% options
opts.bworder = 5;
opts.hfnfactor = 0.5;
opts.hfntype = 'n'; % 'n' or 'u'
opts.hp = 0.016;
opts.lfmaxfreq = 0.01;
opts.lfnfactor = 2;
opts.lp = 0.18;
opts.numcases = 100000;
opts.numdpts = 200;
opts.tr = 2000;

% output
t = zeros(1, 4);

% create signal
hs = hrf('twogamna', 0.001);
hf = 1 / std(hs);

% sampling
h = ceil(opts.tr * rand(1, opts.numcases));
h = h + opts.tr .* ceil(0.8 * opts.numdpts .* (rand(1, opts.numcases) - 1.25));
h = (0:opts.tr:(opts.tr * (opts.numdpts - 1)))' * ones(1, opts.numcases) + ones(opts.numdpts, 1) * h;
hx = (h > 0 & h <= numel(hs));
h(hx) = hf .* hs(h(hx));
h(~hx) = 0;

% create random values (high-frequency noise)
if opts.hfntype == 'n'
    r = opts.hfnfactor .* randn(opts.numdpts, opts.numcases);
else
    r = (opts.hfnfactor * 3.4641016) .* rand(opts.numdpts, opts.numcases);
end

% create additional low-frequency noise
fq = 2 * pi * opts.lfmaxfreq;
mfq = (opts.numdpts - 0.5) * fq;
s = opts.lfnfactor .* ztrans( ...
    sin((0:fq:mfq)' * rand(1,opts.numcases) + ones(opts.numdpts, 1) * (150 .* randn(1,opts.numcases))));

% combine
rs = h + r + s;

% transpose h and r (for later)
h = h';
r = r';

% filter with butterworth
[b, a] = butterworth([opts.hp, opts.lp], opts.bworder);
tic;
bf = filter(b, a, rs);
t(1) = toc;
[~, chbf]  = cov_nd(h, bf');
[~, crbf]  = cov_nd(r, bf');

% also create "reverse-filtered" to improve removal of LF
bfr = filter(b, a, rs(end:-1:1, :));
[~, chbfr] = cov_nd(h, bfr(end:-1:1, :)');
[~, crbfr] = cov_nd(r, bfr(end:-1:1, :)');

% and mixed version
midp = round(opts.numdpts / 2) + 1;
bfr = [bfr(end:-1:midp, :); bf(midp:end, :)];
[~, chbfx] = cov_nd(h, bfr');
[~, crbfx] = cov_nd(r, bfr');

% use tempfilter (smoothing) for approach
tic;
tfs = tempfilter(rs, struct('temphp', round(2/opts.hp), 'templp', round(1/opts.lp)));
t(2) = toc;
[~, chtfs]  = cov_nd(h, tfs');
[~, crtfs]  = cov_nd(r, tfs');

% then try DCT
tic;
tfd = tempfilter(rs, struct('tempdct', 1/opts.hp, 'templp', round(1/opts.lp)));
t(3) = toc;
[~, chtfd]  = cov_nd(h, tfd');
[~, crtfd]  = cov_nd(r, tfd');

% and straight sin/cos regression
tic;
tff = tempfilter(rs, struct('tempsc', floor((opts.numdpts + 1) * opts.hp), 'templp', round(1/opts.lp)));
t(4) = toc;
[~, chtff]  = cov_nd(h, tff');
[~, crtff]  = cov_nd(r, tff');

% compute medians and show histograms
fh = (1 / opts.numcases) .* ([ ...
    histcount(chbf, -1, 0.9995, 0.005); ...
    histcount(crbf, -1, 0.9995, 0.005); ...
    histcount(chbfr, -1, 0.9995, 0.005); ...
    histcount(crbfr, -1, 0.9995, 0.005); ...
    histcount(chbfx, -1, 0.9995, 0.005); ...
    histcount(crbfx, -1, 0.9995, 0.005); ...
    histcount(chtfs, -1, 0.9995, 0.005); ...
    histcount(crtfs, -1, 0.9995, 0.005); ...
    histcount(chtfd, -1, 0.9995, 0.005); ...
    histcount(crtfd, -1, 0.9995, 0.005); ...
    histcount(chtff, -1, 0.9995, 0.005); ...
    histcount(crtff, -1, 0.9995, 0.005)])';
