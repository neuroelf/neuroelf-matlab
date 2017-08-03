% test interpolation with (down-) weighting (of outliers)

% settings
ns = 100;       % number of time points
nn = 5;         % number of noise samples
na = 5;         % noise amplitude
sk = 2.5;       % smoothing kernel
im = 'lanczos8';% interpolation method

% create slightly smoothed test signal
r = flexinterpn(randn(ns, 1), [inf; 1; 1; ns], smoothkern(sk), 1);

% copy signal and add noise
rn = r;
ri = ceil(ns * rand(nn, 1));
n = ztrans(randn(nn, 1));
rn(ri) = rn(ri) + na .* n;

% interpolate this signal twice (at .5 shifts)
rni  = flexinterpn_method(rn,  [inf; 0.5; 1; ns+1], im);
rni2 = flexinterpn_method(rni, [inf; 1.5; 1; ns+1], im);

% find weights (using a DCT-based filtering matrix)
[rf, rx, rw] = tempfilter(rn, struct('tempdct', 30, 'trobust', true));

% interpolate data multiplied by weights, then divide by interpolated weights
rnwi  = flexinterpn_method(rn   .* rw,   [inf; 0.5; 1; ns+1], im);
rnww  = flexinterpn_method(        rw,   [inf; 0.5; 1; ns+1], im);
rnwi  = rnwi ./ rnww;
rnwi2 = flexinterpn_method(rnwi .* rnww, [inf; 1.5; 1; ns+1], im);
rnww2 = flexinterpn_method(        rnww, [inf; 1.5; 1; ns+1], im);
rnwi2 = rnwi2 ./ rnww2;

% plot data
figure;
plot([r, rn + 0.1, rni2+0.2, rnwi2+0.3, r - rnwi2, rw-3, rnww2-2]);
legend('Signal w/o noise', 'Signal w/ noise', 'Signal w/ noise, 2 ip''s', ...
    'Signal w/ noise+weights, 2 ip''s', 'Deviation', ...
    'Original weights', 'Interpolated weights');
