
% test filters

% create data array (10000 time courses)
r = zeros(312,10000);

% fill with sines
for c = 1:10000
    r(:,c) = 4 .* sin((c/40000) .* (0:311) + (2*pi*rand(1)))';
end

% store as sines array
rsin = r;

% generate random data
rr = randn(size(r));

% add together
r = rsin + exp(-1) .* rr;

% build design matrices
[d, Xd] = tempfilter(zeros(312, 1), struct('tempdct', 52));
[d, Xl] = tempfilter(zeros(312, 1), struct('temppoly', 14, 'orthpoly', false));
[d, Xp] = tempfilter(zeros(312, 1), struct('temppoly', 14, 'orthpoly', true));
[d, Xs] = tempfilter(zeros(312, 1), struct('tempsc', 6));
Xd(:,end+1) = 1;
Xl(:,end+1) = 1;
Xp(:,end+1) = 1;
Xs(:,end+1) = 1;

% filter with temphp (same as 12 sin+cos)
rf = tempfilter(r, struct('temp', true, 'tempdt', false, 'temphp', 52));
rf = rf - ones(312, 1) * mean(rf);
rfo = tempfilter(r, struct('temp', true, 'tempdt', false, 'temphp', 52, 'temphpr', true));
rfo = rfo - ones(312, 1) * mean(rfo);

% then filter using design matrices
bd = (Xd' * Xd) \ (Xd' * r);
bl = (Xl' * Xl) \ (Xl' * r);
bp = (Xp' * Xp) \ (Xp' * r);
bs = (Xs' * Xs) \ (Xs' * r);

% get new data
rd = r - Xd * bd;
rl = r - Xl * bl;
rp = r - Xp * bp;
rs = r - Xs * bs;

% compute correlation and Z-score
[cv, crr] = cov_nd(rr', r');
[cv, crf] = cov_nd(rr', rf');
[cv, crfo] = cov_nd(rr', rfo');
[cv, crd] = cov_nd(rr', rd');
[cv, crl] = cov_nd(rr', rl');
[cv, crp] = cov_nd(rr', rp');
[cv, crs] = cov_nd(rr', rs');

% plot
figure;
subplot(4, 2, 1); plot(fisherr2z(crr, false, 312)); title('Raw');
subplot(4, 2, 2); plot(fisherr2z([crf, crfo], false, 312)); title('HP-filtered');
subplot(4, 2, 3); plot(fisherr2z(crs, false, 312)); title('Fourier (s+c)');
subplot(4, 2, 4); plot(fisherr2z(crd, false, 312)); title('DCT');
subplot(4, 2, 5); plot(fisherr2z(crl, false, 312)); title('Legendre');
subplot(4, 2, 6); plot(fisherr2z(crp, false, 312)); title('Orthogonal polynomials');
subplot(4, 2, [7,8]); p = plot(flexinterpn(fisherr2z([crr,crf,crd,crl,crp,crs], false, 312), ...
    [inf,inf;1,1;1,1;10000,6], {smoothkern(50), [0;1;0]}, {1,1})); title('smoothed curves');
legend(p, 'Raw', 'HP-filtered', 'DCT', 'Legendre', 'Orth. poly', 'SIN+COS');

% now take the mixed data r as a model for the noise and add HRF signal
prt = xff('new:prt');
prt.AddCond('block', ...
    [ 10000,  25000; ...
      52000,  68000; ...
      91000, 105000; ...
     134000, 147000; ...
     172000, 189000; ...
     210000, 224000; ...
     249000, 265000; ...
     293000, 307000; ...
     330000, 345000; ...
     372000, 388000; ...
     410000, 424000; ...
     453000, 466000; ...
     491000, 507000; ...
     528000, 543000; ...
     567000, 583000; ...
     610000, 622000]);
sdm = prt.CreateSDM(struct('prtr', 2000, 'nvol', 312, 'rcond', []));
sdmd = repmat(sdm.SDMMatrix(:, 1), 1, 10000);
rsig = r + exp(-1) .* sdmd;

% and recompute filtered data
rsigf = tempfilter(rsig, struct('temp', true, 'tempdt', false, 'temphp', 52));
rsigf = rsigf - ones(312, 1) * mean(rsigf);
rsigfo = tempfilter(rsig, struct('temp', true, 'tempdt', false, 'temphp', 52, 'temphpr', true));
rsigfo = rsigfo - ones(312, 1) * mean(rsigfo);
bd = (Xd' * Xd) \ (Xd' * rsig);
bl = (Xl' * Xl) \ (Xl' * rsig);
bp = (Xp' * Xp) \ (Xp' * rsig);
bs = (Xs' * Xs) \ (Xs' * rsig);

% get new data
rsigd = rsig - Xd * bd;
rsigl = rsig - Xl * bl;
rsigp = rsig - Xp * bp;
rsigs = rsig - Xs * bs;

% compute correlation and Z-score
[cv, crr] = cov_nd(sdmd', rsig');
[cv, crf] = cov_nd(sdmd', rsigf');
[cv, crfo] = cov_nd(sdmd', rsigfo');
[cv, crd] = cov_nd(sdmd', rsigd');
[cv, crl] = cov_nd(sdmd', rsigl');
[cv, crp] = cov_nd(sdmd', rsigp');
[cv, crs] = cov_nd(sdmd', rsigs');

% plot
figure;
subplot(4, 2, 1); plot(fisherr2z(crr, false, 312)); title('Raw');
subplot(4, 2, 2); plot(fisherr2z([crf, crfo], false, 312)); title('HP-filtered');
subplot(4, 2, 3); plot(fisherr2z(crs, false, 312)); title('Fourier (s+c)');
subplot(4, 2, 4); plot(fisherr2z(crd, false, 312)); title('DCT');
subplot(4, 2, 5); plot(fisherr2z(crl, false, 312)); title('Legendre');
subplot(4, 2, 6); plot(fisherr2z(crp, false, 312)); title('Orthogonal polynomials');
subplot(4, 2, [7,8]); p = plot(flexinterpn(fisherr2z([crr,crf,crd,crl,crp,crs], false, 312), ...
    [inf,inf;1,1;1,1;10000,6], {smoothkern(50), [0;1;0]}, {1,1})); title('smoothed curves');
legend(p, 'Raw', 'HP-filtered', 'DCT', 'Legendre', 'Orth. poly', 'SIN+COS');

% and then compute Z scores with filtered designs
rsdmf = tempfilter(sdmd(:, 1), struct('temp', true, 'tempdt', false, 'temphp', 52));
rsdmf = rsdmf - ones(312, 1) * mean(rsdmf);
rsdmf = repmat(rsdmf, 1, 10000);
rsdmfo = tempfilter(sdmd(:, 1), struct('temp', true, 'tempdt', false, 'temphp', 52, 'temphpr', true));
rsdmfo = rsdmfo - ones(312, 1) * mean(rsdmfo);
rsdmfo = repmat(rsdmfo, 1, 10000);
bd = (Xd' * Xd) \ (Xd' * sdmd);
bl = (Xl' * Xl) \ (Xl' * sdmd);
bp = (Xp' * Xp) \ (Xp' * sdmd);
bs = (Xs' * Xs) \ (Xs' * sdmd);
rsdmd = sdmd - Xd * bd;
rsdml = sdmd - Xl * bl;
rsdmp = sdmd - Xp * bp;
rsdms = sdmd - Xs * bs;
[cv, crr] = cov_nd(sdmd', rsig');
[cv, crf] = cov_nd(rsdmf', rsigf');
[cv, crfo] = cov_nd(rsdmfo', rsigfo');
[cv, crd] = cov_nd(rsdmd', rsigd');
[cv, crl] = cov_nd(rsdml', rsigl');
[cv, crp] = cov_nd(rsdmp', rsigp');
[cv, crs] = cov_nd(rsdms', rsigs');

% plot
figure;
subplot(4, 2, 1); plot(fisherr2z(crr, false, 312)); title('Raw');
subplot(4, 2, 2); plot(fisherr2z([crf, crfo], false, 312)); title('HP-filtered');
subplot(4, 2, 3); plot(fisherr2z(crs, false, 312)); title('Fourier (s+c)');
subplot(4, 2, 4); plot(fisherr2z(crd, false, 312)); title('DCT');
subplot(4, 2, 5); plot(fisherr2z(crl, false, 312)); title('Legendre');
subplot(4, 2, 6); plot(fisherr2z(crp, false, 312)); title('Orthogonal polynomials');
subplot(4, 2, [7,8]); p = plot(flexinterpn(fisherr2z([crr,crf,crd,crl,crp,crs], false, 312), ...
    [inf,inf;1,1;1,1;10000,6], {smoothkern(20), [0;1;0]}, {1,1})); title('smoothed curves');
legend(p, 'Raw', 'HP-filtered', 'DCT', 'Legendre', 'Orth. poly', 'SIN+COS');
