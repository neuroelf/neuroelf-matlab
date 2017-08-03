% compute relative power of frequency of PRT onsets

% use NeuroElf
n = neuroelf;

% load VTC
vtc = xff('*.vtc');
vmp = vtc.MeanStd;
vtcd = double(vtc.VTCData(:, :, :, :));
vsz = size(vtcd);
tr = vtc.TR;

% regress out noise
[~, fX] = n.tempfilter(randn(vsz(1), 1), struct('tempdct', 0.3 * vsz(1)));
X = [ones(vsz(1), 1), n.ztrans(vtc.RunTimeVars.MotionParameters), fX];
[b, ixx, ptc] = n.calcbetas(X, vtcd, 1);
vtcd = vtcd - ptc;
vtc.ClearObject;

% load PRT
prt = xff('*.prt');

% find condition with onsets to consider
onsn = zeros(numel(prt.Cond), 1);
for pc = 1:numel(onsn)
    onsn(pc) = size(prt.Cond(pc).OnOffsets, 1);
end
[~, p] = max(onsn);
cond = prt.Cond(p);
prt.ClearObject;

% onsets, number of onsets, and first onset
oo = cond.OnOffsets;
no = size(oo, 1);
o1 = oo(1);

% compute average distance between 
od = mean(diff(oo(:, 1)));
ol = o1 + no * od;

% resample VTCdata from first onset to after last onset
ss = (ol/tr)/512;
vtcd = n.flexinterpn_method(vtcd, ...
    [Inf, Inf, Inf, Inf; 1 + o1/tr, 1, 1, 1; ss, 1, 1, 1; 1 + (o1+ol)/tr - ss, vsz(2:4)], ...
    'cubic');

% compute power along first dim
ptc = (1/512) .* abs(fft(vtcd, [], 1));

% and power of experiment (N repetitions)
pptc = squeeze(ptc(1 + no, :, :, :) ./ sum(ptc(2:4*no, :, :, :)));

% also compute average time course of experiment
avtcd = squeeze(mean(reshape(vtcd, [32, 16, vsz(2:4)]), 2));

% and correlation with that
[~, cvtc] = n.cov_nd(permute(repmat(avtcd, [16, 1, 1, 1]), [2, 3, 4, 1]), permute(vtcd, [2, 3, 4, 1]));

% store in VMP
vmp.Map(1).Type = 2;
vmp.Map(1).LowerThreshold = 0.025;
vmp.Map(1).UpperThreshold = 0.05;
vmp.Map(1).Name = 'relative POWER';
vmp.Map(1).VMPData = single(pptc);
vmp.Map(2).Type = 2;
vmp.Map(2).LowerThreshold = 0.25;
vmp.Map(2).UpperThreshold = 0.5;
vmp.Map(2).Name = 'correlation with average signal';
vmp.Map(2).VMPData = single(cvtc);
vmp.SaveAs;
