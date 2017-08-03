function [H,f,s,c,tau,w] = LSTransform(t,h,TH,Tr,ofac,hifac)
%LSTRANSFORM  Fills up missing samples in time-series with estimates.
%   [H, F, S, C, TAU, W] = LSTRANSFORM(T, H, TH, TR, OFAC, HIFAC) takes
%   the signal in H and samples new data points given the frequecy content
%   of H at points TH. TR is the sampling resolution (this argument is
%   currently not used), OFAC is the oversampling factor, and HIFAC is the
%   highest allowed frequency.
%
%   H must be NrOfTimePoints-by-NrOfSamples matrix.
%
%   Example:
%
%   rdata = randn(512, 200);
%   allidx = (1:size(rdata, 1))';
%   badidx = ceil(numel(allidx) .* rand(10, 1));
%   goodidx = allidx;
%   goodidx(badidx) = [];
%   repdata = LSTRANSFORM(goodidx, rdata(goodidx, :), allidx, 1, 4, 1);

%Input t is a column vector listing the time points for which observations
%are present.  Input h is a matrix with observations in columns and the
%number of rows equals the number the time points.  For our purposes number
%of voxels = number of columns.  Ofac = oversampling frequency (generally
%>=4), hifac = highest frequency allowed.  hifac = 1 means 1*nyquist limit
%is highest frequency sampled.  
%Lasted edited:  Anish Mitra, October 25 2012

% double precision
D = double(h);

% number of time points
N = size(D, 1);

% total time span
t = t(:);
T = max(t) - min(t);

% calculate sampling frequencies
f = (1 / (T * ofac) : 1 / (T * ofac) : hifac * N / (2 * T))';

% angular frequencies and constant offsets
w = 2 * pi * f;
wt = w * t';
tau = atan2(sum(sin(2 .* wt), 2), sum(cos(2 * wt), 2)) ./ (2 * w);
wtau = wt - repmat(w .* tau, 1, length(t));

% spectral power sin and cosine terms
cterm = cos(wtau);
sterm = sin(wtau);

% compute numerator and denominator for cosines
numerator = cterm * D;
denominator = sum(cterm .* cterm,2);
c = diag(1 ./ denominator) * numerator;

% repeat the above for Sine term
numerator = sterm * D;
denominator = sum(sterm .* sterm,2);
s = diag(1 ./ denominator) * numerator;

% the inverse function to re-construct the original time series
prod = TH(:) * w';
H = sin(prod) * s + cos(prod) * c;

% normalize the reconstructed spectrum, needed when ofac > 1
H = H * diag(std(h) ./ std(H));
