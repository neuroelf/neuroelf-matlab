function [num, den] = butterworth(ffreq, forder)
% butterworth  - generate butterworth bandpass filter coefficients
%
% FORMAT:       [num, den] = butterworth(ffreq, forder)
%
% Input fields:
%
%       ffreq       1x2 double frequencies ([high-pass, low-pass])
%       forder      1x1 double Nth order filter
%
% Output fields:
%
%       num         numerator for filter
%       den         denuminator for filter
%
% Note: the cut-off frequencies, ffreq, must be 0.0 < ffreq < 1.0, with
% 1.0 corresponding to half the sample rate.
%
% See also butter.

%   Author(s): J.N. Little, 1-14-87
%   	   J.N. Little, 1-14-88, revised
%   	   L. Shure, 4-29-88, revised
%   	   T. Krauss, 3-24-93, revised

%   References:
%     [1] T. W. Parks and C. S. Burrus, Digital Filter Design,
%         John Wiley & Sons, 1987, chapter 7, section 7.3.3.

% argument check
if nargin < 2 || ...
   ~isa(ffreq, 'double') || ...
   ~any(numel(ffreq) == [1, 2]) || ...
    any(isnan(ffreq)) || ...
   ~isa(forder, 'double') || ...
    numel(forder) ~= 1 || ...
    isinf(forder) || ...
    isnan(forder) || ...
    forder < 1 || ...
    forder > 500
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
ffreq(1) = min(0.5, ffreq(1));
if numel(ffreq) == 1
    ffreq(2) = 1 - sqrt(eps);
else
    ffreq(2) = min(1 - sqrt(eps), max(2 * ffreq(1), ffreq(2)));
end
ffreq(1) = max(sqrt(eps) / (ffreq(2)^2), ffreq(1));
Wn = ffreq(:)';
n = round(forder);

% step 1: get analog, pre-warped frequencies
fs = 2;
u = 2 * fs * tan(pi * Wn / fs);

% step 2: convert to low-pass prototype estimate
Bw = u(2) - u(1);
Wn = sqrt(u(1) * u(2));	% center frequency

% step 3: Get N-th order Butterworth analog lowpass prototype
z = [];
p = exp(1i * (pi * (1:2:n-1) / (2 * n) + pi / 2));
p = [p; conj(p)];
p = p(:);
if rem(n, 2) == 1   % n is odd
    p = [p; -1];
end
k = real(prod(-p));

% Transform to state-space
[a, b, c, d] = zp2ss(z, p, k);

% step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn
a = lp2bp(a, c, Wn, Bw);

% step 5: Use Bilinear transformation to find discrete equivalent:
a = bilinear(a, fs);

% convert to outputs
den = poly(a);
num = buttnum(n, Wn, den);

% lowpass to bandpass analog filter transformation
function at = lp2bp(a,c,wo,bw)
ma = size(c, 2);
q = wo / bw;
at = wo * [a / q eye(ma); -eye(ma) zeros(ma)];

% state-space version of bilinear transformation
function ad = bilinear(z, fp)
t = 1 / fp;
t1 = eye(size(z)) + z * t / 2;
t2 = eye(size(z)) - z * t / 2;
ad = t2 \ t1;

% This internal function returns more exact numerator vectors
% for the num/den case.
% Wn input is two element band edge vector
function b = buttnum(n,Wn,den)
Wn = 2*atan2(Wn,4);
r = [ones(n,1); -ones(n,1)];
w = Wn;
b = poly(r);
% now normalize so |H(w)| == 1:
kern = exp(-1i*w*(0:length(b)-1));
b = real(b*(kern*den(:))/(kern*b(:)));
