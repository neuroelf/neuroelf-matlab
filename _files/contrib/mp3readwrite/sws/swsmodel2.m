function [F,M] = swsmodel2(D,R,P,H,O)
% [F,M] = swsmodel(D,R,H)  Sine wave speech analysis
%       D is a speech example sampled at R samples per second.
%       Return a sinusoid model of up to P/2 components with each sinusoid
%       defined by a row of F (frequencies in Hz) and M (linear magnitude). 
%       Each column of F and M corresponds to H samples at R.
%       Rows of F are sorted with lowest frequency first; 
%       sinusoids cannot cross.
%
%       Relies on lpcfit.m and lpca2frq.m to form LPC model (H and O define
%       window size and overlap, respectively) and convert it into frequencies.
% 2001-03-12 dpwe@ee.columbia.edu

if nargin < 2
  R = 8000; 
end
if nargin < 3
  H = 0.02*R;
end

% Target sampling rate
MyR = P*1000;

% Resample to P kHz, so LPC only picks main formants
if R ~= MyR
  % D = resample(D, round(MyR/1000) ,round(R/1000));
  D = interp1((1:numel(D))', D, (1:(R/MyR):numel(D))');
end

% Step size in units of my sampling rate
HH = 2*round(H/R * MyR/2);

% Form Pth-order LPC model (3 or 4 pole pairs)
lpca = lpcfit(D,P,2*HH,HH);

% Convert poles to sorted freqs and magnitudes
% If only 3 nonzero freqs are found, 4th row will have mag/frq zero
[fa, ma] = lpca2frq(lpca);

% Convert frqs into Hz
F = fa'*MyR/(2*pi);

M = ma';

% Check for overflow values of amplitude and adjust to min/max range
i=find(M>1.0);
if ~isempty(i)
	M(i)=1.0;
end

i=find(M<0.0);
if ~isempty(i)
	M(i)=0.0;
end




