function [F,M] = swsmodel(D,R,H)
% [F,M] = swsmodel(D,R,H)  Sine wave speech analysis
%       D is a speech example sampled at R samples per second.
%       Return a sinusoid model of up to 4 components with each sinusoid
%       defined by a row of F (frequencies in Hz) and M (linear magnitude). 
%       Each column of F and M corresponds to H samples at R.
%       Rows of F are sorted with lowest frequency first; 
%       sinusoids cannot cross.
%
%       Relies on lpcfit.m and lpca2frq.m to form LPC model and convert it 
%       into frequencies.
% 2001-03-12 dpwe@ee.columbia.edu $Header: $

if nargin < 2
  R = 8000; 
end
if nargin < 3
  H = 128;
end

% flip?
if size(D, 2) == length(D)
    D = D';
end

% Target sampling rate
MyR = 8000;

% Resample to 8 kHz, so LPC only picks main formants
if R ~= MyR
  D = resample(D, round(MyR/1000) ,round(R/1000));
end

% Step size in units of my sampling rate
HH = 2*round(H/R * MyR/2);

% Form 8th-order LPC model (3 or 4 pole pairs)
lpca = lpcfit(D,8,2*HH,HH);

% Convert poles to sorted freqs and magnitudes
% If only 3 nonzero freqs are found, 4th row will have mag/frq zero
[fa, ma] = lpca2frq(lpca);

% Convert frqs into Hz
F = fa'*MyR/(2*pi);

M = ma';
