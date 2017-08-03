function s = spectral(m)
%SPECTRAL Black-purple-blue-green-yellow-red-white color map.
%
%         map = spectral(num_colors)
%
% SPECTRAL(M) returns an M-by-3 matrix containing a "spectral" colormap.
% SPECTRAL, by itself, is the same length as the current colormap.
%
% For example, to reset the colormap of the current figure:
%
%           colormap(spectral)
%
% See also HSV, GRAY, PINK, HOT, COOL, BONE, COPPER, FLAG,
%          COLORMAP, RGBPLOT.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
base = [
  0.0000 0.0000 0.0000
  0.4667 0.0000 0.5333
  0.5333 0.0000 0.6000
  0.0000 0.0000 0.6667
  0.0000 0.0000 0.8667
  0.0000 0.4667 0.8667
  0.0000 0.6000 0.8667
  0.0000 0.6667 0.6667
  0.0000 0.6667 0.5333
  0.0000 0.6000 0.0000
  0.0000 0.7333 0.0000
  0.0000 0.8667 0.0000
  0.0000 1.0000 0.0000
  0.7333 1.0000 0.0000
  0.9333 0.9333 0.0000
  1.0000 0.8000 0.0000
  1.0000 0.6000 0.0000
  1.0000 0.0000 0.0000
  0.8667 0.0000 0.0000
  0.8000 0.0000 0.0000
  0.8000 0.8000 0.8000
];
n = length(base);
X0 = linspace (1, n, m);
s = interp1(1:n,base,X0);

return
