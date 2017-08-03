function [a,e] = lpcfit(x,p,w,h)
% [a,e] = lpcfit(x,p,w,h)  Fit LPC to short-time segments
%    x is a stretch of signal.  Using w point windows every 
%    h points, fit order p LPC models.  Return the successive 
%    all-pole coefficients as rows of a, and the residual excitation 
%    in e.
% 2001-02-25 dpwe@ee.columbia.edu

% neuroelf
n = neuroelf;

if nargin < 2
  p = 12;
end
if nargin < 3
  w = 256;
end
if nargin < 4
  h = w/2;
end

% hanning window
hw = n.filtwin(w, 'hann');

npts = length(x);

nhops = 1 + floor((npts-w)/h);

a = zeros(nhops, p+1);
e = zeros(1, npts);

% Pre-emphasis
pre = [1 -0.9];
x = filter(pre,1,x);

for hop = 1:nhops
  % Extract segment of signal
  xx = x((hop - 1)*h + [1:w]);
  if all(xx == 0)
      continue;
  end
  % Apply hanning window
  wxx = xx .* hw;
  % Form autocorrelation (calculates *way* too many points)
  % rxx = xcorr(wxx,wxx);
  rxx = xcorr(wxx,wxx,p+1);
  % extract just the points we need (middle p+1 points)
  % rxx = rxx(w+[0:p]);
  rxx = rxx(2+p+[0:p]);
  % Setup the normal equations
  R = toeplitz(rxx(1:p));
  % Solve for a (horribly inefficient to use full inv())
  an = inv(R)*rxx(2:(p+1));
  % Calculate residual by filtering xx (after windowing?)
  % rs = filter([1 -an'],1,xx);
  rs = filter([1 -an'],1,wxx);
  G = sqrt(mean(rs.^2));
  % Save them; first column is residual MSE (gain)
  a(hop,:) = [G -an'];
  e((hop - 1)*h + [(1+h/2):(w-h/2)]) = rs([(1+h/2):(w-h/2)])/G;  
end
