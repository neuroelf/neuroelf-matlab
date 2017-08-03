function [f,m] = lpca2frq(a)
% [f,m] = lpca2frq(a)  Convert LPC analysis frames into resonant frequencies
%    Each row of a defines an LPC analysis e.g. from lpcfit.  Convert
%    this into poles, and return the frequencies in rows of f.
% 2001-02-25 dpwe@ee.columbia.edu

[nhops,p] = size(a);

f = zeros(nhops, floor((p - 1)/2));
m = f;

for hop = 1:nhops
  aa = a(hop,:);
  G = aa(1);
  an = aa(2:p);
  rts = roots([1 an]);
  frqs = angle(rts');
  mags = G./(1 - abs(rts'));
  [dummy, ix] = sort(frqs);
  keep = frqs(ix) > 0;
  ix = ix(keep);
  f(hop,1:length(ix)) = frqs(ix);
  m(hop,1:length(ix)) = mags(ix);
end
