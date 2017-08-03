function y = rampsound(x,sr);

y = x;
rdur = 0.01; %10 msec ramp
nsamp = length(x);
n = fix(sr*rdur);
mvec = (1:n)./n;

y(1:n) = y(1:n) .* mvec;

start2 = (nsamp-n)+1;
mvec2 = ( (nsamp:-1:start2) - (start2-1)) ./n;
y(start2:nsamp) = y(start2:nsamp) .* mvec2;
