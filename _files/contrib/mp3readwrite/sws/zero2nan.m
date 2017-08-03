function y = zero2nan(x);


t = find((x==0));
x(t) = NaN;
y = x;
