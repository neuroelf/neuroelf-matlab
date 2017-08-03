function [b, ixx, ixxx] = wregress(X, w, y)

X = diag(w) * X;
ixx = inv(X' * X);
ixxx = ixx * X';
b = ixxx * y;
