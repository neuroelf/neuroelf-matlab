function [rxx, lags] = xcorr(x, y, maxlags)

if nargin < 2
    y = x;
end
x = x(:);
y = y(:);
if numel(x) < numel(y)
    x(numel(y)) = 0;
elseif numel(y) < numel(x)
    y(numel(x)) = 0;
end
nx = numel(x);
if nargin < 3
    maxlags = nx - 1;
end
lags = -maxlags:maxlags;

% values
rxx = zeros(numel(lags), 1);
for lc = 2:numel(lags)-1
    l = lags(lc);
    if l <= 0
        cp = corrcoef(x(1-l:end), y(1:end+l));
    else
        cp = corrcoef(y(1+l:end), x(1:end-l));
    end
    rxx(lc) = cp(2);
end
