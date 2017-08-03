function rxx = sxcorr(x, maxlags)

x = x(:);
nx = numel(x);
if nargin < 2
    maxlags = nx - 2;
end

% values
rxx = zeros(maxlags, 1);
for l = 1:maxlags
    cp = corrcoef(x(1+l:end), x(1:end-l));
    rxx(l) = cp(2);
end
