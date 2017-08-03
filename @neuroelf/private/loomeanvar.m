function [m, v] = loomeanvar(y, d)
%LOOMEANVAR  Compute leave-one-out sample mean and variance.
%   [M, V] = LOOMEANVAR(Y) computes the mean M and variance V of the data
%   in sample Y using the jack-knife samples, such that for each data
%   point in Y one mean and variance are computed. Both M and V have the
%   same size as Y, but are forced to double precision. By default, the
%   mean and variance are computed along the first non-singleton dimension.
%
%   [M, V] = LOOMEANVAR(Y, D) computes the mean and variance along
%   dimension D.
%
%   Note: if Y contains Inf or NaN values, the returned mean M and V are
%   estimated by assuming those values are at the M.

% argument check
if nargin < 1 || ~isnumeric(y) || numel(y) < 2
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end

% ensure data is double
y = double(y);

% no dimension given
if nargin < 2 || ~isa(d, 'double') || numel(d) ~= 1 || isinf(d) || isnan(d) || d < 1 || d > ndims(y)
    d = findfirst(size(y) > 1);
else
    d = fix(real(d));
end
ne = size(y, d);
re = ones(1, ndims(d));
re(d) = ne;

% get elements that are valid
ve = ~(isinf(y) | isnan(y));
ge = sum(ve, d);

% replace bad elements
y(~ve) = 0;

% compute mean and sums-of-squares 
m = repmat(sum(y, d) ./ ge, re);
v = repmat(sum(y .* y, d), re);

% correct mean if all elements are valid
if all(ge(:) == ne)
    m = (ne / (ne - 1)) .* (m - (1 / ne) .* y);
    v = (1 / (ne - 2)) .* ((v - y .* y) - (ne - 1) .* m .* m);

% if not all elements are valid
else

    % first, create factors and replicate
    mge = repmat(1 ./ ge, re);
    rge = repmat(ge ./ (ge - 1), re);
    vge = repmat(ge, re) - double(ve);

    % then compute
    mr = rge .* (m - mge .* y);

    % but only replace where valid (invalid samples don't affect the mean)
    m(ve) = mr(ve);
    v = (1 ./ (vge - 1)) .* ((v - y .* y) - vge .* m .* m);
end
