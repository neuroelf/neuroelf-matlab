function r = rms(x, d)
%RMS  Root of mean of squared values
%   Y = RMS(X) computes the root-mean-squared along the first non-singleton
%   dimension in X.
%
%   Y = RMS(X, DIM) operates along the dimension DIM.
%
%   See also MEAN, STD.

% arguments
if nargin < 1 || ~isnumeric(x)
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
if nargin < 2 || ~isa(d, 'double') || numel(d) ~= 1 || isinf(d) || isnan(d) || d < 1 || d ~= fix(d)
    d = findfirst(size(x) > 1);
    if isempty(d)
        d = 1;
    end
end

% compute
if isreal(x)
    r = sqrt(mean(x .* x, d));
else
    r = sqrt(mean(x .* conj(x), d));
end
