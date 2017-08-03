function g = gammapdf(x, h, l)
% gammapdf  - gamma distribution Probability Density Function (PDF)
%
% FORMAT:       g = gammapdf(x, h, l)
%
% Input fields:
%
%       x           points to sample PDF at
%       h           gamma shape parameter (h > 0), scalar 1x1!
%       l           gamma scale parameter (l > 0), scalar 1x1!
%
% Output fields:
%
%       g           PDF of gamma distribution with shape h, scale l
%
% Note: this function has been copied from the SPM2 package published
%       by the Wellcome Department, go here for details:
%       http://www.fil.ion.ucl.ac.uk/spm/

% Definition:
%-----------------------------------------------------------------------
% The PDF of the Gamma distribution with shape parameter h and scale l
% is defined for h>0 & l>0 and for x in [0,Inf) by: (See Evans et al.,
% Ch18, but note that this reference uses the alternative
% parameterisation of the Gamma with scale parameter c=1/l)
%
%           l^h * x^(h-1) exp(-lx)
%    f(x) = ---------------------
%                gamma(h)
%
% Variate relationships: (Evans et al., Ch18 & Ch8)
%-----------------------------------------------------------------------
% For natural (strictly +ve integer) shape h this is an Erlang distribution.
%
% The Standard Gamma distribution has a single parameter, the shape h.
% The scale taken as l=1.
%
% The Chi-squared distribution with v degrees of freedom is equivalent
% to the Gamma distribution with scale parameter 1/2 and shape parameter v/2.
%
% Algorithm:
%-----------------------------------------------------------------------
% Direct computation using logs to avoid roundoff errors.
%
% References:
%-----------------------------------------------------------------------
% Evans M, Hastings N, Peacock B (1993)
%       "Statistical Distributions"
%        2nd Ed. Wiley, New York
%
% Abramowitz M, Stegun IA, (1964)
%       "Handbook of Mathematical Functions"
%        US Government Printing Office
%
% Press WH, Teukolsky SA, Vetterling AT, Flannery BP (1992)
%       "Numerical Recipes in C"
%        Cambridge

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% argument check
if nargin < 3 || ...
   ~isnumeric(x) || ...
   ~isa(h, 'double') || ...
   ~isa(l, 'double') || ...
    numel(h) ~= 1 || ...
    numel(l) ~= 1 || ...
    any(isinf([h, l]) | isnan([h, l])) || ...
    h <= 0 || ...
    l <= 0
    error( ...
        'neuroelf:BadArgument', ...
        'Insufficient or bad arguments.' ...
    )
end
if ~isa(x, 'double')
    x = double(x(:));
else
    x = x(:);
end

% initialise result to zeros
g = zsz(x);
if isempty(x)
    return;
end

% special cases for x == 0
if h < 1
    g(x == 0) = Inf;
elseif h == 1
    g(x == 0) = l;
end

% compute where x > 0
q = find(x > 0);
if isempty(q)
    return;
end
g(q) = exp((h - 1) .* log(x(q)) + h .* log(l) - l .* x(q) - gammaln(h));
