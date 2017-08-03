function d = sdist(type, v, df1, df2, pin)
% sdist  - statistics distributions
%
% FORMAT:       d = sdist(type, variate, df1 [, df2 [, pin]])
%
% Input fields:
%
%       type        either of 'betacdf', 'betainv', 'betapdf', 'fcdf', 'finv',
%                   'fpdf', 'normcdf', 'normpdf', 'norminv', 'tcdf', 'tinv', 'tpdf'
%       variate     variate (or p-value)
%       df1, df2    d.f. (or other required parameter)
%       pin         logical flag, inverse from 0.9 to 0.1 logic
%
% Output fields:
%
%       d           determined value of requested distribution function

% Version:  v0.9a
% Build:    13041414
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% argument check
if nargin < 3 || ...
   ~ischar(type) || ...
    size(type, 2) ~= numel(type) || ...
   ~any(strcmp(type, {'betacdf', 'betainv', 'betapdf', 'fcdf', 'finv', ...
        'fpdf', 'normcdf', 'normpdf', 'norminv', 'tcdf', 'tinv', 'tpdf'})) || ...
   (any('bfn' == type(1)) && ...
    nargin < 4) || ...
   ~isnumeric(v) || ...
   ~isnumeric(df1) || ...
   (numel(df1) ~= numel(v) && ...
    numel(df1) ~= 1 && ...
    numel(v) ~= 1) || ...
    isempty(df1) || ...
    any(isinf(df1(:)) | isnan(df1(:))) || ...
   (nargin > 3 && ...
    (~isnumeric(df2) || ...
     isempty(df2) || ...
     any(isinf(df2(:)) | isnan(df2(:))) || ...
     ~isequal(size(df1), size(df2))))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% get output size
osize = size(v);
if numel(v) == 1 && ...
    numel(df1) > 1
    osize = size(df1);
end
onum = prod(osize);

% and make sure arguments are actual numbers
if ~isa(v, 'double')
    v = double(v(:));
else
    v = v(:);
end
if ~isa(df1, 'double')
    df1 = double(df1(:));
else
    df1 = df1(:);
end
if nargin > 3
    if ~isa(df2, 'double')
        df2 = double(df2(:));
    else
        df2 = df2(:);
    end
end

% handle empty input
if isempty(v)
    d = [];
    return;
end

% what to do, start with more common tasks!
switch (type)

    % Student's t-distribution CDF
    case {'tcdf'}
    % Note: this function has been copied from the SPM2 package published
    %       by the Wellcome Department, go here for details:
    %       http://www.fil.ion.ucl.ac.uk/spm/
    %
    % Definition:
    % --------------------------------------------------------------------
    % The CDF F(x) of the Student's t-distribution with v degrees of
    % freedom is the probability that a realisation of a t random variable
    % X has value less than x; F(x)=Pr{X<x} for X~G(h,c). Student's
    % t-distribution is defined for real x and positive integer v (See
    % Evans et al., Ch37).
    %
    % This implementation is not restricted to whole (positive integer) df
    % v, rather it will compute for any df v>0.
    %
    % Variate relationships: (Evans et al., Ch37 & 7)
    % --------------------------------------------------------------------
    % The Student's t distribution with 1 degree of freedom is the Standard
    % Cauchy distribution, which has a simple closed form CDF.
    %
    % Algorithm:
    % --------------------------------------------------------------------
    % The CDF of the Student's t-distribution with v degrees of freedom
    % is related to the incomplete beta function by:
    %       Pr(|X|<x) = betainc(v/(v+x^2),v/2,1/2)
    % so
    %              {     betainc(v/(v+x^2),v/2,1/2) / 2      for x<0
    %       F(x) = |   0.5                                   for x=0
    %              { 1 - betainc(v/(v+x^2),v/2,1/2) / 2      for x>0
    %
    % See Abramowitz & Stegun, 26.5.27 & 26.7.1; Press et al., Sec6.4 for
    % definitions of the incomplete beta function. The relationship is
    % easily verified by substituting for v/(v+x^2) in the integral of the
    % incomplete beta function.
    %
    % MatLab's implementation of the incomplete beta function is used.
    %
    %
    % References:
    % --------------------------------------------------------------------
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

        % check/adapt input size
        if numel(v) ~= numel(df1) && ...
            numel(v) == 1
            v = v .* ones(onum, 1);
        end

        % initialize output
        d = zeros(onum, 1);

        % special case: f is 0.5 when x == 0
        d(v == 0) = 0.5;

        % special case: Standard Cauchy distribution when v == 1
        v1 = (df1 == 1);
        d(v1) = 0.5 + atan(v(v1)) / pi;

        % compute where defined and no special case
        dc = find((~v1) & (v ~= 0));
        if isempty(dc)
            d = reshape(d, osize);
            return;
        end

        % positive
        xpos = dc(v(dc) > 0);
        if ~isempty(xpos)
            if numel(df1) == 1
                d(xpos) = 1 - 0.5 .* betainc(df1 ./ ...
                    (df1 + v(xpos) .^ 2), 0.5 .* df1, 0.5);
            else
                d(xpos) = 1 - 0.5 .* betainc(df1(xpos) ./ ...
                    (df1(xpos) + v(xpos) .^ 2), 0.5 .* df1(xpos), 0.5);
            end
        end

        % negative
        xpos = dc(v(dc) < 0);
        if ~isempty(xpos)
            if numel(df1) == 1
                d(xpos) = 0.5 .* betainc(df1 ./ ...
                    (df1 + v(xpos) .^ 2), 0.5 .* df1, 0.5);
            else
                d(xpos) = 0.5 .* betainc(df1(xpos) ./ ...
                    (df1(xpos) + v(xpos) .^ 2), 0.5 .* df1(xpos), 0.5);
            end
        end
        d = reshape(d, osize);

    % inverse CDF
    case {'tinv'}

        % calculus
        if numel(v) > 1 && ...
            numel(df1) == 1
            d = reshape((sign(v - 0.5) .* ...
                sqrt(df1 ./ sd_betainv(2 * min(v, 1 - v), 0.5 .* df1, ...
                0.5) - df1)), osize);
        else
            d = reshape((sign(v - 0.5) .* ...
                sqrt(df1 ./ sd_betainv(2 * min(v, 1 - v), 0.5 .* df1, ...
                0.5 .* ones(onum, 1)) - df1)), osize);
        end

    % t PDF
    case {'tpdf'}

        % compute and reshape
        d = reshape(((1 + v .^ 2 ./ df1) .^ (-0.5 .* (df1 + 1))) ./ ...
            (sqrt(df1) .* beta(0.5 .* df1, 0.5)), osize);

    % Fisher-Snedecor F distribution CDF
    case {'fcdf'}
    % Note: this function has been copied from the SPM2 package published
    %       by the Wellcome Department, go here for details:
    %       http://www.fil.ion.ucl.ac.uk/spm/
    %
    % Definition:
    % --------------------------------------------------------------------
    % The CDF F(x) of the F distribution with degrees of freedom v & w,
    % defined for positive integer degrees of freedom v & w, is the
    % probability that a realisation of an F random variable X has value
    % less than x F(x)=Pr{X<x} for X~F(v,w). The F-distribution is defined
    % for v>0 & w>0, and for x in [0,Inf) (See Evans et al., Ch16).
    %
    % Variate relationships: (Evans et al., Ch16 & 37)
    % --------------------------------------------------------------------
    % The square of a Student's t variate with w degrees of freedom is
    % distributed as an F-distribution with [1,w] degrees of freedom.
    %
    % For X an F-variate with v,w degrees of freedom, w/(w+v*X^2) has
    % distributed related to a Beta random variable with shape parameters
    % w/2 & v/2.
    %
    % Algorithm:
    % --------------------------------------------------------------------
    % Using the relationship with the Beta distribution: The CDF of the
    % F-distribution with v,w degrees of freedom is related to the
    % incomplete beta function by:
    %       Pr(X<x) = 1 - betainc(w/(w+v*x^2),w/2,v/2)
    % See Abramowitz & Stegun, 26.6.2; Press et al., Sec6.4 for
    % definitions of the incomplete beta function. The relationship is
    % easily verified by substituting for w/(w+v*x^2) in the integral of
    % the incomplete beta function.
    %
    % MatLab's implementation of the incomplete beta function is used.
    %
    %
    % References:
    % --------------------------------------------------------------------
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

        % check/adapt input size
        if numel(v) ~= numel(df1) && ...
            numel(v) == 1
            v = v .* ones(onum, 1);
        end

        % calculate F
        if nargin > 4 && ...
            islogical(pin) && ...
            numel(pin) == 1 && ...
            pin
            d = reshape(betainc( ...
                df2 ./ (df2 + df1 .* v), 0.5 .* df2, 0.5 .* df1), osize);
        else
            d = reshape(1 - betainc( ...
                df2 ./ (df2 + df1 .* v), 0.5 .* df2, 0.5 .* df1), osize);
        end

    % inverse CDF
    case {'finv'}

        % create output
        d = zeros(osize);
        d(v < 0 | v > 1 | isnan(v) | df1 <= 0 | df2 <= 0) = NaN;
        d(v == 1 & df1 > 0 & df2 > 0) = Inf;

        % regular cases
        k = find(v > 0 & v < 1 & df1 > 0 & df1 > 0);
        if ~isempty(k)
            if nargin > 4 && ...
                islogical(pin) && ...
                numel(pin) == 1 && ...
                pin
                if numel(df1) == 1
                    d(k) = ((1 ./ sd_betainv(v(k), ...
                        0.5 .* df2, 0.5 .* df1) - 1) * df2 ./ df1);
                else
                    d(k) = ((1 ./ sd_betainv(v(k), ...
                        0.5 .* df2(k), 0.5 .* df1(k)) - 1) .* df2(k) ./ df1(k));
                end
            else
                if numel(df1) == 1
                    d(k) = ((1 ./ sd_betainv(1 - v(k), ...
                        0.5 .* df2, 0.5 .* df1) - 1) * df2 ./ df1);
                else
                    d(k) = ((1 ./ sd_betainv(1 - v(k), ...
                        0.5 .* df2(k), 0.5 .* df1(k)) - 1) .* df2(k) ./ df1(k));
                end
            end
        end

        % reshape
        d = reshape(d, osize);

    % F PDF
    case {'fpdf'}
    % Note: this function has been copied from the SPM2 package published
    %       by the Wellcome Department, go here for details:
    %       http://www.fil.ion.ucl.ac.uk/spm/
    %
    % Definition:
    % --------------------------------------------------------------------
    % The PDF of the F-distribution with degrees of freedom v & w, defined
    % for positive integer degrees of freedom v>0 & w>0, and for x in
    % [0,Inf) by: (See Evans et al., Ch16)
    %
    %             gamma((v+w)/2)  * (v/w)^(v/2) x^(v/2-1)
    %    f(x) = --------------------------------------------
    %           gamma(v/2)*gamma(w/2) * (1+(v/w)x)^((v+w)/2)
    %
    % Variate relationships: (Evans et al., Ch16 & 37)
    % --------------------------------------------------------------------
    % The square of a Student's t variate with w degrees of freedom is
    % distributed as an F-distribution with [1,w] degrees of freedom.
    %
    % For X an F-variate with v,w degrees of freedom, w/(w+v*X^2) has
    % distributed related to a Beta random variable with shape parameters
    % w/2 & v/2.
    %
    % Algorithm:
    % --------------------------------------------------------------------
    % Direct computation using the beta function for
    %       gamma(v/2)*gamma(w/2) / gamma((v+w)/2)  =  beta(v/2,w/2)
    %
    % References:
    % --------------------------------------------------------------------
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

        % compute and reshape
        d = reshape((df1 ./ df2) .^ (0.5 .* df1) .* v .^ (0.5 .* df1 - 1) ./ ...
            (1 + (df1 ./ df2) .* v) .^ (0.5 .* (df1 + df2)) ./ ...
            beta(0.5 .* df1, 0.5 .* df2), osize);

    % normal distribution inverse
    case {'norminv'}

        % create output
        d = zeros(osize);
        d(v < 0 | v > 1 | isnan(v) | df2 <= 0) = NaN;

        % compute variate for valid positions
        dv = ~isnan(d);
        if numel(df1) == 1
            d(dv) = df2 .* (-sqrt(2) .* erfcinv(2 .* v(dv))) + df1;
        else
            d(dv) = df2(dv) .* (-sqrt(2) .* erfcinv(2 .* v(dv))) + df1(dv);
        end

    % normal distribution
    case {'normpdf'}

        % compute
        d = reshape((1 ./ (df2 .* sqrt(2 * pi))) .* exp(-0.5 * ((v - df1) ./ df2) .^ 2), osize);

    % normal CDF
    case {'normcdf'}

        % pass on to normcdfc
        d = reshape(normcdfc((1 ./ df2) .* (v - df1)), osize);

    % beta CDF
    case {'betacdf'}

        % just pass on
        d = reshape(sd_betacdf(v, df1, df2), osize);

    % inverse CDF
    case {'betainv'}

        % just pass on
        d = reshape(sd_betainv(v, df1, df2), osize);

    % beta PDF
    case {'betapdf'}

        % just pass on
        d = reshape(sd_betapdf(v, df1, df2), osize);
end



% sub-functions



% beta CDF
function bcdf = sd_betacdf(p, a, b)

% build output
bcdf = zeros(size(p));
bcdf(~(a > 0) | ~(b > 0) | isnan(p)) = NaN;
bcdf(p >= 1 & a > 0 & b > 0) = 1;

% regular cases
k = find(p > 0 & p < 1 & a > 0 & b > 0);
if ~isempty(k)
    if numel(a) == 1 && ...
        numel(b) == 1
        bcdf(k) = betainc(p(k), a, b);
    elseif numel(a) == 1
        bcdf(k) = betainc(p(k), a, b(k));
    elseif numel(b) == 1
        bcdf(k) = betainc(p(k), a(k), b);
    else
        bcdf(k) = betainc(p(k), a(k), b(k));
    end
end


% inverse beta CDF
function binv = sd_betainv(p, a, b)

% build output
binv = zeros(size(p));

% special cases
binv(p < 0 | p > 1 | ~(a > 0) | ~(b > 0) | isnan(p)) = NaN;
binv(p == 1 & a > 0 & b > 0) = 1;

% find normal case
k = find(p > 0 & p < 1 & a > 0 & b > 0);
if ~isempty(k)
    if numel(a) > 1
        a = a(k);
    end
    if numel(b) > 1
        b = b(k);
    end
    if numel(a) ~= 1 || ...
        numel(b) ~= 1
    	y = a ./ (a + b);
    else
        y = a ./ (a + b) * ones(numel(k), 1);
    end
    p = p(k);

    % prepare loop
    y(y < eps) = sqrt(eps);
    y(y > (1 - eps)) = 1 - sqrt(eps);

    y_old = y;
    for lc = 1:10000
        h = (sd_betacdf(y_old, a, b) - p) ./ sd_betapdf(y_old, a, b);
        y_new = y_old - h;
        ind = find(y_new <= eps);
        if ~isempty(ind)
            y_new(ind) = y_old(ind) ./ 10;
        end
        ind = find(y_new >= (1 - eps));
        if ~isempty(ind)
            y_new(ind) = 1 - (1 - y_old(ind)) ./ 10;
        end
        h = y_old - y_new;
        h0 = (h == 0);
        if any(h0)
            binv(k(h0)) = y_new(h0);
            if numel(a) > 1
                a(h0) = [];
            end
            if numel(b) > 1
                b(h0) = [];
            end
            h(h0) = [];
            k(h0) = [];
            p(h0) = [];
            y_new(h0) = [];
        end
        if isempty(h) || ...
           (max(abs(h)) < sqrt(eps))
            break;
        end
        y_old = y_new;
    end
    if ~isempty(k)
        binv(k) = y_new;
    end
end


% beta PDF
function bpdf = sd_betapdf(x, a, b)

% build output
bpdf = zeros(size(x));
bpdf(~(a > 0) | ~(b > 0) | isnan(x)) = NaN;

% regular cases
k = find(x > 0 & x < 1 & a > 0 & b > 0);
if ~isempty(k)
    if numel(a) == 1
        bpdf(k) = exp((a - 1) .* log(x(k)) + ...
		    (b - 1) .* log(1 - x(k))) ./ beta(a, b);
    else
        bpdf(k) = exp((a(k) - 1) .* log(x(k)) + ...
		    (b(k) - 1) .* log(1 - x(k))) ./ beta(a(k), b(k));
    end
end
