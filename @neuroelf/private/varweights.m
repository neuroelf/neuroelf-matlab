function vw = varweights(r, d)

if nargin < 1 || ...
   (~isa(r, 'single') && ...
    ~isa(r, 'double')) || ...
    numel(r) < 3
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isa(d, 'double') || ...
    numel(d) ~= 1 || ...
    isinf(d) || ...
    isnan(d)
    d = ndims(r);
    if size(r, d) < 3
        d = findfirst(size(r) > 2, -1);
        if isempty(d)
            error( ...
                'neuroelf:BadArgument', ...
                'Weights require at least N=3.' ...
            );
        end
    end
end

% n
n = size(r, d);

% set inf/nan to 0
if ~isa(r, 'double')
    r = double(r);
end
r(isinf(r) | isnan(r)) = 0;

% compute squared term
r = r .* r;

% compute relative contribution (assuming 0 is missing or bad values)
rc = ((1 / n) .* sum(r > 0, d)) .^ 4;
rc = rc(:);
rcg = (rc > 0);
if all(rcg)
    rcg = ':';
end
rc = rc(rcg);

% create vector for variance weights
vw = zeros(n, 1);

% create accessor
a = repmat({':'}, 1, ndims(r));

% for each N
for c = 1:n
    
    % access data
    a{d} = c;
    m = r(a{:});
    
    % linearize and remove invalid data
    m = m(:);
    m = m(rcg);
    
    % min/max/mean variance (for this N)
    mmm = minmaxmean(m);
    
    % log variances?
    if mmm(2) > (4 * mmm(3))
        m = log(1 + m);
        lv = true;
        mmm = minmaxmean(m);
    else
        lv = false;
    end
    
    % weighted histcount (between min and max)
    hstep = 1e-4 * (mmm(2) - mmm(1));
    hc = histcount(m, mmm(1), mmm(2)+hstep, hstep, 1, rc);
    
    % trim histogram (in case of rounding errors)
    hc(end) = 0;
    
    % and find weighted median as (linear interpolated) at CDF == 0.5
    shc = sum(hc);
    hc = (1 / shc) .* cumsum(hc);
    hx = findfirst(hc > 0.5);
    hx = (hx - 2) + (0.5 - hc(hx-1)) / (hc(hx) - hc(hx-1));
    vw(c) = mmm(1) + hx * hstep;
    if lv
        vw(c) = exp(vw(c)) - 1;
    end
end

% 1 / sqrt(var)
vw = 1 ./ sqrt(vw);
