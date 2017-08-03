function r = fwhm2resel(f)

persistent k;
if isempty(k)
    k = load([neuroelf_path('spm') '/smoothkernmax.mat']);
end

if nargin < 1 || ...
   (~isa(f, 'double') && ...
    ~isa(f, 'single') && ...
    ~istransio(f))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
f = double(f);
r = f;
if isempty(r)
    return;
end
f = f(:);

% which entries to use
nu = isinf(f) | isnan(f) | f <= 0;

% invalid to 0
r(nu) = 0;

% two sub-categories
u1 = ~nu & f <= 10;
u2 = ~nu & f > 10;

% values up to 10 from data
r(u1) = flexinterpn_method(k.ki, 1 + 1000 * f(u1), 'lanczos3');

% values > 10 as an estimation
r(u2) = k.ki(end) ./ (f(u2) / 10);
