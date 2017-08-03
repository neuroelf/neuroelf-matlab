function [r, cv] = svcvolume(v, fwhm, t)
% svcvolume  - compute the actually covered volume for a SVC

% argument check
if nargin < 2 || ...
   (~isa(v, 'double') && ...
    ~isa(v, 'single') && ...
    ~islogical(v)) || ...
    ndims(v) ~= 3 || ...
   (~isa(fwhm, 'double') && ...
    ~isa(fwhm, 'single')) || ...
    any(isinf(fwhm(:)) | isnan(fwhm(:)) | fwhm(:) <= 0)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% force type
v = double(double(v) > 0);
sv = size(v);
fwhm = double(fwhm);

% single kernel given
if numel(fwhm) == 1

    % smooth data, and threshold
    sk = smoothkern(fwhm, 1e-8, false, 'lanczos8');
    skm = max(sk);
    cv = flexinterpn(v, [Inf, Inf, Inf; ones(2, 3); size(v)], {sk, [0; 1; 0], [0; 1; 0]}, {1, 1, 1}, 0);
    cv = flexinterpn(cv, [Inf, Inf, Inf; ones(2, 3); size(v)], {[0; 1; 0], sk, [0; 1; 0]}, {1, 1, 1}, 0);
    cv = flexinterpn(cv, [Inf, Inf, Inf; ones(2, 3); size(v)], {[0; 1; 0], [0; 1; 0], sk}, {1, 1, 1}, 0);
    cv = (1 / (skm ^ 3)) .* cv;
    cv(cv > 1) = 1;

    % and generate idealized RPV image
    rpv = (skm ^ 3) .* ones(sv);

% 3D kernel given
elseif numel(fwhm) == 3

    % also smooth data and threshold
    sk1 = smoothkern(fwhm(1), 1e-6, false, 'lanczos8');
    sk2 = smoothkern(fwhm(2), 1e-6, false, 'lanczos8');
    sk3 = smoothkern(fwhm(3), 1e-6, false, 'lanczos8');
    skm = max(sk1) * max(sk2) * max(sk3);
    cv = (flexinterpn(v, [Inf, Inf, Inf; ones(2, 3); size(v)], ...
        {sk1, sk2, sk3}, {1, 1, 1}, 0) > (0.5 * skm));

    % and generate idealized RPV image
    rpv = skm .* ones(sv);

% invalid data given
elseif ~isequal(size(fwhm), size(v))
    error( ...
        'neuroelf:BadArgument', ...
        'FWHM image must match volume in size.' ...
    );

% 3D FWHM (or RPV) image given
else
    cv = false(sv);
    vi = find(v(:));

    % depending on type -> FWHM
    if nargin < 3 || ...
        any(fwhm(:) > 1.1) || ...
       (ischar(t) && ...
        ~isempty(t) && ...
        lower(t(1)) == 'f')

        % compute rpv
        rpv = fwhm2resel .^ 3;

    % for RPV image
    else

        % estimate smoothness image first
        rpv = fwhm;
        fwhm = 1 ./ (rpv .^ (1/3));
    end

    % find voxels
    [vf, vs] = sort(fwhm(vi));
    [vx, vy, vz] = ind2sub(sv, vi);

end

% compute resels of covered volume
r = sum(cv(:) .* rpv(:));
