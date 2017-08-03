function cv = simcluster(d, do, conn)

% argument check
if nargin < 2 || ...
   ~isnumeric(d) || ...
    isempty(d) || ...
    ndims(d) ~= 4 || ...
   ~isa(do, 'double') || ...
    size(do, 1) ~= size(d, 1) || ...
    size(do, 2) ~= size(d, 2) || ...
    size(do, 3) ~= size(d, 3) || ...
    any(isinf(do(:)) | isnan(do(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
   ~isa(conn, 'double') || ...
    numel(conn) ~= 1 || ...
    isinf(conn) || ...
    isnan(conn) || ...
   ~any(conn == [1, 2, 3])
    conn = 2;
end

% output
cv = zeros(size(do));

% search order
[sov, so] = sort(abs(do(:)), 'descend');

% remove empty values
so(sov == 0) = [];
sov(sov == 0) = [];

% volume size
vsz = size(d);
vsz(4) = [];

% reserve space for partial correlations in directions as needed
switch (conn)

    % only major directions
    case {1}
        cdiag = zeros([vsz, 6]);

    % major directions and 2-way diagonals
    case {2}
        cdiag = zeros([vsz, 18]);

    % all diagonals
    case {3}
        cdiag = zeros([vsz, 26]);
end

% compute main elements
[ccv, cdiag(1:end-1, :, :, 1)] = cov_nd(d(1:end-1, :, :, :), d(2:end, :, :, :));
cdiag(2:end, :, :, 2) = cdiag(1:end-1, :, :, 1);
[ccv, cdiag(:, 1:end-1, :, 3)] = cov_nd(d(:, 1:end-1, :, :), d(:, 2:end, :, :));
cdiag(:, 2:end, :, 4) = cdiag(:, 1:end-1, :, 3);
[ccv, cdiag(:, 1:end-1, :, 3)] = cov_nd(d(:, 1:end-1, :, :), d(:, 2:end, :, :));
cdiag(:, :, 2:end, 6) = cdiag(:, :, 1:end-1, 5);

% compute main diagonals
if conn > 1

    % compute all diagonals
    if conn > 2
    end
end

% remove bad entries
cdiag(isnan(cdiag) | isinf(cdiag)) = 0;

% remove values outside mask
cdiag(repmat(do == 0, [1, 1, 1, size(cdiag, 4)])) = 0;

% reshape
cdiag = reshape(cdiag, prod(vsz), size(cdiag, 4));

% mark voxels
cv(so) = 1:numel(so);

% search
for sc = 1:numel(so)

    % find maximum correlation
    mcpos = maxpos(cdiag(so(sc), :));

    %
end
