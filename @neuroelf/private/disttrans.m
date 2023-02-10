function x = disttrans(x, opts)
%DISTTRANS  Distance transform binary data (2D/3D)
%
% D = DISTTRANS(M) creates the distance transform of binary mask M
%
% D = DISTTRANS(M, OPTS) takes options into consideration, where OPTS is
%     a 1x1 struct with optional fields and {default} values:
%
%   'conn'    - connectivity, either 'face', 'edge', or {'vertex'}
%               (for 2D data, 'face' and 'edge' produce the same result)
%
% D = DISTTRANS(X [, OPTS]) if X is a numeric array, M will be produced as
%     M = X > 0

% check inputs
if nargin < 1 || (~isnumeric(x) && ~islogical(x))
    error('neuroelf:badArgument', 'Bad or missing argument.')
end
if isempty(x)
    return;
end
if ~islogical(x)
    x = x > 0;
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'conn') || ~ischar(opts.conn) || isempty(opts.conn) || ...
   ~any(lower(opts.conn(1)) == 'efv')
    conn = 'v';
else
    conn = lower(opts.conn(1));
end

% take care of size
xso = size(x);
x = squeeze(x);
xs = size(x);
if numel(xs) > 3
    error('neuroelf:badArgument', 'Function only supports up to 3D input.');
end

% special cases
if all(x(:))
    x = (1/0) .* ones(xso);
    return;
elseif ~any(x(:))
    x = zeros(xso);
    return;
end

% dilate
if numel(xs) < 3
    x = shiftdim(x, numel(xs) - 3);
end
d = zeros(3, 3);
d(2, 2, :) = 1;
d(2, :, 2) = 1;
d(:, 2, 2) = 1;
if conn == 'e'
    d(2, :, :) = 1;
    d(:, 2, :) = 1;
    d(:, :, 2) = 1;
elseif conn == 'v'
    d(:) = 1;
end
xd = dilate3d(x, d) & ~x;
[dx, dy, dz] = ind2sub(size(x), find(xd(:)));
xi = find(x(:));
[x, y, z] = ind2sub(size(x), xi);
xd = xi;
for c = 1:numel(xi)
    xd(c) = min(sqrt((dx - x(c)) .^ 2 + (dy - y(c)) .^ 2 + (dz - z(c)) .^ 2));
end

% set in volume
x = zeros(xso);
x(xi) = xd;
