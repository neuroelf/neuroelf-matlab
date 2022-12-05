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

% cluster positive side
connf = find('fev' == conn);
if ndims(x) < 3
    connf = min(5, connf + 3);
end

% iterate over clusters
d = zeros(xs);
nd = numel(xs);
[cs, cv, cl, cc] = clustercoordsc(x, connf);
ar = struct('type', '()', 'subs', {repmat({':'}, [1, nd])});
for c = 1:numel(cs)
    
    % cutout
    nc = size(cc{c}, 1);
    mnc = min(cc{c});
    mxc = max(cc{c});
    for dc = 1:nd
        ar.subs{dc} = mnc(dc):mxc(dc);
    end
    cx = subsref(cv, ar) == c;
    
    % embed into zeros, create and find boundary pixels
    cxx = false(size(cx) + 2);
    if nd < 3
        cxx(2:end-1, 2:end-1) = cx;
        cxx(1:end-2, 2:end-1) = cxx(1:end-2, 2:end-1) | cx;
        cxx(3:end, 2:end-1) = cxx(3:end, 2:end-1) | cx;
        cxx(2:end-1, 1:end-2) = cxx(2:end-1, 1:end-2) | cx;
        cxx(2:end-1, 3:end) = cxx(2:end-1, 3:end) | cx;
        cxx(2:end-1, 2:end-1) = cxx(2:end-1, 2:end-1) & ~cx;
        [cxx, cxy] = ind2sub(size(cxx), find(cxx(:)));
        cxx = [cxx, cxy] - 1;
    else
        cxx(2:end-1, 2:end-1, 2:end-1) = cx;
        cxx(1:end-2, 2:end-1, 2:end-1) = cxx(1:end-2, 2:end-1, 2:end-1) | cx;
        cxx(3:end, 2:end-1, 2:end-1) = cxx(3:end, 2:end-1, 2:end-1) | cx;
        cxx(2:end-1, 1:end-2, 2:end-1) = cxx(2:end-1, 1:end-2, 2:end-1) | cx;
        cxx(2:end-1, 3:end, 2:end-1) = cxx(2:end-1, 3:end, 2:end-1) | cx;
        cxx(2:end-1, 2:end-1, 1:end-2) = cxx(2:end-1, 2:end-1, 1:end-2) | cx;
        cxx(2:end-1, 2:end-1, 3:end) = cxx(2:end-1, 2:end-1, 3:end) | cx;
        cxx(2:end-1, 2:end-1, 2:end-1) = cxx(2:end-1, 2:end-1, 2:end-1) & ~cx;
        [cxx, cxy, cxz] = find(cxx);
        cxx = [cxx, cxy, cxz] - 1;
    end
    nbc = size(cxx, 1);
    
    % compute distances
    ds = zeros(nc, nbc);
    for dc = 1:nd
        ds = ds + ((cc{c}(:, dc) - (mnc(dc) - 1)) * ones(1, nbc) - ones(nc, 1) * cxx(:, dc)') .^ 2;
    end
    ds = sqrt(min(ds, [], 2));
    
    % set in partial volume
    dx = subsref(d, ar);
    dx(cx) = ds;
    
    % set in output
    d = subsasgn(d, ar, dx);
end

% resize again
x = reshape(d, xso);
