function x = disttrans(x, opts)
%DISTTRANS  Distance transform binary data
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

% special cases
if all(x(:))
    x = (1/0) .* ones(xs);
    return;
elseif ~any(x(:))
    x = zeros(xs);
    return;
end

% take negative mask (for growth limitation)
nx = ~x;
d = zeros(xs);
nd = numel(xs);
ar = struct('type', '()', 'subs', {repmat({':'}, [1, nd])});
am = ar.subs;
an = ar.subs;
for dc = 1:nd
    am{dc} = 1:(xs(dc)-1);
    an{dc} = 2:xs(dc);
end

% dilate negative mask
dnx = nx;
for dc = 1:nd
    ars = ar;
    art = ar;
    ars.subs{dc} = am{dc};
    art.subs{dc} = an{dc};
    dnx = subsasgn(dnx, art, subsref(dnx, art) | subsref(nx, ars));
    ars.subs{dc} = an{dc};
    art.subs{dc} = am{dc};
    dnx = subsasgn(dnx, art, subsref(dnx, art) | subsref(nx, ars));
end
if conn == 'e' ||
elseif conn == 'v' && nd > 3
end

% only border of negative mask!
dnx = dnx & ~nx;


% process until everything is marked up
while ~all(nx(:))
    
end

% resize again
x = reshape(d, xso);
