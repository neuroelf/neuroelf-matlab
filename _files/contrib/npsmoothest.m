function [rrad, rr] = npsmoothest(r)

% use neuroelf
using(neuroelf, 'all');

% variance thresholds to check
vthresh = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99];

% data check
if nargin < 1 || ...
   ~isnumeric(r) || ...
    ndims(r) ~= 4 || ...
    size(r, 1) < 13 || ...
    size(r, 2) < 13 || ...
    size(r, 3) < 13
    r = randn(66,52,56,150);
    
    % for random data, use set of smoothing kernels
    sv = 1:0.25:5;
    
% for given data
else
    
    % make sure it's double (and resolve transio)
    r = double(r);
    
    % no smoothing
    sv = 0;
end

% size
rsz = size(r);

% number of voxels
nsz = prod(rsz(1:3));

% middle subscript (index) of volume
cxyz = round(0.5 .* rsz(1:3));

% maximum radius
dx = min(12, cxyz(1) - 1);
dy = min(12, cxyz(2) - 1);
dz = min(12, cxyz(3) - 1);

% full box (of grid)
[gx, gy, gz] = ndgrid(-dx:dx, -dy:dy, -dz:dz);
gx = gx(:);
gy = gy(:);
gz = gz(:);

% distance from center (sorted and trimmed)
dxyz = sqrt(gx .* gx + gy .* gy + gz .* gz);
[dxyz, sx] = sort(dxyz);
rx = (dxyz > 13);
dxyz(rx) = [];
sx(rx) = [];

% unique radii and positions
rrad = unique(dxyz);
rrad(1) = [];
urad = 1 + find(diff(dxyz(2:end)));
urad(end+1) = numel(dxyz);

% generate output
if numel(sv) == 1
    rr = NaN .* zeros(numel(rrad), numel(vthresh));
else
    rr = NaN .* zeros(numel(rrad), numel(sv));
end

% generate indices (relative to center)
iidx = sub2ind(rsz(1:3), gx(sx) + cxyz(1), gy(sx) + cxyz(2), gz(sx) + cxyz(3));

% iterate over smoothing values
for sk = sv
    
    % smooth data and reshape
    if sk > 0
        rs = reshape(smoothdata3(r, [sk, sk ,sk]), nsz, rsz(4));
        
    % or simply reshape
    else
        rs = reshape(r, nsz, rsz(4));
    end
    
    % iterate over radii
    
    % access
end
