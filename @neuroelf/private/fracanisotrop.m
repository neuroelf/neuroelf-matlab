function fa = fracanisotrop(ev)
% fracanisotrop  - calculating the Fractional Anisotropy over voxels
%
% FORMAT:       fa = fracanisotrop(ev)
%
% Input fields:
%
%       ev          3xV (or Vx3) eigenvalues
%
% Output fields:
%
%       fa          Vx1 fractional anisotropy

% argument check
if nargin < 1 || ...
   (~isa(ev, 'double') && ...
    ~isa(ev, 'single')) || ...
    ndims(ev) > 2 || ...
   ~any(size(ev) == 3)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% make sure computation is valid
if size(ev, 1) == 3
    mask = all(~isnan(ev), 1);
    l1 = ev(1, mask);
    l2 = ev(2, mask);
    l3 = ev(3, mask);
else
    mask = all(~isnan(ev), 2);
    l1 = ev(mask, 1);
    l2 = ev(mask, 2);
    l3 = ev(mask, 3);
end

% computation
d1 = l1 - l2;
d1 = d1 .* d1;
d2 = l2 - l3;
d2 = d2 .* d2;
d1 = d1 + d2;
d2 = l3 - l1;
d2 = d2 .* d2;
d1 = d1 + d2;
d1 = sqrt(d1);
l1 = l1 .* l1;
l2 = l2 .* l2;
l3 = l3 .* l3;
l1 = l1 + l2;
l1 = sqrt(2 .* (l1 + l3));
fa = NaN(numel(mask), 1);
fa(mask) = d1(:) ./ l1(:);
