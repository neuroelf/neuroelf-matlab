function r = lmrainbow(n, r)
%LMRAINBOW Returns a set of luminance matched RGB rainbow colors.
%   C = LMRAINBOW returns the 7x3 RGB color codes (0-255) matched for
%   luminosity.
%
%   C = LMRAINBOW(N) returns a version with interpolated colors added.
%
%   C = LMRAINBOW(N, C) returns a version of interpolated colors C.

% colors
if nargin < 2 || ~isnumeric(r) || size(r, 2) ~= 3 || ndims(r) ~= 2
    r = [255, 0, 0; 224, 112, 8; 192, 192, 64; 16, 160, 16; 16, 160, 160; 0, 0, 255; 64, 0, 112; 128, 0, 224];
else
    r = min(255, max(0, double(r)));
end

% steps given
if nargin > 0 && numel(n) == 1 && isa(n, 'double') && ~isinf(n) && ~isnan(n) && n >= 1
    n = floor(n);
    ri = r;
    rf = (1:(1/n):size(ri, 1))';
    rt = ceil(rf);
    rp = rt - rf;
    rf = rt - 1;
    rf(1) = 1;
    r = rp(:, [1, 1, 1]) .* ri(rf, :) + (1 - rp(:, [1, 1, 1])) .* ri(rt, :);
end

% compute hue, saturation, and luma
[mx, mxp] = max(r, [], 2);
mxp(mx == 0) = 0;
mn = min(r, [], 2);
c = mx - mn;
h = (mxp == 1) .* mod((r(:, 2) - r(:, 3)) ./ c, 6) + ...
    (mxp == 2) .* (2 + (r(:, 3) - r(:, 1)) ./ c) + ...
    (mxp == 3) .* (4 + (r(:, 1) - r(:, 2)) ./ c);
hf = floor(h);
hf(hf == 6) = 0;

% (roughly) match luminance
l = (1 / 255) .* (0.3 * r(:, 1) + 0.59 * r(:, 2) + 0.11 * r(:, 3));
l = max(l) .* ones(size(l));

% recomposite
c = (1 / 255) .* c;
x = c .* (1 - abs(mod(h, 2) - 1));
b = x .* (hf == 2 | hf == 5) + c .* (hf == 3 | hf == 4);
g = x .* (hf == 0 | hf == 3) + c .* (hf == 1 | hf == 2);
r = x .* (hf == 4 | hf == 1) + c .* (hf == 5 | hf == 0);
m = l - (0.3 .* r + 0.59 .* g + 0.11 .* b);
r = min(255, max(0, 255 .* [r + m, g + m, b + m]));

% uint8
r = uint8(round(r));
