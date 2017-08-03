function ia = image_analyze(im)
% image_analyze  - perform some analyses on an image
%
% FORMAT:       ia = image_analyze(im)
%
% Input fields:
%
%       im          Y-by-X-by-color image data
%
% Output fields:
%
%       ia          image analysis

% first, compute HSV equivalent
imhsv = hsvconv(im, 2);

% compute separate histograms
h = imhsv(:, :, 1);
h = [h(:); 1 + h(:); 2 + h(:)];
s = imhsv(:, :, 2);
s = [s(:); s(:); s(:)];
v = imhsv(:, :, 3);
v = [v(:); v(:); v(:)];
hh = histcount(h, 0, 3, 0.01);
hh = flexinterpn(hh(:), [Inf; 101; 1; 201], smoothkern(2.5), 1);
hs = histcount(s, 0, 1, 0.01);
hs = flexinterpn(hs(:), [Inf; 1; 1; 101], smoothkern(2.5), 1);
hv = histcount(v, 0, 1, 0.01);
hv = flexinterpn(hv(:), [Inf; 1; 1; 101], smoothkern(2.5), 1);

% joint histograms
hhs = histcount(h, 0, 3, 0.01, s, 0, 1, 0.01);
hhs = flexinterpn(hhs, [Inf, Inf; 101, 1; 1, 1; 201, 101], smoothkern(2.5), 1);
hhv = histcount(h, 0, 3, 0.01, v, 0, 1, 0.01);
hhv = flexinterpn(hhv, [Inf, Inf; 101, 1; 1, 1; 201, 101], smoothkern(2.5), 1);
hsv = histcount(s, 0, 1, 0.01, v, 0, 1, 0.01);
hsv = flexinterpn(hsv, [Inf, Inf; 1, 1; 1, 1; 101, 101], smoothkern(2.5), 1);
