function [Q, O, D] = spmm2q(M)
% spmm2q  - convert from rotation matrix to quaternion form
%
% FORMAT:       [Q, O, D] = spmm2q(M)
%
% Input fields:
%
%       M           4x4 transformation matrix
%
% Output fields:
%
%       Q           1x3 quaternion B, C, and D values in NIftI header
%       O           1x3 quaternion X, Y, and Z offset in NIftI header
%       D           1x4 voxel sizes and orientation flag in ImgDim.PixSpacing
%
% See also: http://skal.planet-d.net/demo/matrixfaq.htm
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%
% taken from SPM8/@nifti/private/encode_qform0.m and SPM8/@nifti/private/M2Q.m

%
% $Id: M2Q.m 1143 2008-02-07 19:33:33Z spm $

% argument check
if nargin ~= 1 || ...
   ~isa(M, 'double') || ...
    ndims(M) > 2 || ...
    size(M, 1) ~= size(M, 2) || ...
    any(size(M) < 3 | size(M) > 4) || ...
    any(isinf(M(:)) | isnan(M(:))) || ...
   (size(M, 1) == 4 && ...
    any(M(4, :) ~= [0, 0, 0, 1]))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid M matrix supplied.' ...
    );
end

% 1-based to 0-based voxels
M = M * [eye(4, 3), ones(4, 1)];

% offset
O = M(1:3, 4);

% computation
R = M(1:3, 1:3);
vx = sqrt(sum(M(1:3, 1:3) .^ 2));
vx(vx == 0) = 1;
R = R * diag(1 ./ vx);
[U, S, V] = svd(R);
R = U * V';
if det(R) > 0
    D = [1, vx];
else
    R = R * diag([1, 1, -1]);
    D = [-1, vx];
end
d = diag(R);
t = sum(d) + 1;
if t > 0.5
    s = sqrt(t) * 2;
    Q = [(R(3, 2) - R(2, 3)) / s, ...
         (R(1, 3) - R(3, 1)) / s, ...
         (R(2, 1) - R(1, 2)) / s, ...
         0.25 * s]';
else
    t = find(d == max(d));
    t = t(1);

    switch (t)

        case {1}
            s = 2 * sqrt(1 + R(1, 1) - R(2, 2) - R(3, 3));
            Q = [0.25 * s, ...
                 (R(1, 2) + R(2, 1)) / s, ...
                 (R(3, 1) + R(1, 3)) / s, ...
                 (R(3, 2) - R(2, 3)) / s]';

        case {2}
            s = 2 * sqrt(1 + R(2, 2) - R(1, 1) - R(3, 3));
            Q = [(R(1, 2) + R(2, 1)) / s, ...
                 0.25 * s, ...
                 (R(2, 3) + R(3, 2)) / s, ...
                 (R(1, 3) - R(3, 1)) / s]';

        case {3}
            s = 2 * sqrt(1 + R(3, 3) - R(1, 1) - R(2, 2));
            Q = [(R(3, 1) + R(1, 3)) / s, ...
                 (R(2, 3) + R(3, 2)) / s, ...
                 0.25 * s, ...
                 (R(2, 1) - R(1, 2)) / s]';
    end
end

% w must be +ve
if Q(4) < 0
    Q = -Q;
end
Q(4) = [];
