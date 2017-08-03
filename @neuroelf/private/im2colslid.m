function B = im2colslid(A, w)

m=w(1);
n=w(2);

indices = (1:m)' * ones(1, n) + ones(m, 1) * (0:size(A,1):(size(A,1) * n - 1));
access = (0:(size(A, 1)-m))' * ones(1, size(A, 2) + 1 - n) + ...
    ones(size(A, 1) + 1 -m, 1) * (0:size(A, 1):(size(A, 1) * (size(A, 2) - n)));
indices = reshape(indices, numel(indices), 1) * ones(1, numel(access)) + ...
    ones(numel(indices), 1) * reshape(access, 1, numel(access));
B = reshape(A(reshape(indices, numel(indices), 1)), size(indices));
