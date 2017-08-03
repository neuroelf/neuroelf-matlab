% load srf and smp
srf = xff('*.srf');
smp = xff('*.smp');

% get triangle vertices
tv = srf.TriangleVertex;

% sort according to value in map
tvv = reshape(smp.Map(1).SMPData(tv(:)), size(tv));
[tvm, tvmp] = min(tvv, [], 2);
tv(tvmp == 2, :) = tv(tvmp == 2, [2, 3, 1]);
tv(tvmp == 3, :) = tv(tvmp == 3, [3, 1, 2]);
tvv = reshape(smp.Map(1).SMPData(tv(:)), size(tv));

% get vertex coordinates and then coordinates for each triangle vertex
vc = srf.VertexCoordinate;
tc = cat(3, vc(tv(:, 1), :), vc(tv(:, 2), :), vc(tv(:, 3), :));

% for each triangle, subtract the first coordinate (which is then [0,0,0])
tc = tc - repmat(tc(:, :, 1), [1, 1, 3]);
tvv = tvv - repmat(tvv(:, 1), [1, 3]);

% compute lengths
la = sqrt(sum(tc(:, :, 2) .^ 2, 2));
lb = sqrt(sum(tc(:, :, 3) .^ 2, 2));

% re-scale values
tvv(:, 2:3) = tvv(:, 2:3) ./ [la, lb];

% compute absolute max. gradient per triangle
tvv = max(tvv(:, 2:3), [], 2);

% distribute across vertices
map = zeros(size(vc, 1), 1);
for tc = 1:size(tv, 1)
    for vc = 1:3
        map(tv(tc, vc)) = map(tv(tc, vc)) + tvv(tc);
    end
end

% correctly scale (according to number of vertices/neighbors/triangles)
nei = mesh_neighborsarray(srf.Neighbors);
map = map ./ sum(nei ~= 0, 2);

% store in new SMP
smp = smp.CopyObject;
smp.NrOfMaps = 1;
smp.Map = smp.Map(1);
smp.Map.SMPData = single(map);
smap = sort(map);
smp.Map.LowerThreshold = double(smap(ceil(0.95 * numel(smap))));
smp.Map.UpperThreshold = double(smap(end));
