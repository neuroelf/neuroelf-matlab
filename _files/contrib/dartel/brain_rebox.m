function brain_rebox(file)

% neuroelf library
using(neuroelf, {'findfirst'});

% load file
h = xff(file);

% load VoxelData
h.LoadVoxelData;

% get data
vd = h.VoxelData;
vd(isnan(vd)) = 0;

% get all-0 planes
na = (vd ~= 0);
n1 = any(any(na, 2), 3);
n2 = any(any(na, 1), 3);
n3 = any(any(na, 1), 2);

% find first and last index
n11 = max(1, findfirst(n1) - 1);
n12 = min(numel(n1), findfirst(n1, -1) + 1);
n21 = max(1, findfirst(n2) - 1);
n22 = min(numel(n2), findfirst(n2, -1) + 1);
n31 = max(1, findfirst(n3) - 1);
n32 = min(numel(n3), findfirst(n3, -1) + 1);

% shrink data
svd = size(vd);
vd = vd(n11:n12, n21:n22, n31:n32);
if isequal(svd, size(vd))
    h.ClearObject;
    return;
end

% compute new offsets
no = h.CoordinateFrame.Trf * [n11; n21; n31; 1];

% set fields
h.ImgDim.Dim(2:4) = size(vd);
h.DataHist.NIftI1.QuatOffsetX = no(1);
h.DataHist.NIftI1.QuatOffsetY = no(2);
h.DataHist.NIftI1.QuatOffsetZ = no(3);
h.DataHist.NIftI1.AffineTransX(end) = no(1);
h.DataHist.NIftI1.AffineTransY(end) = no(2);
h.DataHist.NIftI1.AffineTransZ(end) = no(3);
h.VoxelData = vd;

% save and clear
h.Save;
h.ClearObject;
