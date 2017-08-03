% load surface
if exist('usrf', 'var') == 0 || ...
    numel(usrf) ~= 1 || ...
   ~isxff(usrf, 'srf')
    x = xff;
    srfs = x.Documents('srf');
    if isempty(srfs)
        srf = [];
    else
        for c = 1:numel(srfs)
            srf = x.Document(srfs{c});
            if srf.NrOfVertices >= 40000
                break;
            else
                srf = [];
            end
        end
    end
    if isempty(srf)
        srf = xff([neuroelf_path('colin') '/colin_LH_SPH_ICBMnorm.srf']);
        srf.Browse;
    end
else
    srf = usrf;
end

% generate faces+vertices struct
FV = struct;

% make the coordinates in range -1 .. 1
FV.vertices = (1 / 128) .* (128 - srf.VertexCoordinate);

% and faces
FV.faces = srf.TriangleVertex;%(:, [3, 2, 1]); %(end:-1:1, :);

% get iso normals form the volume data;
FV.normals = patchnormals(FV);

% Make the buffers RGB, and depth
Img = zeros(512, 512, 6); 
Img(:, :, 5) = 1;

% Set the ModelViewMatrix
FV.modelviewmatrix = [0, 1, 0, 0; -1, 0, 0, 0; 0, 0, -1, 0; 0, 0, 0, 1];
  
% Set color to blue
FV.color = srf.ConvexRGBA(1:3);

% enable light
FV.culling = -1;
%FV.depthfunction = 2;
FV.enableshading = 1;
FV.lightposition = [0, 0, -2, 1];
FV.material = [0.3 0.8 0.0 Inf 1.0];
%FV.stencilfunction = 2;

% Render the patch
oImg = renderpatch(Img,FV);
figure;
imshow(oImg(:, :, 1:3));
set(gca, 'YDir', 'normal');
