function mesh_tet(input_file, output_file_base)

%MESH_TET finds lengths of edges of a tetrahedral mesh
%
% MESH_TET( INPUT_FILE, OUTPUT_FILE_BASE);
%
% INPUT_FILE is the name of a single 4D image file with multiple  
% frames, or a matrix of image file names, each with a single 3D frame,
% either ANALYZE (.img) or MINC (.mnc) format. Extra blanks are ignored. 
% File separator can be / or \ on Windows. Gzipped files are gunzipped on     
% unix. The first 3 dimesnions are space, the fourth is the coordinates. 
%
% Each 3D cube of 8 neighbouring points (voxels) is divided into 5 
% tetrahedra in a checkerboard pattern e.g. 3 x 4 x 2 array:
%
% z=1:                  z=2:                      y=1:
%           y                     y                     z
%     1   2   3   4         1   2   3   4             1   2
%   1 +---+---+---+       1 -------------           1 +---+
%     |  /|\  |  /|         |\  |  /|\  |             |  /|
%     | / | \ | / |         | \ | / | \ |             | / |
%     |/  |  \|/  |         |  \|/  |  \|             |/  |
% x 2 +---+---+---+     x 2 +---+---+---+         x 2 +---+
%     |\  |  /|\  |         |  /|\  |  /|             |\  |
%     | \ | / | \ |         | / | \ | / |             | \ |
%     |  \|/  |  \|         |/  |  \|/  |             |  \|
%   3 +---+---+---+       3 +---+---+---+           3 +---+    
%
%  
% x=1:                  x=2:                  x=3:                  
%           y                     y                     y           
%     1   2   3   4         1   2   3   4         1   2   3   4     
%   1 +---+---+---+       1 +---+---+---+       1 +---+---+---+     
%     |  /|\  |  /|         |\  |  /|\  |         |  /|\  |  /|     
% z   | / | \ | / |     z   | \ | / | \ |     z   | / | \ | / |     
%     |/  |  \|/  |         |  \|/  |  \|         |/  |  \|/  |     
%   2 +---+---+---+       2 +---+---+---+       2 +---+---+---=     
%
% Squared edge lengths using the normalised coordinates are stored in
% a 4D array in OUTPUT_FILE_BASE_tet.img or .mnc with 6 frames for edges
% in directions x, y, xy, z, xz, yz respectively. In cubes of 8 voxels, 
% the squared lengths of the 6 edges closest to x=0, y=0, z=Inf are 
% stored in the voxel of that cube that is closest to x=0, y=0, z=Inf. 

numfiles=size(input_file,1)
d=fmris_read_image(input_file,0,0);
d.dim
n=d.dim(4);
numslices=d.dim(3);
J=d.dim(2);
I=d.dim(1);

[base,ext]=fileparts2(input_file(1,:));
lam.file_name=[deblank(output_file_base(1,:)) '_tet' ext];
lam.dim=[I J numslices 6];
lam.parent_file=input_file(1,:);

% Set up:

i=kron(ones(1,J),1:I);
j=kron(1:J,ones(1,I));

IJ=I*J;
ex=find(i<I);
ex1=[ex; ex+IJ]';
ex2=[find(i>1); find(i>1)+IJ]';

ey=find(j<J);
ey1=[ey; ey+IJ]';
ey2=[find(j>1); find(j>1)+IJ]';

ez=1:(IJ);
ez1=ez;
ez2=ez+IJ;

exye=find((rem(i+j,2)==0)&(i<I)&(j<J));
exyo=find((rem(i+j,2)==1)&(i<I)&(j<J));
exy=[exye exyo];
exy1=[exye     exyo+1; exye+1+IJ exyo+IJ]';
exy2=[exye+1+I exyo+I; exye+I+IJ exyo+1+I+IJ]';

exze=find((rem(i+j,2)==0)&(i<I));
exzo=find((rem(i+j,2)==1)&(i<I));
exz =[exze exzo];
exz1=[exze exzo+1];
exz2=[exze+1+IJ exzo+IJ];

eyze=find((rem(i+j,2)==0)&(j<J));
eyzo=find((rem(i+j,2)==1)&(j<J));
eyz =[eyze eyzo];
eyz1=[eyze eyzo+I];
eyz2=[eyze+I+IJ eyzo+IJ];

% START:

lams=zeros(1,IJ*9);
u=zeros(2*IJ,n);
lam.data=zeros(I,J,6);
flip=1;
for slice=1:numslices
   slice
   flip=3-flip;
   if numfiles==1
      d=fmris_read_image(input_file,slice,1:n);
      u((1:IJ)+(flip-1)*IJ,:)=reshape(d.data,IJ,n);
   else
      for i=1:n
         d=fmris_read_image(input_file(i,:),slice,1);
         u((1:IJ)+(flip-1)*IJ,i)=reshape(d.data,IJ,1);
      end
   end
   start=(flip-1)*6*IJ;
   lams(start+      ex)=sum((u( ex1(:,flip),:)-u( ex2(:,flip),:)).^2,2);
   lams(start+  IJ+ ey)=sum((u( ey1(:,flip),:)-u( ey2(:,flip),:)).^2,2);
   lams(start+2*IJ+exy)=sum((u(exy1(:,flip),:)-u(exy2(:,flip),:)).^2,2);
   lam.data(:,:,1:3)=reshape(lams((flip-1)*6*IJ+(1:(3*IJ))),I,J,3);
   if slice>1
      lams(3*IJ+ ez)=sum((u(ez1 ,:)-u(ez2 ,:)).^2,2);
      lams(4*IJ+exz)=sum((u(exz1,:)-u(exz2,:)).^2,2);
      lams(5*IJ+eyz)=sum((u(eyz1,:)-u(eyz2,:)).^2,2);
      lam.data(:,:,4:6)=reshape(lams(3*IJ+(1:(3*IJ))),I,J,3);
   end
   fmris_write_image(lam,slice,1:6);
end   

return

