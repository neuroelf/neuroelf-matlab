function mask_mesh(input_file, output_file_base, mask_file, mask_thresh, ...
   normalize, write_mask)

%MASK_MESH makes a mesh and mask for input to mesh_tet.m 
%
% MASK_MESH( INPUT_FILE, OUTPUT_FILE_BASE, MASK_FILE [, MASK_THRESH 
%       [, NORMALIZE ]]);
%
% INPUT_FILE is the name of a single 4D image file with multiple  
% frames, or a matrix of image file names, each with a single 3D frame,
% either ANALYZE (.img) or MINC (.mnc) format. Extra blanks are ignored. 
% File separator can be / or \ on Windows. Gzipped files are gunzipped on     
% unix. The first 3 dimesnions are space, the fourth is the coordinates. 
%
% MASK_FILE: A continuous volume to define the mask. If it has multiple
% frames, then the first frame is used. 
%
% MASK_THRESH defines the search volume as the first frame of MASK_FILE 
% > MASK_THRESH. If MASK_THRESH is a vector [a b], a<=b, then mask 
% is a < MASK_FILE <= b. If empty (default), calls fmri_mask_thresh. 
%
% NORMALIZE: If =1, then interpolate the norms of INPUT_FILE as well as
% the vectors, so that if INPUT_FILE is normalized, then so is the mesh.
% Use =1 for wresid files from fmrilm.m or multistat.m. Default is 0.
%
% What it does is to 'move' the cubic mesh so it touches the boundary of 
% the mask and 'fills' the mask. To do this, it finds edges of the 
% tetrahedral mesh that cross the boundary, i.e. one vertex in, one vertex
% out. Then it linearly interpolates the input_file to the point at which
% MASK_FILE = MASK_THRESH, and replaces the values at the outer vertex
% by the interpolated values. The interpolated data is written to 
% OUTPUT_FILE_BASE_mesh.mnc or .img. Thus the mask is enlarged to include
% the outer vertex, which mow 'touches' the boundary. Such an outer vertex 
% may have more than one such 'boundary crossing edge' joined to it - 
% it doesn't matter which edge is used for the interpolated values, so the  
% first one encountered is used. Also writes OUTPUT_FILE_BASE_mask.mnc 
% or .img with 1 for voxels where the interpolated mesh is inside the mask, 
% and 0 outside.

if nargin < 4
    mask_thresh=[];
end
if isempty(mask_thresh)
    mask_thresh=fmri_mask_thresh(mask_file);
end
mask_thresh1=mask_thresh(1)
if length(mask_thresh)>=2
    mask_thresh2=mask_thresh(2)
else
    mask_thresh2=Inf
end

if nargin < 5
   normalize=0
end
if normalize==0
   f=1;
end

if nargin < 6
   write_mask=0
end

numfiles=size(input_file,1)
d=fmris_read_image(input_file,0,0);
d.dim
n=d.dim(4);
numslices=d.dim(3);
J=d.dim(2);
I=d.dim(1);

[base,ext]=fileparts2(input_file(1,:));

dd.file_name=[deblank(output_file_base(1,:)) '_mesh' ext];
dd.dim=[I J numslices n];
dd.parent_file=input_file(1,:);

mm.file_name=[deblank(output_file_base(1,:)) '_mask' ext];
mm.dim=[I J numslices 1];
mm.parent_file=input_file(1,:);

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

% edges:

edges_start1=[ex1(:,1)' ey1(:,1)' exy1(:,1)' ex2(:,1)' ey2(:,1)' exy2(:,1)'];
edges_start2=[ex2(:,1)' ey2(:,1)' exy2(:,1)' ex1(:,1)' ey1(:,1)' exy1(:,1)'];
 
edge1=[ex1(:,1)' ey1(:,1)' exy1(:,1)' ...
       ez1 exz1 eyz1 ...
       ex1(:,2)' ey1(:,2)' exy1(:,2)'];
edge2=[ex2(:,1)' ey2(:,1)' exy2(:,1)' ...
       ez2 exz2 eyz2 ...
       ex2(:,2)' ey2(:,2)' exy2(:,2)'];
    
edges1=[edge1 edge2]; 
edges2=[edge2 edge1]; 

% START:

u=zeros(2*IJ,n);
v=zeros(2*IJ,n);
mask=zeros(2*IJ,1);
nask=zeros(2*IJ,1);
flip=1;
for slice=1:numslices
   slice
   flip=3-flip;
   if numfiles==1
      d=fmris_read_image(input_file,slice,1:n);
      u((1:IJ)+(flip-1)*IJ,:)=reshape(d.data,IJ,n);
      v((1:IJ)+(flip-1)*IJ,:)=reshape(d.data,IJ,n);
   else
      for i=1:n
         d=fmris_read_image(input_file(i,:),slice,1);
         u((1:IJ)+(flip-1)*IJ,i)=reshape(d.data,IJ,1);
         v((1:IJ)+(flip-1)*IJ,i)=reshape(d.data,IJ,1);
      end
   end
   
   % Find voxels that are inside (m1) and outside (m0) the mask boundary:
   m=fmris_read_image(mask_file,slice,1);
   mask((1:IJ)+(flip-1)*IJ)=reshape(m.data,IJ,1);
   nask((1:IJ)+(flip-1)*IJ)= ...
       reshape(m.data>mask_thresh1 & m.data<=mask_thresh2 ,IJ,1);
   if slice==1
      surf=find(~nask(edges_start1) & ...
          mask(edges_start2)>mask_thresh1 & mask(edges_start2)<=mask_thresh2);
      m0=edges_start1(surf);
      m1=edges_start2(surf);
   else
      surf=find(~nask(edges1) & ...
          mask(edges2)>mask_thresh1 & mask(edges2)<=mask_thresh2);
      m0=edges1(surf);
      m1=edges2(surf);
   end
   
   % Linearly interpolate the data (u) and its norm, so that if the data
   % is normalized, then so is the interpolated data:
   if ~isempty(surf)
      dm=mask(m1)-mask(m0);
      mt=mask_thresh1*ones(length(surf),1);
      mt(mask(m0)> mask_thresh2)=mask_thresh2;
      w=(mt-mask(m0))./(dm+(dm==0)).*(dm~=0);
      u0=u(m0,:);
      u1=u(m1,:);
      if normalize==1
         su0=sum(u0.^2,2);
         su1=sum(u1.^2,2);
         su01=sum(u0.*u1,2);
         svt=su0.*(1-w).^2+2*su01.*(1-w).*w+su1.*w.^2;
         f=sqrt((su0.*(1-w).^2+su1.*w.^2)./((1-w).^2+w.^2) ...
            ./(svt+(svt<=0)).*(svt>0));
      end
      for i=1:n
         v(m0,i)=(u0(:,i).*(1-w)+u1(:,i).*w).*f;
      end
      nask(m0)=1;
   end
   
   if slice>1
      dd.data=reshape(v((1:IJ)+(2-flip)*IJ,:),I,J,n);
      fmris_write_image(dd,slice-1,1:n);
      mm.data=reshape(nask((1:IJ)+(2-flip)*IJ),I,J);
      fmris_write_image(mm,slice-1,1);
   end
   if slice==numslices
      dd.data=reshape(v((1:IJ)+(flip-1)*IJ,:),I,J,n);
      fmris_write_image(dd,slice,1:n);
      mm.data=reshape(nask((1:IJ)+(flip-1)*IJ),I,J);
      fmris_write_image(mm,slice,1);
   end
end   

return






