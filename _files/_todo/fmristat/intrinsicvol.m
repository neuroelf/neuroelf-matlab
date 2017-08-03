function invol=intrinsicvol(tet_file, mask_file, mask_thresh)

%INTRINSICVOL intrinsic volumes for a non-isotropic field.
%
% INTRINSICVOL( TET_FILE, MASK_FILE, MASK_THRESH)
%
% TET_FILE: Tetrehedral mesh file created by mesh_tet - 
% see help of mesh_tet.
%
% MASK_FILE: A continuous volume to define the mask. If it has multiple
% frames, then the first frame is used. 
%
% MASK_THRESH: Mask is all voxels where MASK_FILE > MASK_THRESH.
%
% INVOL: Intrinsic volumes for dimensions 0:3. Note that 
%        resels = invol ./ sqrt(4*log(2)).^(0:3)
%
% For fMRI data, get whitened residuals from fmrilm by which_stats(7)=1,
% which will create a file base_wresid.ext (ext = mnc or img). Then:
% 
% mask_mesh(base_wresid.ext, base, mask_file)
% mesh_tet(base_mesh.ext, base)
% intrinsicvol(base_tet.ext, base_mask.ext, 0.5)


d=fmris_read_image(tet_file,0,0);
d.dim
numslices=d.dim(3);
J=d.dim(2);
I=d.dim(1);
IJ=I*J;

% Set up:

i=kron(ones(1,J),1:I);
j=kron(1:J,ones(1,I));

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

edge1=[ex1(:,1)' ey1(:,1)' exy1(:,1)' ...
       ez1 exz1 eyz1 ...
       ex1(:,2)' ey1(:,2)' exy1(:,2)'];
edge2=[ex2(:,1)' ey2(:,1)' exy2(:,1)' ...
       ez2 exz2 eyz2 ...
       ex2(:,2)' ey2(:,2)' exy2(:,2)'];
 
node=[ex ey+IJ exy+2*IJ ...
      ez+3*IJ exz+4*IJ eyz+5*IJ ...
      ex+6*IJ ey+7*IJ exy+8*IJ];
spnode=sparse([edge1 edge2]+2*IJ*([edge2 edge1]-1), ...
   ones(1,2*(9*IJ-5*(I+J)+2)), [node node]);

% triangles:

tri1=[exy1(:,1)' exy1(:,1)' exz1 exz1 eyz1 eyz1 exy1(:)' exy1(:)' exy1(:,2)' exy1(:,2)'];
tri2=[exy2(:,1)' exy2(:,1)' exz2 exz2 eyz2 eyz2 exy2(:)' exy2(:)' exy2(:,2)' exy2(:,2)'];
tri3=[exye+1 exyo exye+I exyo+1+I ...
    exze+1 exzo exze+IJ exzo+1+IJ ...
    eyze+I eyzo eyze+IJ eyzo+I+IJ ...
    exye+1+IJ exyo+IJ exye exyo+1 ...
    exye+I+IJ exyo+1+I+IJ exye+1+I exyo+I ...
    exye+IJ exyo+1+IJ exye+1+I+IJ exyo+I+IJ];

% tetrahedra:

tet1=[exye exyo         exye exyo+1+I+IJ      exye exyo+1+I+IJ ...
      exye exyo+1+I+IJ  exye+1+I exyo+1+I+IJ];
tet2=[exye+1 exyo+1     exye+I exyo+IJ        exye+1+IJ exyo+IJ ...
      exye+1+IJ exyo+IJ exye+1+IJ exyo+1 ];
tet3=[exye+1+I exyo+I   exye+1+I exyo+1+IJ    exye+I+IJ exyo+I+IJ ...
      exye+I+IJ exyo+1  exye+I+IJ exyo+I ];
tet4=[exye+1+IJ exyo+IJ exye+I+IJ exyo+1      exye+1+I exyo+I ...
      exye+IJ exyo+I    exye+1+I+IJ exyo+1+I ];

% START:

lams=zeros(1,IJ*9);
d=fmris_read_image(tet_file,1,1:3);
lams(6*IJ+1:9*IJ)=reshape(d.data,1,3*IJ);
mask=zeros(IJ,2);
m=fmris_read_image(mask_file,1,1);
mask(:,2)=(reshape(m.data,IJ,1)>mask_thresh);
flip=2;
mink=zeros(1,4);

for slice=1:numslices
   slice
   flip=3-flip;
   if slice<numslices
      d=fmris_read_image(tet_file,slice+1,1:3);
      lams((flip-1)*6*IJ+1:((flip-1)*6+3)*IJ)=reshape(d.data,1,3*IJ);
      d=fmris_read_image(tet_file,slice+1,4:6);
      lams(3*IJ+1:6*IJ)=reshape(d.data,1,3*IJ);
      m=fmris_read_image(mask_file,slice+1,1);
      mask(:,flip)=(reshape(m.data,IJ,1)>mask_thresh);
   else
      mask(:,flip)=zeros(IJ,1);
   end
   
   % Intrinsic volume for points:
   
   pc=IJ*(2-flip);
   pf=(1+pc):(IJ+pc);
   
   mink00=sum(mask(pf));
   
   % Intrinsic volume for edges:
   
   ec=(3*IJ-2*I-2*J+1)*(2-flip);
   ef=(1+ec):(6*IJ-3*I-3*J+1+ec);
   edge=find(mask(edge1(ef)) & mask(edge2(ef)))+ec;
   
   mink10=length(edge);
   mink11=sum(psqrt(lams(node(edge))));
   
   % Intrinsic volume for tringles:
   
   tc=2*(IJ-I-J+1)*(2-flip);
   tf=(1+tc):(10*IJ-8*I-8*J+6+tc);
   tri=find(mask(tri1(tf)) & mask(tri2(tf)) & mask(tri3(tf)))+tc;
   l12=lams(spnode(tri1(tri)+2*IJ*(tri2(tri)-1)));
   l13=lams(spnode(tri1(tri)+2*IJ*(tri3(tri)-1)));
   l23=lams(spnode(tri2(tri)+2*IJ*(tri3(tri)-1)));
   
   mink20=length(tri);
   mink21=0.5*sum(psqrt([l12 l13 l23]));
   mink22=sum(psqrt(4*l12.*l13-(l12+l13-l23).^2))/4;
   
   % Intrinsic volume for tetrahedra:
   
   tet=find(mask(tet1) & mask(tet2) & mask(tet3) & mask(tet4));
   l12=lams(spnode(tet1(tet)+2*IJ*(tet2(tet)-1)));
   l13=lams(spnode(tet1(tet)+2*IJ*(tet3(tet)-1)));
   l14=lams(spnode(tet1(tet)+2*IJ*(tet4(tet)-1)));
   l23=lams(spnode(tet2(tet)+2*IJ*(tet3(tet)-1)));
   l24=lams(spnode(tet2(tet)+2*IJ*(tet4(tet)-1)));
   l34=lams(spnode(tet3(tet)+2*IJ*(tet4(tet)-1)));
   
   d12=4*l12.*l34-(l13+l24-l23-l14).^2;
   d13=4*l13.*l24-(l12+l34-l23-l14).^2;
   d14=4*l14.*l23-(l12+l34-l24-l13).^2;
   
   a1=4*l24.*l34-(l24+l34-l23).^2;
   a2=4*l14.*l34-(l14+l34-l13).^2;
   a3=4*l14.*l24-(l14+l24-l12).^2;
   a4=4*l13.*l23-(l13+l23-l12).^2;
   
   h=(a1<=0)|(a2<=0);
   delta12=sum(psqrt(l34).*pacos((d12-a1-a2)./psqrt(a1.*a2+h)/2.*(1-h)+h));
   h=(a1<=0)|(a3<=0);
   delta13=sum(psqrt(l24).*pacos((d13-a1-a3)./psqrt(a1.*a3+h)/2.*(1-h)+h));
   h=(a1<=0)|(a4<=0);
   delta14=sum(psqrt(l23).*pacos((d14-a1-a4)./psqrt(a1.*a4+h)/2.*(1-h)+h));
   h=(a2<=0)|(a3<=0);
   delta23=sum(psqrt(l14).*pacos((d14-a2-a3)./psqrt(a2.*a3+h)/2.*(1-h)+h));
   h=(a2<=0)|(a4<=0);
   delta24=sum(psqrt(l13).*pacos((d13-a2-a4)./psqrt(a2.*a4+h)/2.*(1-h)+h));
   h=(a3<=0)|(a4<=0);
   delta34=sum(psqrt(l12).*pacos((d12-a3-a4)./psqrt(a3.*a4+h)/2.*(1-h)+h));
   
   mink30=length(tet);
   mink31=(delta12+delta13+delta14+delta23+delta24+delta34)/(2*pi);
   mink32=sum(psqrt(a1)+psqrt(a2)+psqrt(a3)+psqrt(a4))/8;
   mink33=sum(psqrt((4*a1.*a2-(a1+a2-d12).^2)./ ...
      (l34+(l34<=0)).*(l34>0)))/48;
   
   % Intrinsic volume for mask:
   
   mink(1)=mink(1)+mink00-mink10+mink20-mink30;
   mink(2)=mink(2)+mink11-mink21+mink31;
   mink(3)=mink(3)+mink22-mink32;
   mink(4)=mink(4)+mink33;
   
end   

invol=mink

return

function y=psqrt(x)
y=sqrt(max(x,0));
return

function y=pacos(x)
y=acos(min(abs(x),1).*sign(x));
return

%End