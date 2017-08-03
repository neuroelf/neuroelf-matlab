function tet_curv(tet_file, output_file_base, dfw)

%TET_CURV Calculates 1D Lipschitz-Killing curvature density.
%
% TET_CURV( TET_FILE, OUTPUT_FILE_BASE [,DFW])
%
% TET_FILE: Tetrehedral mesh file from mesh_tet, see help of mesh_tet.
%
% OUTPUT_FILE_BASE_curv.ext is the 1D L-K curvature density, mm^(-1).
%
% DFW=[DFW1 DFW2]: If TET_FILE is created from multistat.m, then 
% DFW1 and DFW2 are the numerator and denominator df of the wresid's. 
% Set DFW=[DF_RESID DF] from output of multistat.m. This is only used to 
% correct the bias in estimating curvature. It is not needed for wresid 
% files from fmrilm, where DFW1=DFW2=DF, since curvature is unbiased
% if DFW1=DFW2. Default is DFW=[] which skips the bias correction.
%
% The integral of the curvature density^3 is the contribution from the 
% interior of the search region to the 1D Lipschitz-Killing curvature. 

if nargin<3
   dfw=[];
end
biasr=ones(1,4);
if ~isempty(dfw)
   df_resid=dfw(1);
   df=dfw(2);
   alphar=1/2;
   dr=df_resid/df;
   for D=1:3
      dv=df_resid-dr-(0:D-1);
      biasr(D+1)=exp(sum(gammaln(dv/2+alphar)-gammaln(dv/2)) ...
         +gammaln(df/2-D*alphar)-gammaln(df/2))*dr^(-D*alphar);
   end
end

d=fmris_read_image(tet_file,0,0);
d.dim
numslices=d.dim(3);
J=d.dim(2);
I=d.dim(1);

[base,ext]=fileparts2(tet_file(1,:));
curv.file_name=[deblank(output_file_base(1,:)) '_curv' ext];
curv.dim=[I J numslices 1];
curv.parent_file=tet_file(1,:);

f=(1/abs(prod(d.vox))/biasr(2))^(1/3);

% Set up:

IJ=I*J;
i=kron(ones(1,J),1:I);
j=kron(1:J,ones(1,I));
nxy=(((i>1)+(i<I)).*((j>1)+(j<J)))';

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

% tetrahedra:

tet1e=[exye      exye+I+IJ exye+1+IJ exye+IJ   exye+1+I    ];
tet2e=[exye+1    exye+1+I  exye      exye+1+IJ exye+1+IJ   ];
tet3e=[exye+1+I  exye      exye+I+IJ exye+I+IJ exye+1+I+IJ ];
tet4e=[exye+1+IJ exye+I    exye+1+I  exye      exye+I+IJ   ];

tet1o=[exyo+I  exyo+1      exyo+I+IJ   exyo+IJ     exyo+1+I+IJ  ];
tet2o=[exyo    exyo+IJ     exyo+IJ     exyo+I      exyo+1       ];
tet3o=[exyo+1  exyo+1+IJ   exyo+1+I+IJ exyo+1      exyo+1+I     ];
tet4o=[exyo+IJ exyo+1+I+IJ exyo+I      exyo+1+I+IJ exyo+I       ];

% edges of tetrahedra:

e12e=spnode(tet1e+2*IJ*(tet2e-1))';
e13e=spnode(tet1e+2*IJ*(tet3e-1))';
e14e=spnode(tet1e+2*IJ*(tet4e-1))';
e23e=spnode(tet2e+2*IJ*(tet3e-1))';
e24e=spnode(tet2e+2*IJ*(tet4e-1))';
e34e=spnode(tet3e+2*IJ*(tet4e-1))';

e12o=spnode(tet1o+2*IJ*(tet2o-1))';
e13o=spnode(tet1o+2*IJ*(tet3o-1))';
e14o=spnode(tet1o+2*IJ*(tet4o-1))';
e23o=spnode(tet2o+2*IJ*(tet3o-1))';
e24o=spnode(tet2o+2*IJ*(tet4o-1))';
e34o=spnode(tet3o+2*IJ*(tet4o-1))';

e12=[e12e e12o];
e13=[e13e e13o];
e14=[e14e e14o];
e23=[e23e e23o];
e24=[e24e e24o];
e34=[e34e e34o];

ee=1:length(exye)*5;
eo=(length(exye)*5+1):(length(exye)+length(exyo))*5;

store=[1:6*IJ; 6*IJ+1:9*IJ 3*IJ+1:6*IJ]';

ptet=[((j>1)+(j<J)).*(i<I)/4 ...
      ((i>1)+(i<I)).*(j<J)/4 ...
      (i<I).*(j<J)/2 ...
      ((i>1)+(i<I)).*((j>1)+(j<J))/4 ...
      ((j>1)+(j<J)).*(i<I)/2 ...
      ((i>1)+(i<I)).*(j<J)/2 ...
      ((j>1)+(j<J)).*(i<I)/4 ...
      ((i>1)+(i<I)).*(j<J)/4 ...
      (i<I).*(j<J)/2 ];

im=1:I-1;
jm=1:J-1;

% START:

lams=zeros(1,IJ*9);
b=zeros(I,J,2);
flip=1;
for slice=1:numslices
   slice
   flip=3-flip;
   d=fmris_read_image(tet_file,slice,1:6);
   lams(store(:,flip))=reshape(d.data,1,6*IJ);
   
   if slice>1
      
      l12=lams(e12);
      l13=lams(e13);
      l14=lams(e14);
      l23=lams(e23);
      l24=lams(e24);
      l34=lams(e34);
      
      d12=4*l12.*l34-(l13+l24-l23-l14).^2;
      d13=4*l13.*l24-(l12+l34-l23-l14).^2;
      d14=4*l14.*l23-(l12+l34-l24-l13).^2;
      
      a1=4*l24.*l34-(l24+l34-l23).^2;
      a2=4*l14.*l34-(l14+l34-l13).^2;
      a3=4*l14.*l24-(l14+l24-l12).^2;
      a4=4*l13.*l23-(l13+l23-l12).^2;
      
      angle=zeros(1,IJ*9);
      
      h=(a1<=0)|(a2<=0);
      g=acos(-(d12-a1-a2)./sqrt(a1.*a2+h)/2.*(1-h)+h);
      angle(e34e)=angle(e34e)+g(ee);
      angle(e34o)=angle(e34o)+g(eo);
      
      h=(a1<=0)|(a3<=0);
      g=acos(-(d13-a1-a3)./sqrt(a1.*a3+h)/2.*(1-h)+h);
      angle(e24e)=angle(e24e)+g(ee);
      angle(e24o)=angle(e24o)+g(eo);
      
      h=(a1<=0)|(a4<=0);
      g=acos(-(d14-a1-a4)./sqrt(a1.*a4+h)/2.*(1-h)+h);
      angle(e23e)=angle(e23e)+g(ee);
      angle(e23o)=angle(e23o)+g(eo);
      
      h=(a2<=0)|(a3<=0);
      g=acos(-(d14-a2-a3)./sqrt(a2.*a3+h)/2.*(1-h)+h);
      angle(e14e)=angle(e14e)+g(ee);
      angle(e14o)=angle(e14o)+g(eo);
      
      h=(a2<=0)|(a4<=0);
      g=acos(-(d13-a2-a4)./sqrt(a2.*a4+h)/2.*(1-h)+h);
      angle(e13e)=angle(e13e)+g(ee);
      angle(e13o)=angle(e13o)+g(eo);
      
      h=(a3<=0)|(a4<=0);
      g=acos(-(d12-a3-a4)./sqrt(a3.*a4+h)/2.*(1-h)+h);
      angle(e12e)=angle(e12e)+g(ee);
      angle(e12o)=angle(e12o)+g(eo);
      
      a=reshape((ptet-angle/(2*pi)).*sqrt(lams),I,J,9);

      b(:,:,flip)=zeros(I,J,1);
      b=b+convn(a(im, :,[1 7]),ones(2,1,1))/2 ...
         +convn(a( :,jm,[2 8]),ones(1,2,1))/2 ...
         +convn(a(im,jm,[3 9]),ones(2,2,1))/4 ...
         +convn(a( :, :,4),ones(1,1,2))/2 ...
         +convn(a(im, :,5),ones(2,1,2))/4 ...
         +convn(a( :,jm,6),ones(1,2,2))/4;
      
      curv.data=abs(b(:,:,3-flip)).^(1/3).*sign(b(:,:,3-flip))*f;
      fmris_write_image(curv,slice-1,1);
   end
end

curv.data=abs(b(:,:,flip)).^(1/3).*sign(b(:,:,flip))*f;
fmris_write_image(curv,numslices,1);

return   

%End

