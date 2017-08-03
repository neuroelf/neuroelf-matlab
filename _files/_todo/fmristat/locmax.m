function lm = locmax(input_file,input_thresh,mask_file,mask_thresh, ...
    fwhm,flip,cluster_file);

%LOCMAX
%
% LM = LOCMAX( INPUT_FILE, INPUT_THRESH [, MASK_FILE [, MASK_THRESH ...
%              [,FWHM [, FLIP [, CLUSTER_FILE]]]]] )
%
% Finds local maxima and clusters of INPUT_FILE above INPUT_THRESH.
% Clusters are voxels connected in any of the 2*D directions. A local maximum 
% is a voxel which is greater than or equal to all its 2*D neighbours,
% and strictly greater than at least one of them.
%
% MASK_FILE is a mask file. If empty, it is ignored. Default is [].
%
% MASK_THRESH defines the search volume as the first frame of MASK_FILE 
% > MASK_THRESH. If MASK_THRESH is a vector [a b], a<=b, then mask 
% is a < MASK_FILE <= b. If empty (default), calls fmri_mask_threshold. 
%
% FWHM: If a file is provided, cluster volume is measured in resels.
% If FWHM is a scalar, this is taken as the FWHM for the whole volume.
% Default is 1, i.e. cluster volume is measured in mm^D.
%
% FLIP: INPUT_FILE is multiplied by FLIP before processing. For T statistic
% images, FLIP = -1 will look for negative peaks nd clusters. Default is 1.
% 
% CLUSTER_FILE: If nonempty, a volume of cluster indexes is written to 
% CLUSTER_FILE. Default is [].
%
% LM is a matrix with 7 columns. 
% Col 1:    values of local maxima, sorted in descending order.
% Cols 2-4: i,j,k coords of local maxima in voxels, starting at 0, as in 'register'.
% Col 5:    index of cluster to which they belong, in descending order of 
%           cluster resels (1=largest cluster).
% Col 6:    volume of cluster in mm^D.
% Col 7:    resels of cluster.

%############################################################################
% COPYRIGHT:   Copyright 2002 K.J. Worsley, 
%              Department of Mathematics and Statistics,
%              McConnell Brain Imaging Center, 
%              Montreal Neurological Institute,
%              McGill University, Montreal, Quebec, Canada. 
%              worsley@math.mcgill.ca
%
%              Permission to use, copy, modify, and distribute this
%              software and its documentation for any purpose and without
%              fee is hereby granted, provided that this copyright
%              notice appears in all copies. The author and McGill University
%              make no representations about the suitability of this
%              software for any purpose.  It is provided "as is" without
%              express or implied warranty.
%############################################################################

% Defaults:

if nargin<3
   mask_file=[]
end
if nargin<4
   mask_thresh=[]
end
if nargin < 5
   fwhm=1
end
if nargin < 6
   flip=1
end
if nargin < 7
   cluster_file=[]
end
    
d=fmris_read_image(input_file,1,1);
numslices=d.dim(3)
numys=d.dim(2)
numxs=d.dim(1)
numpix=numxs*numys;
D=2+(numslices>1)

i2=2:numxs;
j2=2:numys;
i1=i2-1;
j1=j2-1;
excurset=zeros(numxs,numys);
vox=[];
voxid=[];
edge1=[];
edge2=[];
n=0;

xm1=[1 1:(numxs-1)];
xp1=[2:numxs numxs];
ym1=[1 1:(numys-1)];
yp1=[2:numys numys];
zp1=[2:numslices numslices];
X=d.data*flip;

if ~isempty(mask_file)
   if isempty(mask_thresh)
      mask_thresh=fmri_mask_thresh(mask_file);
   end
   mask_thresh1=mask_thresh(1);
   if length(mask_thresh)>=2
      mask_thresh2=mask_thresh(2);
   else
      mask_thresh2=Inf;
   end
   d=fmris_read_image(mask_file,1,1);
   mask=d.data;
   X(mask<=mask_thresh1 | mask>mask_thresh2)=-Inf;
end

[base,ext]=fileparts2(input_file(1,:));
isminc=strcmp(lower(ext),'.mnc');

X2=X;
imat=(1:numxs)'*ones(1,numys);
if d.vox(1)<0 & isminc
   imat=numxs+1-imat;
end
jmat=ones(numxs,1)*(1:numys);
if d.vox(2)<0 & isminc
   jmat=numys+1-jmat;
end
lm=[];
lmvox=[];
resels=[];

if exist('bwlabeln')==2
   XX=zeros(d.dim(1:3));
end

for slice=1:numslices
   X1=X;
   X=X2;
   d=fmris_read_image(input_file,zp1(slice),1);
   X2=d.data*flip;
   if ~isempty(mask_file)
      d=fmris_read_image(mask_file,zp1(slice),1);
      mask=d.data;
      X2(mask<=mask_thresh1 | mask>mask_thresh2)=-Inf;
   end
   X3=X(xm1,:);
   X4=X(xp1,:);
   X5=X(:,ym1);
   X6=X(:,yp1);
   
   % Excursion set:
   excursetm=excurset;
   excurset=X>input_thresh;
   
   % Local maxima:
   islm=(X>=X1)&(X>=X2)&(X>=X3)&(X>=X4)&(X>=X5)&(X>=X6) ...
      &((X>X1)+(X>X2)+(X>X3)+(X>X4)+(X>X5)+(X>X6)>0)&excurset;
   lmvox=[lmvox; find(islm)+(slice-1)*numpix];
   i=imat(islm);
   j=jmat(islm);
   k=slice*ones(length(i),1);
   Xlm=X(islm);
   lm=[lm; [Xlm i j k]];

   if isstr(fwhm)
      d=fmris_read_image(fwhm,0,0);
      if d.dim(4)>=2
         d=fmris_read_image(fwhm,slice,2);
      else
         d=fmris_read_image(fwhm,slice,1);
         d.data=abs(prod(d.vox(1:3)))./(d.data+(d.data<=0)).^3.*(d.data>0);
      end
      resels=[resels; d.data(excurset)];
   end   
   
   % Clusters:
   voxidm=voxid;
   voxid=reshape(cumsum(excurset(:)),numxs,numys)+n;
   n=n+sum(sum(excurset));
   vox=[vox; find(excurset)+numpix*(slice-1)];
   
   if exist('bwlabeln')~=2
      edge=find([zeros(1,numys); excurset(i1,:).*excurset(i2,:)]);
      edge1=[edge1; voxid(edge-1)];
      edge2=[edge2; voxid(edge)];
      edge=find([zeros(numxs,1)  excurset(:,j1).*excurset(:,j2)]);
      edge1=[edge1; voxid(edge-numxs)];
      edge2=[edge2; voxid(edge)];
      if slice>1
         edge=find(excurset.*excursetm);
         edge1=[edge1; voxidm(edge)];
         edge2=[edge2; voxid(edge)];
      end
   else
      XX(:,:,slice)=excurset;
   end

end

if n<1
   lm=[];
   return
end

if exist('bwlabeln')~=2
   % Find cluster id's in nf (from Numerical Recipes in C, page 346):
   nf=1:n;
   for l=1:length(edge1)
      j=edge1(l);
      k=edge2(l);
      while nf(j)~=j j=nf(j); end
      while nf(k)~=k k=nf(k); end
      if j~=k nf(j)=k; end
   end
   for j=1:n
      while nf(j)~=nf(nf(j)) nf(j)=nf(nf(j)); end
   end   
else
   L=bwlabeln(XX,6);
   nf=L(find(XX))';
   clear L XX;
end

% find the unique cluster id's corresponding to the local maxima:
ivox=find(ismember(vox,lmvox));
clmid=nf(ivox);
[uclmid,iclmid,jclmid]=unique(clmid);

% find their volumes:
ucid=unique(nf);
ucvol=hist(nf,ucid)*prod(abs(d.vox));   
if isstr(fwhm)
    ucrsl=zeros(1,length(ucid));
    for i=1:length(ucid)
        ucrsl(i)=sum(resels(nf==ucid(i)));
    end
else
    ucrsl=ucvol/fwhm^D;
end
f=find(ismember(ucid,uclmid));
uclmvol=ucvol(f);
uclmrsl=ucrsl(f);

% and their ranks (in ascending order):
[sortuclmrsl,iuclmrsl]=sort(uclmrsl);
len=length(iuclmrsl);
rankrsl=zeros(1,len);
rankrsl(iuclmrsl)=1:len;

% add these to lm as extra columns:
lm=[lm rankrsl(jclmid)' uclmvol(jclmid)' uclmrsl(jclmid)'];

% sort by cluster resels, then local maxima, then flip:
lm=flipud(sortrows(lm,[5 1]));
lm(:,5)=len+1-lm(:,5);

% -1 for register:
lm(:,2:4)=lm(:,2:4)-1;

if ~isempty(cluster_file)
   % volume of cluster indexes:
   clus.file_name=cluster_file;
   clus.dim=[numxs numys numslices 1];
   clus.parent_file=input_file;
   ucllab=zeros(1,length(ucid));
   ucllab(f)=len+1-rankrsl;
   cllab=interp1([0 ucid],[0 ucllab],nf);
   for slice=1:numslices
      clus.data=zeros(numxs,numys);
      inslice=find((slice-1)*numpix<vox & vox<=slice*numpix);
      clus.data(vox(inslice)-(slice-1)*numpix)=cllab(inslice);
      fmris_write_image(clus,slice,1);
   end
end

return




