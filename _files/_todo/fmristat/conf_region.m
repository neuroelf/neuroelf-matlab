function [conf_region_file, conf_regions] = conf_region(input_file, ...
    input_thresh,mask_file,mask_thresh,pvalue,flip,is_tstat);

%CONF_REGION finds confidence regions for the location of local maxima
%
% [CONF_REGION_FILE, CONF_REGIONS] = 
%        CONF_REGION( INPUT_FILE [, INPUT_THRESH [, MASK_FILE ...
%              [, MASK_THRESH [, PVALUE [, FLIP [, IS_TSTAT ]]]]]] )
%
% Finds confidence regions for the spatial location of the center of all local
% maxima in INPUT_FILE above INPUT_THRESH. INPUT_FILE must be a T statistic 
% image with high (>40) df or an F statistic with high (>40) denominator df.
% For valid results, INPUT_THRESH should be >= D for T stats, or >= D^2 for  
% F stats in D dimensions (default is INPUT_THRESH=D). Output is a confidence  
% region image (INPUT_FILE minus extension plus _100*(1-PVALUE)cr.mnc or .img) 
% with values equal to the local maximum value in each confidence region. 
% All other voxels are set to 0.
%
% Details: The confidence regions are the set of connected voxels with
% value > sqrt( (local max)^2 - chisq_pvalue), where the probability
% of a chisq random variable with D df exceeding chisq_pvalue is PVALUE
% (chisq_pvalue = 7.81 for D=3, PVALUE = 0.05). For more details, see
% Ma, L., Worsley, K.J. and Evans, A.C. (1999). Variability of spatial 
% location of activation in fMRI and PET CBF images. NeuroImage, 9:S178.
%
% MASK_FILE is a mask file. If empty, it is ignored. Default is [].
%
% MASK_THRESH defines the search volume as the first frame of MASK_FILE 
% > MASK_THRESH. If MASK_THRESH is a vector [a b], a<=b, then mask 
% is a < MASK_FILE <= b. If empty (default), calls fmri_mask_threshold. 
%
% PVALUE: Confidence = (1 - PVALUE) * 100%. Default is 0.05.
%
% FLIP: INPUT_FILE is multiplied by FLIP before processing. For T statistic
% images, FLIP = -1 will look for negative peaks and clusters. Default is 1.
%
% IS_TSTAT=1 if INPUT_FILE is a T-stat image, 0 if an F-stat image. Default is 1.
% 
% CONF_REGION_FILE is the file name of the confidence region image, equal to
% INPUT_FILE minus extension plus _100*(1-PVALUE)cr.mnc or .img. 
%
% CONF_REGIONS is a matrix with 5 columns:
% Col 1:    values of local maxima, sorted in descending order.
% Cols 2-4: x,y,z coords of local maxima in C notation i.e. starting at 0
%           so that they agree with 'register'.
% Col 5:    Volume of confidence region in mm^D.

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

d=fmris_read_image(input_file,0,0);
numslices=d.dim(3)
numys=d.dim(2)
numxs=d.dim(1)
numpix=numxs*numys;
D=2+(numslices>1)

if nargin<2
   input_thresh=D
end
if nargin<3
   mask_file=[]
end
if nargin<4
   mask_thresh=[]
end
if nargin < 5
   pvalue=0.05
end
if nargin < 6
   flip=1
end
if nargin < 7
   is_tstat=1
end
a=3-is_tstat;

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
end

i2=2:numxs;
j2=2:numys;
i1=i2-1;
j1=j2-1;
nbr=[-numpix +numpix -1 +1 -numxs +numxs 0]';

vox=[];
voxvalue=[];
edge1=[];
edge2=[];
X2=repmat(-Inf,numxs,numys);
X=X2;
excurset2=(X2>-Inf);
voxid=zeros(numxs,numys,3);
n=0;

% chisq:
u=1:0.1:20;
df=repmat(D,1,length(u));
pchisq=gammainc(u/2,df/2);
chisq=interp1(pchisq,u,1-pvalue,'pchip');
thresh=sqrt((input_thresh^a-chisq).*(input_thresh^a>chisq))-0.01

for slice=0:numslices
   X1=X;
   X=X2;
   if slice<numslices
      d=fmris_read_image(input_file,slice+1,1);
      X2=d.data*flip;
      if ~isempty(mask_file)
         d=fmris_read_image(mask_file,slice+1,1);
         mask=d.data;
         X2(mask<=mask_thresh1 | mask>mask_thresh2)=-Inf;
      end
      X2(X2<=thresh)=-Inf;
   else
      X2=repmat(-Inf,numxs,numys);
   end
   
   excurset=excurset2;
   excurset2=(X2>-Inf);
   voxid(:,:,1:2)=voxid(:,:,2:3);
   voxid(:,:,3)=reshape(cumsum(excurset2(:)),[numxs,numys,1])+n;
   n=n+sum(sum(excurset2));
   
   if slice>0
      vox=[vox; find(excurset)+numpix*(slice-1)];
      voxvalue=[voxvalue; X(excurset(:))];
      
      X3=[repmat(-Inf,1,numys); X(i1,:)];
      X4=[X(i2,:); repmat(-Inf,1,numys)];
      X5=[repmat(-Inf,numxs,1)  X(:,j1)];
      X6=[X(:,j2)  repmat(-Inf,numxs,1)];
      
      % find max nbr and its id:
      [Xmax, Imax]=max(reshape([X1 X2 X3 X4 X5 X6 X],[numxs numys 7]),[],3);
      edge=find(Imax<7 & excurset);
      edge1=[edge1; voxid(edge+numpix)];
      edge2=[edge2; voxid(edge+numpix+nbr(Imax(edge)))];
   end
end

if n<1
   return
end

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

% make a volume of confidence regions labelled by their maximum value:
d.data=zeros(numxs,numys,numslices);
lm=[];
for i=unique(nf)
   vi=vox(nf==i);
   vv=voxvalue(nf==i);
   [maxvv, maxvox]=max(vv);
   if maxvv>=input_thresh
%      d.data(vi((maxvv^a-vv.^a)<=chisq))=i;
      d.data(vi((maxvv^a-vv.^a)<=chisq))=maxvv;
      lm=[lm; maxvv vi(maxvox) length(vi) i];
   end
end

% write out local maxima and the volumes of their conf regions:

[base,ext]=fileparts2(input_file(1,:));
isminc=strcmp(lower(ext),'.mnc');

lm=flipud(sortrows(lm,1));
[x,y,z]=ind2sub([numxs,numys,numslices],lm(:,2));
x=x-1;
if d.vox(1)<0 & isminc
   x=numxs-1-x;
end
y=y-1;
if d.vox(2)<0 & isminc
   y=numys-1-y;
end
z=z-1;
conf_regions=[lm(:,1) x y z lm(:,3)*prod(abs(d.vox(1:D)))];
num2str(conf_regions)

% write out results:

% [j,i]=sort(lm(:,4));
% d.data=interp1([0; j],[0; i],d.data); 

if strcmp(ext,'.BRIK')
   addon=['_' num2str(100*(1-pvalue)) 'cr'];
   base1=base(1:(length(base)-5));
   base2=base(length(base)-5+(1:5));
   file_w=[base1 addon ext];
   file_h=[base1 addon base2 '.HEAD'];
   file_b=[base1 addon base2 ext];
   if exist(file_h);  delete(file_h); end
   if exist(file_b);  delete(file_b); end
   d.file_name=file_w;
else
   conf_region_file=[base '_' num2str(100*(1-pvalue)) 'cr' ext]
   if exist(conf_region_file,'file'); delete(conf_region_file); end;
   d.file_name=conf_region_file;
end
d.dim=[numxs numys numslices 1];
d.parent_file=input_file;
fmris_write_image(d);

return

