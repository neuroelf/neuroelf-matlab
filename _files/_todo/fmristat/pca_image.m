function [V, D]=pca_image(input_file, exclude, p, mask_file, ...
    mask_thresh, output_file_base, iscov, X_remove, X_interest)

%PCA_IMAGE PCA and PLS of a multivariate linear model for image data 
% ignoring the temporal correlation, as described in:
% Worsley, K.J., Poline, J.B., Friston, K.J. and Evans, A.C. (1998). 
% Characterizing the response of PET and fMRI data using 
% Multivariate Linear Models (MLM). NeuroImage, 6:305-319,
%
% The frame components are normalized so that their sum of
% squares is 1, and the spatial components are divided by their sd
% over voxels. Produces nice figures of time and space components.
% 
% [V, D]=PCA_IMAGE(INPUT_FILE [, EXCLUDE [,P [, MASK_FILE [, MASK_THRESH 
%    [, OUTPUT_FILE_BASE [, ISCOV [, X_REMOVE [, X_INTEREST ]]]]]]]])
%
% INPUT_FILE is a matrix of image file names, padded with blanks if  
% necessary, in AFNI (.BRIK), ANALYZE (.img) or MINC format (.mnc).
%
% EXCLUDE is a matrix of frames to be excluded, one row for each 
% row of INPUT_FILE, or a vector of indicators for the excluded
% frames (1=exclude, 0=keep) if the length of EXCLUDE equals the 
% total number of frames in all INPUT_FILEs. Empty by default.
%
% P is number of desired components, not more than numX, the number
% of frames or files. Default is 1.
%
% MASK_FILE is a mask file. If empty, it is ignored. Default is [].
%
% MASK_THRESH defines the search volume as the first frame of MASK_FILE 
% > MASK_THRESH. If MASK_THRESH is a vector [a b], a<=b, then mask 
% is a < MASK_FILE <= b. If empty (default), calls fmri_mask_thresh. 
%
% OUTPUT_FILE_BASE_pca.mnc or .img is output file of the first P 
% spatial components as frames, none if empty (default), in which 
% case a figure of the temporal and spatial components is produced.
%
% ISCOV: 1 for PCA on covariances, 0 for PCA on correlations (default).
%
% X_REMOVE is a design matrix with a row for each frame of 
% the covariates to be removed from the data before doing the PCA. 
% Default is ones(numX,1), i.e. removing the mean image over time or 
% files; use [ones(numX,1) (1:numX)'] to remove a linear drift as well.
%
% X_INTEREST is a design matrix with a row for each frame of 
% the covariates of interest. The PCA is done on the effects for these 
% covariates, i.e. a Parital Least-Squares (PLS) analysis. 
% Default is eye(numX), i.e. a PCA on all the frames or  
% files; use repmat(eye(12),numX/12,1) for a PLS on the average time 
% course of each block of 12 consecutive frames. 
%
% V is numframes x |P| matrix of the first |P| temporal or file components.
%
% D is a diagonal matrix of all the eigen values of the sum of squares matrix. 

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

if nargin < 2
   exclude=[]
end
if nargin < 3
   p=1
end
if nargin < 4
   mask_file=[];
end
if nargin < 5
   mask_thresh=[];
end
if nargin < 6
   output_file_base=[];
end
if nargin < 7
   iscov=[];
end
if isempty(iscov)
   iscov=0
end

numfiles=size(input_file,1)
numframes=zeros(1,numfiles);
for i=1:numfiles
   d=fmris_read_image(input_file(i,:),0,0);
   numframes(i)=d.dim(4);
end
if numfiles<=5
   numframes
end
numslices=d.dim(3)
numys=d.dim(2)
numxs=d.dim(1)
numpix=numxs*numys;
[base,ext]=fileparts2(input_file(1,:));
isminc=strcmp(lower(ext),'.mnc');
isanal=strcmp(lower(ext),'.img');
isafni=strcmp(lower(ext),'.brik');

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
   d=fmris_read_image(mask_file,1:numslices,1);
   mask= (d.data>mask_thresh1 & d.data<=mask_thresh2);
   xf=find(max(max(mask,[],2),[],3));
   xr=min(xf):max(xf);
   yf=find(max(max(mask,[],1),[],3));
   yr=min(yf):max(yf);
   mask=reshape(mask,numpix,numslices);
   N=sum(sum(mask));
else
   xr=1:numxs;
   yr=1:numys;
   N=prod(d.dim(1:3));
end

numX=sum(numframes);
if length(exclude)==numX
   keep=find(exclude==0);
else
   keep=setxor(1:numX,exclude);
end
n=length(keep)

if nargin < 8
   X_remove=ones(numX,1);
end
if nargin < 9
   X_interest=eye(numX);
end

X=X_interest(keep,:);
if ~isempty(X_remove)
   Z=X_remove(keep,:);
   XZ=X-Z*(pinv(Z)*X);
else
   XZ=X;
end
[UX,SX,VX]=svd(XZ);
n1=rank(XZ)
UX=UX(:,1:n1);

A=zeros(n1);
for slice=1:numslices
   slice
   Y=zeros(numpix,n);
   sumin=0;
   sumframes=0;
   for i=1:numfiles
      keepin=intersect(1:numframes(i),keep-sumframes);
      d=fmris_read_image(input_file(i,:),slice,keepin);
      Y(:,sumin+(1:length(keepin)))=reshape(d.data,numpix,length(keepin));
      sumin=sumin+length(keepin);
      sumframes=sumframes+numframes(i);
   end
   if ~isempty(mask_file)
      Y=Y(mask(:,slice),:);
   end
   YX=Y*UX;
   if ~iscov
      if ~isempty(X_remove)
         S=sum((Y-(Y*pinv(Z)')*Z').^2,2);
      else
         S=sum(Y.^2,2);
      end
      Smhalf=(S>0)./sqrt(S+(S<=0));
      for i=1:n1
         YX(:,i)=YX(:,i).*Smhalf;
      end
   end
   A=A+YX'*YX;
end
[Vs, D]=eig(A);
[ds,is]=sort(-diag(D));
ds=-ds;
sd=sqrt(ds(1:p)/N)
pcntvar=ds(1:p)/sum(ds)*100

VX=UX*Vs(:,is(1:p));
out_svd.data=zeros(numxs,numys,numslices,p);
for slice=1:numslices
   slice
   Y=zeros(numpix,n);
   sumin=0;
   sumframes=0;
   for i=1:numfiles
      keepin=intersect(1:numframes(i),keep-sumframes);
      d=fmris_read_image(input_file(i,:),slice,keepin);
      Y(:,sumin+(1:length(keepin)))=reshape(d.data,numpix,length(keepin));
      sumin=sumin+length(keepin);
      sumframes=sumframes+numframes(i);
   end
   if ~isempty(mask_file)
      Y=Y.*repmat(mask(:,slice),1,n);
   end
   u=Y*VX;
   if ~iscov
      if ~isempty(X_remove)
         S=sum((Y-(Y*pinv(Z)')*Z').^2,2);
      else
         S=sum(Y.^2,2);
      end
      Smhalf=(S>0)./sqrt(S+(S<=0));
      for i=1:p
         u(:,i)=u(:,i).*Smhalf;
      end
   end
   out_svd.data(:,:,slice,:)=reshape(u,[numxs numys p]);
end

% change the signs so that max=max(abs) of each image:
umins=squeeze(min(min(min(out_svd.data,[],1),[],2),[],3));
umaxs=squeeze(max(max(max(out_svd.data,[],1),[],2),[],3));
umaxa=max(abs(umins),abs(umaxs));
ok=(umaxs>-umins);
for i=1:p
   if ~ok(i)
      out_svd.data(:,:,:,i)=-out_svd.data(:,:,:,i)/umaxa(i);
      VX(:,i)=-VX(:,i);
   else
      out_svd.data(:,:,:,i)=out_svd.data(:,:,:,i)/umaxa(i);
   end
end

if ~isempty(mask_file)
   mask=reshape(double(mask),[numxs,numys,numslices]);
   mask(mask==0)=NaN;
   for i=1:p
      out_svd.data(:,:,:,i)=out_svd.data(:,:,:,i).*mask;
   end
end

V=zeros(numX,p);
V(keep,1:p)=VX;

LW=1;
subplot(2,1,1);
a=max(abs(VX),[],1);
plot([0 numX],repmat(-(1:p),2,1),'k'); 
hold on; 
plot(keep,VX./repmat(a,n,1)*0.5+repmat(-(1:p),n,1),'LineWidth',LW); 
hold off;
xlim([0 numX]*1.17);
ylabel('-Component'); 
if numfiles==1
   xlabel('Frame');
   title('Temporal components (sd, % variance explained)');
else
   xlabel('File');
   title('File components (sd, % variance explained)');
end
for i=1:p
   text(numX,-i,[num2str(round(sd(i)*100)/100), ', ', ...
         num2str(round(pcntvar(i)*10)/10) '%']);
end

subplot(2,1,2);
if isanal 
   xr=flipdim(xr,2);
end
if isanal %| isafni
   yr=flipdim(yr,2);
end
nrow=round(sqrt(numslices/3.25/p));
nrow=max(nrow,1);
ncol=ceil(numslices/nrow);
numxr=length(xr);
numyr=length(yr);
bigmat=zeros(p*nrow*numyr,ncol*numxr);
r=0;
for k=1:p
   c=0; 
   for i=1:numslices
      bigmat((1:numyr)+r*numyr,(1:numxr)+c*numxr)=flipud(out_svd.data(xr,yr,i,k)');
      c=c+1;
      if c>=ncol & i<numslices
         r=r+1;
         c=0;
      end
   end
   r=r+1;
end
xtick=((1:ncol*numxr)-1)/(ncol*numxr-1)*ncol-0.5;
ytick=((1:p*nrow*numyr)-1)/(p*nrow*numyr-1)*p+0.5;
zlim=[-1 1]; % [min(bigmat(:)) max(bigmat(:))];
imagesc(xtick,ytick,bigmat,zlim);
xlabel('Slice (0 based)'); ylabel('Component');
title('Spatial components');
colormap(spectral); colorbar;

if ~isempty(output_file_base)
   out_svd.parent_file=deblank(input_file(1,:));
   out_svd.dim=[numxs numys numslices p];
   out_svd.file_name=[deblank(output_file_base) '_pca' ext];
   fmris_write_image(out_svd,1:numslices,1:p); 
end

return

   
   
