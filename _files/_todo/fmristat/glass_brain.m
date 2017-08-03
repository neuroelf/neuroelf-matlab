function [search_volume, num_voxels]= ...
   glass_brain(input_file,input_thresh,mask_file,mask_thresh,flip);

%GLASS_BRAIN makes an SPM stlye maximum intensity projection.
%
% [SEARCH_VOLUME, NUM_VOXELS] = 
% GLASS_BRAIN( INPUT_FILE, INPUT_THRESH [, MASK_FILE, MASK_THRESH [, FLIP]])
%
% Makes an SPM stlye maximum projection of INPUT_FILE above INPUT_THRESH.
% Coordinates are in voxels, starting at 0, as in 'register'. Since matlab 
% can't reverse axis labels, they are negative if the step sizes are negative 
% - ignore the minus signs.
%
% MASK_FILE is a mask file. If empty, it is ignored. Default is [].
%
% MASK_THRESH defines the search volume as the first frame of MASK_FILE 
% > MASK_THRESH. If MASK_THRESH is a vector [a b], a<=b, then mask 
% is a < MASK_FILE <= b. If empty (default), calls fmri_mask_thresh. 
%
% FLIP: INPUT_FILE is multiplied by FLIP before processing. For T statistic
% images, FLIP = -1 will look for negative peaks and clusters. Default is 1.
%
% SEARCH_VOLUME is the volume of the masked region in mm^D. 
%
% NUM_VOXELS is the number of voxels (3D) or pixels (2D) in SEARCH_VOLUME.

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
if nargin < 4
   mask_thresh=[]
end
if nargin < 5
   flip=1
end

d=fmris_read_image(input_file,0,0);
numslices=d.dim(3);

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
   mask=d.data;
   mask= (mask>mask_thresh1 & mask<=mask_thresh2);
   num_voxels=sum(sum(sum(mask)));
   
   xf=find(max(max(mask,[],2),[],3));
   xr=min(xf):max(xf);
   yf=find(max(max(mask,[],1),[],3));
   yr=min(yf):max(yf);
   zf=find(max(max(mask,[],1),[],2));
   zr=min(zf):max(zf);
   mask=mask(xr,yr,zr);
   
   d=fmris_read_image(input_file);
   X=d.data(xr,yr,zr)*flip;
   maxX=max(X(:));
   X(X<=input_thresh)=input_thresh-(maxX-input_thresh)/10;
   X(~mask)=input_thresh-(maxX-input_thresh)/5;
else
   d=fmris_read_image(input_file);
   X=d.data*flip;
   X(X<=input_thresh)=input_thresh-1;
   xr=1:d.dim(1);
   yr=1:d.dim(2);
   zr=1:d.dim(3);
   num_voxels=prod(d.dim(1:3));
end

voxel_volume=abs(prod(d.vox(1:3)));
search_volume=num_voxels*voxel_volume
num_voxels

sx=abs(d.vox(1))*length(xr);
sy=abs(d.vox(2))*length(yr);
sz=abs(d.vox(3))*length(zr);

[path,name,ext]=fileparts(deblank(input_file(1,:)));
if strcmp(ext,'.gz')
   [path,name,ext]=fileparts([path name]);
end   
isminc=strcmp(lower(ext),'.mnc');
isanal=strcmp(lower(ext),'.img');

xp=xr-1; signx=' '; if  isminc & d.vox(1)<0;  xp=-(d.dim(1)-xp-1); signx='-'; end
yp=yr-1; signy=' '; if  isminc & d.vox(2)<0;  yp=-(d.dim(2)-yp-1); signy='-'; end
zp=zr-1;
if isanal; xp=-xp; signx='-'; yp=-yp;  signy='-';  end 

h=0.75/max(sx+sz,(sy+sx)*0.75);
clf;
axes('position',[.1 .15+sx*h sy*h*0.75 sz*h])
imagesc(yp,zp,squeeze(max(X,[],1))'); grid on; 
if d.vox(3)>0 
   axis xy; 
   ylabel(['z=' num2str(min(zp)) ':' num2str(max(zp))]);
else
   ylabel(['z=' num2str(max(zp)) ':' num2str(min(zp))]);
end 
title([strrep(input_file,'_',' ') ' > ' num2str(input_thresh)]);

axes('position',[.15+sy*h*0.75 .15+sx*h sx*h*0.75 sz*h])
imagesc(xp,zp,squeeze(max(X,[],2))'); grid on; 
if d.vox(3)>0 
    axis xy; 
end 
xlabel([signx 'x=' num2str(min(xp)) ':' num2str(max(xp))]);
if flip==1
    title('Positive values');
else
    title('Negative values');
end

axes('position',[.1 .1 sy*h*7/8 sx*h])
imagesc(yp,xp,squeeze(max(X,[],3))); grid on; colorbar;
ylabel([signx 'x=' num2str(min(xp)) ':' num2str(max(xp))]);
xlabel([signy 'y=' num2str(min(yp)) ':' num2str(max(yp))]);
colormap(spectral);
   
return







