function blob_brain(input_file,input_thresh,colour_file,colour_thresh,blob_file);

%BLOB_BRAIN makes an isosurface 3D plot of blobs.
%
%    BLOB_BRAIN(INPUT_FILE, INPUT_THRESH [, COLOUR_FILE [, COLOUR_THRESH
%         [, BLOB_FILE]]])
%
% Makes 3D blob plot of INPUT_FILE above INPUT_THRESH coloured by values of
% COLOUR_FILE (default is COLOUR_FILE=INPUT_FILE).
%
% COLOUR_THRESH: if it is supplied, then only plots blobs with 
% COLOUR_FILE > COLOUR_THRESH. If COLOUR_THRESH is a vector [a b], a<=b, 
% then values are a < COLOUR_FILE <= b. If COLOUR_THRESH is non-empty
% (default is empty), then to make good plots of coloured clusters, 
% the colours are first 'spread' to the neighbours by taking the 
% max of the colour at each voxel and its 6 3D neigbours. 
%
% BLOB_FILE: If supplied, outputs mask in frame 1 and colours in frame 2.
%
% Coordinates are in voxels, starting at 0, as in 'register'. 
% Since matlab can't reverse axis labels, they are negative if the step 
% sizes are negative - ignore the minus signs.

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

if nargin < 3
   colour_file=input_file;
end
if nargin < 4
   colour_thresh=[];
end
if nargin < 5
   blob_file=[];
end

[base,ext]=fileparts2(deblank(input_file(1,:)));
isminc=strcmp(lower(ext),'.mnc');
d=fmris_read_image(input_file,0,0);

if isempty(colour_thresh)
   d=fmris_read_image(input_file,1:d.dim(3),1);
   mask=d.data>input_thresh;
   xf=find(max(max(mask,[],2),[],3));
   yf=find(max(max(mask,[],1),[],3));
   zf=find(max(max(mask,[],1),[],2));
   xr=max(min(xf)-1,1):min(max(xf)+1,d.dim(1));
   yr=max(min(yf)-1,1):min(max(yf)+1,d.dim(2));
   zr=max(min(zf)-1,1):min(max(zf)+1,d.dim(3));
   mask=d.data(xr,yr,zr);
   d=fmris_read_image(colour_file,1:d.dim(3),1);
   c=d.data(xr,yr,zr);
else
   d=fmris_read_image(colour_file,1:d.dim(3),1);
   colour_thresh1=colour_thresh(1);
   if length(colour_thresh)>=2
      colour_thresh2=colour_thresh(2);
   else
      colour_thresh2=Inf;
   end
   mask=d.data>colour_thresh1 & d.data<=colour_thresh2;
   xf=find(max(max(mask,[],2),[],3));
   yf=find(max(max(mask,[],1),[],3));
   zf=find(max(max(mask,[],1),[],2));
   xr=max(min(xf)-1,1):min(max(xf)+1,d.dim(1));
   yr=max(min(yf)-1,1):min(max(yf)+1,d.dim(2));
   zr=max(min(zf)-1,1):min(max(zf)+1,d.dim(3));
   
   mask=d.data(xr,yr,zr);
   c=mask;
   c1=mask([1 1:(length(xr)-1)],:,:);
   c=c1.*(c1>mask)+c.*(c1<=mask);
   c1=mask([2:length(xr) length(xr)],:,:);
   c=c1.*(c1>mask)+c.*(c1<=mask);
   c1=mask(:,[1 1:(length(yr)-1)],:);
   c=c1.*(c1>mask)+c.*(c1<=mask);
   c1=mask(:,[2:length(yr) length(yr)],:,:);
   c=c1.*(c1>mask)+c.*(c1<=mask);
   c1=mask(:,:,[1 1:(length(zr)-1)]);
   c=c1.*(c1>mask)+c.*(c1<=mask);
   c1=mask(:,:,[2:length(zr) length(zr)]);
   c=c1.*(c1>mask)+c.*(c1<=mask);
   mask=mask>colour_thresh1 & mask<=colour_thresh2;
   
   d=fmris_read_image(input_file,1:d.dim(3),1);
   mask=d.data(xr,yr,zr).*mask;
end

if ~isempty(blob_file)
   d=fmris_read_image(input_file,0,0);
   d.dim(4)=2;
   d.file_name=blob_file;
   d.data=zeros(d.dim(1:3));
   d.data(xr,yr,zr)=mask;
   fmris_write_image(d,1:d.dim(3),1);
   d.data=zeros(d.dim(1:3));
   d.data(xr,yr,zr)=c;
   fmris_write_image(d,1:d.dim(3),2);
end

xr=xr-1;
yr=yr-1;
zr=zr-1;
if d.vox(1)<0 & isminc
   xr=d.dim(1)-1-xr;
end
if d.vox(2)<0 & isminc
   yr=d.dim(2)-1-yr;
end
yr=-yr;
if d.vox(3)<0 
   zr=-zr;
end

[x,y,z]=meshgrid(yr,xr,zr);
p = patch(isosurface(x,y,z,mask,input_thresh,c));
isonormals(x,y,z,mask*sign(-d.vox(3)),p);
set(p,'FaceColor','interp','EdgeColor','none')
view(3); daspect(1./abs(d.vox(1:3))); axis tight;
camlight('right'); camlight('left'); 
lighting gouraud; 
if ~isempty(colour_thresh)
   if colour_thresh1>-Inf & colour_thresh2<Inf
      set(gca,'CLim',[colour_thresh1 colour_thresh2]);
   end
end
colorbar;
ylabel('x');
if d.vox(3)>0; 
   zlabel('z'); 
else 
   zlabel('-z'); 
end
xlabel('-y'); 
return

