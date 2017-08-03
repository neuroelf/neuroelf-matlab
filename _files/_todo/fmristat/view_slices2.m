function view_slices(input_file, mask_file, mask_thresh, slices, frames, zlim, layout, expr)

%VIEW_SLICES displays the slices in a 3D volume.
%
% VIEW_SLICES( INPUT_FILE [, MASK_FILE [, MASK_THRESH [, SLICES 
%           [, FRAMES [, ZLIM [, LAYOUT [, EXPR ]]]]]]])
%
% INPUT_FILE is the name of a single 4D image file with one or more  
% volumes, or a matrix of image file names, each with a single 3D volume.
%
% MASK_FILE is a mask file. If empty, it is ignored. Default is [].
%
% MASK_THRESH defines the search volume as the first frame of MASK_FILE 
% > MASK_THRESH. If MASK_THRESH is a vector [a b], a<=b, then mask 
% is a < MASK_FILE <= b. If empty (default), calls fmri_mask_thresh. 
%
% SLICES is a vector of desired slices, starting at 0, as in 'register',
% displayed at top left. Default is 0:(numslices-1), i.e. all the slices. 
%
% FRAMES is a vector of frames, displayed in the top right. Default is 1.
% If INPUT_FILE is a matrix, then FRAMES refers to the rows.
%
% If SLICES and FRAMES is a scalar, then axis values are given in voxels,
% starting at 0, as in 'register'. Since matlab can't reverse axis labels, 
% they are negative if the step sizes are negative - ignore the minus signs.
%
% ZLIM is the limits [min max] for the image, ignored if empty (default).
%
% LAYOUT is a matrix whose entries are the image numbers (0 is empty).
%
% EXPR: matlab expression applied to input in d.data, 
% e.g. 'd.data=sqrt(d.data);'. 

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
   mask_file=[]
end
if nargin < 3
   mask_thresh=[]
end
if nargin < 5
   frames=1
end
if nargin < 6
   zlim=[]
end
if nargin < 7
   layout=[];
end
if nargin < 8
   expr=[];
end

layout=fliplr(layout');

numfiles=size(input_file,1);
d=fmris_read_image(deblank(input_file(1,:)),0,0);
[base,ext]=fileparts2(input_file(1,:));
isminc=strcmp(lower(ext),'.mnc');
isanal=strcmp(lower(ext),'.img');
isafni=strcmp(lower(ext),'.brik');
isnifti=strcmp(lower(ext),'.nii');

if nargin < 4
   slices=(1:d.dim(3))-1
end

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
   m=fmris_read_image(mask_file,slices+1,1);
   m.data= (m.data>mask_thresh1 & m.data<=mask_thresh2);
   xf=find(max(max(m.data,[],2),[],3));
   xr=min(xf):max(xf);
   yf=find(max(max(m.data,[],1),[],3));
   yr=min(yf):max(yf);
else
   xr=1:d.dim(1);
   yr=1:d.dim(2);
end
xp=xr-1; signx=' '; if isminc & d.vox(1)<0; xp=-(d.dim(1)-xp-1); signx='-'; end
yp=yr-1; signy=' '; if isminc & d.vox(2)<0; yp=-(d.dim(2)-yp-1); signy='-'; end
if isanal %| isafni
   xp=fliplr(-xp);
   signx='-';
end

mask_value=NaN;
nx=length(xr);
ny=length(yr);
nslices=length(slices)*length(frames);
Slices=repmat(slices,1,length(frames));
Frames=repmat(frames',1,length(slices))';
Frames=Frames(:)';
if nslices>1
   if isempty(layout)
      mx=round(sqrt(nx*ny*nslices)/nx);
      my=ceil(nslices/mx);
      layout=zeros(mx,my);
      layout(1:nslices)=1:nslices;
      layout=fliplr(layout);
   else
      mx=size(layout,1);
      my=size(layout,2);
   end
   M=zeros(mx*nx,my*ny).*NaN;
   for i=1:mx
      for j=1:my
         islice=layout(i,j);
         if islice>0
            slice=Slices(islice);
            frame=Frames(islice);
            if numfiles==1
               d=fmris_read_image(input_file,slice+1,frame);
            else
               d=fmris_read_image(deblank(input_file(frame,:)),slice+1,1);
            end
            if ~isempty(expr)
               eval(expr);
            end
            v=d.data(xr,yr);
            if ~isempty(mask_file)
               v(m.data(xr,yr,find(slices==slice))==0)=mask_value;
            end
            if isanal
               v=flipdim(v,2);
            end
            M((i-1)*nx+(1:nx),(j-1)*ny+(1:ny))=v;
         end
      end
   end
   xa=(1:size(M,1))*abs(d.vox(1));
   ya=(1:size(M,2))*abs(d.vox(2));
   sx=abs(d.vox(1))*length(xa)*7;
   sy=abs(d.vox(2))*length(ya)*8;
   pos=get(gca,'Position');
   h=pos(4)/max(sx,sy);
   set(gca,'Position',[pos(1)+(pos(3)-sx*h)/2 pos(2)+(pos(4)-sy*h)/2 sx*h sy*h]);
   if isempty(zlim)
      M(isnan(M))=min(min(M));
      imagesc(xa,ya,M');
   else
      M(isnan(M))=min(min(min(M)),min(zlim));
      imagesc(xa,ya,M',zlim);
   end
   set(gca,'XTick',[]); set(gca,'YTick',[]); 
   if length(slices)>1
      %title([strrep(deblank(input_file(1,:)),'_',' ') ', frame ' num2str(frames)]);  
   end
   if length(frames)>1
      %title([strrep(deblank(input_file(1,:)),'_',' ') ', slice ' num2str(slices)]);  
   end
   %  ', max=' num2str(max(M(:))) ', num>3=' num2str(sum(M(:)>3))]);
   axis off;
   for i=1:mx
      for j=1:my
         islice=layout(i,j);
         if islice>0
            slice=Slices(islice);
            frame=Frames(islice);
            if length(slices)>1
               %text((i-1+0.05)*nx*abs(d.vox(1)),(j-1+0.9)*ny*abs(d.vox(2)), ...
               %num2str(slice),'Color','w');
            end
            if length(frames)>1
               %text((i-1+0.9)*nx*abs(d.vox(1)),(j-1+0.9)*ny*abs(d.vox(2)), ...
               %num2str(frame),'Color','w');
            end
         end
      end
   end
else
   d=fmris_read_image(input_file,slices+1,frames);
   if ~isempty(expr)
      eval(expr);
   end
   v=d.data(xr,yr);
   if ~isempty(mask_file)
      v(m.data(xr,yr,1)==0)=mask_value;
   end
   
   if isanal %| isafni
      v=flipdim(v,2);
   end
   
   sx=abs(d.vox(1))*length(xr)*7;
   sy=abs(d.vox(2))*length(yr)*8;
   pos=get(gca,'Position');
   h=(pos(3)+pos(4))/2/max(sx,sy);
   set(gca,'Position',[pos(1)+(pos(3)-sx*h)/2 pos(2)+(pos(4)-sy*h)/2 sx*h sy*h]);
   if isempty(zlim)
      v(isnan(v))=min(min(v));
      imagesc(xp,yp,v');
   else
      v(isnan(v))=min(min(min(v)),min(zlim));
      imagesc(xp,yp,v',zlim);
   end
   grid on; 
   title([strrep(input_file,'_',' ') ', slice ' num2str(slices) ', frame ' num2str(frames)]);
end
xlabel([signx 'x=' num2str(min(xp)) ':' num2str(max(xp)) ' voxels']);
ylabel([signy 'y=' num2str(min(yp)) ':' num2str(max(yp)) ' voxels']);
colorbar; colormap(spectral);
if isminc | isafni | isnifti
   axis xy; 
end

return





   
