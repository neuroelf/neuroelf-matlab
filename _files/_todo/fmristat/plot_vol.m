function [m,xl,ycx,yl,xcy]=plot_vol(input_file_x,input_file_y,xlims,ylims,...
   mask_file,mask_thresh,nx,ny)

%PLOT_VOL scatter plot of one volume versus another.
%
% [M,X,YCX,Y,XCY] = PLOT_VOL( INPUT_FILE_X, INPUT_FILE_Y, XLIMS, YLIMS, ...
%                  [MASK_FILE, MASK_THRESH [, NX NY]])
%
% Plots values of INPUT_FILE_Y in the range [YLIMS(1) YLIMS(2)] versus 
% values of INPUT_FILE_X in the range [XLIMS(1) XLIMS(2)] where MASK_FILE 
% is above MASK_THRESH (ignored if empty; empty by default). If MASK_THRESH 
% is a vector [a b], a<=b, then mask is a<MASK_FILE<=b. Output is a NX by NY
% (default 64 by 64) matrix M of counts. INPUT_FILE_X and/or 
% INPUT_FILE_Y can be 3D arrays instead of .img or .mnc file names.
% X=x values, YCX=mean(Y) conditional on X;
% Y=Y values, XCY=mean(X) conditional on Y.

if nargin<5
   mask_file=[]
   mask_thresh=[]
end
if nargin<7
   nx=64;
   ny=64;
end

if ~isempty(mask_thresh)
   mask_thresh1=mask_thresh(1);
   if length(mask_thresh)>=2
      mask_thresh2=mask_thresh(2);
   else
      mask_thresh2=Inf;
   end
end

if ischar(input_file_x)
   x=fmris_read_image(input_file_x,0,0);
else
   x.dim=size(input_file_x);
end
h=zeros(1,(nx+2)*(ny+2));
mask=find(ones(x.dim(1),x.dim(2)));
for slice=1:x.dim(3)
   if ischar(input_file_x)
      x=fmris_read_image(input_file_x,slice,1);
   else
      x.data=input_file_x(:,:,slice);
   end
   if ischar(input_file_y)
      y=fmris_read_image(input_file_y,slice,1);
   else
      y.data=input_file_y(:,:,slice);
   end
   if ~isempty(mask_file)
      if ischar(mask_file)
         d=fmris_read_image(mask_file,slice,1);
      else
         d.data=mask_file(:,:,slice);
      end
      mask=find(d.data>mask_thresh1 & d.data<=mask_thresh2);
   end
   xx=round((x.data(mask)-xlims(1))/(xlims(2)-xlims(1))*nx+0.5);
   xx=max(min(xx,nx+1),0)+1;
   yy=round((y.data(mask)-ylims(1))/(ylims(2)-ylims(1))*ny+0.5);
   yy=max(min(yy,ny+1),0)+1;
   h=h+hist(xx+(yy-1)*(nx+2),1:((nx+2)*(ny+2)));
end
mm=reshape(h,nx+2,ny+2);
m=mm(2:(nx+1),2:(ny+1));
xl=((1:nx)-0.5)/nx*(xlims(2)-xlims(1))+xlims(1);
yl=((1:ny)-0.5)/ny*(ylims(2)-ylims(1))+ylims(1);

subplot(2,2,3)
imagesc(xl,yl,m'); axis xy; title('counts');
if ischar(input_file_x)
   xlabel(input_file_x); 
end
if ischar(input_file_y)
   ylabel(input_file_y);
end

subplot(2,2,4)
mt=repmat(max(m,[],1),nx,1);
mdx=m./(mt+(mt<=0)).*(mt>0);
imagesc(xl,yl,mdx'); axis xy; title('divided by max in x');

subplot(2,2,1)
mt=repmat(max(m,[],2),1,ny);
mdy=m./(mt+(mt<=0)).*(mt>0);
imagesc(xl,yl,mdy'); axis xy; title('divided by max in y');

subplot(2,2,2)
ycx=m*yl'./sum(m,2);
xcy=xl*m./sum(m,1);
plot(xl,ycx,xcy,yl); 
xlim([min(xl) max(xl)]) ; ylim([min(yl) max(yl)]);
title('means')

return











