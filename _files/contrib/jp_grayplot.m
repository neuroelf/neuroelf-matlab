function grayplot(varargin)

% This is a limited script just intended to show you how to use Matlab to
% make slice images of .nii.gz files, and also how to make grayplots. You
% run it right in the current folder with the test data provided.
%
%
% IMPORTANT
% You need FSL installed because this script makes system calls to 'fslhd'
%
% You also need to put the NifTI matlab library on your matlab path
%   (the library is included in this folder and is added right below here)
%
% This script was developed on a mac system, it wasn't tested on Windows

pdir=pwd;
%addpath(genpath([pwd]));

% load the images in RAS orientation
fname.bold=[ pdir '/data/bold1_at_EPI.nii.gz'];
fname.mpr=[ pdir '/data/mprage_at_EPI.nii.gz'];
fname.ribbonmask=[ pdir '/data/ribbonmask_at_EPI.nii.gz'];

img.bold=nii_give_RAS(fname.bold);
img.mpr=nii_give_RAS(fname.mpr);
img.ribbonmask=nii_give_RAS(fname.ribbonmask);

% prepare the output directory
odir = [pdir '/output'];
system(['mkdir -p ' odir]);



%%%
%%% First show how to make slices of the mprage (or any image)
%%%

% how big are the arrays
dd=size(img.bold);

% divide the dimensions into 11 slices (0.05-1 in 0.1 steps)
% slc will be 3 columns, corresponding to slices in X, Y, and Z dimensions
slc=ceil([.05:.1:1]'*dd(1:3));

% grab slices of each orientation
im=img.mpr;
t=1; % 4th dimension index to sample from
orin=0; % horizontal concatenation of images
arr.side=multslc(im,'side',slc(:,1),t,orin);
arr.front=multslc(im,'front',slc(:,2),t,orin);
arr.top=multslc(im,'top',slc(:,3),t,orin);

% now plot them and save an image
close all;
h=figure;
set(h,'position',[100 100 1920 1080]);
set(h,'visible','on');

subplot(3,1,1);
imagesc(arr.side);
colormap(gray);
axis('equal','off');

subplot(3,1,2);
imagesc(arr.front);
colormap(gray);
axis('equal','off');

subplot(3,1,3);
imagesc(arr.top);
colormap(gray);
axis('equal','off');

oname = [odir '/00_OUTPUT_MPRAGE_IMAGESC.png'];
set(h,'paperpositionmode','auto');
print(gcf,'-dpng',oname);

% NOTE
% axis('equal') is appropriate because these are isotropic voxels
% if anositropic you need to set 'daspect' appropriately



%%%
%%% NOW PUT AN OVERLAY ON
%%%

% here we will convert the array to grayscale and use image, not imagesc
% this fixes the bounds of the image. Imagesc does it automatically

% I already know the histogram of this image's intensities is ~0-450
glimz=[0 400];

% recreate the above image with a fixed scale set by glimz

rgb.mprside=convert2rgb(arr.side,gray(256),glimz);
rgb.mprfront=convert2rgb(arr.front,gray(256),glimz);
rgb.mprtop=convert2rgb(arr.top,gray(256),glimz);

close all;
h=figure;
set(h,'position',[100 100 1920 1080]);
set(h,'visible','on');

subplot(3,1,1);
image(rgb.mprside);
axis('equal','off');

subplot(3,1,2);
image(rgb.mprfront);
axis('equal','off');

subplot(3,1,3);
image(rgb.mprtop);
axis('equal','off');

oname = [odir '/00_OUTPUT_MPRAGE_IMAGE.png'];
set(h,'paperpositionmode','auto');
print(gcf,'-dpng',oname);

% now get the same slices of the ribbon image
im2=img.ribbonmask;
arr2.side=multslc(im2,'side',slc(:,1),t,orin);
arr2.front=multslc(im2,'front',slc(:,2),t,orin);
arr2.top=multslc(im2,'top',slc(:,3),t,orin);

% overlay the ribbon on the mpr underlay with a red color
rrgb=[0 0 0;1 0 0];
rlimz=[0 1];
rgb.ribside=convert2rgb(arr2.side,rrgb,rlimz);
rgb.ribfront=convert2rgb(arr2.front,rrgb,rlimz);
rgb.ribtop=convert2rgb(arr2.top,rrgb,rlimz);

o.side=olay(rgb.mprside,rgb.ribside,arr2.side);
o.front=olay(rgb.mprfront,rgb.ribfront,arr2.front);
o.top=olay(rgb.mprtop,rgb.ribtop,arr2.top);

close all;
h=figure;
set(h,'position',[100 100 1920 1080]);
set(h,'visible','on');

subplot(3,1,1);
image(o.side);
axis('equal','off');

subplot(3,1,2);
image(o.front);
axis('equal','off');

subplot(3,1,3);
image(o.top);
axis('equal','off');

oname = [odir '/00_OUTPUT_MPRAGE_IMAGE_RIBBON.png'];
set(h,'paperpositionmode','auto');
print(gcf,'-dpng',oname);



%%%
%%% NOW TAKE EPI PICTURES IN TIME
%%%

im=img.bold;
slc=ceil([.3 .5 .7]'*dd(1:3));
blimz=[0 2000]; % again, I know this image's intensity histogram

% just the first 5 volumes
for t=1:10
    arr.top=multslc(im,'top',slc(:,3),t,orin);
    o.top=convert2rgb(arr.top,gray(256),blimz);
    
    close all;
    h=figure;
    set(h,'position',[100 100 1920 1080]);
    set(h,'visible','on');

    
    image(o.top);
    axis('equal','off');
    
    % sprintf can pad numbers with zeros, more useful than num2str
    oname = [odir '/01_OUTPUT_BOLD_t' sprintf('%3.3d',t) '.png'];
    set(h,'paperpositionmode','auto');
    print(gcf,'-dpng',oname);
    
end


%%%
%%% AND NOW A GRAY PLOT (or 3)
%%%

im=reshape(img.bold,[dd(1)*dd(2)*dd(3) dd(4)]);
fd=load([pdir '/data/fd.txt']);

close all;
h=figure;
set(h,'position',[100 100 1920 1080]);
set(h,'visible','on');

% plot FD as a reference
subplot(5,1,1);
plot([1:dd(4)],fd,'r');
xlim([1 dd(4)]);
ylim([0 2]);

% plot the gray ribbon timeseries
subplot(5,1,[2:5]);
imagesc(im(~~img.ribbonmask(:),:));
colormap(gray);

oname = [odir '/02_grayplot.png'];
set(h,'paperpositionmode','auto');
print(gcf,'-dpng',oname);



% mean intensity of each voxel dominates this image, not helpful
% remove mean and trend terms

% mean and trend regressors
r0=ones(1,dd(4));
r1=linspace(0,1,dd(4));
% also set up a temporal mask ignoring the first 4 volumes
tmask=r0;
tmask(1:4)=0;

% remove those terms
im2=myregress([r0;r1],im,~~tmask);

close all;
h=figure;
set(h,'position',[100 100 1920 1080]);
set(h,'visible','on');

% plot FD as a reference
subplot(5,1,1);
plot([1:dd(4)],fd,'r');
xlim([1 dd(4)]);
ylim([0 2]);

% plot the gray ribbon timeseries
% due to the first 4 TRs (high intensity) you gotta set the scale
% note that these images are not "mode 1000" scaled
% if they were I would set it to -20 to 20 scale
subplot(5,1,[2:5]);
imagesc(im2(~~img.ribbonmask(:),:),[-30 30]);
colormap(gray);

oname = [odir '/02_grayplot_dmdt.png'];
set(h,'paperpositionmode','auto');
print(gcf,'-dpng',oname);


% let's add a DVARS and global signal trace for good measure
dv=rms(diff(im2(~~img.ribbonmask(:),:),1,2),1);
dv=[0 dv];
gs=nanmean(im2(~~img.ribbonmask(:),:),1);

close all;
h=figure;
set(h,'position',[100 100 1920 1080]);
set(h,'visible','on');

% plot FD as a reference
subplot(5,1,1);
plot([1:dd(4)],fd,'r');
xlim([1 dd(4)]);
ylim([0 2]);

% add the dv plot
ax1=gca;
ax1_pos=get(ax1,'Position');

ax2 = axes('Position',ax1_pos,...
    'YAxisLocation','right',...
    'Color','none');
hold on;
plot(ax2,[1:dd(4)],dv,'b');
set(ax2,'ylim',[0 50],'ytick',[],'xlim',[1 dd(4)],'xtick',[]);
set(ax2,'color','none','YAxisLocation','right','ycolor',[0 0 0],'box','off');
hold off;

% add the gs plot
ax3 = axes('Position',ax1_pos,...
    'YAxisLocation','right',...
    'Color','none');
hold on;
plot(ax3,[1:dd(4)],gs,'k');
set(ax3,'ylim',[-50 50],'ytick',[],'xlim',[1 dd(4)],'xtick',[]);
set(ax3,'color','none','YAxisLocation','right','ycolor',[0 0 0],'box','off');
hold off;

subplot(5,1,[2:5]);
imagesc(im2(~~img.ribbonmask(:),:),[-30 30]);
colormap(gray);

oname = [odir '/02_grayplot_dmdt_dvgs.png'];
set(h,'paperpositionmode','auto');
print(gcf,'-dpng',oname);

close all;

%%%
%%% TIPS
%%%
%
% This is a barebones script, intended to give you the skeleton for making
% images of whatever kind you like in Matlab. Here's what I wanted to
% illustrate:
%
% 1. Load an image into matlab in a known orientation
% 2. Illustrate slices with automatic color scaling (imagesc)
% 3. Control the color scaling, and overlay slices (image)
% 4. Show how to make slices in time
% 5. Show how to make gray plots
% 6. Show subplot to divide images
% 7. Show how to make large figures
% 8. Show quick regression code with temporal masking
% 9. Show how to overlay traces on the same plot
% 10. Show a DVARS and GS calculation
%
% I use subaxis instead of subplot usually, and often some other
% community-made functions from Matlabcentral. Here I kept it simple so
% that a stock install with few libraries are needed.
%
% Also if making lots of images you might want to turn 'visible' to 'off'
% This makes figures the same way, but as a background process.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img dims o]=nii_give_RAS(im,varargin)

% NIFTI qform and sform codes are relative to RAS
%
% does not matter which code used for showing the arrays
% but we do want to flip the arrays to be RAS or close to it
%
% also we need to scale the image, sometimes uint8 and a scaling factor are used to save space instead of single or double (float)
%
% If a filename or structure is given as the sole input, it returns the
% RAS-oriented 4D array of the input as a double
%
% If an additional array is provided, this array is converted from RAS to
% the orientation of the first input but is otherwise not touched


if ischar(im)
    % if a filename given
    tmp=load_untouch_nii(im);
else
    % if a cell array given (from load_untouch_nii)
    tmp=im;
    im=tmp.fileprefix;
end

% if flipping RAS to the original orientation
if ~isempty(varargin)
    tmp.img=varargin{1,1};
end

tmp.img=double(tmp.img);

if isempty(varargin)
if tmp.hdr.dime.scl_slope~=0
    fprintf('Rescaling the image\n');
    tmp.img=double(tmp.img)*tmp.hdr.dime.scl_slope+tmp.hdr.dime.scl_inter;
end
end

dims=tmp.hdr.dime.pixdim(2:4);

[jk o.sx]=system(['fslhd ' im ' | grep sform_xorient | awk ''{print $2}''' ]); o.sx(end)=[];
[jk o.sy]=system(['fslhd ' im ' | grep sform_yorient | awk ''{print $2}''' ]); o.sy(end)=[];
[jk o.sz]=system(['fslhd ' im ' | grep sform_zorient | awk ''{print $2}''' ]); o.sz(end)=[];
[jk o.scode]=system(['fslhd ' im ' | grep sform_code | awk ''{print $2}''' ]); o.scode(end)=[]; o.scode=str2num(o.scode);
[jk o.sRx]=system(['fslhd ' im ' | grep sto_xyz:1 | awk ''{print $2 " " $3 " " $4}''' ]); o.sRx(end)=[]; o.sR(1,1:3)=str2num(o.sRx);
[jk o.sRy]=system(['fslhd ' im ' | grep sto_xyz:2 | awk ''{print $2 " " $3 " " $4}''' ]); o.sRy(end)=[]; o.sR(2,1:3)=str2num(o.sRy);
[jk o.sRz]=system(['fslhd ' im ' | grep sto_xyz:3 | awk ''{print $2 " " $3 " " $4}''' ]); o.sRz(end)=[]; o.sR(3,1:3)=str2num(o.sRz);


[jk o.qx]=system(['fslhd ' im ' | grep qform_xorient | awk ''{print $2}''' ]); o.qx(end)=[];
[jk o.qy]=system(['fslhd ' im ' | grep qform_yorient | awk ''{print $2}''' ]); o.qy(end)=[];
[jk o.qz]=system(['fslhd ' im ' | grep qform_zorient | awk ''{print $2}''' ]); o.qz(end)=[];
[jk o.qcode]=system(['fslhd ' im ' | grep qform_code | awk ''{print $2}''' ]); o.qcode(end)=[]; o.qcode=str2num(o.qcode);
[jk o.qRx]=system(['fslhd ' im ' | grep qto_xyz:1 | awk ''{print $2 " " $3 " " $4}''' ]); o.qRx(end)=[]; o.qR(1,1:3)=str2num(o.qRx);
[jk o.qRy]=system(['fslhd ' im ' | grep qto_xyz:2 | awk ''{print $2 " " $3 " " $4}''' ]); o.qRy(end)=[]; o.qR(2,1:3)=str2num(o.qRy);
[jk o.qRz]=system(['fslhd ' im ' | grep qto_xyz:3 | awk ''{print $2 " " $3 " " $4}''' ]); o.qRz(end)=[]; o.qR(3,1:3)=str2num(o.qRz);

if o.scode>0
    if ~isempty(varargin)
        [img perm]=perm_and_flip(tmp,o.sR,1);
    else
        [img perm]=perm_and_flip(tmp,o.sR);
    end
elseif o.qcode>0
    if ~isempty(varargin)
        [img perm]=perm_and_flip(tmp,o.qR,1);
    else
        [img perm]=perm_and_flip(tmp,o.qR);
    end
else
    fprintf('ERROR WARNING in nii_give_RAS.m: qformcode and sformcode of %s are both 0 so conversion to RAS is unknown',niifile);
end

if numel(perm)==4
    perm(4)=[];
end
dims=dims(perm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img perm] = perm_and_flip(tmp,R,varargin)

% define the principal directions
signR=sign(R);
absR=abs(R);

[s ss]=sort(absR(:));
s3(ss)=[1:9];
s3=reshape(s3,[3 3]);
topR=zeros(3,3);
[a b]=find(s3==max(s3(:)));
s3(a,:)=[0 0 0];
s3(:,b)=[0 0 0]';
topR(a,b)=1;
[a b]=find(s3==max(s3(:)));
s3(a,:)=[0 0 0];
s3(:,b)=[0 0 0]';
topR(a,b)=1;
[a b]=find(s3==max(s3(:)));
s3(a,:)=[0 0 0];
s3(:,b)=[0 0 0]';
topR(a,b)=1;



flipR=topR.*signR;
% this is to go from RAS to original orientation
if ~isempty(varargin)
    flipR=inv(flipR);
end

for i=1:3
    [jk perm(i) flipr(i) ]=find(flipR(i,:));
end

% have to preserve the time dimensions explicitly if present
d=size(tmp.img);
if numel(d)>3
    perm(4)=4;
end


img=permute(tmp.img,perm);

% now flip where needed to get RAS
for i=1:3
    if flipr(i)<0
        img=flipdim(img,i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img]=showbrik(V,orien,slc,t,varargin)

%d1ATL: left to right; front down; bottom left
%d2ATL: back to front; left top; bottom left
%d3ATL: bottom to top; left top; front right
%d1EPI: left to right; front down; bottom left
%d2EPI: back to front; left top; bottom left
%d3EPI: bottom to top; left top; front right

% presumes RAS orientation, just saying.

switch orien
    case 'side' % rot 90 ccw: L=back;R=front; T=top; B=down
        img=(rot90(squeeze(V(slc,:,:,t))));
    case 'front' % rot 90 ccw: L=left; R=right; T=top; B=down
        img=(rot90(squeeze(V(:,slc,:,t))));
    case 'top' % rot 90 ccw: L=left; R=right; T=front; B=back
        img=(rot90(squeeze(V(:,:,slc,t))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I
function [ming]=multslc(V,orien,slcz,t,varargin)

ming=[];

for z=1:numel(slcz)
    [img]=showbrik(V,orien,slcz(z),t);
    if isempty(varargin) | varargin{1,1}==0
        ming=[ming img];
    else
        ming=[ming; img];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rgbval] = convert2rgb(tmp,scl,varargin)

% tmp is an image
% scl is like jet or copper or self-made
% varargin is limits for scaling the image



%if varargin then we're scaling image
if ~isempty(varargin)
    minz=varargin{1,1}(1);
    maxz=varargin{1,1}(2);
    
    % add an extra column with limits
    tmp=[tmp zeros(size(tmp,1),1)];
    tmp(end,end)=minz;
    tmp(end-1,end)=maxz;
    
    tmp(tmp<minz)=minz;
    tmp(tmp>maxz)=maxz;
    tmp=tmp-min(tmp(:));
    tmp=tmp/max(tmp(:));
    tmp=uint8(tmp*256);
    
    % delete the extra column
    tmp(:,end)=[];
else
    tmp=tmp-min(tmp(:));
    tmp=tmp/max(tmp(:));
    tmp=uint8(tmp*256);
end
rgbval=ind2rgb(tmp,scl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nimg] = olay(bimg,timg,msk)

% overlay one image on another

msk=~~msk;
msk=repmat(msk,[1 1 size(bimg,3)]);
nimg=bimg;
nimg(msk)=timg(msk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [resid pred b] = removepoly(tc,ord,varargin)

% remove legendre polynomials up to ORD
% 0 - mean
% 1 - trend
% 2 - bow
% 3 - sinusoid
% 4 - w-shape
% 5 - w-sinusoid
% etc etc
%
% be careful, these are high-pass filters also

% presumes tc is vox x time

% tmask is 1 when used, 0 when censored

d=size(tc);

for i=0:ord
    b=legendre(i,linspace(-1,1,d(2)));
    if i==0
        r=b(1,:);
    else
        r=[r; b(1,:)];
    end
end


if isempty(varargin)
    % use all timepoints
    b=r'\tc';
    pred=r'*b;
    resid=tc-pred';
    
else
    % use only specified timepoints
    tmask=~~varargin{1,1};
    b=r(:,tmask)'\tc(:,tmask)';
    pred=r'*b;
    resid=tc-pred';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [resid pred b]=myregress(r,tc,varargin)

% presumes vox x time input variable structure

if isempty(varargin)
    % use all timepoints
    b=r'\tc';
    pred=r'*b;
    resid=tc-pred';
    
else
    % use only specified timepoints
    tmask=varargin{1,1};
    b=r(:,tmask)'\tc(:,tmask)';
    pred=r'*b;
    resid=tc-pred';
end

