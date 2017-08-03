function resels=mask_resels(fwhm_info, dfw, mask_file, mask_thresh, ...
    output_file_base)

%MASK_RESELS finds (non-)isotropic resels for a smooth mask
%
% RESELS = MASK_RESELS( FWHM_INFO, DFW, MASK_FILE [, MASK_THRESH
%            [, OUTPUT_FILE_BASE]]);
%
% FWHM_INFO: for isotropic data, the scalar FWHM (mm); for non-
% isotropic data, the name of a file(s) of whitened normalized residuals
% e.g. wresid file from fmrilm.m or multistat.m with which_stats(7)=1.  
% It can be a single 4D image file with multiple frames, 
% or a matrix of image file names, each with a single 3D frame, either 
% ANALYZE (.img) or MINC (.mnc) format. Extra blanks are ignored. File 
% separator can be / or \ on Windows. Gzipped files are gunzipped on unix.     
% The first 3 dimesnions are space, the fourth is the whitened residuals.
%
% DFW=[DFW1 DFW2]: If FWHM_INFO is a wresid file from multistat.m, then 
% DFW1 and DFW2 are the numerator and denominator df of the wresid's. 
% Set DFW=[DF_RESID DF] from output of multistat.m. This is only used to 
% correct the bias in estimating the resels. It is not needed for wresid 
% files from fmrilm, where DFW1=DFW2=DF, since the resels are unbiased
% if DFW1=DFW2, in which case set DFW=[] to skip the bias correction. 
% 
% MASK_FILE: A continuous volume to define the mask. If it has multiple
% frames, then the first frame is used. 
%
% MASK_THRESH defines the search volume as the first frame of MASK_FILE 
% > MASK_THRESH. If MASK_THRESH is a vector [a b], a<=b, then mask 
% is a < MASK_FILE <= b. If empty (default), calls fmri_mask_thresh. 
%
% RESELS is the vector of resels for dimensions 0:3 suitable as the
% first parameter of stat_threshold.m.
%
% OUTPUT_FILE_BASE: If given, calls tet_curv.m to calculate 1D
% Lipschitz Killing curvature density in OUTPUT_FILE_BASE_curv.ext. 
% The integral of the curvature density^3 is the contribution from the 
% interior of the search region to the 1D Lipschitz-Killing curvature. 
%
% If FWHM_INFO is a scalar, it makes a 4D image file of Euclidean
% coordinates of the mask voxels. It passes this file or the wresid file 
% to mask_mesh.m which 'moves' the coordinate mesh out to the boundary
% of the mask. Then it calls mesh_tet.m to find edge lengths of a 
% tetrahedral mesh that fills the mask. This and the extended mask are 
% passed to intrinsicvol.m to find the intrinsic volumes, which are 
% finally converted to resels. The program needs scratch disk space, 
% which it finds in either the directory of mask_file (isotropic case) 
% or wresid_file (non-isotropic case) using these files minus extension 
% as a base. It cleans up as it goes, producing spurious warning 
% messages if it tries to delete a file that was not already created. 

if nargin<4 
   mask_thresh=fmri_mask_thresh(mask_file);
end
if isempty(mask_thresh)
   mask_thresh=fmri_mask_thresh(mask_file);
end

if isstr(fwhm_info)
   fwhm=sqrt(4*log(2));
   [base,ext]=fileparts2(fwhm_info(1,:));
else
   fwhm=fwhm_info;
   d=fmris_read_image(mask_file,0,0);
   [base,ext]=fileparts2(mask_file(1,:));
   base=[base '_coord'];
   if exist([base ext],'file'); delete([base ext]); end
   d.file_name=[base ext];
   d.dim(4)=3;
   d.parent_file=mask_file;
   [x,y]=ndgrid((0:d.dim(1)-1)*d.vox(1)+d.origin(1), ...
      (0:d.dim(2)-1)*d.vox(2)+d.origin(2));
   d.data=zeros(d.dim(1),d.dim(2),3);
   d.data(:,:,1)=x;
   d.data(:,:,2)=y;
   for k=1:d.dim(3)
      d.data(:,:,3)=ones(d.dim(1),d.dim(2))*((k-1)*d.vox(3)+d.origin(3));
      fmris_write_image(d,k,1:3);
   end
end

if exist([base '_mesh' ext],'file');  delete([base '_mesh' ext]); end
if exist([base '_mask' ext],'file');  delete([base '_mask' ext]); end
mask_mesh([base ext], base, mask_file, mask_thresh, isstr(fwhm_info));
if ~isstr(fwhm_info) & exist([base ext],'file'); delete([base ext]); end

if exist([base '_tet' ext],'file');  delete([base '_tet' ext]); end
mesh_tet([base '_mesh' ext], base);
if exist([base '_mesh' ext],'file');  delete([base '_mesh' ext]); end

invol=intrinsicvol([base '_tet' ext], [base '_mask' ext], 0.5);

if nargin>4 & ~isempty(dfw)
   tet_curv([base '_tet' ext], output_file_base, dfw);
end

if exist([base '_mask' ext],'file');  delete([base '_mask' ext]); end
if exist([base '_tet' ext],'file');  delete([base '_tet' ext]); end

resels=invol./fwhm.^(0:3);

if ~isempty(dfw)
   df_resid=dfw(1);
   df=dfw(2);
   alphar=1/2;
   dr=df_resid/df;
   biasr=ones(1,4);
   for D=1:3
      dv=df_resid-dr-(0:D-1);
      biasr(D+1)=exp(sum(gammaln(dv/2+alphar)-gammaln(dv/2)) ...
         +gammaln(df/2-D*alphar)-gammaln(df/2))*dr^(-D*alphar);
   end
   biasr
   resels=resels./biasr;
end
resels

return
