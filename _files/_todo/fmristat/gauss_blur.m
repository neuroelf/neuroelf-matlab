function blur_file=gauss_blur(input_file,fwhm,output_file_base)

%GAUSS_BLUR A Gaussian smoother that works well for fMRI data. 
%
%  BLUR_FILE = GAUSS_BLUR( INPUT_FILE [, FWHM, [, OUTPUT_FILE_BASE ]] )
%
% Gaussian bluring of FWHM (default 10mm) of all frames in INPUT_FILE to give
% BLUR_FILE = OUPUT_FILE_BASE_blurFWHM.mnc or .img whose default is INPUT_FILE 
% minus extension. If FWHM is a scalar, uses the same value in all directions;
% if FWHM is a vector, uses the components in the x, y and z directions.
% If INPUT_FILE is a 4-D array, BLUR_FILE is the blurred 4-D array.
%
% Method: fMRI data often has a few slices that do not cover the entire brain 
% in the z direction, but the slices often extend well beyond the brain in  
% the x and y directions. Uses the MATLAB function CONV2 to convolve with a
% Gaussian filter in the x and y directions, so in effect the slices
% are padded with zeros before smoothing. Then uses straight matrix 
% mutiplication in the z direction, using a filter truncated by the slices 
% and rescaled to sum to 1, so that means are preserved in the z direction. 

%############################################################################
% COPYRIGHT:   Copyright 2001 K.J. Worsley, 
%              Department of Mathematics and Statistics,
%              McConnell Brain Imaging Center, 
%              Montreal Neurological Institute,
%              McGill University, Montreal, Quebec, Canada. 
%              worsley@math.mcgill.ca
%
%              Permission to use, copy, modify, and distribute this
%              software and its documentation for any purpose and without
%              fee is hereby granted, provided that the above copyright
%              notice appear in all copies.  The author and McGill University
%              make no representations about the suitability of this
%              software for any purpose.  It is provided "as is" without
%              express or implied warranty.
%############################################################################

% Defaults:

if nargin < 2
   fwhm=10
end

if isstr(input_file)
   [base,ext]=fileparts2(input_file(1,:));
   if nargin < 3 
      output_file_base=base;
   end
   d=fmris_read_image(input_file,0,0);
   DimSizes=d.dim;
   numframes=DimSizes(4);
   numslices=DimSizes(3);
   numys=DimSizes(2);
   numxs=DimSizes(1);
   numpix=numxs*numys;
   Steps=d.vox;
   
   if length(fwhm)==1
      blur_file=[deblank(output_file_base) '_blur' num2str(round(fwhm)) ext]
      fwhm=repmat(fwhm,1,3);
   else
      blur_file=[deblank(output_file_base) '_blur' num2str(round(fwhm(1))) '-' ...
            num2str(round(fwhm(2))) '-' num2str(round(fwhm(3))) ext]
   end
   if exist(blur_file,'file'); delete(blur_file); end
   
   out_handle=d;
   out_handle.file_name=blur_file;
   out_handle.parent_file=input_file;
   
   fwhm_x=fwhm(1)/abs(Steps(1));
   ker_x=exp(-(-ceil(fwhm_x):ceil(fwhm_x)).^2*4*log(2)/fwhm_x^2);
   ker_x=ker_x/sum(ker_x);
   fwhm_y=fwhm(2)/abs(Steps(2));
   ker_y=exp(-(-ceil(fwhm_y):ceil(fwhm_y)).^2*4*log(2)/fwhm_y^2);
   ker_y=ker_y/sum(ker_y);
   fwhm_z=fwhm(3)/abs(Steps(3));
   ker_z=exp(-(0:(numslices-1)).^2*4*log(2)/fwhm_z^2);
   K=toeplitz(ker_z);
   K=K./(ones(numslices)*K);
   
   for frame=1:numframes
      frame
      out_handle.data=zeros(numpix,numslices);
      for slice=1:numslices
         d=fmris_read_image(input_file,slice,frame);
         out_handle.data(:,slice)=reshape(conv2(ker_x,ker_y,d.data,'same'),numpix,1);   
      end
      out_handle.data=reshape(out_handle.data*K,[numxs numys numslices]);
      fmris_write_image(out_handle,1:numslices,frame);   
   end
   
else
   numxs=size(input_file,1);
   numys=size(input_file,2);
   numslices=size(input_file,3);
   numframes=size(input_file,4);
   numpix=numxs*numys;
   
   if length(fwhm)==1
      fwhm=repmat(fwhm,1,3);
   end
   fwhm_x=fwhm(1)
   ker_x=exp(-(-ceil(fwhm_x):ceil(fwhm_x)).^2*4*log(2)/fwhm_x^2);
   ker_x=ker_x/sum(ker_x);
   fwhm_y=fwhm(2)
   ker_y=exp(-(-ceil(fwhm_y):ceil(fwhm_y)).^2*4*log(2)/fwhm_y^2);
   ker_y=ker_y/sum(ker_y);
   fwhm_z=fwhm(3)
   ker_z=exp(-(0:(numslices-1)).^2*4*log(2)/fwhm_z^2);
   K=toeplitz(ker_z);
   K=K./(ones(numslices)*K);
   
   blur_file=zeros(numxs,numys,numslices,numframes);
   for frame=1:numframes
      frame
      out_handle.data=zeros(numpix,numslices);
      for slice=1:numslices
         temp(:,slice)=reshape(conv2(ker_x,ker_y,input_file(:,:,slice,frame),'same'),numpix,1);   
      end
      blur_file(:,:,:,frame)=reshape(temp*K,[numxs numys numslices 1]);
   end
   
end

sdif=sqrt(sum(ker_x.^2)*sum(ker_y.^2)*mean(sum(K.^2,1)));

return

