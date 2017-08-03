function [volume, num_voxels]=mask_vol(mask_file,mask_thresh);

%MASK_VOL finds the volume (mm^3) and number of voxels in a mask.
%
% [VOLUME, NUM_VOXELS] = MASK_VOL( MASK_FILE [, MASK_THRESH ])
%
% MASK_THRESH defines the search volume as the first frame of MASK_FILE 
% > MASK_THRESH. If MASK_THRESH is a vector [a b], a<=b, then mask 
% is a < MASK_FILE <= b. If empty (default), calls fmri_mask_thresh. 
%
% VOLUME is the volume (mm^3 in 3D) or area (mm^2 in 2D) in the mask.
%
% NUM_VOXELS is the number of voxels (3D) or pixels (2D) in the mask.

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
   mask_thresh=[]
end
if isempty(mask_thresh)
   mask_thresh=fmri_mask_thresh(mask_file);
end
mask_thresh1=mask_thresh(1);
if length(mask_thresh)>=2
   mask_thresh2=mask_thresh(2);
else
   mask_thresh2=Inf;
end
d=fmris_read_image(mask_file,0,0);
numslices=d.dim(3);
d=fmris_read_image(mask_file,1:numslices,1);
num_voxels=sum(d.data(:)>mask_thresh1 & d.data(:)<=mask_thresh2);
if numslices==1
   D=2;
else
   D=3;
end
volume=num_voxels*abs(prod(d.vox(1:D)));

return







