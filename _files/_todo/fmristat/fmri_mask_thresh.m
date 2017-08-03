function mask_thresh = fmri_mask_thresh(fmri_file);

%FMRI_MASK_THRESH finds a threshold that masks the brain from fmri data.
%
% MASK_THRESH = FMRI_MASK_THRESH( FMRI_FILE );
%
% MASK_THRESH = the first 'dip' in the histogram of the first frame of
% FMRI_FILE by searching up to the largest local max; 
% uses max(FMRI_FILE)/4 if unsuccessful.

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

nbin=100;
dbin=10;
dm=fmris_read_image(fmri_file(1,:),0,0);
numslices=dm.dim(3);
dm=fmris_read_image(fmri_file(1,:),1:numslices,1);
[freq, mask]=hist(dm.data(:),nbin);
fm=freq(1:nbin-2*dbin);
f0=freq(1+dbin:nbin-dbin);
fp=freq(1+2*dbin:nbin);
h=(abs(f0-fm)+abs(f0-fp)).*(f0>fp).*(f0>fm);
if any(h)
    mh=min(find(h==max(h))+dbin);
    mask_thresh=max(mask(find(freq==min(freq(1:mh)))))
else
    mask_thresh=max(mask)/4
end

return
