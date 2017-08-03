function values=extract(voxels,input_file_base,extension) 

%EXTRACT
%
% Extracts values of all frames in files at specified voxels 
% (all files must have the same dimensions).
%
% VALUES = EXTRACT ( VOXELS, INPUT_FILE_BASE [, EXTENSION] )
%
% VOXELS is a 3 column matrix of voxel locations in 'register' format,
% i.e. starting at zero.
%
% INPUT_FILE_BASE is a matrix of file bases
%
% EXTENSION is a file extension to be added to INPUT_FILE_BASE.
% Default is [].
%
% VALUES is a #VOXELS x #FILEs x #FRAMES array of extracted values.

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
   extension=[];
end

numvoxels=size(voxels,1);
numfiles=size(input_file_base,1);

input_file=[deblank(input_file_base(1,:)) extension];
d=fmris_read_image(input_file,0,0);
numframes=d.dim(4);
values=zeros(numvoxels,numfiles,numframes);

[base,ext]=fileparts2(input_file(1,:));
isminc=strcmp(lower(ext),'.mnc');

for file=1:numfiles
   input_file=[deblank(input_file_base(file,:)) extension];
   d=fmris_read_image(input_file,0,0);
   if file==1
      i=voxels(:,1);
      if d.vox(1)<0 & isminc
         i=d.dim(1)-1-i;
      end
      j=voxels(:,2);
      if d.vox(2)<0 & isminc
         j=d.dim(2)-1-j;
      end
      k=voxels(:,3);
   end
   for slice=unique(k)'
      d=fmris_read_image(input_file,slice+1,1:numframes);
      for kk=find(k==slice)'
         values(kk,file,:)=squeeze(d.data(i(kk)+1,j(kk)+1,1,:));
      end
   end
end

return
