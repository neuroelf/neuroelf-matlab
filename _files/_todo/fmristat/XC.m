function X_confounds=fmri_interact(X_cache,confounds);

%FMRI_INTERACT forms interactions for FMRILM
%
% X_CONFOUNDS = FMRI_INTERACT(X_CACHE,CONFOUNDS); 
%
% X_CONFOUNDS is a 3D array. For each slice (the last dimension), 
% the first two dimensions is a matrix whose columns are:
% [ the columns of CONFOUNDS, followed by 
% the product of each column of CONFOUNDS with each
% column of X_CACHE (CONFOUNDS running fastest) ]. 

[numframes,numX,numbases,numslices]=size(X_cache.X);
numC=size(confounds,2);
if size(size(confounds))==2
   sliceC=ones(1,numslices);
else
   sliceC=1:numslices;
end
X_confounds=zeros(numframes,(numX+1)*numC,numslices);
for slice=1:numslices
   X_confounds(:,1:numC,slice)=confounds(:,:,sliceC(slice)); 
   for i=1:numX
      for j=1:numC
         X_confounds(:,i*numC+j,slice)=X_cache.X(:,i,1,slice).* ...
            confounds(:,j,sliceC(slice));
      end
   end
end

return
   
