function X_interp=fmri_interp(times,X,frametimes,slicetimes,method);

%FMRI_INTERP resamples X at frametimes offset by slicetimes
%
% X_INTERP = FMRI_INTERP(TIMES, X, FRAMETIMES [, SLICETIMES, 
%                  [, METHOD]]); 
%
% Uses matlab INTERP1 to interpolate the columns of X.
% TIMES is a matrix of times at which X is measured.
% If TIMES has one column, it is used for all columns of X.
% X_INTERP is a 3D array whose first two dimensions are the columns 
% of X resampled at FRAMETIMES + SLICETIMES(SLICE), and whose
% third dimension is SLICE, SLICE=1:length(SLICETIMES). 
% Default SLICETIMES is 0, in which case FMRI_INTERP is INTERP1.
% Uses method = METHOD (see help interp1); default is 'spline'.

if nargin < 4
   slicetimes=0;
end
if nargin < 5
   method='spline';
end

numframes=length(frametimes);
numslices=length(slicetimes);
numxs=size(X,2);
X_interp=zeros(numframes,numxs,numslices);
for slice=1:numslices
   if length(size(times))==1
      X_interp(:,:,slice)=interp1(times,X, ...
         frametimes+slicetimes(slice),method);
   else
      for col=1:size(times,2)
         X_interp(:,:,slice)=interp1(times(:,col),X(:,col), ...
            frametimes+slicetimes(slice),method);
      end
   end
end

return
   
