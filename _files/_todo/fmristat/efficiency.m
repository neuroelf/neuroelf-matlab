function [sd_ef, Y]=efficiency(X_cache,contrast,exclude,rho,n_temporal, ...
   confounds, contrast_is_delay, mag_t, del_ef, num_hrf_bases, basis_type)

%EFFICIENCY - the relative sd of effects for the efficiency of a design  
%
% [SD_EF, Y] = EFFICIENCY( X_CACHE, CONTRAST [, EXCLUDE [, RHO  
%     [, N_TEMPORAL [, CONFOUNDS [, CONTRAST_IS_DELAY [,MAG_T  
%     [, DEL_EF [, NUM_HRF_BASES [, BASIS_TYPE ]]]]]]]]]
% 
% X_CACHE: A structure usually supplied by FMRIDESIGN.
% X_CACHE.TR: TR, average time between frames (secs), only used for 
% calculating the number of temporal drift terms (see N_TEMPORAL below).  
% X_CACHE.X: A cache of the design matrices (hrf*X) stored as a 4D array. 
% Dim 1: frames; Dim 2: response variables; Dim 3: 4 values, corresponding 
% to the stimuli convolved with: hrf, derivative of hrf, first and second 
% spectral basis functions over the range in SHIFT; Dim 4: slices. 
% X_CACHE.W: A 3D array of coefficients of the basis functions in X_CACHE.X.
% Dim 1: frames; Dim 2: response variables; Dim 3: 5 values: 
% coefficients of the hrf and its derivative, coefficients of the first and 
% second spectral basis functions, shift values.
%
% CONTRAST is a matrix whose rows are contrasts for the
% response variables, and columns are the contrast. Extra columns can 
% be added to estimate contrasts in the polynomial terms and the confounds 
% (in that order - see N_POLY and CONFOUNDS below). 
% 
% EXCLUDE is a list of frames that should be excluded from the
% analysis. This must be used with Siemens EPI scans to remove the
% first few frames, which do not represent steady-state images.
% If NUMLAGS=1, the excluded frames can be arbitrary, otherwise they 
% must be from the beginning and/or end. Default is [1].
% 
% RHO is a row vector of temporal autocorrelations for an AR model.
% The length of RHO determines the order of the AR model.
% Default is 0, i.e. indendent errors. 
% 
% N_TEMPORAL: number of cubic spline temporal trends to be removed per 6 
% minutes of scanner time (so it is backwards compatible). Temporal  
% trends are modeled by cubic splines, so for a 6 minute run, N_TEMPORAL
% <=3 will model a polynomial trend of degree N_TEMPORAL in frame times, 
% and N_TEMPORAL>3 will add (N_TEMPORAL-3) equally spaced knots.
% N_TEMPORAL=0 will model just the constant level and no temporal trends.
% N_TEMPORAL=-1 will not remove anything, in which case the design matrix 
% is completely determined by X_CACHE.X. Default is 3.
%
% CONFOUNDS: A matrix or array of extra columns for the design matrix
% that are not convolved with the HRF, e.g. movement artifacts. 
% If a matrix, the same columns are used for every slice; if an array,
% the first two dimensions are the matrix, the third is the slice.
% For functional connectivity with a single voxel, use fmri_interp
% to resample the reference data at different slice times. 
% Default is [], i.e. no confounds.
% 
% CONTRAST_IS_DELAY is a logical vector indicating if the row of CONTRAST
% is a contrast in the delays (1) or magnitudes (0). Delays are shifts
% of the time origin of the HRF, measured in seconds. Note that you cannot
% estimate delays of the polynomial terms or confounds. Note that
% F statistics are not yet available with this option. If the length
% of CONTRAST_IS_DELAY is less then the number of contrasts, it is padded
% with zeros. Default is 0, i.e. all contrasts refer to magnitudes.
%
% MAG_T is a row vector of t stats of the magnitudes of each 
% response involved in any contrast in the delays (in order), and
%
% DEL_EF are the corresponding delays (temporal shifts) in the hrf. 
% They are only used if you want the sd of contrasts in delays. 
% Default is MAG_T=5, DEL_EF=0 for each such response.
%
% NUM_HRF_BASES is a row vector indicating the number of basis functions
% for the hrf for each response, either 1 or 2 at the moment. At least  
% one basis functions is needed to estimate the magnitude, but two basis 
% functions are needed to estimate the delay. If empty (default), then 
% NUM_HRF_BASES is 2 for each response where CONTRAST is non-zero and  
% CONTRAST_IS_DELAY = 1, otherwise NUM_HRF_BASES = 1. By setting 
% NUM_HRF_BASES = 2 you can allow for an unknown delay, without  
% actually estimating it.  Example:   
%     CONTRAST=[1 -1 0 0; 1 0 -1 0];    CONTRAST_IS_DELAY=[0 1];  
% The first contrast is the magnitude of response_1 - response_2;
% the second contrast is the delay of response_1 - response_3. 
% The default setting of NUM_HRF_BASES is [2 1 2 1]. By setting        
% NUM_HRF_BASES=[2 2 2 1] unknown delays are allowed for the second 
% response but not actually estimated.
%
% BASIS_TYPE selects the basis functions for the hrf used for delay
% estimation, or whenever NUM_HRF_BASES = 2. These are convolved with
% the stimulus to give the responses in Dim 3 of X_CACHE.X:
% 'taylor' - use hrf and its first derivative (components 1 and 2), or 
% 'spectral' - use first two spectral bases (components 3 and 4 of Dim 3).
% Ignored if NUM_HRF_BASES = 1, in which case it always uses component 1,  
% i.e. the hrf is convolved with the stimulus. Default is 'spectral'. 
%
% NUMLAGS is the order (p) of the autoregressive model. Default is 1.
%
% SD_EF is the sd of each effect in the contrast. For magnitudes, 
% SD_EF is relative to the sd of the fMRI data (~6); for delays,
% SD_EF is in seconds and does not depend on the sd of the data.
%
% Y is a matrix whose rows are the linear combinations of the 
% fMRI frames that make up the contrast. There is one row for
% each magnitude contrast (there is no such result for delays).

%############################################################################
% COPYRIGHT:   Copyright 2002 K.J. Worsley
%              Department of Mathematics and Statistics,
%              McConnell Brain Imaging Center, 
%              Montreal Neurological Institute,
%              McGill University, Montreal, Quebec, Canada. 
%              worsley@math.mcgill.ca, liao@math.mcgill.ca
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

if nargin < 3
   exclude=[1];
end
if nargin < 4
   rho= 0;
end
if nargin < 5
   n_temporal=3;
end
if nargin < 6
   confounds=[];
end
if nargin < 7
   contrast_is_delay=0;
end
if nargin < 10
   num_hrf_bases=[];
end
if nargin < 11
   basis_type='spectral';
end

switch lower(basis_type)
case 'taylor',    
   basis1=1;
   basis2=2;
case 'spectral',    
   basis1=3;
   basis2=4;
otherwise, 
   disp('Unknown basis_type.'); 
   return
end

numlags=length(rho);

if ~isempty(X_cache)
   numframes=size(X_cache.X,1);
   numslices=size(X_cache.X,4);
   numresponses=size(X_cache.X,2);
else
   numframes=size(confounds,1);
   numslices=size(confounds,3);
   numresponses=0;
end

allpts = 1:numframes;
allpts(exclude) = zeros(1,length(exclude));
keep = allpts( find( allpts ) );
n=length(keep);

% Create polynomial trends:

n_spline=round(n_temporal*X_cache.TR*n/360);
if n_spline>=0 
   trend=((2*keep-(max(keep)+min(keep)))./(max(keep)-min(keep)))';
   if n_spline<=3
      trend=(trend*ones(1,n_spline+1)).^(ones(n,1)*(0:n_spline));
   else
      trend=(trend*ones(1,4)).^(ones(n,1)*(0:3));
      knot=(1:(n_spline-3))/(n_spline-2)*(max(keep)-min(keep))+min(keep);
      for k=1:length(knot)
         cut=keep'-knot(k);
         trend=[trend (cut>0).*(cut./max(cut)).^3];
      end
   end
else
   trend=[];
end 

% Add confounds:

numtrends=n_spline+1+size(confounds,2);
Trend=zeros(n,numtrends,numslices);
for slice=1:numslices
   if isempty(confounds)
      Trend(:,:,slice)=trend;
   else  
      if length(size(confounds))==2
         Trend(:,:,slice)=[trend confounds(keep,:)];
      else
         Trend(:,:,slice)=[trend confounds(keep,:,slice)];
      end
   end
end

% Make full contrasts:

numcontrasts=size(contrast,1);
contrast_is_delay=[contrast_is_delay ...
      zeros(1,numcontrasts-length(contrast_is_delay))];
if isempty(num_hrf_bases)
   num_hrf_bases=ones(1,numresponses);
end
is_in_delay_contrast=ones(1,numresponses);
for k=find(contrast_is_delay)
   num_hrf_bases(find(contrast(k,:)~=0))=2;
   is_in_delay_contrast(find(contrast(k,:)~=0))=1;
end
num_in_delay_contrast=sum(is_in_delay_contrast);
contrasts=[contrast zeros(numcontrasts,numresponses+numtrends-size(contrast,2))];

% Check for estimability:

tolerance=0.0000001;
for slice=1:numslices
   if ~isempty(X_cache)
      X=[squeeze(X_cache.X(keep,:,1,slice)) Trend(:,:,slice)];
   else
      X=Trend(:,:,slice);
   end
   NullSpaceX=null(X);
   Cmhalf=diag(1./sqrt(diag(contrasts*contrasts')));
   ContrastNullSpaceX=Cmhalf*contrasts*NullSpaceX;
   nonest=sum(abs(ContrastNullSpaceX)>tolerance,2);
   if sum(nonest)>0
      fprintf(['Error: the following contrasts are nonestimable in slice ' ...
            num2str(slice) ':']);
      RowsNonEstContrast=find(nonest>0)
      NonEstContrast=contrasts(RowsNonEstContrast,:)
      NullSpaceX
      ContrastNullSpaceX
      return
   end
end

indk1=((keep(2:n)-keep(1:n-1))==1);
k1=find(indk1)+1;

X_type=[ones(1,sum(num_hrf_bases==1))*1 ones(1,sum(num_hrf_bases==2))*2 ...
      ones(1,sum(num_hrf_bases==2))*3 ones(1,numtrends)*4 ];
find_X_is_u1=find(X_type==2);   
find_X_is_u2=find(X_type==3);         
find_X_is_u12=[find_X_is_u1 find_X_is_u2];         
find_X_is_mag=[find(X_type==1) find(X_type==2) find(X_type==4)];

find_contrast_is_mag=find(~contrast_is_delay);
find_response_is_mag=[find(num_hrf_bases==1) find(num_hrf_bases==2) ...
      numresponses+(1:numtrends)];
contr_mag=contrasts(find_contrast_is_mag,find_response_is_mag);

if any(contrast_is_delay)

   find_contrast_is_delay=find(contrast_is_delay);
   find_response_is_delay=find(num_hrf_bases==2);
   contr_delay=contrast(find_contrast_is_delay,find_response_is_delay);
   numcontr_delay=length(find_contrast_is_delay);
   contr_delay2=repmat(contr_delay,1,2)';
   contr_delay_is_1_col=zeros(numcontr_delay,1);
   find_delay_is_1_col=ones(numcontr_delay,1);
   for i=1:numcontr_delay
      pos=find(contr_delay(i,:)~=0);
      if length(pos)==1
         contr_delay_is_1_col(i)=1;
         find_delay_is_1_col(i)=pos;
      end
   end
   contr_delay_is_1_col;
   find_delay_is_1_col;
   
   % Fit a tangent function to the basis coefficients, W:
   cv=zeros(numresponses,1);
   dd0v=zeros(numresponses,1);
   for k=1:numresponses
      delta=X_cache.W(:,k,5);
      R=X_cache.W(:,k,basis2)./X_cache.W(:,k,basis1);
      ddelta=gradient(delta)./gradient(R);
      dd0=ddelta(delta==0);
      c=max(delta)/(pi/2);
      deltahat=atan(R/c*dd0)*c;
      for niter=1:5
         c=pinv(deltahat/c-cos(deltahat/c).^2.*R/c*dd0)*(delta-deltahat)+c;
         deltahat=atan(R/c*dd0)*c;
      end
      cv(k)=c;
      dd0v(k)=dd0;
   end
   C=cv(find_response_is_delay);
   C*pi/2;
   Dd0=dd0v(find_response_is_delay);
   
   if nargin < 8
      mag_t=ones(1,num_in_delay_contrast)*5;
   end
   if nargin < 9
      del_ef=zeros(1,num_in_delay_contrast);
   end
   
   T0=mag_t';
   delay=del_ef';
   T0_2=T0.^2;
   c0=T0_2./(T0_2+1);
   rhat=tan(delay./C).*C./Dd0./c0;
   drs=cos(delay./C).^2.*Dd0;
   
end

if numlags==1
   factor=1./sqrt(1-rho^2);
else
   Coradj_pix=rho;
   [Ainvt posdef]=chol(toeplitz([1 Coradj_pix]));
   nl=size(Ainvt,1);
   A=inv(Ainvt');
   B=ones(n-nl,1)*A(nl,:);
   Vmhalf=spdiags(B,1:nl,n-nl,n);
   %Vhalf=inv([A zeros(nl,n-nl); Vmhalf]);
   Anl=A(nl,nl);
   a=A(nl,(nl-1):-1:1)';
end

sd_ef=zeros(numcontrasts,numslices);
Y=zeros(size(contr_mag,1),numframes,numslices);
for slice=1:numslices
   if ~isempty(X_cache)
      X=[squeeze(X_cache.X(keep,num_hrf_bases==1,1,slice)) ... 
            squeeze(X_cache.X(keep,num_hrf_bases==2,basis1,slice)) ...
            squeeze(X_cache.X(keep,num_hrf_bases==2,basis2,slice)) ...
            Trend(:,:,slice)];
   else
      X=Trend(:,:,slice);
   end
   Xstar=X;
   
   if numlags==1
      Xstar(k1,:)=(X(k1,:)-rho*X(k1-1,:))*factor;
   else
      Xstar(1:nl,:)=A*X(1:nl,:);
      Xstar((nl+1):n,:)=Vmhalf*X;
   end
   pinvXstar=pinv(Xstar);
   V=pinvXstar*pinvXstar';
   if any(~contrast_is_delay)
      mag_sd=sqrt(diag(contr_mag*V(find_X_is_mag,find_X_is_mag)*contr_mag'));
      sd_ef(find_contrast_is_mag,slice)=mag_sd;
      Ystar=contr_mag*pinvXstar(find_X_is_mag,:);
      if numlags==1
         Y(:,keep,slice)=Ystar;
         for i=k1'
            Y(:,keep(i),slice)=Ystar(:,i)/factor+rho*Y(:,keep(i-1),slice);
         end
      else
         %Y(:,keep,slice)=Ystar*Vhalf';
         Y(:,keep(1:nl),slice)=Ystar(:,1:nl)*Ainvt;
         for i=(nl+1):n
            Y(:,keep(i),slice)=(Ystar(:,i)-Y(:,keep(i-(1:(nl-1))),slice)*a)/Anl;
         end
      end
   end
   if any(contrast_is_delay)
      sdef0=sqrt(diag(V(find_X_is_u1,find_X_is_u1)));
      gdot1=c0./(T0_2+1).*(1-T0_2).*rhat.*drs./T0./sdef0;
      gdot2=c0.*drs./T0./sdef0;
      gdot=[gdot1; gdot2]; 
      gdotc=kron(gdot,ones(1,numcontr_delay)).*contr_delay2;
      del_sd=sqrt(sum((V(find_X_is_u12,find_X_is_u12)*gdotc).*gdotc,1));
      sd_ef(find_contrast_is_delay,slice)=del_sd';
   end
end

return      
      