function [fdr_thresh, bon_thresh] = fdr_threshold( input_file, input_thresh, ... 
   mask_file, mask_thresh, df, fdr, flip, nconj, nvar)

%FALSE_DISCOVERY_RATE 
%
% [FDR_THRESH, BON_THRESH] = FDR_THRESHOLD( INPUT_FILE, [INPUT_THRESH, 
%   [MASK_FILE, [MASK_THRESH, [DF, [FDR [, FLIP [, NCONJ [, NVAR ]]]]]]]] )
%
% Calculates the threshold (FDR_THRESH) for a t or F image for 
% controlling the false discovery rate (FDR). The FDR is the expected 
% proportion of false positives among the voxels above FDR_THRESH. 
% This threshold is higher than that for controlling the expected proportion 
% of false positives in the whole search volume, and usually lower than the 
% Bonferroni threshold (printed out as BON_THRESH) which controls the 
% probability that *any* voxels are above the threshold. References: 
% Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery 
% rate: a practical and powerful approach to multiple testing. Journal of 
% the Royal Statistical Society, Series B, 57:289-300, and
% Genovese et al. (2002). Thresholding statistical maps in functional
% neuroimaging using the false discovery rate. Neuroimage, 15:722-786.
% NOTE: the procedure may set FDR_THRESH=Inf, i.e. no discoveries,
% particularly if INPUT_FILE is pure noise (no activation).
%
% INPUT_FILE is the input t or F statistic image.
%
% INPUT_THRESH is the lower limit for searching the t or F statistic image.
% If empty, it is the threshold corresponding to a P-value of FDR at a single 
% point. Default is []. 
%
% MASK_FILE is a mask file. If empty, it is ignored. Default is [].
%
% MASK_THRESH defines the search volume as the first frame of MASK_FILE 
% > MASK_THRESH. If MASK_THRESH is a vector [a b], a<=b, then mask 
% is a < MASK_FILE <= b. If empty (default), calls fmri_mask_thresh. 
%
% DF: If a scalar, then it is the df of the t statistic image;
% if DF >=1000 then DF is set it Inf, so that it calculates 
% thresholds for a Gaussian image (if DF is very large the t-dbn 
% is almost identical to the Gaussian dbn).
% If DF=[DF1, DF2] then these are the df's of the F statistic image.
% If DF2<0 then [DF1, -DF2] are the df's of the Hotelling's T^2 statistic.
% If DF2 >= 1000 then DF2 is set to Inf. Default is Inf.
%
% FDR is the desired false discovery rate. 
% If the first element is greater than 1, then they are treated as 
% thresholds and Q-values are returned. Default is 0.05.
%
% FLIP: INPUT_FILE is multiplied by FLIP before processing. For T statistic
% images, FLIP = -1 will look for negative peaks and clusters. Default is 1.
%
% NCONJ is the number of conjunctions. If NCONJ > 1, calculates Q-values and 
% thresholds for peaks of the minimum of NCONJ independent 
% SPM's - see Friston, K.J., Holmes, A.P., Price, C.J., Buchel, C.,
% Worsley, K.J. (1999). Multi-subject fMRI studies and conjunction analyses.
% NeuroImage, 10:385-396. Default is NCONJ = 1 (no conjunctions). 
%
% NVAR is the number of variables for multivariate equivalents of T and F 
% statistics, found by maximizing T^2 or F over all linear combinations of 
% variables, i.e. Hotelling's T^2 for DF1=1, Roy's maximum root for DF1>1. 
% Default is 1, i.e. univariate statistics.

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
if nargin<2;  input_thresh=[]; end
if nargin<3;  mask_file=[]; end
if nargin<4;  mask_thresh=[]; end
if nargin<5;  df=[];  end
if nargin<6;  fdr=[];  end
if nargin<7;  flip=[];  end
if nargin<8;  nconj=[];  end
if nargin<9;  nvar=[];  end

if isempty(df);  df=Inf;  end
if isempty(fdr);  fdr=0.05;  end
if isempty(flip);  flip=1;  end
if isempty(nconj);  nconj=1;  end
if isempty(nvar);  nvar=1;  end

if ~isempty(mask_file)
   if isempty(mask_thresh)
      mask_thresh=fmri_mask_thresh(mask_file);
   end
   mask_thresh1=mask_thresh(1);
    if length(mask_thresh)>=2
        mask_thresh2=mask_thresh(2);
    else
        mask_thresh2=Inf;
    end
end

if length(df)==1
    df=[df 0];
end
if size(df,1)==1
    df=[df; Inf Inf]
end
if size(df,2)==1
    df=[df [0; df(2,1)]]
end

% is_tstat=1 if it is a t statistic
is_tstat=(df(1,2)==0);

if is_tstat
   df1=1;
   df2=df(1,1);
else
   df1=df(1,1);
   df2=df(1,2);
end
if df2 >= 1000
   df2=Inf;
end

% Values of the F statistic or T^2/df1 based on squares of t values:
t=[100:-0.1:10.1 10:-0.01:0].^2;

% Find the upper tail probs of the F distribution by 
% cumulating the F density using the mid-point rule:
n=length(t);
n1=1:(n-1);
tt=(t(n1)+t(n1+1))/2;
if df2==Inf
   u=df1*tt;
   b=exp(-u/2-gammaln(df1/2)-df1/2*log(2)+(df1/2-1)*log(u));
else  
   u=df1*tt/df2;
   b=exp(-(df1+df2)/2*log(1+u)+(df1/2-1)*log(u)-betaln(df1/2,df2/2))*df1/df2;
end
D=0;
tau=zeros(D+nvar,n);
tau(1,:)=[0 -cumsum(b.*diff(t))];

% Find the EC densities:
y=df1*t;
for N=1:D+nvar-1
   s1=0;
   for i=0:(N-1)
      j=0:min(N-1-i,i);
      if df2==Inf
         s2=sum(exp(nchoosekln(df1-1+j-j,N-1-i-j) ...
            -gammaln(j+1)-gammaln(i-j+1)-j*log(2)));
      else
         s2=sum(exp(nchoosekln(df1-1+j-j,N-1-i-j) ...
            +nchoosekln((df1+df2-N)/2+j-1,j) ...
            +nchoosekln(df2-1+j-j,i-j)-i*log(df2)));
      end
      if s2>0
         s1=s1+(-1)^(N-1-i)*y.^(i+(df1-N)/2)*s2;
      end
   end
   if df2==Inf
      cons=-gammaln(df1/2)-N/2*log(2*pi)-(df1-2)/2*log(2)+gammaln(N);
      tau(N+1,:)=exp(cons-y/2).*s1;
   else   
      cons=-gammaln(df1/2)-N/2*log(2*pi)-(df1-2)/2*log(2)+gammaln(N) ...
         +gammaln((df1+df2-N)/2)-gammaln(df2/2)-(df1-N)/2*log(df2/2);
      tau(N+1,:)=exp(cons-(df1+df2-2)/2*log(1+y/df2)).*s1;
   end
end

% For multivariate statistics, add a sphere to the search region:
j=(nvar-1):-2:0;
invol_sphere=zeros(1,nvar);
invol_sphere(j+1)=exp(j*log(2)+j/2*log(pi) ...
   +gammaln((nvar+1)/2)-gammaln((nvar+1-j)/2)-gammaln(j+1));
rho=toeplitz([invol_sphere(1) zeros(1,D)]',[invol_sphere zeros(1,D)])*tau;

pt=rho(1,:);
if is_tstat
   t=[sqrt(t(n1)) -fliplr(sqrt(t))];
   pt=[pt(n1)/2 1-fliplr(pt)/2];
end
pt=pt.^nconj;

if isempty(input_thresh)
   input_thresh=stat_threshold(0,1,0,df,fdr,[],[],nconj,nvar);
end

d=fmris_read_image(input_file,0,0);
d.dim
numslices=d.dim(3);
numys=d.dim(2);
numxs=d.dim(1);
numpix=numxs*numys;
Steps=d.vox
voxel_volume=abs(prod(Steps))

Y=[];
num_voxels=0;
for slice=1:numslices
   d=fmris_read_image(input_file,slice,1);
   img=reshape(d.data,numpix,1)*flip;
   if ~isempty(mask_file)
      dm=fmris_read_image(mask_file,slice,1);
      mask=reshape(dm.data,numpix,1);
      mask= (mask>mask_thresh1 & mask<=mask_thresh2);
      Y=[Y; img(find(img>input_thresh & mask))];
      num_voxels=num_voxels+sum(mask);
   else
      Y=[Y; img(find(img>input_thresh))];
      num_voxels=num_voxels+length(img);
   end
end

% Bonferroni: 
search_volume=num_voxels*voxel_volume
pval_bon=num_voxels*pt;

% False discovery rate
P_val=interp1(t,pt,Y,'linear',0);
[P_sort, index]=sort(P_val);

nfdr=length(fdr);
fdr_thresh=zeros(1,nfdr);
if fdr(1)<1
   bon_thresh=minterp1(pval_bon,t,fdr);
   % False discovery rate threshold:
   q=P_sort./(1:length(P_sort))'*num_voxels;
   for i=1:nfdr
      r=max(find(q <= fdr(i)));
      if isempty(r)
         fdr_thresh(i)=Inf;
      else
         fdr_thresh(i)=Y(index(r));
      end
   end
   if nfdr<=1
      bon_thresh
      fdr_thresh
   end
else
   P0=interp1(t,pt,fdr);
   bon_thresh=P0*num_voxels;
   % Q-value:
   for i=1:nfdr
      K=find(P_sort>=P0(i));
      fdr_thresh(i)=min(P_sort(K)./K*num_voxels);
   end
   if nfdr<=1
      bon_thresh
      q_value=fdr_thresh
   end
end

return

function x=nchoosekln(n,k);
i=find(n>=0 & k>=0 & n>=k);
x=-Inf+n+k;
if ~isempty(i)
   x(i)=-log(n(i)+1)-betaln(k(i)+1,n(i)-k(i)+1);
end
return

function iy=minterp1(x,y,ix);
% interpolates only the monotonically increasing values of x at ix
n=length(x);
mx=x(1);
my=y(1);
xx=x(1);
for i=2:n
   if x(i)>xx
      xx=x(i);
      mx=[mx xx];
      my=[my y(i)];
   end
end
iy=interp1(mx,my,ix);
return













