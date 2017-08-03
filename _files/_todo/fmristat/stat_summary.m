function [summary_clusters, summary_peaks, thresh_peaks]= ...
   stat_summary(input_file, fwhm, ...
   df, mask_file, mask_thresh, input_thresh, flip, nconj, nvar);

%STAT_SUMMARY produces SPM-style summary analyses of T or F statistic images
%
% [SUMMARY_CLUSTERS SUMMARY_PEAKS THRESH_PEAKS] = 
%     STAT_SUMMARY( INPUT_FILE [, FWHM [, DF [, MASK_FILE [, 
%     MASK_THRESH [, INPUT_THRESH [, FLIP [, NCONJ [, NVAR ]]]]]]]])
%
% Produces an SPM-style glass brain and summary analysis of a T or F statistic
% image. If an FWHM image file is provided, P-values for local maxima and 
% cluster sizes are based on non-isotropic random field theory, or Bonferroni,
% or discrete local maxima (DLM), whichever is smaller.
%
% The random field theory is based on the assumption that the search region 
% is a sphere (in isotropic space), which is a very tight lower bound for 
% any non-spherical region, unless you supply all the resels. 
%
% The DLM P-value is an upper bound, like Bonferroni, but it is more accurate
% over a greater range of FWHM. 
%
% It also produces a volume of clusters labelled by their index
% (as in SUMMARY_CLUSTERS) with '_cluster' before the
% extension, handy for identifying the clusters in 'register'.
%
% It also produces a volume of corrected P-values with '_Pval' before the 
% extension. Values are multiplied by -1 (so bigger values are more 
% significant), with -1.1 outside the mask (see MASK_FILE below).
%
% INPUT_FILE: Finds its local maxima and clusters above INPUT_THRESH.
% Clusters are voxels connected in any of the 2*D directions. A local maximum 
% is a voxel which is greater than or equal to all its 2*D neighbours,
% and strictly greater than at least one of them. If empty, just prints out
% average FWHM (see below).
%
% FWHM is the fwhm in mm of a smoothing kernel applied to the data, either 
% as a fwhm file from fmrilm, multistat or glim_image, or as a scalar. 
% If FWHM is a vector, these are treated as resels of the mask.
% If empty (default), looks for a fwhm attribute in INPUT_FILE, and if
% it can't find it, sets FWHM=0.
%
% DF=[DF1 DF2; DFW1 DFW2] is a 2 x 2 matrix of degrees of freedom.
% If DF2 is 0, then DF1 is the df of the T statistic image.
% If DF1=Inf then it calculates thresholds for the Gaussian image. 
% If DF2>0 then DF1 and DF2 are the df's of the F statistic image.
% DFW1 and DFW2 are the numerator and denominator df of the FWHM image. 
% If DF=[DF1 DF2] (row) then DFW1=DFW2=Inf, i.e. FWHM is fixed.
% If DF=[DF1; DFW1] (column) then DF2=0 and DFW2=DFW1, i.e. t statistic. 
% If DF=DF1 (scalar) then DF2=0 and DFW1=DFW2=Inf, i.e. t stat, fixed FWHM.
% If any component of DF >= 1000 then it is set to Inf. Default is Inf. 
% If FWHM is estimated by FMRILM, set DFW1=DFW2=DF outputted by FMRILM;
% if FWHM is estimated by MULTISAT, set DFW1=DF_RESID, DFW2=DF outputted 
% by MULTISTAT. If empty (default), finds DF from INPUT_FILE, and if
% it can't find it, sets DF=Inf.
%
% MASK_FILE is a mask file. If empty (default), it is ignored.
%
% MASK_THRESH defines the search volume as the first frame of MASK_FILE 
% > MASK_THRESH. If MASK_THRESH is a vector [a b], a<=b, then mask 
% is a < MASK_FILE <= b. If empty (default), calls fmri_mask_thresh. 
%
% INPUT_THRESH: If <= 1 then the second element is taken as a probability and 
% the threshold is chosen so that the uncorrected P-value is this probability.
% If INPUT_THRESH is a scalar, the second element is set equal to the first.
% The default is 0.001, i.e. the threshold satisfies P=0.001 (uncorrected).
%
% FLIP: INPUT_FILE is multiplied by FLIP before processing. For T statistic
% images, FLIP = -1 will look for negative peaks and clusters. Default is 1.
%
% NCONJ: use is now discouraged. If NCONJ > 1, calculates P-values
% for peaks (but not clusters) of the minimum of NCONJ independent 
% SPM's - see Friston, K.J., Holmes, A.P., Price, C.J., Buchel, C.,
% Worsley, K.J. (1999). Multi-subject fMRI studies and conjunction analyses.
% NeuroImage, 10:385-396. If empty (default), finds NCONJ from INPUT_FILE, 
% and if it can't find it, sets NCONJ=1 (no conjunctions).
%
% NVAR is the number of variables for multivariate equivalents of T and F 
% statistics, found by maximizing T^2 or F over all linear combinations of 
% variables, i.e. Hotelling's T^2 for DF1=1, Roy's maximum root for DF1>1. 
% Default is 1, i.e. univariate statistics.
%
% SUMMARY_CLUSTERS is a matrix with 6 columns:
% Col 1: index of cluster, in descending order of cluster resels. 
% Col 2: volume of cluster in mm^3.
% Col 3: resels of cluster.
% Col 4: P-value of cluster extent.
% Col 5: P-value if the cluster was chosen in advance, e.g. nearest to an ROI.
%
% SUMMARY_PEAKS is a matrix with 11 columns. 
% Col 1: index of cluster. 
% Col 2: values of local maxima, sorted in descending order withihn cluster.
% Col 3: P-value of local maxima.
% Col 4: P-value if the peak was chosen in advance, e.g. nearest to an ROI.
% Col 5: Q-value or false discovery rate ~ probability that voxel is not signal.
% Cols 6-8: i,j,k coords of local maxima in voxels, starting at 0, as in 'register'.
% Cols 9-11: x,y,z coords of local maxima in world coordinates (mm).
%
% THRESH_PEAKS is the P=0.05 threshold for peaks.

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

summary_clusters=[];
summary_peaks=[];
peaks_thresh=[];
      
if nargin<2; fwhm=[]; end
if nargin<3; df=[]; end
if nargin<4; mask_file=[]; end
if nargin<5; mask_thresh=[]; end
if nargin<6; input_thresh=[]; end
if nargin<7; flip=[]; end
if nargin<8;  nconj=[];  end
if nargin<9;  nvar=[];  end

if isempty(flip)
   flip=1
end
if isempty(nvar)
   nvar=1
end

if ~isempty(input_file)
   d=fmris_read_image(input_file,0,0);
   if isempty(df) 
      if isfield(d,'df') 
         df=d.df
      else
         df=Inf
      end
   end
   if isempty(nconj) 
      if isfield(d,'nconj') 
         nconj=d.nconj
      else
         nconj=1
      end
   end
   if isempty(fwhm) 
      if isfield(d,'fwhm')
         fwhm=d.fwhm
      else
         fwhm=0
      end
   end
end

if ischar(fwhm)
   d=fmris_read_image(fwhm,0,0);
   if isfield(d,'df')
      df(2,1:2)=d.df
   end
end

if isempty(input_thresh)
   input_thresh=0.001;
end
if length(input_thresh)==1
   input_thresh=[input_thresh input_thresh];
end
if input_thresh(1) > 1
   input_thresh=input_thresh(1)
else
   x=stat_threshold(1,1,0,df,input_thresh,0.001,0.05,nconj,nvar);
   input_thresh=x(2)
end

if ~isempty(mask_file) & isempty(mask_thresh)
   mask_thresh=fmri_mask_thresh(mask_file);
end

if ~isempty(input_file)
   [base,ext]=fileparts2(input_file);
   if strcmp(ext,'.BRIK')
      base1=base(1:(length(base)-5));
      base2=base(length(base)-5+(1:5));
      file_w=[base1 '_cluster' ext];
      file_h=[base1 '_cluster' base2 '.HEAD'];
      cluster_file=[base1 '_cluster' base2 ext];
      if exist(file_h);  delete(file_h); end
   else
      cluster_file=[base '_cluster' ext];
      file_w=cluster_file;
   end
   if exist(cluster_file);  delete(cluster_file); end
   
   % Test that you can write files:
   %fid=fopen(cluster_file,'w');
   fid=-1
   if fid==-1
      canwritefile=0;
      ['Can''t write files, so no 3D cluster blobs or Pval file']
      lm=locmax(input_file,input_thresh,mask_file,mask_thresh,fwhm,flip);
   else
      canwritefile=1;
      fclose(fid);
      delete(cluster_file);
      lm=locmax(input_file,input_thresh,mask_file,mask_thresh,fwhm,flip,file_w);
   end

   if isempty(lm)
      return
   end
   
   [search_volume, num_voxels]= ...
      glass_brain(input_file,input_thresh,mask_file,mask_thresh,flip);
   colormap(spectral);
   
   d=fmris_read_image(input_file,0,0);
else
   d=fmrs_read_image(fwhm_file,0,0);
end

if ~isempty(mask_file)
   mask_thresh1=mask_thresh(1);
   if length(mask_thresh)>=2
      mask_thresh2=mask_thresh(2);
   else
      mask_thresh2=Inf;
   end
   d=fmris_read_image(mask_file,1:d.dim(3),1);
   mask=d.data;
   mask= (mask>mask_thresh1 & mask<=mask_thresh2);
else
   mask=ones(d.dim(1),d.dim(2),d.dim(3));
end

if ~isempty(input_file)
   d=fmris_read_image(input_file);
   d.data=d.data*flip;
   ivox=find(mask & d.data>input_thresh);
   z=d.data(ivox);
end
nlm1=size(lm,1)+1;

if isstr(fwhm)
   d=fmris_read_image(fwhm,0,0);
   if d.dim(4)>=2
      d=fmris_read_image(fwhm,1:d.dim(3),2);
   else
      d=fmris_read_image(fwhm,1:d.dim(3),1);
      d.data=abs(prod(d.vox(1:3)))./(d.data+(d.data<=0)).^3.*(d.data>0);
   end
   search_resels=sum(sum(sum(mask.*d.data)))
   num_voxels=sum(sum(sum(mask)))
   search_volume=num_voxels*abs(prod(d.vox(1:3)))
   average_fwhm=(search_volume/search_resels)^(1/3)
   if isempty(input_file)
      return
   end
   
   [p_peak, p_cluster, p_peak1, p_cluster1]= ...
      stat_threshold(search_resels, num_voxels, 1, ... 
      df, [10; lm(:,1); min(z); z], input_thresh, [10; lm(:,7)], nconj, nvar);
   p_vol=p_peak(1+nlm1+(1:length(z)));
   
   if d.dim(4)>2
      fwhms=fwhm;
   else
      fwhms=average_fwhm./d.vox(1:3);
   end

   pdlm=dlm([10; lm(:,1); min(z)],fwhms,df,mask_file,mask_thresh,nconj,nvar)';
   [lmz,i,j]=unique([lm(:,1); min(z)]);
   logpratio=log(pdlm(i+1)./p_peak(i+1));
   if length(lmz)>=2
      logpratioz=interp1(lmz, logpratio, z);
      pdlmz=exp(logpratioz).*p_vol;
      p_vol=min(p_vol,pdlmz);
   end
   p_peak=min(p_peak(1:nlm1),pdlm(1:nlm1));
elseif length(fwhm)==1
   [p_peak, p_cluster, p_peak1, p_cluster1]= ...
      stat_threshold(search_volume, num_voxels, fwhm, ... 
      df, [10; lm(:,1); z], input_thresh, [10; lm(:,6)], nconj, nvar);
   p_vol=p_peak(nlm1+(1:length(z)));
else
   [p_peak, p_cluster, p_peak1, p_cluster1]= ...
      stat_threshold(fwhm, num_voxels, 1, ... 
      df, [10; lm(:,1); z], input_thresh, [10; lm(:,6)], nconj, nvar);   
   p_vol=p_peak(nlm1+(1:length(z)));
end
p_peak=p_peak(2:nlm1);
p_cluster=p_cluster(2:nlm1);
p_peak1=p_peak1(2:nlm1);
p_cluster1=p_cluster1(2:nlm1);

if canwritefile
   if strcmp(ext,'.BRIK')
      file_w=[base1 '_Pval' ext];
      file_h=[base1 '_Pval' base2 '.HEAD'];
      pval_file=[base1 '_Pval' base2 ext];
      if exist(file_h);  delete(file_h); end
   else
      pval_file=[base '_Pval' ext];
   end
   if exist(pval_file);  delete(pval_file); end
   out.dim=d.dim;
   out.dim(4)=1;
   out.data=zeros(d.dim(1:3))-1.1;
   out.data(find(mask))=-1;
   out.data(ivox)=-min(p_vol,1);
   out.parent_file=input_file;
   out.file_name=pval_file
   fmris_write_image(out,1:out.dim(3),1);
end

q_value = fdr_threshold( input_file, input_thresh, ... 
   mask_file, mask_thresh, df, lm(:,1), flip, nconj, nvar);

if isnan(p_cluster(1))
   summary_clusters=unique(lm(:,5:7),'rows');
   n=size(summary_clusters,1);
   ['clus    vol  resel']
   [repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,1))) repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,2))) repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,3)*100)/100) ]
else
   summary_clusters=unique([lm(:,5:7) p_cluster p_cluster1],'rows');
   n=size(summary_clusters,1);
   ['clus    vol  resel  Pval   (one)']
   [repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,1))) repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,2))) repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,3)*100)/100) repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,4)*1000)/1000) repmat(' (',n,1)  ... 
         num2str(round(summary_clusters(:,5)*1000)/1000) repmat(')  ',n,1)]
   n_clus=sum(summary_clusters(:,4)<=0.05);
   if n_clus>0 & canwritefile 
      p=get(gca,'Position');
      axes('position',[p(1)+p(3)*7/6+0.05 p(2) 0.9-p(1)-p(3)*7/6 p(4)]);
      if d.dim(3)>1
          blob_brain(input_file,input_thresh*flip,cluster_file,[0.5 n_clus+0.5]);
      else
          view_slices(cluster_file,mask_file,mask_thresh,0,1,[0.5 n_clus+0.5]);
      end
      title('Clusters, P<0.05, index=');
   end
end

isminc=strcmp(lower(ext),'.mnc');
if isminc
   h=openimage(input_file);
   coord=voxeltoworld(h,lm(:,2:4)','xyzorder zerobase noflip')';
   closeimage(h);
else
   d=fmris_read_image(input_file,0,0);
   coord=(lm(:,2:4)+1-ones(size(lm,1),1)*d.origin).*(ones(size(lm,1),1)*d.vox(1:3));
end
summary_peaks=[lm(:,5) lm(:,1) p_peak p_peak1 q_value' lm(:,2:4) coord ];
summary_peaks=flipud(sortrows(summary_peaks,2));
n=size(lm,1);
['clus   peak   Pval   (one)   Qval    (i   j   k)  (   x      y      z )']
[repmat('  ',n,1) ...
   num2str(round(summary_peaks(:,1))) repmat('  ',n,1) ...
   num2str(round(summary_peaks(:,2)*100)/100) repmat('  ',n,1)  ...  
   num2str(round(summary_peaks(:,3)*1000)/1000) repmat(' (',n,1)  ...  
   num2str(round(summary_peaks(:,4)*1000)/1000) repmat(')  ',n,1)  ... 
   num2str(round(summary_peaks(:,5)*1000)/1000) repmat('  (',n,1)  ... 
   num2str(round(summary_peaks(:,6:8))) repmat(')  (',n,1)  ... 
   num2str(round(summary_peaks(:,9)*10)/10)  repmat('  ',n,1) ...
   num2str(round(summary_peaks(:,10)*10)/10) repmat('  ',n,1) ...
   num2str(round(summary_peaks(:,11)*10)/10) repmat(')',n,1)]

if min(p_peak)<=0.05 & max(p_peak)>=0.05
   thresh_peaks=interp1(log(p_peak),lm(:,1),log(0.05),'spline')
else
   thresh_peaks=nan
end

return

