function [summary_cluster, summary_peaks]=conjunction(input_files, output_file_base);

%CONJUNCTION (minimum) of images
%
%    CONJUNCTION( INPUT_FILES [, OUTPUT_FILE_BASE [, INPUT_THRESH ]])
%
% INPUT_FILES is a matrix of image files, one per row, padded with blanks.
% Minimum (conjunction) is written to OUTPUT_FILE_BASE_conj.ext, which
% if empty (default), is the first row of INPUT_FILES minus ext.
% Summary statistics are produced as in stat_summary.m, with a default
% input_thresh of -0.99 (no p-values for clusters are available). 

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

[base,ext]=fileparts2(input_files(1,:));
if nargin<2 
   output_file_base=base;
end
if nargin<3
   input_thresh=-0.99;
end

numfiles=size(input_files,1);
d=fmris_read_image(input_files(1,:),0,0);
for ifile=1:numfiles
   d=fmris_read_image(input_files(ifile,:),1:d.dim(3),1);
   if ifile==1
      out.data=d.data;
   else
      out.data=min(out.data,d.data);
   end
end

out.dim=d.dim;
out.dim(4)=1;
out.parent_file=deblank(input_files(1,:));
out.file_name=[output_file_base '_conj' ext];
fmris_write_image(out,1:out.dim(3),1);

   [base,ext]=fileparts2(out.file_name);
   if strcmp(ext,'.BRIK')
      base1=base(1:(length(base)-5));
      base2=base(length(base)-5+(1:5));
      file_w=[base1 '_cluster' ext];
      file_h=[base1 '_cluster' base2 '.HEAD'];
      file_b=[base1 '_cluster' base2 ext];
      if exist(file_h);  delete(file_h); end
      if exist(file_b);  delete(file_b); end
      lm=locmax(out.file_name,input_thresh,out.file_name,-1.05,1,1,file_w);
      cluster_file=file_b;
   else
      cluster_file=[base '_cluster' ext];
      if exist(cluster_file);  delete(cluster_file); end
      lm=locmax(out.file_name,input_thresh,out.file_name,-1.05,1,1,cluster_file);
   end
   if isempty(lm)
      summary_clusters=[];
      summary_peaks=[];
      return
   end
   glass_brain(out.file_name,input_thresh,out.file_name,-1.05,1);
   colormap(spectral);
   
   summary_clusters=unique(lm(:,5:6),'rows');
   n=size(summary_clusters,1);
   ['clus    vol']
   [repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,1))) repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,2))) ]

isminc=strcmp(lower(ext),'.mnc');
if isminc
   h=openimage(out.file_name);
   coord=voxeltoworld(h,lm(:,2:4)','xyzorder zerobase noflip')';
   closeimage(h);
else
   d=fmris_read_image(out.file_name,0,0);
   coord=lm(:,2:4).*(ones(size(lm,1),1)*d.vox)+ones(size(lm,1),1)*d.origin;
end
summary_peaks=[lm(:,5) -lm(:,1) lm(:,2:4) coord ];
n=size(lm,1);
['clus peak=Pval (i   j   k)  (   x      y      z )']
[repmat('  ',n,1) ...
   num2str(round(summary_peaks(:,1))) repmat('  ',n,1) ...
   num2str(round(summary_peaks(:,2)*1000)/1000) repmat('  (',n,1)  ...  
   num2str(round(summary_peaks(:,3:5))) repmat(')  (',n,1)  ... 
   num2str(round(summary_peaks(:,6)*10)/10) repmat('  ',n,1) ...
   num2str(round(summary_peaks(:,7)*10)/10) repmat('  ',n,1) ...
   num2str(round(summary_peaks(:,8)*10)/10) repmat(')',n,1)]

return
