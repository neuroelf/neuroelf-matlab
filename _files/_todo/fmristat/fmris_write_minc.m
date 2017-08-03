function	[d]=fmris_write_minc(d,Z,T);

if ~isfield(d,'dim') 
   d2=fmris_read_image(d.parent_file,0,0);
   d.dim=d2.dim;
end

if nargin==1
   Z=1:d.dim(3);
   T=1:d.dim(4);
end

if isfield(d,'precision')
   precision=d.precision;
else
   precision='float';
end

if T(1)==1 & Z(1)==1
   dim=d.dim;
   if dim(4)==1
      dim(4)=0;
   end
   if exist(d.file_name,'file')
      delete(d.file_name);
   end
   newh=newimage(d.file_name,fliplr(dim),d.parent_file,precision);
   if isfield(d,'df')
      miwriteatt(d.file_name,'df','df',d.df);
   end
   if isfield(d,'nconj')
      miwriteatt(d.file_name,'nconj','nconj',d.nconj);
   end
   if isfield(d,'fwhm')
      miwriteatt(d.file_name,'fwhm','fwhm',d.fwhm);
   end
else
   newh=openimage(d.file_name,'w');
end

d.data=squeeze(reshape(d.data,d.dim(1)*d.dim(2),length(Z),length(T)));

if d.dim(4)<=1
   putimages(newh,d.data,Z);
elseif length(T)==1|length(Z)==1
   putimages(newh,d.data,Z,T);
else
   for i=1:length(T)
      putimages(newh,squeeze(d.data(:,:,i)),Z,T(i));
   end
end

closeimage(newh);

return
