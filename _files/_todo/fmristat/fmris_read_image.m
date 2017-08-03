function [d] = fmris_read_image(file,Z,T)

% (c) John Aston & Roger Gunn

file = deblank(file);
if ~(file(1)=='/' | file(2)==':')
   if file(1:2)=='./'
      file=[pwd file(3:length(file))];
   else
      file=[pwd '/' file];
   end
end
   
[path,name,ext]=fileparts(file);
if strcmp(ext,'.gz')
   if isunix
      unix(['gunzip ' file]);
   else
      ['Can''t gunzip on non-unix system']
      return
   end
   [path,name,ext]=fileparts(name)
end 

switch lower(ext)
case '.mnc'
   if nargin == 3
      d = fmris_read_minc(file,Z,T);
   else
      d = fmris_read_minc(file);
   end
case '.img'
   if nargin == 3
      d = fmris_read_analyze(file,Z,T);
   else
      d = fmris_read_analyze(file);
   end
   if isfield(d,'data') d.data(isnan(d.data))=0; end;
case '.brik'
   if nargin == 3
      d = fmris_read_afni(file,Z,T);
   else
      d = fmris_read_afni(file);
   end
case '.nii'
   if nargin == 3
      d = fmris_read_nifti(file,Z,T);
   else
      d = fmris_read_nifti(file);
   end
otherwise
   ['Unknown file extension']
end

d.parent_file=file;

return