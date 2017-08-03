function [d]=fmris_read_nifti(file,Z,T)

d.file_name=file;
machineformats='nlbdgcas';
for i=1:length(machineformats)
  fid=fopen(file,'r',machineformats(i));
  if fid == -1;
     fprintf('Invalid file %s',file)
  end
  sizeof_hdr=fread(fid,1,'int');
  if sizeof_hdr==348;
     %Must have found a formt that works, so break and continue using this format
     break;
  else
    %Need to close if file is not native format
    %else if the file format is 's', then 7 stranded files are left orphaned and never closed.
    fclose(fid);
  end
end
data_type=fread(fid,10,'char');
db_name=fread(fid,18,'char');
extents=fread(fid,1,'int');
session_error=fread(fid,1,'short');
regular=char(fread(fid,1,'char')');
dim_info=char(fread(fid,1,'char')');
dim=fread(fid,8,'short');
intent_p =fread(fid,2,'float');
intent_q =fread(fid,2,'uint16');
intent_code =fread(fid,1,'short');
datatype=fread(fid,1,'short');
bitpix=fread(fid,1,'short');
slice_start=fread(fid,1,'short');
pixdim=fread(fid,8,'float');
vox_offset=fread(fid,1,'float');
scl_slope =fread(fid,1,'float');
scl_inter =fread(fid,1,'float');
slice_end=fread(fid,1,'short');
slice_code =char(fread(fid,1,'char')');
xyzt_units =char(fread(fid,1,'char')');
cal_max=fread(fid,1,'float');
cal_min=fread(fid,1,'float');
slice_duration=fread(fid,1,'float');
toffset=fread(fid,1,'float');
glmax=fread(fid,1,'int');
glmin=fread(fid,1,'int');
descrip=char(fread(fid,80,'char')');
aux_file=char(fread(fid,24,'char')');
qform_code =fread(fid,1,'short');
sform_code =fread(fid,1,'short');
quatern_b =fread(fid,1,'float');
quatern_c =fread(fid,1,'float');
quatern_d =fread(fid,1,'float');
qoffset_x =fread(fid,1,'float');
qoffset_y =fread(fid,1,'float');
qoffset_z =fread(fid,1,'float');
srow_x =fread(fid,4,'float');
srow_y =fread(fid,4,'float');
srow_z =fread(fid,4,'float');
intent_name=char(fread(fid,16,'char')');
magic =char(fread(fid,4,'char')');

d.machineformat=machineformats(i);
d.dim=ones(1,4);
if dim(1)>4
   fclose(fid);
   return
end
d.dim(1:dim(1))=dim((1:dim(1))+1);
d.vox=zeros(1,3);
d.vox(1:dim(1))=pixdim((1:dim(1))+1);
%Attempt to fill out more information for a complete nifti description
  d.vox_offset  = vox_offset;
  d.scale       = scl_slope;
  d.intercept   = scl_inter;
  d.global      = [glmin glmax];
  d.calib       = [cal_min cal_max];
  if qform_code>0;
    d.origin    = [qoffset_x qoffset_y qoffset_z];
  else
    d.origin    = [srow_x(4) srow_y(4) srow_z(4)];
  end
  d.descrip     = descrip;
d.parent_file=file;

if intent_p(1)>0; d.df(1)=intent_p(1); end;
if intent_p(2)>0; d.df(2)=intent_p(2); end;

intent_q=intent_q/2^16*100;
if intent_q(1)>0; d.fwhm(1)=intent_q(1); end;
if intent_q(2)>0; d.fwhm(2)=intent_q(2); end;

if nargin<2
   Z=1:d.dim(3);
   T=1:d.dim(4);
end

if Z(1)==0 & T(1)==0
   fclose(fid);
   return
end

types=lower(['UINT8  ';'INT16  ';'INT32  ';'FLOAT32';'FLOAT64']);
codes=[2; 4; 8; 16;  64];
d.precision=[];
for i=1:5
   if datatype==codes(i)
      d.precision=deblank(types(i,:));
   end
end
if isempty(d.precision)
   unknown_datatype=datatype
   fclose(fid);
   return
end

d.byte=bitpix/8;
if all(Z>0)&all(Z<=d.dim(3))&all(T>0)&all(T<=d.dim(4))
   d.data=zeros(d.dim(1),d.dim(2),length(Z),length(T));
   for t=1:length(T)
      for z=1:length(Z)
         position=d.byte*((T(t)-1)*prod(d.dim(1:3))+(Z(z)-1)*prod(d.dim(1:2)))+vox_offset;
         status=fseek(fid,position,'bof');
         d.data(:,:,z,t)=fread(fid,[d.dim(1) d.dim(2)],d.precision);
      end
   end
  fclose(fid);
else
   Z
   T
   fclose(fid);
   return
end

if scl_slope~=0
   d.data=d.data*scl_slope+scl_inter;
end


return;
