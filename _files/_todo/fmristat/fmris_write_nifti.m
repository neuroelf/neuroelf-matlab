function [d]=fmris_write_nifti(d,Z,T)

if nargin<2
   Z=1;
   T=1;
end

if Z(1)==1 & T(1)==1
   fidout=fopen(d.file_name,'w');
   isniftiparent=0;
   if isfield(d,'parent_file')
      if exist(d.parent_file)
         pf=deblank(d.parent_file);
         ext=pf((length(pf)-2):length(pf));
         isniftiparent=(ext=='nii');
      end
   end
   if isniftiparent
      machineformats='nlbdgcas';
      for i=1:length(machineformats)
         fid=fopen(d.parent_file,'r',machineformats(i));
         sizeof_hdr=fread(fid,1,'int');
         if sizeof_hdr==348;
            %Must have found a formt that works, so break and continue using this format
            break
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
      intent_p =fread(fid,3,'float');
      intent_code =fread(fid,1,'short');
      datatype=fread(fid,1,'short');
      bitpix=fread(fid,1,'short');
      slice_start=fread(fid,1,'short');
      pixdim=fread(fid,8,'float');
      vox_offset=fread(fid,1,'float');
      scl_slope=fread(fid,1,'float');
      scl_inter=fread(fid,1,'float');
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
      magic =fread(fid,4,'char');
      fclose(fid);
   else
      data_type='          ';
      db_name='                  ';
      extents=0;
      session_error=0;
      regular='r';
      dim_info=' ';
      slice_start=0;
      pixdim=ones(1,8);
      slice_code =' ';
      xyzt_units =' ';
      cal_max=25500;
      cal_min=3;
      slice_duration=0;
      toffset=0;
      aux_file='                        ';
      qform_code =1;
      sform_code =0;
      quatern_b =0;
      quatern_c =1;
      quatern_d =0;
      qoffset_x =0;
      qoffset_y =0;
      qoffset_z =0;
      srow_x =[0 0 0 0];
      srow_y =[0 0 0 0];
      srow_z =[0 0 0 0];
      intent_name='                ';
   end
   
   sizeof_hdr=348;
   datatype=16;
   bitpix=32;
   vox_offset=352;
   scl_slope =0;
   scl_inter =0;
   slice_end=0;
   glmax=0;
   glmin=0;
   descrip=['FMRISTAT' repmat(' ',1,72)];
   magic =[double('n+1') 0]';
   
   dim=ones(1,8);      
   dim(1)=max(find(d.dim>1));
   dim((1:dim(1))+1)=d.dim(1:dim(1));
   if isfield(d,'vox')
      pixdim=ones(1,8);
      pixdim(1)=-1;
      pixdim((1:dim(1))+1)=d.vox(1:dim(1));
   end
   
   intent_p=zeros(1,2);
   if isfield(d,'df') 
      intent_p(1:length(d.df))=d.df; 
      intent_code=length(d.df)+2;
   else
      intent_code=0;
   end
   
   intent_q=zeros(1,2);
   if isfield(d,'fwhm'); 
      intent_q(1:length(d.fwhm))=d.fwhm; 
   end;
   intent_q=round(intent_q/100*2^16);
   intent_q=(intent_q.*(intent_q>=0)-2^16).*(intent_q<2^16)+2^16;

   descrip=['FMRISTAT' repmat(' ',1,72)];
   
   fwrite(fidout,sizeof_hdr,'int');
   fwrite(fidout,data_type,'char');
   fwrite(fidout,db_name,'char');
   fwrite(fidout,extents,'int');
   fwrite(fidout,session_error,'short');
   fwrite(fidout,regular,'char');
   fwrite(fidout,dim_info,'char');
   fwrite(fidout,dim,'short');
   fwrite(fidout,intent_p,'float');
   fwrite(fidout,intent_q,'short');
   fwrite(fidout,intent_code,'short');
   fwrite(fidout,datatype,'short');
   fwrite(fidout,bitpix,'short');
   fwrite(fidout,slice_start,'short');
   fwrite(fidout,pixdim,'float');
   fwrite(fidout,vox_offset,'float');
   fwrite(fidout,scl_slope ,'float');
   fwrite(fidout,scl_inter ,'float');
   fwrite(fidout,slice_end,'short');
   fwrite(fidout,slice_code,'char');
   fwrite(fidout,xyzt_units,'char');
   fwrite(fidout,cal_max,'float');
   fwrite(fidout,cal_min,'float');
   fwrite(fidout,slice_duration,'float');
   fwrite(fidout,toffset,'float');
   fwrite(fidout,glmax,'int');
   fwrite(fidout,glmin,'int');
   fwrite(fidout,descrip,'char');
   fwrite(fidout,aux_file,'char');
   fwrite(fidout,qform_code,'short');
   fwrite(fidout,sform_code,'short');
   fwrite(fidout,quatern_b,'float');
   fwrite(fidout,quatern_c,'float');
   fwrite(fidout,quatern_d,'float');
   fwrite(fidout,qoffset_x,'float');
   fwrite(fidout,qoffset_y,'float');
   fwrite(fidout,qoffset_z,'float');
   fwrite(fidout,srow_x,'float');
   fwrite(fidout,srow_y,'float');
   fwrite(fidout,srow_z,'float');
   fwrite(fidout,intent_name,'char');
   fwrite(fidout,magic,'char');
   
   fwrite(fidout,0,'float');
else
   fidout=fopen(d.file_name,'r+');
end

if nargin<2
   Z=1:d.dim(3);
   T=1:d.dim(4);
end

vox_offset=352;
if ~isfield(d,'precision'); d.precision='float32'; end;
if ~isfield(d,'byte'); d.byte=4; end;

if Z(1)==1 & T(1)==1 & ~(all(Z==1:d.dim(3)) & all(T==1:d.dim(4)))
   for t=1:d.dim(4)
      fwrite(fidout,zeros(1,prod(d.dim(1:3))),d.precision);
   end
end

if all(Z>0)&all(Z<=d.dim(3))&all(T>0)&all(T<=d.dim(4))
   for t=1:length(T)
      for z=1:length(Z)
         position=d.byte*((T(t)-1)*prod(d.dim(1:3))+(Z(z)-1)*prod(d.dim(1:2)))+vox_offset;
         status=fseek(fidout,position,'bof');
         if length(size(d.data))==4
            fwrite(fidout,d.data(:,:,z,t),d.precision);
         elseif length(T)==1 
            fwrite(fidout,d.data(:,:,z),d.precision);
         elseif length(Z)==1  
            fwrite(fidout,d.data(:,:,t),d.precision);
         else
            'Slices and/or frames do not match data'
         end
      end
   end
else
   'Slices and/or frames out of range:'
   Z
   T
end

fclose(fidout);

return;
