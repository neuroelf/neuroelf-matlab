function	[d]=fmris_write_afni(d,Z,T);

existinfo=0;
if ~isfield(d,'dim') & isfield(d,'parent_file')
   [path,name,ext]=fileparts(deblank(d.parent_file));
   [err, Infoparent] = BrikInfo([path '/' name '.HEAD']);
   d.dim=[Infoparent.DATASET_DIMENSIONS(1:3) Infoparent.DATASET_RANK(2)];
   existinfo=1;
end

if nargin==1 
   Z=1:d.dim(3);
   T=1:d.dim(4);
end

Opt.Prefix=d.file_name(1:(length(d.file_name)-5));
Opt.Slices=Z;
Opt.Frames=T;
Opt.NoCheck=1;

if Z(1)==1 & T(1)==1
   if ~existinfo
      [path,name,ext]=fileparts(deblank(d.parent_file));
      [err, Infoparent] = BrikInfo([path '/' name '.HEAD']);
   end
   [path,name,ext]=fileparts(deblank(d.file_name));
   Info.LABEL_1=name;
   Info.DATASET_NAME=['./' name];
   if isfield(d,'origin')
      Info.ORIGIN=d.origin;
   else
      Info.ORIGIN=Infoparent.ORIGIN;
   end
   if isfield(d,'vox')
      Info.DELTA=d.vox;
   else
      Info.DELTA=Infoparent.DELTA;
   end
   Info.SCENE_DATA=Infoparent.SCENE_DATA;
   Info.ORIENT_SPECIFIC=Infoparent.ORIENT_SPECIFIC;
   Info.TYPESTRING=Infoparent.TYPESTRING;
   Opt.NoCheck=0;
end

Info.DATASET_DIMENSIONS=[d.dim(1:3) 0 0];
Info.DATASET_RANK=[3 d.dim(4) 0 0 0 0 0 0];
Info.BRICK_TYPES=repmat(3,1,d.dim(4));
Info.TypeName='float';
Info.TypeBytes=4;
Info.BYTEORDER_STRING='MSB_FIRST';
Info.MachineFormat='ieee-be';

if isfield(d,'df')
   if ~isempty(d.df)
      Info.WORSLEY_DF=d.df;
   end
end

if isfield(d,'nconj')
   if ~isempty(d.nconj)
      Info.WORSLEY_NCONJ=d.nconj;
   end
end

if isfield(d,'fwhm')
   if ~isempty(d.fwhm)
      Info.WORSLEY_FWHM=d.fwhm;
   end
end

[err, ErrMessage, Info] = WriteBrik(d.data, Info, Opt);

return


