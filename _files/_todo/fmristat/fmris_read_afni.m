function [d]=fmris_read_afni(file,Z,T)

[err, Info]=BrikInfo(file);

d.dim 	= [Info.DATASET_DIMENSIONS(1:3) Info.DATASET_RANK(2)];

if nargin == 1 
   Z = 1:d.dim(3);
   T = 1:d.dim(4);
end

if (T(1)~=0)&(Z(1)~=0)
   Opt.Slices=Z;
   Opt.Frames=T;
   [err, d.data, Info, ErrMessage] = BrikLoad(file, Opt);
   d.calib		= [min(min(min(min(d.data)))) max(max(max(max(d.data))))];
end

d.vox 					= Info.DELTA;
d.vox_units				= '';
d.vox_offset			= 0;
d.precision				= '';
d.calib_units			= '';
d.origin 				= Info.ORIGIN;
d.descrip				= '';

if isfield(Info,'WORSLEY_DF')
   df=Info.WORSLEY_DF;
else
   df=[];
end
if ~isempty(df)
   d.df=df;
end

if isfield(Info,'WORSLEY_NCONJ')
   nconj=Info.WORSLEY_NCONJ;
else
   nconj=[];
end
if ~isempty(nconj)
   d.nconj=nconj;
end

if isfield(Info,'WORSLEY_FWHM')
   fwhm=Info.WORSLEY_FWHM;
else
   fwhm=[];
end
if ~isempty(fwhm)
   d.fwhm=fwhm;
end

return;
