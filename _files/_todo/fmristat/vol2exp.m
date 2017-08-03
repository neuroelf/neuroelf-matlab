function vol2exp(d, sample, output_file_base, frames, expr)

%VOL2EXP writes input for IRIS Explorer latin and multin modules.
%
% VOL2EXP( INPUT [, SAMPLE [, OUTPUT_FILE_BASE [, FRAMES [, EXPR]]]] )
%
% INPUT is the input .mnc or .img file or structure from 
% fmris_read_image.
%
% SAMPLE is a 3 x 3 matrix whose rows are the start, step and stop
% indices, and whose columns are the x y and z directions. Default is 
% [1 1 d.dim(1); 1 1 d.dim(2); 1 1 d.dim(3)] i.e. all the data.
%
% OUPUT_FILE_BASE: writes INPUT to OUPUT_FILE_BASE.exp. 
% Default is INPUT minus extension if INPUT is file name.
%
% FRAMES: frames, default=1.
%
% EXPR: matlab expression applied to input in d.data, 
% e.g. 'd.data=sqrt(d.data);'. 

if nargin < 3
   output_file_base=[];
end
if nargin < 4
   frames=1;
end
if nargin < 5
   expr=[];
end

if ischar(d)
   input_file=d;
   d=fmris_read_image(input_file,0,0);
   isfile=1;
else
   input_file=d.file_name;
   isfile=0;
   d.dim=size(d.data);
   if length(d.dim)==2
      d.dim(3)=1;
   end
end

[base,ext]=fileparts2(input_file(1,:));
isminc=strcmp(lower(ext),'.mnc');
if isempty(output_file_base)
   output_file_base=base;
end

if nargin < 2
   sample=[];
end
if isempty(sample)
   sample=[1 1 d.dim(1); 1 1 d.dim(2); 1 1 d.dim(3)]
end

% For .mnc files, emma flips the slices if the voxel sizes are <0, 
% so we have to flip the voxel sizes and origin:

for i=1:2
   if d.vox(i)<0 & isminc
      d.origin(i)=d.origin(i)+d.vox(i).*(d.dim(i)-1);
      d.vox(i)=-d.vox(i);
   end
end

x=sample(1,1):sample(1,2):sample(1,3);
y=sample(2,1):sample(2,2):sample(2,3);
z=sample(3,1):sample(3,2):sample(3,3);

sizes=[length(x) length(y) length(z)]
m=[min(x) min(y) min(z); max(x) max(y) max(z)];
bbox=(m-1).*[d.vox; d.vox]+[d.origin; d.origin]

fid=fopen([output_file_base '.exp'],'w');
% first 3 are dimension sizes as 4-byte integers:
fwrite(fid,sizes,'int');
if length(frames)>1
    length(frames)
   %fwrite(fid,length(frames),'int');
end
% next 6 are bounding box as 4-byte floats: 
fwrite(fid,bbox,'float');
% rest are the data as 4-byte floats:
if isfile
   for k=z;
      d=fmris_read_image(input_file,k,frames);
      if ~isempty(expr)
         eval(expr);
      end
      if length(frames)==1
         fwrite(fid,d.data(x,y),'float');
      else
         fwrite(fid,permute(d.data(x,y,:),[3,1,2]),'float');
      end
   end
else
   if ~isempty(expr)
      eval(expr);
   end
   if length(frames)==1
      fwrite(fid,d.data(x,y,z),'float');
   else
      fwrite(fid,permute(d.data(x,y,z,frames),[4 1 2 3]),'float');
   end
end
fclose(fid);

return

