function math_vol(varargin);

%MATH_VOL evaluates simple expressions on volumes
%
% MATH_VOL( [INFILE1, ..., INFILEN,] EXPR ...
%   [, OUTFILE1, ..., OUTFILEN] [, START [, FINISH]])
%
% Evaluates EXPR on INFILEs to produce OUTFILEs. In EXPR,
% INFILEi, frame j is a{i,j}; OUTFILEi, frame j is b{i,j}.
% If there is only one file then i can be dropped; 
% If there is only one frame then j can be dropped. 
%
% START and FINISH are other matlab expressions evaluated 
% before and after looping over voxels. 
%
% EXPR and START must contain ';' somewhere.
%
% EXPR, START and FINISH can contain any matlab expressions. 
% x, y and z are pre-defined as the x, y and z coordinates 
% in voxels (1 based). 
% Note a,b,x,y are actually 2D slice matrices, z is a scalar. 
% Volume dimensions are in d.dim(1:3), sizes in d.vox(1:3),
% origins on d.origin(1:3). They should not alter variables 
% named x,y,z,d,exprid,exprids,la,expr.
%
% E.g. to calculate a t-statistic from an ef and sd file:
%
% file1='c:/keith/results/ha_100326_hmw_mag_ef.mnc';
% file2='c:/keith/results/ha_100326_hmw_mag_sd.mnc';
% file3='c:/keith/results/ha_100326_hmw_mag_t_test.mnc';
% expr='b=a{1}./(a{2}+(a{2}<=0)).*(a{2}>0);'
% mathvol(file1,file2,expr,file3);
%
% To find upper and lower error bars:
%
% file4='c:/keith/results/ha_100326_hmw_mag_ef_up.mnc';
% file5='c:/keith/results/ha_100326_hmw_mag_ef_dn.mnc';
% expr='b{1}=a{1}+a{2}; b{2}=a{1}-a{2};'
% mathvol(file1,file2,expr,file4,file5);
%
% or to put them as 2 frames in 1 file:
%
% file6='c:/keith/results/ha_100326_hmw_mag_ef_up_dn.mnc';
% mathvol(file1,file2,expr,file6);
%
% To count the number of ef voxels above 0.5:
%
% mathvol(file1, 'n=n+sum(a(:)>0.5);', 'n=0;', 'n');

exprids=[];
for i=1:nargin
   if isstr(varargin{i})
      if findstr(varargin{i},';')
         exprids=[exprids i];
      end
   end
end
exprid=exprids(1);
expr=varargin{exprid};
nargins=nargin;
if length(exprids)>1
   nargins=exprids(2)-1;
   eval(varargin{exprids(2)});
end

first_file=varargin{1+(exprid==1)};
d=fmris_read_image(first_file,0,0);

for i=1:(exprid-1)
   d=fmris_read_image(varargin{i},0,0);
   la(i)=d.dim(4);
end
for i=(exprid+1):nargins
   if exist(varargin{i},'file')
      delete(varargin{i});
   end
end

[x, y]=ndgrid(1:d.dim(1),1:d.dim(2));
for z=1:d.dim(3)
   z
   for i=1:(exprid-1)
      for j=1:la(i)
         d=fmris_read_image(varargin{i},z,j);
         if exprid==2 
            if la(i)==1
               a=d.data;
            else
               a{j}=d.data;
            end
         else
            if la(i)==1
               a{i}=d.data;
            else
               a{i,j}=d.data;
            end
         end
      end
   end
   
   eval(expr);
   
   if exist('b','var')
      if iscell(b)
         isemptyb=cellfun('isempty',b)';
      else
         isemptyb=0;
      end
      if (nargins-exprid)==1
         isemptyb=isemptyb';
      end
      
      for i=1:size(isemptyb,1)
         fisb=find(isemptyb(i,:)==0);
         lb=max(fisb);
         for j=fisb
            if size(isemptyb,1)==1 
               if lb==1
                  d.data=b;
               else
                  d.data=b{j};
               end
            else
               if lb==1
                  d.data=b{i};
               else
                  d.data=b{i,j};
               end
            end
            d.file_name=varargin{exprid+i};
            d.parent_file=first_file;
            d.dim(4)=lb;
            fmris_write_image(d,z,j);
         end
      end
   end
end

if length(exprids)==2 
   if nargin>exprids(2)
      eval(varargin{nargin});
   end
end

return






       
      
      
      

