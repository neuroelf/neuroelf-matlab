function S=readstr(File);
%S=readstr(File)
%
% reads a file into MATLAB String-MATRIX (equal num of columns)
%
% by I. Bucher

 fid=fopen(File);
 S=[];
	    while 1
 	        line = fgetl(fid);
 		if ~ischar(line), break, end
 		S=str2mat(S,line);
 	    end
 	    fclose(fid);
