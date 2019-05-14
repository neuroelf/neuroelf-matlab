function swxwrite(nOscs, timeslice, F, M);

% nOscs=2;
% 
% timeslice=10:10:50;
% timeslice=timeslice';
% 
% F=[0 1000 1000 1500 0; 0 2000 2500 3500 0];
% M=[0 0.235 0.512 0.150 0;0 0.418 0.400 0.319 0];

% sjf 2/00

%	constants
mask = '*.swx';				% SWX file mask


[file, pathName] = uiputfile(mask, 'Save data to SWX file');
	fileName = [pathName file];
	
if isnumeric(fileName)	%make sure didn't CANCEL
	return;
end;

%	open the file
[fid, msg] = fopen(fileName, 'wb');
if fid == -1							
	error([msg, ' (', fileName, ')']);
end;

%	write SWX data
fprintf(fid, '%d\n', nOscs);

for thistime=1:length(timeslice)
    fprintf(fid, '%7.0f\t', timeslice(thistime));
    for thisosc=1:nOscs
	    fprintf(fid, '%f \t', F(thisosc,thistime));
            fprintf(fid, '%f\t',M(thisosc,thistime));
    end
   fprintf(fid,'\n');
end
fclose(fid);

%	change file type/creator to Excel
filetype(fileName,'TEXT','XCEL');
