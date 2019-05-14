function swiwrite(nOscs, timeslice, F, M);

% nOscs=2;
% 
% timeslice=10:10:50;
% timeslice=timeslice';
% 
% F=[0 1000 1000 1500 0; 0 2000 2500 3500 0];
% M=[0 0.235 0.512 0.150 0;0 0.418 0.400 0.319 0];

% sjf  9/98

%	constants
mask = '*.swi';				% SWI file mask

%	get filename 
[file, pathName] = uiputfile(mask, 'Save data to SWI file');
	fileName = [pathName file];
	
if isnumeric(fileName)	%in case of CANCEL
	return;
end;

%	open the file
[fid, msg] = fopen(fileName, 'wb');
if fid == -1							
	error([msg, ' (', fileName, ')']);
end;

%	write SWI data
fprintf(fid, '%d\n', nOscs);

for thistime=1:length(timeslice)
    fprintf(fid, '%f\n', timeslice(thistime));
    for thisosc=1:nOscs
	 fprintf(fid, '%f \t', F(thisosc,thistime));
	 fprintf(fid, '%f\n',M(thisosc,thistime));
    end
end
fclose(fid);
