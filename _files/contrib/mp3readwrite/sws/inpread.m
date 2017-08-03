function [timeslice,F,M,nOscs,mat_error]=INPread(alldata); 

nOscs=sscanf(alldata(1,:), '%f', 1);
nts = (tall(alldata)-1) ./ (nOscs+1);

%Check to see if added or deleted an extra line
temp=fix(nts);
if temp~=nts
 mat_error=errordlg('Matrix dimension mismatch: Too many/few lines', '' );
 F=[];
 M=[];
 timeslice=[];
 return;
else
 mat_error=1;
end	

F = zeros(nOscs,nts);
M = zeros(nOscs,nts);

thisline=2;
c = 1;


timeslice=[];
while thisline < tall(alldata)

	[time] = sscanf(alldata(thisline,:), '%f',1);
	timeslice=[timeslice; time];
	thisline=thisline+1;
	for osc = 1:nOscs

		temp = sscanf(alldata(thisline,:), '%f', 2);
		F(osc,c) = temp(1,1);
		M(osc,c) = temp(2,1);
		
		
		thisline=thisline+1;
	end
	c = c + 1;
end

