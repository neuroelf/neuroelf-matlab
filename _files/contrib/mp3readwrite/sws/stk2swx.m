function stk2swx;

%Convert Scicon STK file to SWX format for use in MATLAB SWS.
%To use, simply type stk2swx; at the MATLAB command line.

%get filename interactively
[file, pathName] = uigetfile('*', 'Select Scicon format STK file');
FILENAME = [pathName file];
if isnumeric(FILENAME); return; end;		%make sure didn't cancel
[p,n,ext]=fileparts(FILENAME);
ext=lower(ext(2:end));	

%change working directory
cd(pathName);

if strcmp(ext,'stk')|strcmp(ext,'STK')
	%open file for processing
	fid = fopen(FILENAME, 'r');
	if (fid == -1)
	  fprintf(1, 'readswi: unable to read %s\n', FILENAME);
	end
	
	firstline=str2num(fgetl(fid));
	secondline=str2num(fgetl(fid));
	data=fread(fid)';
	fclose(fid);
	
	%parameters
	duration=firstline(1);
	samprate=firstline(2);
	timeinterval=firstline(3);
	synthtype=secondline(1);	%0=parallel(serial); 1=cascade/parallel(parallel)
	TIMESLICES=(timeinterval:timeinterval:duration);
	
	allDone = 0;
	c = 1;
	while ~allDone,
	  if c> length(data-1), 		% if hit eof before finding f0 then break out 
	    error('File is ill-formed.');
	    allDone = 1;
	  else
	    if (data(c) == 102) & (data(c+1) == 48),	%find ASCII f next to ASCII 0
		  allDone = 1;
	    else
		  c = c+1;			%if don't find f and 0 keep going thru data
	    end
	  end
	  
	end
	
	data = data(c:end);		%data with header row
	k = find(data==13);			%find end of header line
	data = data(k+1:end);	%chop off header line
	
	data = str2num(char(data));	%convert data to numeric matrix
	
	%get formant values for F1-F5
	FREQ=zeros(size(data,1),5);
	FREQ(:,1)=data(:,13);
	FREQ(:,2)=data(:,12);
	FREQ(:,3)=data(:,11);
	FREQ(:,4)=data(:,10);
	FREQ(:,5)=data(:,9);
	
	%get amplitude values if synthtype=0; get bandwidth values and convert to amplitude if synthtype=1
	%MAG=zeros(size(data,1),5);
	% if synthtype=0
	% 	MAG(:,1)=data(:,30);
	% 	MAG(:,2)=data(:,29);
	% 	MAG(:,3)=data(:,28);
	% 	MAG(:,4)=data(:,27);
	% 	MAG(:,5)=data(:,26);
	% else
	% 	MAG(:,1)=data(:,21);
	% 	MAG(:,2)=data(:,20);
	% 	MAG(:,3)=data(:,19);
	% 	MAG(:,4)=data(:,18);
	% 	MAG(:,5)=data(:,17);
	% end
	
	%amplitude 1 or 0 for now
	MAG=ones(size(data,1),5);
	
	%not using this matrix any more so get rid of it
	clear data
	
	%set maxfreq and minfreq values
	if samprate>=40000 & samprate<=49000
		maxfreq=20000;
	elseif samprate>=20000 & samprate<=23000
		maxfreq=10000;
	else 
		maxfreq=5000;
	end
	
	minfreq=1;
	
	%error check and clean-up data
	i=find(FREQ>maxfreq | FREQ<minfreq);	%if outside freq range set FREQ and MAG to 0
	FREQ(i)=0;
	MAG(i)=0;
	
	%now check for Fn+1<Fn at same time slice
	a=find(FREQ==0);
	FREQ(a)=[nan];	%change 0 to nan to avoid checking these values
	
	tempdat = zeros(size(FREQ));
	tempdat(:,2:5) = FREQ(:,2:5) - FREQ(:,1:4); %tempd holds diff from prev F to this one
	[i,j] = find(tempdat<0);
	
	FREQ(i,j)=FREQ(i,j-1);
	%MAG(i,j)=MAG(i,j-1);
	
	%change nans back to zero
	FREQ(a)=0;
	
	SWXwrite(5,TIMESLICES',FREQ',MAG');
else
	error('File must have .stk extension');
end
