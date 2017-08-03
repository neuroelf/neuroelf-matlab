function [F,M,timeslice,nOscs]= SWImultiply(F,M,timeslice,nOscs,samprate);


prompts={'Multiply Frequency(F), Amplitude(A), or Time(T):',...
				['SinWave 1-'num2str(nOscs) ' (default=ALL):'],...
				'Multiplier:',...
				'Starting timeslice (default=0):',...
				'Ending timeslice (default=END)'};
				
initvals ={ '', '', '','0','end'};
% set up prompt for starting time slice value so that...
% default value=0 and ending to end


sw_out=inputdlg(prompts, '', 1, initvals);

if isempty(sw_out); return; end;

tbm=(sw_out{1,:});

swave=str2num(sw_out{2,:});
%Default to mult all oscillators
if isempty(swave)
	swave=1:nOscs(1,1);
end


multiplier=str2num(sw_out{3,:});

starttime=str2num(sw_out{4,:});

endtime=str2num(sw_out{5,:});


if isempty(endtime)
	timer=find(timeslice>=starttime);
else
	timer=find(timeslice>=starttime & timeslice<=endtime);
end


if isequal('F',tbm) | isequal ('f', tbm)
	for x=1:size(swave,2)
		F(swave(1,x),timer(1,1):timer(end,1))=F(swave(1,x),timer(1,1):timer(end,1)) .* multiplier;
	end
elseif isequal('A',tbm) | isequal ('a', tbm)
	for x=1:size(swave,2)
		M(swave(1,x),timer(1,1):timer(end,1))=M(swave(1,x),timer(1,1):timer(end,1)) .* multiplier;
	end
elseif isequal('T',tbm) | isequal ('t', tbm)
	timeslice(timer(1,1):timer(end,1),:)=timeslice(timer(1,1):timer(end,1),:) .* multiplier;
else
	errordlg('Can only multiply Frequency, Amplitude, or Time', '');
	return;
end;

% Check for overflow values of amplitude and adjust to min/max range
i=find(M>1.0);		
if ~isempty(i)
	msgbox('Amplitude Overflow: Amplitude Values Reset to Min/Max Range', '');
	M(i)=1.0;
end


i=find(M<0.0);
if ~isempty(i)
	msgbox('Amplitude Overflow: Amplitude Values Reset to Min/Max Range', '');
	M(i)=0.0;
end


% Check for overflow values of frequency and adjust to min/max range
freqmax=samprate/2;

i=find(F>freqmax);
if ~isempty(i)
	msgbox('Frequency Overflow: Frequencies Reset to 1/2 Sampling Rate', '');
	F(i)=freqmax;
end
