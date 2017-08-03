function M= Setamp(F,M,timeslice,nOscs);

state=get(gcf,'userdata');

prompts={['SinWave 1-'num2str(nOscs) ' (default=ALL):'],...
				'Amplitude Value (0-1):',...
				'Starting timeslice (default=0):',...
				'Ending timeslice (default=END)'};
				
initvals ={ '', '','0','end'};
% set up prompt for starting time slice value so that...
% default value=0 and ending to end


sw_out=inputdlg(prompts, '', 1, initvals);

if isempty(sw_out); return; end;


swave=str2num(sw_out{1,:});

%Default to all oscillators
if isempty(swave)
	swave=1:nOscs(1,1);
end

%value to which freq will be set
value=str2num(sw_out{2,:});

starttime=str2num(sw_out{3,:});
endtime=str2num(sw_out{4,:});


if isempty(endtime)
	timer=find(timeslice>=starttime);
else
	timer=find(timeslice>=starttime & timeslice<=endtime);
end

%set Amp to Value
for x=1:size(swave,2)
		M(swave(1,x),timer(1,1):timer(end,1))=value;
end



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

for temp5=1:size(i,1)
	M(i,j)=0.0;
end

