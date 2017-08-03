function F= Setfreq(F,M,timeslice,nOscs);

state=get(gcf,'userdata');

prompts={['SinWave 1-'num2str(nOscs) ' (default=ALL):'],...
				'Frequency Value:',...
				'Starting timeslice (default=0):',...
				'Ending timeslice (default=END)'};
				
initvals ={ '', '','0','end'};
% set up prompt for starting time slice value so that...
% default value=0 and ending to end


sw_out=inputdlg(prompts, '', 1, initvals);

if isempty(sw_out); return; end;


swave=str2num(sw_out{1});

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

%set Freqs to Value
for x=1:size(swave,2)
		F(swave(1,x),timer(1,1):timer(end,1))=value;
end



% Check for overflow values of frequency and adjust to min/max range
freqmax=state.SRATE/2;

i=find(F>freqmax);
if ~isempty(i)
	msgbox('Frequency Overflow: Frequencies Reset to 1/2 Sampling Rate', '');
	F(i)=freqmax;
end

