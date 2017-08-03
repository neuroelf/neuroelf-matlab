function [alldata, nOscs, timeslice]= SWIinput;

state=get(gcf,'userdata');

prompts={'Enter the number of oscillators:',...
				'Enter the number of time slices:',...
				'Enter the value in msec for time slice increments:',...
				'Enter the value in msec for the first time slice:'};
				
initvals ={ ' ', ' ', ' ','0'};
% set up prompt for first time slice value so that...
% default value=0


sw_output=inputdlg(prompts, '', 1, initvals);

if isempty(sw_output); return; end;

nOscs=str2num(sw_output{1,:});

timeslices=str2num(sw_output{2,:});

increment=str2num(sw_output{3,:});

initialvalue=str2num(sw_output{4,:});

if initialvalue==0
	finalvalue=(timeslices*increment)-increment;
else
	finalvalue=timeslices*increment;
end


timeslice=initialvalue:increment:finalvalue;


alldata=pad(num2str(nOscs),30);
for currentdata=1:timeslices
	ts=timeslice(1,currentdata);
	ts = pad(num2str(ts),30);
	alldata = [alldata; ts];
	
	for thisosc = 1:nOscs
		fa = pad('0.00000	0.00000',30);
		alldata=[alldata;fa];
	end
end

timeslice=timeslice';

set(state.DATAWINDOW, 'String', alldata);
set(state.FILENAME,'string',[]);
state.SOUND_OUT=[];
cla;
set(gcf,'userdata',state);
