function display_win(F, M, timeslice);

if isempty(timeslice);
	return;
end;

prompts={'Enter the start time for the display',...
			'Enter the end time for the display'};			

tmp_start=num2str(timeslice(1));
tmp_end=num2str(timeslice(size(timeslice,1)));

initvals={tmp_start, tmp_end};

disp_parameters=inputdlg(prompts, '', 1, initvals);

if isempty(disp_parameters); return; end; %make sure didn't hit cancel

start_val=str2num(disp_parameters{1,:});

end_val=str2num(disp_parameters{2,:});

%get index of times to plot
i=find(timeslice>=start_val & timeslice<=end_val);

plot(timeslice(i,:),F(:,i));

xlabel('Time (msec)');
ylabel('Frequency (kHz)');

hold on;	%Do this in order to provide second line indicating relative amplitude

width=M .*80;
D=F+width;

plot(timeslice(i,:), D(:,i));

%white-out places where M=0
j=find(M~=0);	
newF=F;
newF(j)=[nan];
newline=line(timeslice(i,:),newF(:,i));
set(newline(:),'color','w');

hold off;
