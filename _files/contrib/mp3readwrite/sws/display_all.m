function display_all(F, M, timeslice);

if isempty(timeslice);
	return;
end;


plot(timeslice,F);

xlabel('Time (msec)');
ylabel('Frequency (kHz)');

hold on;	%Do this in order to provide second line indicating relative amplitude

width=M .*80;
D=F+width;

plot(timeslice, D);


%white-out places where M=0
j=find(M~=0);	
newF=F;
newF(j)=[nan];
newline=line(timeslice,newF);
set(newline(:),'color','w');


hold off;
