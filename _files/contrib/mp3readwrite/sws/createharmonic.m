function createharmonic(swave,start_val,end_val,mode,value,ampmult,iterations);

%note to self: if start and end are not whole swave then need
%to pad out newswave to make same length

state=get(gcf,'userdata');

%create harmonic
j=find(state.TIMESLICES>=start_val & state.TIMESLICES<=end_val);	%find indices between start:end

if strcmp(mode,'A')|strcmp(mode,'a')
	newFREQ=zeros(size(state.FREQ,2),1)';	%create a zeros matrix of same length
	newAMP=newFREQ;
	for i=1:iterations
		if i>1,swave=size(state.FREQ,1);end;		%iterate from last one
		newFREQ(j)=state.FREQ(swave,j)+value;		%fill in non-zero values
		newAMP(j)=state.MAG(swave,j).*ampmult;
		state.FREQ=[state.FREQ;newFREQ];			%add new harmonic to matrix
		state.MAG=[state.MAG;newAMP];
		state.NUM_OSCS=state.NUM_OSCS+1;
	end
end

if strcmp(mode,'M')|strcmp(mode,'m')
	newFREQ=zeros(size(state.FREQ,2),1)';
	newAMP=newFREQ;
	for i=1:iterations
		if i>1,swave=size(state.FREQ,1);end;	
		newFREQ(j)=state.FREQ(swave,j).*value;
		newAMP(j)=state.MAG(swave,j).*ampmult;
		state.FREQ=[state.FREQ;newFREQ];		
		state.MAG=[state.MAG;newAMP];
		state.NUM_OSCS=state.NUM_OSCS+1;
	end
end

%check for violations of Nyquist limit and set to max
Nylimit=state.SRATE/2;
i=find(state.FREQ>Nylimit);
state.FREQ(i)=Nylimit;


state.SOUND_OUT=synthtrax(state.FREQ,state.MAG,state.TIMESLICES,state.SRATE);	%re-synth data for PLAY

set(gcf,'userdata',state);

