function displayparameters(state,ext);


%(SWI format)
if strcmp(ext,'swi')|strcmp(ext,'SWI')

	nOscsmatrix=sprintf('%d',state.NUM_OSCS);
	alldata=[];
	
	for thistime=1:length(state.TIMESLICES)
		 timeslicematrix=sprintf('%7.0f', state.TIMESLICES(thistime));
	 paramatrix=[];
		    for thisosc=1:state.NUM_OSCS
			    Fmatrix=sprintf('%8.2f\t', state.FREQ(thisosc,thistime));
		            Mmatrix=sprintf('%7.6f',state.MAG(thisosc,thistime));
			    FMmatrix=strcat(Fmatrix,Mmatrix);
			    paramatrix=strvcat(paramatrix,FMmatrix);
		    end
		nextmatrix=[];
		nextmatrix=strvcat(timeslicematrix,paramatrix);
		alldata=strvcat(alldata,nextmatrix);
	end
	
	state.ALLDATA=strvcat(nOscsmatrix,alldata);
	
%(SWX format)
else
	nOscsmatrix=sprintf('%d', state.NUM_OSCS);

	alldata=[];

	for thistime=1:length(state.TIMESLICES)
	    timeslicematrix=sprintf('%7.0f', state.TIMESLICES(thistime));
	    paramatrix=[];
	    for thisosc=1:state.NUM_OSCS
		    Fmatrix=sprintf('%8.2f\t', state.FREQ(thisosc,thistime));
	            Mmatrix=sprintf('%7.6f\t',state.MAG(thisosc,thistime));
		    FMmatrix=strcat(Fmatrix,Mmatrix);
		    paramatrix=strvcat(paramatrix,FMmatrix);
	    end
	nextmatrix=[];
	nextmatrix=strvcat(timeslicematrix,paramatrix);
	alldata=strvcat(alldata,nextmatrix);
	end

	state.ALLDATA=strvcat(nOscsmatrix,alldata);

end

%won't display the whole thing correctly so just display 1st 100 rows
if size(state.ALLDATA,1)>=100
	set(state.DATAWINDOW,'string',state.ALLDATA(1:100,:));
else 
	set(state.DATAWINDOW,'string',state.ALLDATA(1:end,:));
end

%	deal with slider values to display sections of ALLDATA
set(state.SLIDER,'value',1);
[i,j]=size(state.ALLDATA);
set(state.SLIDER,'min',-i);
set(state.SLIDER,'enable','on');	

set(gcf,'userdata',state);
