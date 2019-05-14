function sws_gui(action)

%SWS_GUI  - SineWaveSynthesis Control window support

% 06/16/98  steve frost
%1/31/00 2/26/02 revised


%	branch by action...

switch(action),

%-----------------------------------------------------------------------------
%	OPEN:	open SWI or SWX file and prepare for synthrax

	case 'OPEN',

	state=get(gcf,'userdata');
		
		[file, pathName] = uigetfile('*', 'Select Haskins format SWI or SWX file');
		state.FILENAME = [pathName file];
		if isnumeric(state.FILENAME); return; end;					%make sure didn't cancel 
		[p,n,ext]=fileparts(state.FILENAME);
		ext=lower(ext(2:end));
		
		if strcmp(ext,'swi')|strcmp(ext,'SWI')|strcmp(ext,'swx')|strcmp(ext,'SWX')
			cla; set(state.DATAWINDOW, 'string', []);			%clear inputwindow, sound buffer, clear axis
			
			[state.FREQ,state.MAG, state.TIMESLICES, state.NUM_OSCS]=swiread(state.FILENAME);	%read in parameters
			state.SOUND_OUT=synthtrax(state.FREQ,state.MAG,state.TIMESLICES,state.SRATE);	%prepare for PLAY
					
			set(state.FILEDISPLAY,'string',file);
			state.INPUT_FORMAT='SYN';
			set(state.LIST,'enable','on')
		else
			errordlg('File must have .swi or .swx extension',' ');
			return;
		end;
		
	set(gcf,'userdata',state);	
	
	if state.OPENSETTINGS(1)			%check prefs for list on open
		displayparameters(state,ext);
	end
	
	if state.OPENSETTINGS(2)			%check prefs for display on open
		display_all(state.FREQ, state.MAG, state.TIMESLICES);
	end
	
%-----------------------------------------------------------------------------
%	EXTRACT: 	extract Frequency,Magnitude, and Time paramters from WAV file

	case 'EXTRACT',

	state=get(gcf,'userdata');
	
		[file, pathName] = uigetfile('*', 'Select .WAV or .AU file');
		AUDIONAME = [pathName file];
		if isnumeric(AUDIONAME); return; end;
		[p,n,ext]=fileparts(AUDIONAME);
		ext=lower(ext(2:end));
		
		if strcmp(ext,'wav')|strcmp(ext,'WAV')
			cla; set(state.DATAWINDOW, 'string', []);			%clear inputwindow, sound buffer, clear axis
			
			%open wav file
			[state.SOUND_OUT,state.SRATE] =wavread(AUDIONAME);
			
			%dialog box for LPC analysis
			prompts={'Enter the number of formants (default = 5)     ',...
			'Enter window size in msec (default = 20% sampling rate)'...
			'Enter overlap in msec (default = 1/2 window size)'};			

			winsize=num2str(samps2msecs(.02*state.SRATE,state.SRATE));
			overlap=num2str(samps2msecs(.010*state.SRATE,state.SRATE));

			initvals={'5', winsize, overlap};

			disp_parameters=inputdlg(prompts, 'LPC Settings', 1, initvals);

			if isempty(disp_parameters); return; end; %make sure didn't hit cancel

			%get values from dialog box
			poles = 2*(str2num(disp_parameters{1,:}));
			winsize=floor(msecs2samps(str2num(disp_parameters{2,:}),state.SRATE));
			overlap=floor(msecs2samps(str2num(disp_parameters{3,:}),state.SRATE));
			
			%generate freq and mag parameters
			[state.FREQ,state.MAG]=swsmodel2(state.SOUND_OUT,state.SRATE,poles,winsize,overlap);
			%calculate timeslices and number of oscillators
			state.TIMESLICES = [0:length(state.FREQ)-1]*length(state.SOUND_OUT)/(state.SRATE*length(state.FREQ));
			state.TIMESLICES = 10*round(state.TIMESLICES*100)'; %maintains even increment timeslices
			state.NUM_OSCS = size(state.FREQ,1);
			
			state.SOUND_OUT=synthtrax(state.FREQ,state.MAG,state.TIMESLICES,state.SRATE);	%prepare for PLAY
					
			set(state.FILEDISPLAY,'string',file);
			state.INPUT_FORMAT='SYN';
			set(state.LIST,'enable','on')
			else
				errordlg('File must have .wav  extension',' ');
			return;
		end;
		
	set(gcf,'userdata',state);
	
	if state.OPENSETTINGS(1)			%check prefs for list on open
		displayparameters(state,ext);
	end
	
	if state.OPENSETTINGS(2)			%check prefs for display on open
		display_all(state.FREQ, state.MAG, state.TIMESLICES);
	end
	
	
%-----------------------------------------------------------------------------
%	SAVE:	save SWI file

	case 'SAVE',

	state=get(gcf,'userdata');
	if isempty(state.SOUND_OUT), return; end;	%make sure some data to save
	
	switch state.PARAMETER_FORMAT,
		case '.swi',
			swiwrite(state.NUM_OSCS,state.TIMESLICES,state.FREQ,state.MAG);
		case '.swx',
			swxwrite(state.NUM_OSCS,state.TIMESLICES,state.FREQ,state.MAG);	 	
	end;
		
%-----------------------------------------------------------------------------
%	PLAY ALL:	play synthtraxed formatted file

	case 'PLAY',

		state=get(gcf,'userdata');
		
		if isempty(state.SOUND_OUT);return;end;	%make sure some data to play
		sound(state.SOUND_OUT,state.SRATE);
		
%-----------------------------------------------------------------------------
%	PLAY DISPLAYED SELECTION:	play synthtraxed formatted file

	case 'PLAYDISPLAY'

		state=get(gcf,'userdata');
		
		if isempty(state.SOUND_OUT);return;end; %	PLAY:	play synthtraxed formatted file
	
		times=get(gca,'xlim');
		%convert from msec to samples
		times=times*(state.SRATE/1000);
		%check for outside limits
		if times(1) <= 1
			start_time=1;
		else
			start_time=fix(times(1));
		end
		if times(2) >= size(state.SOUND_OUT,2)
			end_time=size(state.SOUND_OUT,2);
		else
			end_time=fix(times(2));
		end
	
		sound(state.SOUND_OUT(start_time:end_time),state.SRATE);
		
%-----------------------------------------------------------------------------
%	WAVEOPEN:	open sound file

	case 'WAVEOPEN',

		state=get(gcf,'userdata');
			
		[file, pathName] = uigetfile('*', 'Select .WAV or .AU file');
		AUDIONAME = [pathName file];
		if isnumeric(AUDIONAME); return; end;
		[p,n,ext]=fileparts(AUDIONAME);
		ext=lower(ext(2:end));
		
		if strcmp(ext,'WAV')|strcmp(ext,'wav')
			set(state.DATAWINDOW, 'string', []);	 cla;	set(state.LIST,'enable','off');%clear inputwindow, sound buffer, clear axis
			
			[state.SOUND_OUT, state.SRATE]=wavread(AUDIONAME);
			set(state.FILEDISPLAY,'string',file);
			state.INPUT_FORMAT='WAV';
			
		elseif strcmp(ext,'AU')|strcmp(ext,'au')
			set(state.DATAWINDOW, 'string', []);	 cla; set(state.LIST,'enable','off');%clear inputwindow, sound buffer, clear axis

			[state.SOUND_OUT, state.SRATE]=auread(AUDIONAME);
			set(state.FILEDISPLAY,'string',file);
			state.INPUT_FORMAT='AU';
		else
			errordlg('File must have .wav or .au extension',' ');
			return;
		end;
	
		set(gcf,'userdata',state);
%-----------------------------------------------------------------------------
	%AUDIOSAVE:	Save output of synthtrax in .AU or .WAV format
	
	case 'AUDIOSAVE',
	
		state=get(gcf,'userdata');
		
		if isempty(state.SOUND_OUT), return; end;
		if isempty(state.OUTPUT_FORMAT)
			SOUND_FORMAT='.wav';
		else
			SOUND_FORMAT=state.OUTPUT_FORMAT;
		end;

		[file, pathName] = uiputfile(SOUND_FORMAT, ['Save audio as ' SOUND_FORMAT ' format']);
		fileName = [pathName file];
	
		if isnumeric(fileName), return; end;

		if strcmp( SOUND_FORMAT,'.WAV')|strcmp(SOUND_FORMAT,'.wav')
			wavwrite2(state.SOUND_OUT, state.SRATE, fileName);
		else
			auwrite2(state.SOUND_OUT, state.SRATE, fileName);
		end;
		
%-----------------------------------------------------------------------------
	%INPUTDATA:	interactive SWI matrix creation
	%NOT USED IN THIS VERSION
	
	case 'INPUTDATA',

		[alldata, nOscs, timeslice]= swiinput;
		
%-----------------------------------------------------------------------------
	%UPDATEDATA: Freq,Amplitude,Timeslice Parameters update from Edit box
	%NOT USED IN THIS VERSION
	
	case 'UPDATEDATA',

	state=get(gcf,'userdata');

	alldata=get(state.DATAWINDOW, 'string');
	alldata=num2str(alldata);
	alldata = greps(alldata,',', ' ');
	set(state.DATAWINDOW, 'string', alldata);
	
	[state.TIMESLICES,state.FREQ,state.MAG,state.NUM_OSCS,mat_error]=inpread(alldata);
	if mat_error ~=1; return; end;
		
	%Check to see (on callback) if timeslice info error
	temp = diff(state.TIMESLICES);
	templist = find(~(temp>0));
	timelist = state.TIMESLICES(templist);

	if length(templist) > 0
	  ts_error=errordlg(['Error near time point(s) ' num2str(timelist)], ' ');
	  return;
	end
	
	state.SOUND_OUT=synthtrax(state.FREQ,state.MAG,state.TIMESLICES,state.SRATE);
	
	set(gcf,'userdata',state);
	
%-----------------------------------------------------------------------------
	%MOVETHROUGHLIST:	   Move up/down in parameter list based on current value of slider 

	case 'MOVETHRULIST',

	state=get(gcf,'userdata');
	
	topofdisplay=abs(get(state.SLIDER,'value'));
	topofdisplay=floor(topofdisplay);
	
	if topofdisplay==0,
		topofdisplay=1;
	end;
	
	if topofdisplay+100<=abs(get(state.SLIDER,'min'))
		set(state.DATAWINDOW,'string',state.ALLDATA(topofdisplay:topofdisplay+100,:));
	else 
		set(state.DATAWINDOW,'string',state.ALLDATA(topofdisplay:end,:));
	end
	
%-----------------------------------------------------------------------------
	%LIST FUCTION:	List Parameter in Edit Box
	
	case 'LISTPARAMETERS',
	
	state=get(gcf,'userdata');
	if isempty(state.FREQ), return; end;
		
	
	[p,n,ext]=fileparts(state.FILENAME);
	ext=lower(ext(2:end));
	
	displayparameters(state,ext);
	
%-----------------------------------------------------------------------------
	%DISPLAYDAT: Graphical display of Frequency x Time
	
	case 'DISPLAYDAT',
	
	state=get(gcf,'userdata');
	
	display_win(state.FREQ, state.MAG, state.TIMESLICES);
	
%-----------------------------------------------------------------------------
	%UNDO: undo the last edit change
	
	case 'UNDOLAST',

	state=get(gcf,'userdata');
	set(gcf,'userdata',state.COPY);
	
	set(findobj('tag','undo'), 'enable','off');

%-----------------------------------------------------------------------------
	%MULTDATA:	Multiply data values for Freq, Amp, or Time

	case 'MULTDATA',

	state=get(gcf,'userdata');
	
	if strcmp(state.INPUT_FORMAT,'SYN')
		%make backup for undo and enable undo option
		state.COPY=state;
		set(findobj('tag','undo'),'enable','on');
		
		[state.FREQ,state.MAG,state.TIMESLICES,state.NUM_OSCS]= swimultiply(state.FREQ,state.MAG,state.TIMESLICES,state.NUM_OSCS,state.SRATE);	
		state.SOUND_OUT=synthtrax(state.FREQ,state.MAG,state.TIMESLICES,state.SRATE);
		
		set(state.SWS_FIG,'userdata',state);
	end;
	
%-----------------------------------------------------------------------------
	%ADDDATA:	Add data values for Freq, Amp, or Time

	case 'ADDDATA',
	
	state=get(gcf,'userdata');
	
	if strcmp(state.INPUT_FORMAT,'SYN')
		%make backup for undo and enable undo option
		state.COPY=state;
		set(findobj('tag','undo'),'enable','on');
		
		[state.FREQ,state.MAG,state.TIMESLICES,state.NUM_OSCS]= swiaddition(state.FREQ,state.MAG,state.TIMESLICES,state.NUM_OSCS,state.SRATE);	
		state.SOUND_OUT=synthtrax(state.FREQ,state.MAG,state.TIMESLICES,state.SRATE);
		set(state.SWS_FIG,'userdata',state);
	end;

%-----------------------------------------------------------------------------
	%SETFREQ;	Set Freq of selected sinewaves to fixed value
	
	case 'SETFREQ',
	
	state=get(gcf,'userdata');
	
	if strcmp(state.INPUT_FORMAT,'SYN')
		%make backup for undo and enable undo option
		state.COPY=state;
		set(findobj('tag','undo'),'enable','on');
		
		state.FREQ= setfreq(state.FREQ,state.MAG,state.TIMESLICES,state.NUM_OSCS);
		state.SOUND_OUT=synthtrax(state.FREQ,state.MAG,state.TIMESLICES,state.SRATE);
		set(state.SWS_FIG,'userdata',state);
	end
	
%-----------------------------------------------------------------------------
	%SETAMP;	Set Amplitude of selected sinewaves to fixed value
	
	case 'SETAMP',
	
	state=get(gcf,'userdata');
	
	if strcmp(state.INPUT_FORMAT,'SYN')
		%make backup for undo and enable undo option
		state.COPY=state;
		set(findobj('tag','undo'),'enable','on');
		
		state.MAG= setamp(state.FREQ,state.MAG,state.TIMESLICES,state.NUM_OSCS);
		state.SOUND_OUT=synthtrax(state.FREQ,state.MAG,state.TIMESLICES,state.SRATE);
		set(state.SWS_FIG,'userdata',state);
	end	
	
%-----------------------------------------------------------------------------
	%WAVEFRM:	Waveform display
	
	case 'WAVEFORM',

	state=get(gcf,'userdata');
	
	if isempty(state.SOUND_OUT); return; end;				%make sure data to plot and delete any previous waveform
	if findobj('tag','waveformfigure'), delete(findobj('tag','waveformfigure')), end;
	
	state.WAVEFORM_FIG = figure('Color',[0.8 0.8 0.8], ...
		'MenuBar','none', ...
		'NumberTitle','off', ...
		'Tag','waveformfigure');
	plot(state.SOUND_OUT);
	xlabel('samples');
	ylabel('amplitude');
	
	% if checked map amplitude to -1/+1 range else use min/max range
	if strcmp(get(state.MENUS(1,12),'checked'),'on')
		set(gca,'ylim',[-1 1]);
	end;

	set(state.SWS_FIG,'userdata',state);
	
%-----------------------------------------------------------------------------
	%AUDIOFORMAT:	Select Audio File Format for Saving Audio
	
	case 'AUDIOFORMAT',

	state=get(gcf,'userdata');	
	outputvalue=get(gcbo,'userdata');
	
	if outputvalue==1
		set(state.MENUS(1,3), 'checked', 'on');
		set(state.MENUS(1,4), 'checked', 'off');
		state.OUTPUT_FORMAT='.wav';
	else
		set(state.MENUS(1,3), 'checked', 'off');
		set(state.MENUS(1,4), 'checked', 'on');
		state.OUTPUT_FORMAT='.au';
	end

	set(gcf,'userdata',state);
%-----------------------------------------------------------------------------
	%PARAMETERFORMAT:	Select Output Format (SWI or SWX)
	
	case 'PARAMETERFORMAT',

	state=get(gcf,'userdata');
	parametervalue=get(gcbo,'userdata');
	
	if parametervalue==1
		set(state.MENUS(1,1),'checked','on');
		set(state.MENUS(1,2),'checked','off');
		state.PARAMETER_FORMAT='.swi';
	else
		set(state.MENUS(1,1),'checked','off');
		set(state.MENUS(1,2),'checked','on');
		state.PARAMETER_FORMAT='.swx';
	end

	set(gcf,'userdata',state);

%-----------------------------------------------------------------------------


	%SAMPLERATE:	Select Sampling Rate for Saving Audio
	
	case 'SAMPLERATE',

	state=get(gcf,'userdata');	
	samplingratevalue=get(gcbo,'userdata');
	
	if samplingratevalue==1
		set(state.MENUS(1,5), 'checked', 'on');
		set(state.MENUS(1,6:10), 'checked', 'off');
		state.SRATE=8000;
	elseif samplingratevalue==2
		set(state.MENUS(1,5:10), 'checked', 'off');
		set(state.MENUS(1,6), 'checked', 'on');
		state.SRATE=10000;
	elseif samplingratevalue==3
		set(state.MENUS(1,5:10), 'checked', 'off');
		set(state.MENUS(1,7), 'checked', 'on');
		state.SRATE=11025;
	elseif samplingratevalue==4
		set(state.MENUS(1,5:10), 'checked', 'off');
		set(state.MENUS(1,8), 'checked', 'on');
		state.SRATE=20000;
	elseif samplingratevalue==5
		set(state.MENUS(1,5:10), 'checked', 'off');
		set(state.MENUS(1,9), 'checked', 'on');
		state.SRATE=22050;
	elseif samplingratevalue==6
		set(state.MENUS(1,5:9), 'checked', 'off');
		set(state.MENUS(1,10), 'checked', 'on');
		state.SRATE=44100;
	end;

	set(gcf,'userdata',state);
%-----------------------------------------------------------------------------
	%WAVEFORMSCALE:	set preference for amplitude scaling of waveform
	
	case 'WAVEFORMSCALE',
	
	state=get(gcf,'userdata');
	amplitudescalevalue=get(gcbo,'userdata');
	
	if amplitudescalevalue==1
		set(state.MENUS(1,11),'checked','on');
		set(state.MENUS(1,12),'checked','off');
	elseif amplitudescalevalue==2
		set(state.MENUS(1,11),'checked','off');
		set(state.MENUS(1,12),'checked','on');
	end

	if ~isempty(findobj('tag','waveformfigure'))
		figure(state.WAVEFORM_FIG);
		if strcmp(get(state.MENUS(1,11),'checked'),'on')
			set(gca,'ylimmode','auto');
		else
			set(gca,'ylim',[-1 1]);
		end
	end

%-----------------------------------------------------------------------------
	%OPENSETTINGS:	set preference for list and display on opening .SWI and .SWX files
	
	case 'OPENSETTINGS',
	
	state=get(gcf,'userdata');
	
	config_dialog=dialog('units', 'normalized',...
		'position', [.30 .55 .20 .15],...
		'Color',[0.8 0.8 0.8], ...
		'MenuBar','none', ...
		'Name','Display Options', ...
		'NumberTitle','off', ...
		'Tag','openingsettings',...
		'ButtonDownFcn', '', ...
		'closerequestfcn', 'uiresume');
	
	display_box1=uicontrol(config_dialog,...
		'units', 'normalized',...
		'style', 'checkbox',...
		'position', [.2 .6 .75 .25],...
		'value', state.OPENSETTINGS(2),...
		'string', 'Display');
		
	display_box2=uicontrol(config_dialog,...
		'units', 'normalized',...
		'style', 'checkbox',...
		'position', [.2 .35 .75 .25],...
		'value', state.OPENSETTINGS(1),...
		'string', 'Parameter List');
		

	%OK, cancel buttons
	uicontrol(config_dialog, ...
		'units', 'normalized',...
		'Position',[.6 .01 .25 .17], ...
		'String','OK', ...
		'Callback','set(gcbf,''UserData'',1);uiresume');
	uicontrol(config_dialog, ...
		'units', 'normalized',...
		'Position',[.1 .01 .25 .17], ...
		'String','Cancel', ...
		'Callback','uiresume');
	
	% wait for input	
	uiwait(config_dialog);
	if get(config_dialog, 'userdata'),
		state.OPENSETTINGS(2) = get(display_box1, 'value');
		state.OPENSETTINGS(1)=get(display_box2, 'value');
		set(state.SWS_FIG, 'userdata', state);
	end;	
	delete(config_dialog);	

%-----------------------------------------------------------------------------
	%CREATE HARMONICS
	
	case 'HARMONIC',
	
	state=get(gcf,'userdata');
	
	if isempty(state.FREQ); return; end;
		
	prompts={['Sinewave 1-' num2str(state.NUM_OSCS)],...
				'Start',...
				'End',...
				'Mode ([A]dd or [M]ultiply)',...
				'Value to Add/Multiply',...
				'Amplitude Multiplier',...
				'Iterations '};
				
	tmp_swave='1';
	tmp_start=num2str(state.TIMESLICES(1));
	tmp_end=num2str(state.TIMESLICES(size(state.TIMESLICES,1)));
	tmp_amp='1';
	tmp_iterate='1';
	
	initvals={tmp_swave,tmp_start,tmp_end,'','',tmp_amp,tmp_iterate};
	
	harmonics=inputdlg(prompts,'',1,initvals);
	
	if isempty(harmonics); return;end;	%in case hit cancel
		

	
	swave=str2num(harmonics{1});
	start_val=str2num(harmonics{2});
	end_val=str2num(harmonics{3});
	mode=harmonics{4};
	value=str2num(harmonics{5});
	ampmult=str2num(harmonics{6});
	iterations=str2num(harmonics{7});
	
%checks
if swave==0 | swave>size(state.FREQ,1)
	errordlg('Sinewave not in range',' '); return;
end;

	%make backup for undo and enable undo option
	state.COPY=state;
	set(findobj('tag','undo'),'enable','on');
	set(gcf,'userdata',state);

createharmonic(swave,start_val,end_val,mode,value,ampmult,iterations);

%-----------------------------------------------------------------------------
	%DELETE SINEWAVE
	
	case 'SINEWAVEDELETE',
		
	state=get(gcf,'userdata');
	
	if isempty(state.FREQ); return; end;
	
	prompt={['Enter the  sinewave(s) to delete (1-' num2str(state.NUM_OSCS) '):']};	%ex. 1 5:10 would delete 1, 5-10
	
	parameters=inputdlg(prompt,'',1);
	
	if isempty(parameters); return; end; %make sure didn't hit cancel
		
	%make backup for undo and enable undo option
	state.COPY=state;
	set(findobj('tag','undo'),'enable','on');
	set(gcf,'userdata',state);
	
	waves2delete=str2num(parameters{1});
		
	sinewaves2keep=1:size(state.FREQ,1);
	sinewaves2keep=find(~ismember(sinewaves2keep,waves2delete));
	
	state.FREQ=state.FREQ(sinewaves2keep,:);
	state.MAG=state.MAG(sinewaves2keep,:);
	state.NUM_OSCS=state.NUM_OSCS-length(waves2delete);
	
	state.SOUND_OUT=synthtrax(state.FREQ,state.MAG,state.TIMESLICES,state.SRATE);	%re-synth data for PLAY
	
	set(gcf,'userdata',state);
	
%-----------------------------------------------------------------------------
	%AVERAGE:	report average freq and amplitude for all sinewaves between START and END
	
	case 'AVERAGE',
	
	state=get(gcf,'userdata');
	
	ave_freqandamp(state.FREQ,state.MAG,state.TIMESLICES);
	

%-----------------------------------------------------------------------------

	%ZOOM: Zoom view on display
	
	case 'ZOOMVIEW',
	zoomstate=get(gcbo,'userdata');
	if isequal(zoomstate, 0)
		zoom on;
		set(gcbo, 'userdata', 1);
	elseif isequal(zoomstate, 1);
		zoom off;
		set(gcbo, 'userdata', 0);

	end;

%-----------------------------------------------------------------------------
	%QUIT
	
	case 'BYEBYE',
	close all;return

%-----------------------------------------------------------------------------
%	error

	otherwise,
		error(['sws_gui:  unrecognized action (', action, ')']);
	
end;
