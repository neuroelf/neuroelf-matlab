function varargout = ne_audioplayer(varargin)
% ne_audioplayer  - open an audioplayer window
%
% FORMAT:       ne_audioplayer(SRC, EVT, data)
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       data        either data or filename, or a command:
%                   'close'  - close window
%                   'pause'  - pause playback
%                   'resume' - resume playback
%                   'stop'   - stop playback
%
% No output fields.
%
% Example:
%
%       ne_audioplayer(0, 0, 'somefile.mp3');

% Version:  v1.0
% Build:    15111512
% Date:     Nov-15 2015, 12:59 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% global variables
global ne_gcfg;
fcfg = ne_gcfg.fcfg.AP;
ch = ne_gcfg.h.AP;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% echo
if ne_gcfg.c.echo
    ne_echo('neuroelf_gui(''audioplayer'');');
end

% all commands other than 'close' require player to be open
try
    if nargin < 3 || ...
       ~ischar(varargin{3}) || ...
       ~strcmpi(varargin{3}(:)', 'close')

        % open satellite window
        if ~isstruct(ch) || ...
           ~isfield(ch, 'APFig') || ...
           ~isxfigure(ch.APFig, true)
            hFig = xfigure([neuroelf_path('tfg') '/ne_audioplayer.tfg']);
            hTag = hFig.TagStruct;

            % add line
            axspect = hTag.IM_NE_AUDIO_spect.MLHandle;
            hold(axspect, 'on');
            ispect = get(axspect, 'Children');
            lspect = plot(axspect, [1; 1], [0; 320], 'Color', [1, 1, 1]);

            % store
            ne_gcfg.h.AP = struct( ...
                'APFig', hFig, ...
                'APSpect', ispect, ...
                'APSLine', lspect, ...
                'hTags', hTag);
            ch = ne_gcfg.h.AP;

            % set callbacks
            hTag.UIM_NE_Audio_LoadFile.Callback = {@ne_audioplayer, 'load'};

            % create player with timer function
            p = audioplayer(zeros(1152, 2), 44100);
            p.TimerFcn = @ne_audioplayer_timer;
            p.TimerPeriod = 1152 / 44100;

            % color config
            spcol = hsvconv([((2/3):-1/300:0)', ones(201,2)], 1);
            spcol = double([spcol(2:end, :); spcol(2:end, :)]);

            % fill config
            fcfg = struct( ...
                'bpos',   0, ...
                'dlen',   2, ...
                'dmode',  'spect', ...
                'dopts',  {{}}, ...
                'dpos',   0, ...
                'freq',   0, ...
                'player', p, ...
                'plist',  {{}}, ...
                'plisti', [], ...
                'plists', false, ...
                'ppos',   0, ...
                'spcol',  spcol, ...
                'wfunc',  filtwin(1152, 'hamming'), ...
                'y',      zeros(1152, 2));
            ne_gcfg.fcfg.AP = fcfg;

            % settings
            hFig.Resize = 'off';
        end
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg(sprintf('Error opening audioplayer: %s.', ...
        ne_eo.message), 'NeuroElf - error', 'modal'));
    return;
end

% close window
if nargin > 2 && ...
    ischar(varargin{3}) && ...
    strcmpi(varargin{3}(:)', 'close')

    % window configured
    if isstruct(fcfg) && ...
        isfield(fcfg, 'player') && ...
        isa(fcfg.player, 'audioplayer')
        stop(fcfg.player);
        delete(fcfg.player);
    end
    ne_gcfg.fcfg.AP = [];

    % window is open
    if isstruct(ch) && ...
        isfield(ch, 'APFig') && ...
        isxfigure(ch.APFig, true)
        ch.APFig.Delete;
    end
    ne_gcfg.h.AP = [];
    return;
end

% we now know ch must be valid!
hFig = ch.APFig;
hTag = ch.hTags;

% nothing to do
if nargin < 3
    return;
end

% character argument
if ischar(varargin{3})
    
    % commands
    switch (lower(varargin{3}(:)'))
        
        % display type
        case {'display'}
            
        % set playback to first song
        case {'first'}
            
        % load song/playlist
        case {'load'}
            
            % no input
            if nargin < 4 || ...
               ~ischar(varargin{4}) || ...
                numel(varargin{4}) < 5 || ...
               ~any(varargin{4}(:) == '.') || ...
                isempty(regexpi(varargin{4}(:)', '\.(mp3|nal|swi|sws|wav)$'))

                % get input
                [ifile, ipath] = uigetfile( ...
                    {'*.mp3', 'MP3 audio files (*.mp3)'; ...
                     '*.nal', 'NeuroElf Audio List files (*.nal)'; ...
                     '*.swi;*.sws', 'Sine-Wave information/speech files (*.swi, *.sws)'; ...
                     '*.wav', 'WAV audio files (.wav)'}, ...
                    'Please select an NeuroElf-compatible file...');
                if isequal(ifile, 0) || ...
                    isequal(ipath, 0)
                    return;
                end
                if isempty(ipath)
                    ipath = pwd;
                end
                ifile = fullfile(ipath, ifile);
            else
                ifile = varargin{4}(:)';
            end

            % check input
            if exist(ifile, 'file') ~= 2
                uiwait(warndlg('Requested file does not exist.', ...
                    'NeuroElf - audioplayer', 'modal'));
                return;
            end

            % try to load input
            iext = lower(ifile(end-2:end));
            try
                switch (iext)
                    
                    % MP3
                    case {'mp3'}
                        [ne_gcfg.fcfg.AP.y, ne_gcfg.fcfg.AP.freq] = mp3read(ifile);

                    % play list
                    case {'nal'}

                        % read list
                        

                        % play
                        ne_audioplayer(0, 0, 'first');
                        
                        % return
                        return;

                    % sine-wave information file
                    case {'swi'}
                    case {'sws'}

                    % WAV audio file
                    case {'wav'}
                        [ne_gcfg.fcfg.AP.y, ne_gcfg.fcfg.AP.freq] = audioread(ifile);
                end
            catch ne_eo;
                uiwait(warndlg(['Error reading audio file ' ifile '.' char(10), ...
                    ne_eo.message], 'NeuroElf - audioplayer', 'modal'));
                return;
            end
            
            % unless returned already, now expects y and f
            stop(fcfg.player);
            delete(fcfg.player);
            ne_gcfg.fcfg.AP.player = audioplayer(ne_gcfg.fcfg.AP.y, ne_gcfg.fcfg.AP.freq);
            
            % settings
            fcfg = ne_gcfg.fcfg.AP;
            fcfg.player.TimerFcn = @ne_audioplayer_timer;
            fcfg.player.TimerPeriod = 1152 / 44100;
            ne_gcfg.fcfg.AP.dpos = Inf;
            ne_gcfg.fcfg.AP.wfunc = filtwin( ...
                round(ne_gcfg.fcfg.AP.freq * 1152 / 44100), 'nutall');

            % start playback
            play(fcfg.player);
            
        % set playback to next song
        case {'next'}
            
        % play current song
        case {'play'}
            
            % restart playback (of current item)
            play(fcfg.player);
            
        % pause
        case {'pause'}
            
            % pause playback
            pause(fcfg.player);
            
        % set playback to previous song (or return to beginning of song)
        case {'previous'}
            
        % resume playing
        case {'resume'}
            
            % resume playback
            resume(fcfg.player);
            
        % shuffle setting
        case {'shuffle'}
            
        % stop playback
        case {'stop'}
            
            % stop playback
            stop(fcfg.player);
            
        % set visible property of figure
        case {'visible'}
            
            % input argument
            if nargin > 3
                if islogical(varargin{4}) && ...
                    numel(varargin{4}) == 1
                    if varargin{4}
                        varg = 'on';
                    else
                        varg = 'off';
                    end
                elseif isa(varargin{4}, 'double') && ...
                    numel(varargin{4}) == 1 && ...
                   ~isinf(varargin{4}) && ...
                   ~isnan(varargin{4})
                    if varargin{4} > 0
                        varg = 'on';
                    else
                        varg = 'off';
                    end
                elseif ischar(varargin{4}) && ...
                    any(strcmpi(varargin{4}(:)', {'on', 'off'}))
                    varg = lower(varargin{4}(:)');
                else
                    varg = hFig.Visible;
                    if numel(varg) == 3
                        varg = 'on';
                    else
                        varg = 'off';
                    end
                end
            else
                varg = hFig.Visible;
                if numel(varg) == 3
                    varg = 'on';
                else
                    varg = 'off';
                end
            end
            hFig.Visible = varg;
            
        % filename?
        otherwise
            if exist(varargin{3}(:)', 'file') == 2
                ne_audioplayer(0, 0, 'load', varargin{3}(:)');
            end
    end
end


function varargout = ne_audioplayer_timer(varargin)
global ne_gcfg;
ch = ne_gcfg.h.AP;
hTag = ch.hTags;
fcfg = ne_gcfg.fcfg.AP;
varargout = cell(1, nargout);
if ~isstruct(fcfg)
    return;
end
if strcmpi(ch.APFig.Visible, 'off')
    return;
end
p = fcfg.player;

% spectrogram
if strcmp(fcfg.dmode, 'spect')
    
    % compute current window
    cpos = p.CurrentSample;
    dwin = fcfg.dpos + [0, 0.5 .* fcfg.freq * fcfg.dlen];
    
    % needs update
    if cpos < dwin(1) || ...
        cpos > dwin(2)
        
        % compute power spectrum
        plen = floor(cpos + fcfg.freq * fcfg.dlen);
        ne_gcfg.fcfg.AP.dpos = cpos;
        try
            cmax = min(plen, size(fcfg.y, 1));
            y = mean(fcfg.y(cpos:cmax, :), 2);
            if size(y, 1) < (2 * fcfg.freq)
                y(2 * fcfg.freq) = 0;
            end
            [pw, pf] = pwelch(y, struct( ...
                'average', false, ...
                'freq',    fcfg.freq, ...
                'overlap', 0.764, ...
                'units',   'dB', ...
                'window',  fcfg.wfunc));
            pw = limitrangec(0.04 .* (pw + 35), 0, 1, 0);
            if size(pw, 1) < 320 || ...
                size(pw, 2) < 320
                pw(320, 320) = 0;
            end
            if isnan(pw(1))
                pw = zeros(320, 320, 3);
            else
                pw = threshlutc(pw(320:-1:1, 1:320), fcfg.spcol);
            end
            set(ch.APSpect, 'CData', pw);
        catch ne_eo;
            disp(ne_eo.message);
        end
    end
    
    % update line pos
    set(ch.APSLine, 'XData', (1 + 160 * (cpos - ne_gcfg.fcfg.AP.dpos) / fcfg.freq) .* [1; 1]);
    drawnow;
end
