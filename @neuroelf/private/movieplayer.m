function h = movieplayer(filename, opts)
% movieplayer  - play a movie using mmread
%
% FORMAT:       [handle = ] movieplayer(filename [, opts])
%
% Input fields:
%
%       filename    file name of movie clip to play
%       opts        optional settings
%        .buffer    number of frames to buffer, default: 120
%        .frate     frame-rate override
%        .himage    image handle to play the movie in
%        .hpause    pause button handle
%        .hplay     play button handle
%        .hslide    slider control handle
%        .nffargs   optional cell array with arguments
%        .nffunc    function handle that is called at every next frame
%
% Output fields:
%
%       handle      image handle
%
% Notes:
%
% - the UserData field of the image handle will be set/overwritten!
% - if handles are set, all handles must be given!
% - the first argument to nffunc will be the UserData of the handle!

% Version:  v1.0
% Build:    15040309
% Date:     Apr-03 2015, 9:43 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2015, Jochen Weber
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

% global variable
global ne_ui;

% make sure structure exists and is initializedblabla
if ~isstruct(ne_ui) || ...
    numel(ne_ui) ~= 1
    ne_ui = struct;
end
if ~isfield(ne_ui, 'movieplayer')
    ne_ui.movieplayer = struct;
    ne_ui.movieplayer.b = {};
    ne_ui.movieplayer.f = {};
    ne_ui.movieplayer.h = [];
    ne_ui.movieplayer.t = {};
    ne_ui.movieplayer.tb = {};
end

% argument check
if nargin < 1 || ...
   ~ischar(filename) || ...
    exist(filename(:)', 'file') ~= 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
filename = filename(:)';

% try to read first frame
try
    ocd = pwd;
    firstframe = mmread(filename, 1, [], false, true);
catch ne_eo;
    cd(ocd);
    rethrow(ne_eo);
end

% parse options
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'buffer') || ...
   ~isa(opts.buffer, 'double') || ...
    numel(opts.buffer) ~= 1 || ...
    isinf(opts.buffer) || ...
    isnan(opts.buffer) || ...
    opts.buffer <= 1
    opts.buffer = 120;
else
    opts.buffer = min(600, round(opts.buffer));
end
if ~isfield(opts, 'frate') || ...
   ~isa(opts.frate, 'double') || ...
    numel(opts.frate) ~= 1 || ...
    isinf(opts.frate) || ...
    isnan(opts.frate) || ...
    opts.frate < 0
    opts.frate = firstframe.rate;
else
    opts.frate = min(600, opts.frate);
end
if ~isfield(opts, 'himage') || ...
   (~isa(opts.himage, 'double') && ...
    ~isa(opts.himage, 'matlab.graphics.primitive.Image')) || ...
    numel(opts.himage) ~= 1 || ...
   ~ishandle(opts.himage) || ...
   ~strcmpi(get(opts.himage, 'Type'), 'image')
    opts.himage = [];
end
if ~isfield(opts, 'hpause') || ...
   (~isa(opts.hpause, 'double') && ...
    ~isa(opts.hpause, 'matlab.ui.control.UIControl')) || ...
    numel(opts.hpause) ~= 1 || ...
   ~ishandle(opts.hpause) || ...
   ~strcmpi(get(opts.hpause, 'Type'), 'uicontrol') || ...
   ~strcmpi(get(opts.hpause, 'Style'), 'pushbutton')
    opts.hpause = [];
end
if ~isfield(opts, 'hplay') || ...
   (~isa(opts.hplay, 'double') && ...
    ~isa(opts.hplay, 'matlab.ui.control.UIControl')) || ...
    numel(opts.hplay) ~= 1 || ...
   ~ishandle(opts.hplay) || ...
   ~strcmpi(get(opts.hplay, 'Type'), 'uicontrol') || ...
   ~strcmpi(get(opts.hplay, 'Style'), 'pushbutton')
    opts.hplay = [];
end
if ~isfield(opts, 'hslide') || ...
   (~isa(opts.hslide, 'double') && ...
    ~isa(opts.hslice, 'matlab.ui.control.UIControl')) || ...
    numel(opts.hslide) ~= 1 || ...
   ~ishandle(opts.hslide) || ...
   ~strcmpi(get(opts.hslide, 'Type'), 'uicontrol') || ...
   ~strcmpi(get(opts.hslide, 'Style'), 'slider')
    opts.hslide = [];
end
if isempty(opts.himage) || ...
    isempty(opts.hpause) || ...
    isempty(opts.hplay) || ...
    isempty(opts.hslide)
    try
        h = mp_makefig(filename);
    catch ne_eo;
        rethrow(ne_eo);
    end
    ud = get(h, 'UserData');
    opts.himage = h;
    opts.hpause = ud.hpause;
    opts.hplay = ud.hplay;
    opts.hslide = ud.hslide;
else
    h = opts.himage;
    hi = find(ne_ui.movieplayer.h == h);
    if ~isempty(hi)
        hi = hi(end:-1:1);
        for hc = hi(:)'
            if ~isempty(ne_ui.movieplayer.t{hi})
                try
                    stop(ne_ui.movieplayer.t{hi});
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
                try
                    delete(ne_ui.movieplayer.t{hi});
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            end
            ne_ui.movieplayer.b(hi) = [];
            ne_ui.movieplayer.f(hi) = [];
            ne_ui.movieplayer.h(hi) = [];
            ne_ui.movieplayer.t(hi) = [];
            ne_ui.movieplayer.tb(hi) = [];
        end
    end
end
if ~isfield(opts, 'nffargs') || ...
   ~iscell(opts.nffargs)
    opts.nffargs = {};
else
    opts.nffargs = opts.nffargs(:)';
end
if ~isfield(opts, 'nffunc') || ...
    numel(opts.nffunc) ~= 1 || ...
   ~isa(opts.nffunc, 'function_handle')
    opts.nffunc = [];
end

% initialize image
mp_initimg(h, filename, firstframe, opts);

% start playing
mp_play(0, 0, h);



% sub-functions



% create new player figure
function h = mp_makefig(filename)

% get ROOT properties
rp = get(0);

% check some properties
if rp.ScreenSize(3) < 800 || ...
    rp.ScreenSize(4) < 600
    error( ...
        'neuroelf:BadScreenSize', ...
        'Screen too small to create new movie player UI.' ...
    );
end

% create figure
f = figure;
set(0, 'CurrentFigure', f);
set(f, ...
    'Position', [round(0.5 .* rp.ScreenSize(3:4) - [328, 300]), 656, 524], ...
    'Color', [0.8, 0.8, 0.8], ...
    'NumberTitle', 'off', ...
    'Name', sprintf('Movie player: %s', filename), ...
    'Units', 'pixels');
figure(f);
drawnow;

% add required controls
% - movie video frame axes and image
va = axes;
set(f, 'CurrentAxes', va);
set(va, 'Units', 'pixels');
h = image(uint8(zeros(640, 480, 3)));
set(va, ...
    'Position', [8, 36, 640, 480],  ...
    'Color', 'none', ...
    'Visible', 'off', ...
    'XColor', [0.8, 0.8, 0.8], ...
    'XTick', [], ...
    'YColor', [0.8, 0.8, 0.8], ...
    'YDir', 'reverse', ...
    'YTick', []);
set(h, 'Visible', 'on');

% - play button
iplay = uint8(192 .* ones(16, 16, 3));
for ic = 1:6
    iplay(2+ic, 6:5+ic, :) = 0;
    iplay(15-ic, 6:5+ic, :) = 0;
end
bplay = uicontrol( ...
    'Style', 'pushbutton', ...
    'CData', iplay, ...
    'Position', [8, 8, 20, 20]);

% - pause button
ipause = uint8(192 .* ones(16, 16, 3));
ipause(3:14, [6, 7, 10, 11], :) = 0;
bpause = uicontrol( ...
    'Style', 'pushbutton', ...
    'CData', ipause, ...
    'Position', [34, 8, 20, 20]);

% - slider control
sl = uicontrol( ...
    'Style', 'slider', ...
    'Position', [60, 8, 588, 20], ...
    'Max', 2, ...
    'Min', 1, ...
    'Value', 1, ...
    'SliderStep', [1, 2]);

% set handles
ud = struct( ...
    'hpause', bpause, ...
    'hplay',  bplay, ...
    'hslide', sl);
set(h, 'UserData', ud);

% initialize figure callbacks
set(f, ...
    'WindowButtonDownFcn', {@mp_btdown, h}, ...
    'WindowButtonMotionFcn', {@mp_btmove, h}, ...
    'WindowButtonUpFcn', {@mp_btup, h});


% initialize handles
function mp_initimg(h, filename, firstframe, opts)

% global variable
global ne_ui;

% set UserData
ud = get(h, 'UserData');
if ~isstruct(ud) || ...
    numel(ud) ~= 1
    ud = struct;
end
ud.buffering = false;
ud.file = filename;
ud.opts = opts;
ud.hpause = opts.hpause;
ud.hplay = opts.hplay;
ud.hslide = opts.hslide;
ud.playing = false;
ud.vpos = 1;
if firstframe.nrFramesTotal <= 0
    ud.vframes = round(firstframe.totalDuration * firstframe.rate);
else
    ud.vframes = firstframe.nrFramesTotal;
end
ud.vfrate = opts.frate;
ud.vheight = firstframe.height;
ud.vwidth = firstframe.width;

% set callbacks
set(ud.hpause, 'Callback', {@mp_pause, h});
set(ud.hplay, 'Callback', {@mp_play, h});
set(ud.hslide, ...
    'Callback', {@mp_slide, h}, ...
    'Min', 0, ...
    'Max', 1, ...
    'SliderStep', [1 / ud.vframes, round(5 * firstframe.rate) / ud.vframes]);

% add to global structure
es = struct('cdata', [], 'colormap', []);
nhi = numel(ne_ui.movieplayer.h) + 1;
ne_ui.movieplayer.b(nhi) = {es([])};
ne_ui.movieplayer.f(nhi) = {[]};
ne_ui.movieplayer.h(nhi) = h;
ne_ui.movieplayer.t(nhi) = {[]};
ne_ui.movieplayer.tb(nhi) = {[]};

% update image handle
hp = get(h, 'Parent');
set(h, ...
    'CData', firstframe.frames.cdata, ...
    'DeleteFcn', {@mp_uninitimg, h}, ...
    'UserData', ud);
set(hp, ...
    'XLim', [0.5, 0.5 + ud.vwidth], ...
    'YLim', [0.5, 0.5 + ud.vheight]);
drawnow;

% handle mouse clicks
function varargout = mp_btdown(varargin)
varargout = cell(1, nargout);

function varargout = mp_btmove(varargin)
varargout = cell(1, nargout);

function varargout = mp_btup(varargin)
varargout = cell(1, nargout);


% handle UI controls
function mp_pause(varargin)

% global variable
global ne_ui;

% get UserData
h = varargin{3};
hi = find(ne_ui.movieplayer.h == h);
if numel(hi) ~= 1
    return;
end
ud = get(h, 'UserData');

% stop and delete timer
es = struct('cdata', [], 'colormap', []);
if ud.playing
    ud.playing = false;
    set(h, 'UserData', ud);
    ne_ui.movieplayer.b{hi} = es([]);
    ne_ui.movieplayer.f{hi} = [];
    if ~isempty(ne_ui.movieplayer.t{hi})
        stop(ne_ui.movieplayer.t{hi});
        delete(ne_ui.movieplayer.t{hi});
    end
    ne_ui.movieplayer.t{hi} = [];
end

% enable slider
try
    set(ud.hslide, 'Enable', 'on');
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

% wait while buffering
while ud.buffering
    pause(0.5);
    ud = get(h, 'UserData');
end


function mp_play(varargin)

% global variable
global ne_ui;

% get UserData
h = varargin{3};
hi = find(ne_ui.movieplayer.h == h);
if numel(hi) ~= 1
    return;
end
ud = get(h, 'UserData');

% create timer
if ~ud.playing
    ud.playing = true;
    ne_ui.movieplayer.t{hi} = timer;
    set(h, 'UserData', ud);

    % disable slider
    set(ud.hslide, 'Enable', 'off');

    % then also fill buffer
    mp_buffer(0, 0, h);

    % set up timer
    set(ne_ui.movieplayer.t{hi}, ...
        'ExecutionMode', 'fixedSpacing', ...
        'Period',        max(0.01, 0.001 * floor(1000 / ud.opts.frate)), ...
        'StopFcn',       {@mp_pause, h}, ...
        'TimerFcn',      {@mp_nextframe, h});
    start(ne_ui.movieplayer.t{hi});
end

function mp_slide(varargin)

% global variable
global ne_ui;

% get UserData
h = varargin{3};
hi = find(ne_ui.movieplayer.h == h);
if numel(hi) ~= 1
    return;
end
ud = get(h, 'UserData');
sv = get(ud.hslide, 'Value');

% if playing
wasplaying = false;
if ud.playing

    % pause first
    wasplaying = true;
    mp_pause(0, 0, h);
end

% update vpos and then UserData
ud.playing = false;
ud.vpos = min(ud.vframes, max(1, round(ud.vframes * sv)));
set(h, 'UserData', ud);

% was playing?
if wasplaying

    % continue
    mp_play(0, 0, h);

% not playing
else

    % simply show one frame
    v = mmread(ud.file, ud.vpos, [], false, true);
    set(h, 'CData', v.frames(1).cdata);
    drawnow;
end


% buffer video
function mp_buffer(varargin)

% global variable
global ne_ui;

% get UserData
h = varargin{3};
hi = find(ne_ui.movieplayer.h == h);
if numel(hi) ~= 1
    return;
end
ud = get(h, 'UserData');

% figure out how many and which frames to buffer
bframes = ne_ui.movieplayer.f{hi};
if ~isempty(bframes)
    mpos = min(ud.vframes, ud.vpos + ud.opts.buffer - 1);
    kframes = (bframes > ud.vpos & bframes <= mpos);
    ne_ui.movieplayer.f{hi} = bframes(kframes);
    ne_ui.movieplayer.b{hi} = ne_ui.movieplayer.b{hi}(kframes);
else
    bframes = ud.vpos;
end
fromp = max(bframes) + 1;

% read additional frames
[f, fp] = mp_readframes(ud.file,  ...
    fromp:min(ud.vframes, fromp+(ud.opts.buffer - numel(ne_ui.movieplayer.b{hi}))));
ne_ui.movieplayer.b{hi} = [ne_ui.movieplayer.b{hi}(:)', f(:)'];
ne_ui.movieplayer.f{hi} = [ne_ui.movieplayer.f{hi}(:)', fp(:)'];



% read frames
function [f, fp] = mp_readframes(file, pos)

% read
ocd = pwd;
try
    v = mmread(file, pos, [], false, true);
    f = v.frames;
    fp = 1 + round(v.rate .* v.times);
catch ne_eo;
    cd(ocd);
    es = struct('cdata', [], 'colormap', []);
    f = es([]);
    fp = [];
    warning( ...
        'neuroelf:MMREADError', ...
        'Error reading video frames: %s', ...
        ne_eo.message ...
    );
end


% back-ground buffering
function mp_bbuffer(varargin)

% global variable
global ne_ui;

% find in global storage
h = varargin{3};
hi = find(ne_ui.movieplayer.h == h);
if numel(hi) ~= 1
    return;
end

% buffer
mp_buffer(0, 0, h);

% stop timer then delete and remove
stop(ne_ui.movieplayer.tb{hi});
delete(ne_ui.movieplayer.tb{hi});
ne_ui.movieplayer.tb{hi} = [];

% unset buffering flag
ud = get(h, 'UserData');
ud.buffering = false;
set(h, 'UserData', ud);


% show next frame
function mp_nextframe(varargin)

% global variable
global ne_ui;

% get UserData
h = varargin{3};
hi = find(ne_ui.movieplayer.h == h);
if numel(hi) ~= 1
    return;
end
ud = get(h, 'UserData');

% not playing
if ~ud.playing

    % valid timer
    if ~isempty(ne_ui.movieplayer.t{hi})
        stop(ne_ui.movieplayer.t{hi});
        delete(ne_ui.movieplayer.t{hi});
        ne_ui.movieplayer.t{hi} = [];
        return;
    end
end

% find next frame in buffer
ff = (ne_ui.movieplayer.f{hi} > ud.vpos);
if isempty(ff)
    return;
end
nf = find(ff);

% empty and at the end of the movie
if isempty(nf) && ...
    ud.vpos == ud.vframes
    mp_pause(0, 0, h);
    return;
end

% update position
ud.vpos = ne_ui.movieplayer.f{hi}(nf(1));
set(ud.hslide, 'Value', ud.vpos / ud.vframes);

% too few frames in buffer
if sum(ff) < (ud.opts.buffer / 3) && ...
   ~ud.buffering && ...
   ne_ui.movieplayer.f{hi}(end) < ud.vframes

    % initialize buffering
    ud.buffering = true;
    ne_ui.movieplayer.tb{hi} = timer( ...
        'ExecutionMode', 'singleShot', ...
        'TimerFcn', {@mp_bbuffer, h}, ...
        'StartDelay', 0.1);
    start(ne_ui.movieplayer.tb{hi});
end

% set data
set(h, 'CData', ne_ui.movieplayer.b{hi}(nf(1)).cdata, 'UserData', ud);

% call function if requested
if ~isempty(ud.opts.nffunc)
    try
        feval(ud.opts.nffunc, ud, ud.opts.nffargs{:});
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end
drawnow;

% uninitialize image
function mp_uninitimg(varargin)

% global variable
global ne_ui;

% pass on to pause
h = varargin{3};
mp_pause(0, 0, h);

% clear global storage
hi = find(ne_ui.movieplayer.h == h);
if ~isempty(hi)
    ne_ui.movieplayer.b(hi) = [];
    ne_ui.movieplayer.f(hi) = [];
    ne_ui.movieplayer.h(hi) = [];
    ne_ui.movieplayer.t(hi) = [];
    ne_ui.movieplayer.tb(hi) = [];
end
