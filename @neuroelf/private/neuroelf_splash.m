function timr = neuroelf_splash(ihandle)
% neuroelf_splash  - NeuroElf splash display on an image handle
%
% FORMAT:       tmr = neuroelf_splash(ihandle)
%
% Input fields:
%
%       ihandle     image handle on which the splash is displayed
%
% Output fields:
%
%       tmr         timer object to be started

% Version:  v1.1
% Build:    16060912
% Date:     Jun-09 2016, 12:19 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011 - 2016, Jochen Weber
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

% test input
if nargin < 1 || numel(ihandle) ~= 1 || ...
   (~isa(ihandle, 'double') && ~isa(ihandle, 'matlab.graphics.primitive.Image')) || ...
   ~ishandle(ihandle) || ~strcmpi(get(ihandle, 'Type'), 'image')

    % don't bail out, simply don't do anything
    return;
end

% unset any previous information
if isfield(ne_ui, 'splash') && isstruct(ne_ui.splash) && numel(ne_ui.splash) == 1

    % lingering transimg
    try
        delete(ne_ui.splash.timg);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

% figure out whether we need down-sampling of hi-res images
try
    rssi = true;
    tisize = [624, 224];
    tffile = [tempname '.png'];
    sfig = ancestor(ihandle, 'figure');
    sfigsize = get(sfig, 'Position');
    tfimg = getframe(sfig);
    tfimg = tfimg.cdata;
    if all(round([size(tfimg, 2), size(tfimg, 1)] ./ sfigsize(3:4)) == 2)
        rssi = false;
        tisize = [1248, 448];
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    if exist(tffile, 'file') > 0
        delete(tffile);
    end
end

% setup stuff
ne_ui.splash = struct('alpha', {cell(1, 8)}, 'curve', [], 'ihandle', ihandle, ...
    'iter', 0, 'timer', timer, 'timg', transimg(tisize(1), tisize(2), [1, 1, 1]));

% get handles
talp = ne_ui.splash.alpha;
timg = ne_ui.splash.timg;
timr = ne_ui.splash.timer;

% error handling
try

    % splash images path
    sp = neuroelf_path('splash');

    % load background
    im = imread([sp '/curve.jpg']);
    if ~rssi
        im = image_resize(im, 2 * size(im, 1), 2 * size(im, 2));
    end
    ne_ui.splash.curve = rgb2hsv(im);
    addlayer(timg, im, 1);
    talp{1} = single(1);

    % load and add objects for the right-hand side
    im = imread([sp '/slice.jpg']);
    if ~rssi
        im = image_resize(im, 2 * size(im, 1), 2 * size(im, 2));
    end
    ima = single(min(1, 0.01 .* (765 - sum(single(im), 3))));
    addlayer(timg, im, ima);
    talp{2} = ima;
    im = imread([sp '/shen.jpg']);
    if ~rssi
        im = image_resize(im, 2 * size(im, 1), 2 * size(im, 2));
    end
    addlayer(timg, im, 0);
    talp{3} = single(min(1, 0.01 .* (765 - sum(single(im), 3))));
    im = imread([sp '/surf.jpg']);
    if ~rssi
        im = image_resize(im, 2 * size(im, 1), 2 * size(im, 2));
    end
    addlayer(timg, im, 0);
    talp{4} = single(min(1, 0.01 .* (765 - sum(single(im), 3))));
    im = imread([sp '/render.jpg']);
    if ~rssi
        im = image_resize(im, 2 * size(im, 1), 2 * size(im, 2));
    end
    addlayer(timg, im, 0);
    talp{5} = single(min(1, 0.01 .* (765 - sum(single(im), 3))));

    % load main text, add white middle-layer, then other texts
    [im, imm, ima] = imread([sp '/text1.png']);
    if rssi
        im = image_resize(im, 0.5 * size(im, 1), 0.5 * size(im, 2));
        ima = image_resize(ima, size(im, 1), size(im, 2));
    end
    addlayer(timg, im, ima);
    talp{6} = single((1 / 255) .* double(ima));
    im = imread([sp '/text2.jpg']);
    if rssi
        im = image_resize(im, 0.5 * size(im, 1), 0.5 * size(im, 2));
    end
    addlayer(timg, im, 0);
    talp{7} = single(sqrt(0.5));
    im = imread([sp '/text3.jpg']);
    if rssi
        im = image_resize(im, 0.5 * size(im, 1), 0.5 * size(im, 2));
    end
    addlayer(timg, im, 0);
    talp{8} = single(sqrt(0.5));
    im = imread([sp '/text4.jpg']);
    if rssi
        im = image_resize(im, 0.5 * size(im, 1), 0.5 * size(im, 2));
    end
    addlayer(timg, im, 0);
    talp{9} = single(sqrt(0.5));

    % set alpha information back
    ne_ui.splash.alpha = talp;

    % set handle and render
    sethandle(timg, ihandle);
    render(timg)
    drawnow;

    % make settings to timer object
    set(timr, 'ExecutionMode', 'fixedSpacing', 'Period', 0.05, ...
        'StopFcn', @nesp_cleanup, 'TimerFcn', @nesp_splash, 'UserData', {timg, talp});

catch ne_eo;
    neuroelf_lasterr(ne_eo);

    % try to delete transimg
    try
        delete(timg);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

    % return empty!
    timr = [];
end

% sub functions
function nesp_cleanup(varargin)
global ne_ui;
if ~isfield(ne_ui, 'splash') || ~isstruct(ne_ui.splash) || numel(ne_ui.splash) ~= 1
    return;
end
sp = ne_ui.splash;

% try to set final state
try
    get(sp.ihandle, 'UserData');
    if nargin < 3 || ~islogical(varargin{3}) || numel(varargin{3}) ~= 1 || ~varargin{3}
        setlayeralpha(sp.timg, 2, 0);
        setlayeralpha(sp.timg, 3, 0);
        setlayeralpha(sp.timg, 4, 0);
        setlayeralpha(sp.timg, 5, sp.alpha{5});
        setlayeralpha(sp.timg, 6, 0);
        setlayeralpha(sp.timg, 7, 0);
        setlayeralpha(sp.timg, 8, 0);
        setlayeralpha(sp.timg, 9, sp.alpha{9});
    else
        setlayeralpha(sp.timg, 2, 0);
        setlayeralpha(sp.timg, 3, sp.alpha{3});
        setlayeralpha(sp.timg, 4, 0);
        setlayeralpha(sp.timg, 5, 0);
        setlayeralpha(sp.timg, 6, 0);
        setlayeralpha(sp.timg, 7, sp.alpha{7});
        setlayeralpha(sp.timg, 8, 0);
        setlayeralpha(sp.timg, 9, 0);
    end

    % re-render
    render(sp.timg)
    drawnow;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end
try
    delete(sp.timg);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end
try
    delete(sp.timer);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

% clean up struct
ne_ui.splash = [];

function nesp_splash(varargin)
global ne_ui;
if ~isfield(ne_ui, 'splash') || ~isstruct(ne_ui.splash) || numel(ne_ui.splash) ~= 1 
    try
        stop(varargin{1});
        delete(varargin{1});
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
    return;
else
    try
        get(ne_ui.splash.ihandle, 'UserData');
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        stop(varargin{1});
        delete(varargin{1});
        return;
    end
end
sp = ne_ui.splash;
alp = sp.alpha;

% progress along data
sp.iter = sp.iter + 1;
if sp.iter < 210
    imod = mod(sp.iter, 67);
else
    imod = mod(sp.iter, 50);
end
if imod == 0
    imod = 67;
end
sc = (imod - 50) / 17;

% set back
ne_ui.splash = sp;

% for 1 - 50 and 61 - 110
if imod <= 50

    % recompute hue of curve
    curve = sp.curve;
    curve(:, :, 1) = curve(:, :, 1) + 0.02 .* imod;
    curve(:, :, 1) = curve(:, :, 1) - floor(curve(:, :, 1));
    setlayerpixel(sp.timg, 1, hsvconv(curve));

% for 51 - 67
elseif sp.iter < 75

    % switch text and image
    if sc < 1
        setlayeralpha(sp.timg, 2, (1 - sc) .* alp{2});
        setlayeralpha(sp.timg, 6, (1 - sc) .* alp{6});
    else
        setlayeralpha(sp.timg, 2, 0);
        setlayeralpha(sp.timg, 6, 0);
    end
    setlayeralpha(sp.timg, 3, sc .* alp{3});
    setlayeralpha(sp.timg, 7, sc .* alp{7});

% for 118 - 134
elseif sp.iter < 145

    % switch text and image
    if sc < 1
        setlayeralpha(sp.timg, 3, (1 - sc) .* alp{3});
        setlayeralpha(sp.timg, 7, (1 - sc) .* alp{7});
    else
        setlayeralpha(sp.timg, 3, 0);
        setlayeralpha(sp.timg, 7, 0);
    end
    setlayeralpha(sp.timg, 4, sc .* alp{4});
    setlayeralpha(sp.timg, 8, sc .* alp{8});

% for 185 - 200
elseif sp.iter < 210

    % switch text and image
    if sc < 1
        setlayeralpha(sp.timg, 4, (1 - sc) .* alp{4});
        setlayeralpha(sp.timg, 8, (1 - sc) .* alp{8});
    else
        setlayeralpha(sp.timg, 4, 0);
        setlayeralpha(sp.timg, 8, 0);
    end
    setlayeralpha(sp.timg, 5, sc .* alp{5});
    setlayeralpha(sp.timg, 9, sc .* alp{9});
end

% re-render
render(sp.timg)
drawnow;

% cleanup?
if sp.iter >= 500
    stop(sp.timer);
    return;
end
