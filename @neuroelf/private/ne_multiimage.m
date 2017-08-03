function varargout = ne_multiimage(varargin)
% ne_multiimage  - open a "multi-image" browser
%
% FORMAT:       [hSat, tags, iSat] = ne_multiimage(SRC, EVT)
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%
% Output fields:
%
%       hSat        xfigure handle of satellite window
%       tags        1x1 structure with further handles and settings
%       iSat        string representing the satellite ID
%
% Example:
%
%       [mi_viewer, ~, mi_id] = ne_multiimage(0, 0);

% Version:  v1.0
% Build:    14121715
% Date:     Dec-17 2014, 3:27 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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
global ne_gcfg ne_methods;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% echo
if ne_gcfg.c.echo
    ne_echo('neuroelf_gui(''multiimage'');');
end

% try (required for double figure)
try

    % open satellite window
    hSat = [];
    hMix = [];
    hSat = xfigure([neuroelf_path('tfg') '/ne_imagesat.tfg']);
    hMix = xfigure([neuroelf_path('tfg') '/ne_imagemixer.tfg']);
    tSat = hSat.Tag;
    iSat = tSat(1:8);
    tMix = hMix.Tag;
    iMix = tMix(1:8);

    % create tags
    tSat = hSat.TagStruct;
    tMix = hMix.TagStruct;
    htag = struct;
    htag.Satellite = hSat;
    htag.SatelliteMLH = hSat.MLHandle;
    htag.Mixer = hMix;
    htag.MixerMLH = hMix.MLHandle;
catch ne_eo;
    if ~isempty(hMix)
        hMix.Delete;
    end
    if ~isempty(hSat)
        hSat.Delete;
    end
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg(sprintf('Error opening multi-image viewer/mixer: %s.', ...
        ne_eo.message), 'NeuroElf - error', 'modal'));
    return;
end

% add crosshair lines to images axes objects -> SAG, COR, TRA, Zoom
chax = tSat.(sprintf('AX_%s_Image', iSat)).MLHandle;
set(chax, 'Units', 'pixels');
htag.Image = tSat.(sprintf('IM_%s_Image', iSat)).Children;
htag.ImageAxes = chax;
htag.LineX = line([0; 0.999], [0.5; 0.5], 'Color', [1, 1, 1], 'Parent', chax);
htag.LineY = line([0.5; 0.5], [0.001; 0.999], 'Color', [1, 1, 1], 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');
tSat.(sprintf('IM_%s_Image', iSat)).Visible = 'on';

% callbacks
tSat.(sprintf('UIM_%s_LoadImage', iSat)).Callback = {@ne_multiimage_load, iSat};
hSat.CloseRequestFcn = {@ne_closesatwindow, iSat};
hMix.CloseRequestFcn = {@ne_closesatwindow, iSat};
hSat.KeyPressFcn = {@ne_multiimage_keypress, iSat};
hSat.KeyReleaseFcn = @ne_multiimage_keyrelease;
hSat.WindowButtonDownFcn = @ne_multiimage_btdown;
hSat.WindowButtonMotionFcn = @ne_multiimage_btmove;
hSat.WindowButtonUpFcn = {@ne_multiimage_btup, iSat};

% set configuration
htag.Config = struct( ...
    'flexip',  ne_methods.flexinterpn_method, ...
    'file',    '', ...
    'imobj',   [], ...
    'imrsz',   ne_methods.image_resize, ...
    'interpn', 'nearest', ...
    'isize',   [0, 0], ...
    'mixtag',  iMix, ...
    'mixtags', tMix, ...
    'sattag',  iSat, ...
    'sattags', tSat, ...
    'sattype', 'multiimage', ...
    'size',    [512, 512], ...
    'ti',      transimg(512, 512));
sethandle(htag.Config.ti, htag.Image);

% finally, make satellite visible
hSat.HandleVisibility = 'callback';
hSat.Visible = 'on';
hMix.HandleVisibility = 'callback';
hMix.Visible = 'on';

% allow resize
set(hSat.MLHandle, 'Resize', 'on');
hSat.ResizeFcn = {@ne_multiimage_resize, iSat};

% store
ne_gcfg.cc.(iSat) = htag;

% for output
if nargout > 0
    varargout{1} = hSat;
    if nargout > 1
        varargout{2} = htag;
        if nargout > 2
            varargout{3} = iSat;
        end
    end
end



function varargout = ne_multiimage_load(varargin)
global ne_gcfg;
varargout = cell(1, nargout);
cc = ne_gcfg.cc.(varargin{3});

% clear image if object
if isxff(cc.Config.imobj, true)
    imobj = cc.Config.imobj;
    if ~isempty(imobj.FilenameOnDisk) && ...
        isfield(imobj.RunTimeVars, 'AutoSave') && ...
        islogical(imobj.RunTimeVars.AutoSave) && ...
        numel(imobj.RunTimeVars.AutoSave) == 1 && ...
        imobj.RunTimeVars.AutoSave
        cc.Config.imobj.SaveRunTimeVars;
    end
    cc.Config.imobj.ClearObject;
end

% clear data
cc.Config.file  = '';
cc.Config.imobj = [];

% re-instantiate transio
delete(cc.Config.ti);
cc.Config.ti = transimg(cc.Config.size(1), cc.Config.size(2));
sethandle(cc.Config.ti, cc.Image);
render(cc.Config.ti);

% try to load data
[imfile, impath] = uigetfile( ...
   {'*.bmp;*.jpeg;*.jpg;*.png;*.tif', 'Regular images (*.bmp, *.jpeg, *.jpg, *.png, *.tif)'; ...
    '*.czi', 'Carl-Zeiss images (*.czi)'}, 'Please select an image to load...');
if ~isequal(imfile, 0) && ...
   ~isequal(impath, 0) && ...
    ischar(imfile) && ...
   ~isempty(imfile)
    if isempty(impath)
        impath = pwd;
    end
    imfile = [impath '/' imfile];
    imobj = [];
    try
        if exist(imfile, 'file') == 2
            [impath, imfile, imext] = fileparts(imfile);
            imfile = [impath '/' imfile imext];
            switch lower(imext)
                case {'.bmp', '.jpeg', '.jpg', '.png', '.tif'}
                    imobj = imread(imfile);
                    cc.Config.size = [size(imobj, 2), size(imobj, 1)];
                    delete(cc.Config.ti);
                    cc.Config.ti = transimg(cc.Config.size(1), cc.Config.size(2));
                    sethandle(cc.Config.ti, cc.Image);
                case {'.czi'}
                    imobj = xff(imfile);
            end
            cc.Config.file = imfile;
            cc.Config.imobj = imobj;
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        uiwait(warndlg(sprintf('Error loading image: %s.', ne_eo.message), ...
            'NeuroElf - warning', 'modal'));
    end
end

% update
ne_gcfg.cc.(varargin{3}) = cc;
ne_multiimage_update(varargin{:});



function varargout = ne_multiimage_update(varargin)
global ne_gcfg;
varargout = cell(1, nargout);
cc = ne_gcfg.cc.(varargin{3});

% display image
imobj = cc.Config.imobj;
if isnumeric(imobj)
    simobj = size(imobj);
    if simobj(1) ~= cc.Config.size(2) || ...
        simobj(2) ~= cc.Config.size(1)
        cimobj = class(imobj);
        if ~isa(imobj, 'double')
            cimobj = eval(['@' cimobj]);
        else
            cimobj = [];
        end
        k = cc.Config.flexip('kernels');
        if strcmpi(cc.Config.interpn, 'resample')
            imobj = cc.Config.imrsz(imobj, cc.Config.size(2), cc.Config.size(1));
        elseif ~strcmpi(cc.Config.interpn, 'nearest') && ...
            isfield(k, cc.Config.interpn)
            k = k.(cc.Config.interpn);
            imi = [Inf, Inf; 1, 1; (simobj(1) - 1) / (cc.Config.size(2)-1), (simobj(2) - 1) / (cc.Config.size(1)-1); simobj(1:2)];
            if numel(simobj) == 2
                imk3 = {};
                imks3 = {};
            else
                imi(:, 3) = [Inf; 1; 1; simobj(3)];
                imk3 = {[0; 1; 0]};
                imks3 = {1};
            end
            imobj = cc.Config.flexip(cc.Config.flexip(imobj, ...
                imi, [{k, [0;1;0]}, imk3], [{4096, 1}, imks3], 0), ...
                imi, [{[0;1;0], k}, imk3], [{1, 4096}, imks3], 0);
        else
            imc1 = 1 + round(((simobj(1) - 1) / (cc.Config.size(2)-1)) .* (0:(cc.Config.size(2)-1)));
            imc2 = 1 + round(((simobj(2) - 1) / (cc.Config.size(1)-1)) .* (0:(cc.Config.size(1)-1)));
            imobj = imobj(imc1, imc2, :);
        end
        if ~isempty(cimobj)
            imobj = cimobj(imobj);
        end
    end
    setlayer(cc.Config.ti, 1, imobj);
else
end
display(render(cc.Config.ti));



function varargout = ne_multiimage_btdown(varargin)
varargout = cell(1, nargout);

function varargout = ne_multiimage_btmove(varargin)
varargout = cell(1, nargout);

function varargout = ne_multiimage_btup(varargin)
varargout = cell(1, nargout);



function varargout = ne_multiimage_keypress(varargin)
varargout = cell(1, nargout);

function varargout = ne_multiimage_keyrelease(varargin)
varargout = cell(1, nargout);



function varargout = ne_multiimage_resize(varargin)
global ne_gcfg;
varargout = cell(1, nargout);
cc = ne_gcfg.cc.(varargin{3});
