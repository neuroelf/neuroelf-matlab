function ne_screenshot(varargin)
% ne_screenshot  - create a screenshot file
%
% FORMAT:       ne_screenshot(SRC, EVT, [, object [, filename [, hq [, us]]]])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       object      can be empty, 0 (main fig), or figure/axes handle
%       filename    optional filename (otherwise requested)
%       hq          high-quality flag (set to 'high-q' for oversampling)
%       us          flag if ne_gcfg.c.ini.Surface.SatelliteUpsampling unset
%
% No output fields.
%
% Example:
%
%    ne_screenshot(0, 0, 'BS123456', 'screenshot_1234.png', 'high-q');
%
%    creates the file 'screenshot_1234.png' with high-quality oversampling
%    from satellite window with ID 'BS123456'

% Version:  v1.1
% Build:    16052718
% Date:     May-27 2016, 6:31 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2015, 2016, Jochen Weber
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
global ne_gcfg;

% make sure GUI is loaded
try
    ch = ne_gcfg.h;
    ch.MainFig.MLHandle;
    ch.MainFig.Pointer;

% if not fall back on no handle and standard pointer value
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    ch = struct('MainFig', struct('MLHandle', [], 'Pointer', 'arrow'));
end

% make sure config is loaded
try
    ci = ne_gcfg.c.ini;
    ci.Surface.ScreenshotOversampling;
    ci.Surface.SatelliteUpsampling;

% if not fall back on standard values (anti-alias, no upsampling)
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    ci = struct('Surface', struct('ScreenshotOversampling', 4, ...
        'SatelliteUpsampling', 1));
    if nargin > 5 && isa(varargin{6}, 'double') && ~isempty(varargin{6}) && ...
       ~isinf(varargin{6}(1)) && ~isnan(varargin{6}(1))
        ci.Surface.SatelliteUpsampling = min(4, max(1, varargin{6}(1)));
    end
end
mfp = ch.MainFig.Pointer;

% valid input
if nargin < 3 || ((~ischar(varargin{3}) || isempty(varargin{3}) || ...
   ~strcmp(varargin{3}(:)', makelabel(varargin{3}(:)')) || ...
   ~isfield(ne_gcfg.cc, varargin{3}(:)')) && ...
   ((~isa(varargin{3}, 'double') && ~isa(varargin{3}, 'matlab.ui.Figure') && ...
     ~isa(varargin{3}, 'matlab.graphics.axis.Axes')) || ...
    numel(varargin{3}) ~= 1 || ...
   ~ishandle(varargin{3}) || ...
   ~any(strcmpi(get(varargin{3}, 'Type'), {'figure', 'axes'}))))
    handle = ch.MainFig.MLHandle;
elseif ischar(varargin{3})
    handle = ne_gcfg.cc.(varargin{3}(:)').Satellite.MLHandle;
else
    handle = varargin{3};
end
try
    isfig = strcmpi(get(handle, 'Type'), 'figure');
    if isfig
        figp = get(handle, 'Pointer');
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end

% output filename
if nargin < 4 || ~ischar(varargin{4}) || isempty(varargin{4}) || ...
    isempty(regexpi(varargin{4}(:)', '\.(bmp|eps|fig|jpe?g|png|tiff?)$')) || ...
   (~isfig && strcmpi(varargin{4}(end-3:end), '.fig'))
    tlist = {'*.png', 'Portable Network Graphics (*.png)'; ...
         '*.bmp', 'MS Windows Bitmap (*.bmp)'; ...
         '*.jpg', 'JPEG image (*.jpg)'; ...
         '*.tif', 'TIFF image (*.tif)'; ...
         '*.eps', 'Encapsulated Post-Script (*.eps)'};
    if isfig
        tlist(end+1, :) = {'*.fig', 'Matlab figure file'};
    end
    [savefile, savepath, saveidx] = uiputfile(tlist, ...
        'Save figure screenshot as...');
    if isequal(savefile, 0) || isequal(savepath, 0) || isempty(savefile)
        return;
    end
    if isempty(savepath)
        savepath = pwd;
    end
    [nullpath, savefile, saveext] = fileparts(savefile);
    if isempty(saveext) || ...
       ~any(strcmpi(saveext, {'.bmp', '.eps', '.fig', '.jpeg', '.jpg', '.png', '.tif', '.tiff'}))
        saveext = {'.png', '.bmp', '.jpg', '.tif', '.eps', '.fig'};
        saveext = saveext{saveidx};
    elseif ~isfig && ((~isempty(saveext) && strcmpi(saveext, '.fig')) || saveidx == 6)
        saveext = '.eps';
    end
    savefile = [savepath, '/', savefile, saveext];
else
    savefile = varargin{4}(:)';
end

% parse type
[savepath, savefile, saveext] = fileparts(strrep(savefile, '\', '/'));
if isempty(savepath)
    savepath = pwd;
end
if any(strcmpi(saveext, {'.jpg', '.jpeg'}))
    jpgo = {'Quality', 90};
    saveext = '.jpg';
else
    jpgo = {};
    saveext = lower(saveext);
end
savefile = strrep([savepath, '/', savefile, saveext], '//', '/');

% for eps, path is different
if strcmp(saveext, '.eps')
    aseps = true;
    style = struct( ...
        'Version',         1, ...
        'Format',          'eps', ...
        'Preview',         'none', ...
        'Width',           'auto', ...
        'Height',          'auto', ...
        'Units',           'inches', ...
        'Color',           'rgb', ...
        'Background',      'w', ...
        'FixedFontSize',   10, ...
        'ScaledFontSize',  'auto', ...
        'FontMode',        'scaled', ...
        'FontSizeMin',     8, ...
        'FixedLineWidth',  1, ...
        'ScaledLineWidth', 'auto', ...
        'LineMode',        'none', ...
        'LineWidthMin',    0.5, ...
        'FontName',        'auto', ...
        'FontWeight',      'auto', ...
        'FontAngle',       'auto', ...
        'FontEncoding',    'latin1', ...
        'PSLevel',         2, ...
        'Renderer',        'auto', ...
        'Resolution',      ne_gcfg.c.ini.MainFig.ExportResolution, ...
        'LineStyleMap',    'none', ...
        'ApplyStyle',      0, ...
        'Bounds',          'loose', ...
        'LockAxes',        'on', ...
        'ShowUI',          'on', ...
        'SeparateText',    'off');

% so is for figures
elseif strcmp(saveext, '.fig')

    % try
    try

        % set pointers
        ch.MainFig.Pointer = 'watch';
        set(handle, 'Pointer', 'watch');
        drawnow;

        % use saveas
        saveas(handle, savefile);

        % some additional work-arounds
        try
            
            % first, did the loading create another figure?
            oldfigs = findobj('type', 'figure');
            figc = load(savefile, '-mat');
            newfigs = findobj('type', 'figure');
            if ~isequal(oldfigs, newfigs)
                
                % if it did, find the figure handle
                for nfc = 1:numel(newfigs)
                    if ~any(oldfigs == newfigs(nfc))
                        
                        % delete new figure without triggering xfigure!
                        xfsd = xfigure(0, 'skipdelete', true);
                        delete(newfigs(nfc));
                        xfigure(0, 'skipdelete', xfsd);
                        break;
                    end
                end
            end
            
            % as a special service, we want to remove hidden elements
            figf = fieldnames(figc);
            if numel(figf) < 1 || ...
                numel(figf) > 2 || ...
                isempty(regexpi(figf{1}, '^hgs'))
                ch.MainFig.Pointer = mfp;
                set(handle, 'Pointer', figp);
                drawnow;
                return;
            end

            % reset closerequestfcn and clear deletefcn
            figc.(figf{1}).properties.CloseRequestFcn = 'closereq';
            if isfield(figc.(figf{1}).properties, 'DeleteFcn')
                figc.(figf{1}).properties = ...
                    rmfield(figc.(figf{1}).properties, 'DeleteFcn');
            end
            if isfield(figc.(figf{1}).properties, 'MenuBar')
                figc.(figf{1}).properties = ...
                    rmfield(figc.(figf{1}).properties, 'MenuBar');
            end
            if isfield(figc.(figf{1}).properties, 'Pointer')
                figc.(figf{1}).properties = ...
                    rmfield(figc.(figf{1}).properties, 'Pointer');
            end
            figc.(figf{1}).properties.Renderer = get(handle, 'Renderer');
            figc.(figf{1}).properties.Resize = 'off';
            if isfield(figc.(figf{1}).properties, 'ResizeFcn')
                figc.(figf{1}).properties = ...
                    rmfield(figc.(figf{1}).properties, 'ResizeFcn');
            end
            figc.(figf{1}).properties.Tag = ...
                sprintf('NE_SAVEDFIG_%08X', floor((2 ^ 32 - 1) * rand(1, 1)));
            if isfield(figc.(figf{1}).properties, 'UserData') && ...
                isstruct(figc.(figf{1}).properties.UserData) && ...
                numel(figc.(figf{1}).properties.UserData) == 1 && ...
                isfield(figc.(figf{1}).properties.UserData, 'xfigure_UDTagStruct')
                figc.(figf{1}).properties.UserData = rmfield( ...
                    figc.(figf{1}).properties.UserData, 'xfigure_UDTagStruct');
            end

            % get children
            figch = figc.(figf{1}).children;
            figck = true(size(figch));

            % now we iterate over children
            for cc = 1:numel(figch)

                % for non-uicontrol and non-uimenu components
                if ~any(strcmpi(figch(cc).type, {'uicontrol', 'uimenu'}))

                    % simply unset the DeleteFcn property
                    if isfield(figch(cc).properties, 'DeleteFcn')
                        figch(cc).properties = ...
                            rmfield(figch(cc).properties, 'DeleteFcn');
                    end
                    if isfield(figch(cc).properties, 'HandleVisibility')
                        figch(cc).properties = ...
                            rmfield(figch(cc).properties, 'HandleVisibility');
                    end
                    figch(cc).properties.Tag = sprintf( ...
                        'NE_SAVEDCLD_%08X', floor((2 ^ 32 - 1) * rand(1, 1)));
                    continue;

                % uicontrols and uimenus are to be deleted entirely
                else
                    figck(cc) = false;
                end
            end

            % only keep children we want with updated properties
            figc.(figf{1}).children = figch(figck);

            % assign in current workspace and resave
            eval([figf{1} '=figc.' figf{1} ';']);
            save([savefile(1:end-4) '_fig.mat'], figf{1});

            % now we need to rename the new file!
            try
                renamefile([savefile(1:end-4) '_fig.mat'], savefile);
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
                try
                    delete([savefile(1:end-4) '_fig.mat']);
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                    ch.MainFig.Pointer = mfp;
                    set(handle, 'Pointer', figp);
                    drawnow;
                    return;
                end
            end

        % no error thrown on this!
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

    % error handling
    catch ne_eo;
        uiwait(warndlg(sprintf('Error saving figure file: %s', ne_eo.message), ...
            'NeuroElf - warning', 'modal'));
    end

    % un-set pointers
    set(handle, 'Pointer', figp);
    ne_gcfg.h.MainFig.Pointer = mfp;

    % return
    return;

else
    aseps = false;
end

% set pointers
if isfig
    set(handle, 'Pointer', 'watch');
end
ne_gcfg.h.MainFig.Pointer = 'watch';
drawnow;

% handle type
htype = lower(get(handle, 'Type'));
switch (htype)

    % for figures
    case {'figure'}

        % get figure screenshot
        pause(0.05);
        drawnow;
        pause(0.05);

        % try
        try
            % EPS?
            if aseps

                % keep background color
                bcol = get(handle, 'Color');
                if isa(bcol, 'double')
                    style.Background = bcol;
                end

                % use hgexport with default style
                hgexport(handle, savefile, style, 'Format', 'eps');

            % image format
            else
                
                % simple get a copy?
                if (ci.Surface.ScreenshotOversampling <= 1 && ...
                    ci.Surface.SatelliteUpsampling <= 1) || ...
                    nargin < 5 || isempty(varargin{3}) || ...
                    isempty(varargin{5}) || ~ischar(varargin{5}) || ...
                   ~strcmpi(varargin{5}(:)', 'high-q')

                    % get frame
                    fr = getframe(handle);
                    
                % a bit more effort
                else
                    
                    % compute rescaling factor
                    rsfactor = ceil(get(0, 'ScreenPixelsPerInch') * ...
                        ci.Surface.SatelliteUpsampling * ci.Surface.ScreenshotOversampling);
                    
                    % write into tempfile in hi-res
                    tempfile = [tempname '.png'];
                    ppm = get(handle, 'PaperPositionMode');
                    ihc = get(handle,'InvertHardcopy');
                    set(handle, 'PaperPositionMode', 'auto');
                    set(handle, 'InvertHardcopy', 'off');
                    print(handle, sprintf('-r%d', rsfactor), '-dpng', tempfile);
                    set(handle, 'InvertHardcopy', ihc);
                    set(handle, 'PaperPositionMode', ppm);
                    fr = struct('cdata', imread(tempfile));
                    delete(tempfile);
                    
                    % re-set handles' configs
                    hu = get(handle, 'Units');
                    set(handle, 'Units', 'pixels');
                    rsize = floor(ci.Surface.SatelliteUpsampling .* get(handle, 'Position'));
                    set(handle, 'Units', hu);

                    % remove resampling edges
                    frsize = size(fr.cdata);
                    frback = median(double(cat(1, ...
                        reshape(fr.cdata([4:5, end-4:end-3], :, :), 4 * frsize(2), frsize(3)), ...
                        reshape(fr.cdata(:, [4:5, end-4:end-3], :), 4 * frsize(1), frsize(3)))), 1);
                    for cpc = 1:numel(frback)
                        fr.cdata([1:3, end-2:end], :, cpc) = frback(cpc);
                        fr.cdata(:, [1:3, end-2:end], cpc) = frback(cpc);
                    end

                    % then resize
                    fr.cdata = image_resize(fr.cdata, rsize(4), rsize(3));
                end

                % and save
                imwrite(fr.cdata, savefile, jpgo{:});
            end

        % error handling
        catch ne_eo;
            uiwait(warndlg(sprintf('Error saving image file: %s', ne_eo.message), ...
                'NeuroElf - warning', 'modal'));
        end

    % axes
    case {'axes'}

        % if not EPS
        if ~aseps

            % try
            try
                % get frame
                fr = getframe(handle);

                % and save
                imwrite(fr.cdata, savefile, jpgo{:});

            % error handling
            catch ne_eo;
                uiwait(warndlg(sprintf('Error saving image file: %s', ne_eo.message), ...
                    'NeuroElf - warning', 'modal'));
            end

            % un-set pointers
            if isfig
                set(handle, 'Pointer', figp);
            end
            ne_gcfg.h.MainFig.Pointer = mfp;
            drawnow;

            % return
            return;
        end

        % as EPS! -> get own, parent, and ROOT settings
        sets = get(handle);
        phandle = sets.Parent;
        psets = get(phandle);
        rsets = get(0);

        % set units to pixels and get position within parent
        set(phandle, 'Units', 'pixels');
        set(handle, 'Units', 'pixels');
        set(0, 'Units', 'pixels');
        fpos = get(phandle, 'Position');
        hpos = get(handle, 'Position');

        % create new figure with settings
        tfig = -1;
        try
            tfig = figure;
            set(tfig, 'Colormap', psets.Colormap, 'Resize', 'off', 'Visible', 'off');
            set(tfig, 'Position', [fpos(1:2) hpos(3:4)]);
            set(tfig, 'Units', 'normalized');

            % then change axes parent
            set(handle, 'Parent', tfig);
            set(handle, 'Units', 'normalized');
            set(handle, 'Position', [0, 0, 1, 1]);

            % force draw figure
            drawnow;

            % export
            hgexport(tfig, savefile, style, 'Format', 'eps');
            pause(0.05);
        catch ne_eo;
            uiwait(warndlg(sprintf('An error occurred during ne_screenshot: %s', ...
                ne_eo.message), 'modal'));
        end

        % reset everything
        set(handle, 'Parent', phandle);
        set(handle, 'Units', 'pixels');
        set(handle, 'Position', hpos);
        set(handle, 'Units', sets.Units);
        set(phandle, 'Position', fpos);
        set(0, 'Units', rsets.Units);
        set(phandle, 'Units', psets.Units);
        if double(tfig) > 0 && ...
            ishandle(tfig)
            delete(tfig);
        end
end

% final drawnow
ch.MainFig.Pointer = mfp;
if isfig
    set(handle, 'Pointer', figp);
end
drawnow;
