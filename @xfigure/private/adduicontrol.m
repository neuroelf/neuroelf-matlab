function o = adduicontrol(xo, iStr)
%XFIGURE::ADDUICONTROL  Add a uicontrol element to a figure.
%   UIC = ADDUICONTROL(FIG, UICSTRUCT) adds a uicontrol to figure FIG with
%   the settings in UICSTRUCT.

% Version:  v1.1
% Build:    16042214
% Date:     Apr-22 2016, 2:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% global references and storage
global ne_methods xfigmlup xfigsngl xfigures;

% only valid for figure and uipanel objects
if numel(xo) ~= 1 || xo.T ~= 1
    error('neuroelf:xfigure:invalidObjectType', ...
        'Only figures can be parents of uicontrols.');
end

% test iStr
if ~isstruct(iStr) || numel(iStr) ~= 1
    error('neuroelf:xfigure:badPropertyStruct', ...
        'Uicontrol properties must be of type struct.');
end

% perform some checks on struct
iStr = ne_methods.checkstruct(iStr, xfigsngl.optuic);
if isempty(iStr.Type) || ~isfield(xfigsngl.uictypes, lower(iStr.Type)) || ...
   ~isnumeric(iStr.Position) || isempty(iStr.Position) || numel(iStr.Position) ~= 4
    error('neuroelf:xfigure:badPropertyStruct', ...
        'Bad uicontrol property struct supplied.');
end

% use tag from iStr?
if ~isempty(iStr.Tag) && numel(iStr.Tag) < 28 && isrealvarname(iStr.Tag(:)')
    utag = iStr.Tag(:)';
else
    iStr.Tag = sprintf('UIC_%010.0f', floor(1e10 * rand(1, 1)));
    utag = ['xfigure_' iStr.Tag];
end

% shortcuts to some important settings
iCTyp = lower(iStr.Type);
iCap = iStr.Caption;
iPar = xo.H;
iPos = iStr.Position;
iVis = iStr.Visible;

% create element (object)
o = xfigure('new');

% preset callbacks array
iStr.Callbacks = cell(1, 4);
cbclick = iStr.CallbackClick(:)';

% special type
iRTyp = xfigsngl.uictypes.(iCTyp);
if strcmp(iRTyp, 'BUILTIN')
    iRSpec = true;
else
    iRSpec = false;
end

% background color
hasBGColor = true;
if numel(iStr.ColorBG) == 1
    iStr.ColorBG(1:3) = iStr.ColorBG;
elseif numel(iStr.ColorBG) ~= 3
    if ~any(strcmpi(iCTyp, {'label', 'radiobutton'}))
        hasBGColor = false;
    end
    iStr.ColorBG = get(iPar, 'Color');
end
iCBG = max(0, min(1, iStr.ColorBG));

% foreground color
hasFGColor = true;
if numel(iStr.ColorFG) == 1
    iStr.ColorFG(1:3) = iStr.ColorFG;
elseif numel(iStr.ColorFG) ~= 3
    hasFGColor = false;
    iStr.ColorFG = [0 0 0];
end
iCFG = max(0, min(1, iStr.ColorFG));

% fontsize
if ischar(iStr.FontSize) && ~isempty(iStr.FontSize) && ...
    isfield(xfigsngl.fntsizes, iStr.FontSize(:)')
    iStr.FontSize = xfigsngl.fntsizes.(iStr.FontSize(:)');
elseif ~isnumeric(iStr.FontSize) || numel(iStr.FontSize) ~= 1
    if xfigsngl.fontfactuse
        iStr.FontSize = 10 * xfigsngl.fontfact;
    else
        iStr.FontSize = 10;
    end
elseif xfigsngl.fontfactuse
    iStr.FontSize = iStr.FontSize * xfigsngl.fontfact;
end

% sliderstep
if ~any(strcmp(iCTyp, {'slider', 'xprogress'})) || numel(iStr.SliderStep) < 2
    iStr.SliderStep = [];
else
    iStr.SliderStep = iStr.SliderStep(1:2);
    if strcmp(iCTyp, 'slider')
        if isempty(iStr.Value)
            if ~isempty(iStr.MinMaxTop)
                iStr.Value = iStr.MinMaxTop(1);
            else
                iStr.Value = 1;
            end
        elseif ~isempty(iStr.MinMaxTop)
            iStr.Value = min(iStr.MinMaxTop(2), max(iStr.MinMaxTop(1), iStr.Value(1)));
        end
    end
end

% special type requested ?
xchildren = [];
if iRSpec

    % set default axes struct
    oStr = struct('Parent', iPar, 'Units', iStr.Units, 'Position', iPos, ...
        'ButtonDownFcn', cbclick, 'Color', 'none', 'Tag', utag, ...
        'UserData', iStr.UserData, 'Visible', iVis);

    % what type
    switch iCTyp

        % axes
        case 'xaxes'

            % axes options are in caption
            if ischar(iCap) && ~isempty(iCap)
                [myaxopts, myonum] = ne_methods.splittocell(iCap, ',');

                % remove last if impair number
                if rem(myonum, 2)
                    myaxopts(end) = [];
                    myonum = myonum - 1;
                end

                % basic test on options
                for myon = 1:2:myonum

                    % option name in hyphens
                    if ~isempty(myaxopts{myon}) && myaxopts{myon}(1) == '''' && myaxopts{myon}(end) == ''''
                        myaxopts{myon} = myaxopts{myon}(2:end-1);
                    end

                    % option setting with value evaluation
                    try
                        oStr.(myaxopts{myon}) = eval(myaxopts{myon + 1});
                    catch xfigerror
                        neuroelf_lasterr(xfigerror);
                        break;
                    end
                end
            end

            % create axes object
            hOut = axes(oStr);

        % button bar
        case 'xbutton'

            % create axes object and image
            hOut = axes(oStr);

            % set position information
            set(hOut, 'XLim', 0.5 + [0, iPos(3)], 'YLim', 0.5 + [0, iPos(4)], ...
                'XLimMode', 'manual', 'YLimMode', 'manual');

            % reset axes visibility to off
            set(hOut, 'Visible', 'off');

            % keep track of some settings
            iStr.ButtonList = {};
            iStr.ButtonMap = zeros(0, 8); % 1x4 position, brd, ena, vis, type

        % images
        case 'ximage'

            % no useful argument
            if isempty(iCap)
                error('neuroelf:xfigure:badImageContent', ...
                    'Images either need a filename or a binary content.');
            end

            % filename
            if ischar(iCap)

                % try to read image
                try
                    imgpdata = imread(iCap(:)');
                    if isempty(imgpdata)
                        error('INVALID_IMAGE');
                    end
                catch ne_eo;
                    error( ...
                        'xfigure:BadImageFile', ...
                        'The file (%s) is not readable or corrupt (%s).', ...
                        iCap(:)', ne_eo.message ...
                    );
                end

            % binary image data
            elseif isnumeric(iCap)
                imgpdata = iCap;
                iStr.Caption = '';

                % only 2-D -> grayscale
                if ndims(imgpdata) < 3
                    imgpdata(:, :, 2:3) = imgpdata(:, :, [1, 1]);
                end

            % invalid content
            else
                error( ...
                    'xfigure:BadImageContent', ...
                    'Images either need a filename or a binary content.' ...
                );
            end

            % apply shading if requested (watermarks, etc.)
            if hasBGColor && hasFGColor
                imgpdata(:, :, 1) = ...
                    uint8(floor(double(imgpdata(:, :, 1)) * iCFG(1) + ...
                    255 * iCBG(1) * (1 - iCFG(1))));
                imgpdata(:, :, 2) = ...
                    uint8(floor(double(imgpdata(:, :, 2)) * iCFG(2) + ...
                    255 * iCBG(2) * (1 - iCFG(2))));
                imgpdata(:, :, 3) = ...
                    uint8(floor(double(imgpdata(:, :, 3)) * iCFG(3) + ...
                    255 * iCBG(3) * (1 - iCFG(3))));

            % apply frame
            elseif hasFGColor
                imgodata = imgpdata;
                imsz = size(imgodata);
                imgpdata = repmat(reshape(uint8(round(127 + 128 .* iCBG)), [1, 1, 3]), iPos(4), iPos(3));
                uiCFG = uint8(round(255 .* iCFG));
                ifr = round(0.5 .* (imsz(1:2) - iPos([4, 3])));
                if ifr(1) < 0
                    fis1 = 1;
                    tis1 = imsz(1);
                    fit1 = 1 - ifr(1);
                    tit1 = fit1 + tis1 - 1;
                else
                    fis1 = 1 + ifr(1);
                    tis1 = fis1 + iPos(4) - 1;
                    fit1 = 1;
                    tit1 = iPos(4);
                end
                if ifr(2) < 0
                    fis2 = 1;
                    tis2 = imsz(2);
                    fit2 = 1- ifr(2);
                    tit2 = fit2 + tis2 - 1;
                else
                    fis2 = 1 + ifr(2);
                    tis2 = fis2 + iPos(3) - 1;
                    fit2 = 1;
                    tit2 = iPos(3);
                end
                imgpdata([1, end], :, 1) = uiCFG(1);
                imgpdata(:, [1, end], 1) = uiCFG(1);
                imgpdata([1, end], :, 2) = uiCFG(2);
                imgpdata(:, [1, end], 2) = uiCFG(2);
                imgpdata([1, end], :, 3) = uiCFG(3);
                imgpdata(:, [1, end], 3) = uiCFG(3);
                imgpdata(fit1:tit1, fit2:tit2, :) = imgodata(fis1:tis1, fis2:tis2, :);
            end

            % create axes object and image
            hOut = axes(oStr);
            xchildren = image(imgpdata, 'Parent', hOut);
            
            % set callback?
            iStr.ImageData = imgpdata;
            if isempty(cbclick) && ~isempty(iStr.Callback)
                cbclick = iStr.Callback;
            end
            if ~isempty(cbclick)
                set(xchildren, 'ButtonDownFcn', {@xclick, o});
                iStr.XClick = cbclick;
            end

            % reset axes visibility to off
            set(hOut, 'Visible', 'off');
            set(xchildren, 'Visible', iVis, 'UserData', o);

        % extended labels (with tex func :)
        case {'xlabel', 'xlink'}

            % replace empty labels with 'tex:' label
            if isempty(iCap)
                iCap = 'tex:';
            end

            % create axes object and text child
            hOut = axes(oStr);
            myct = struct('Parent', hOut);
            myct.Units           = 'normalized';
            myct.Position        = [0.0 0.0];
            myct.HorizontalAlign = iStr.HAlign;
            myct.VerticalAlign   = iStr.VAlign;
            myct.FontAngle       = iStr.FontItalic;
            myct.FontName        = iStr.FontName;
            myct.FontSize        = iStr.FontSize;
            myct.FontWeight      = iStr.FontWeight;
            myct.FontUnits       = 'points';
            if ~isempty(iStr.Rotation)
                myct.Rotation = iStr.Rotation(1);
            end
            myct.UserData = o;

            % interpreter
            if numel(iCap) > 3 && strcmpi(iCap(1:4), 'tex:')
                myct.Interpreter = 'tex';
                iCap(1:4) = [];
            else
                myct.Interpreter = 'none';
            end

            % y - wise rotation
            if numel(iCap) > 1 && strcmpi(iCap(1:2), 'y:')
                myct.Rotation = 90;
                iCap(1:2) = [];
            end

            % special "align" treatment
            switch lower(iStr.HAlign)
                case 'center'
                    myct.Position(1) = 0.5;
                case 'right'
                    myct.Position(1) = 1.0;
            end
            switch lower(iStr.VAlign)
                case 'middle'
                    myct.Position(2) = 0.5;
                case 'top'
                    myct.Position(2) = 1.0;
            end

            % do we deal with a link ?
            LinkTarget = [];
            if strcmp(iCTyp, 'xlink')
                CapParts = ne_methods.splittocell(iCap, '|', 1);
                iCap = CapParts{1};
                if numel(CapParts) > 1
                    LinkTarget = CapParts{2};
                end
            end

            % finally, set Position
            xchildren = text(myct.Position(1), myct.Position(2), iCap, myct);
            set(hOut, 'Visible', 'off');
            set(xchildren, 'Visible', iVis);
            iStr.Caption = iCap;
            oStr.Label = myct;

            % set link handler
            if ~isempty(LinkTarget) && ~isempty(xfigsngl.linkhandler)
                set(xchildren, 'ButtonDownFcn', ['!' xfigsngl.linkhandler{1} ...
                    deblank(LinkTarget) xfigsngl.linkhandler{2}]);
                try
                    set(xchildren, 'Color', 'red');
                catch xfigerror
                    neuroelf_lasterr(xfigerror);
                end
            end

            % keep note of children and color
            oStr.xcolor = get(xchildren, 'Color');

        % progress bars
        case {'xprogress'}

            % set correct empty caption
            if isempty (iCap)
                iCap = 'x: ';
            end

            % set progress bar range
            if numel(iStr.MinMaxTop) > 2 && ~any(isnan(iStr.MinMaxTop) || ...
                isinf(iStr.MinMaxTop))
                iStr.MinMaxTop = iStr.MinMaxTop(1:3);
            else
                iStr.MinMaxTop = [0, 1, 0];
            end

            % create axes object
            hOut = axes(oStr);
            set(hOut, 'Units', 'normalized');

            % the problem with images is that they do not have a position
            % property for the parent axes object, so we need TWO axes objects
            % so as to make this work "properly"...

            % with no sliderstep value make "flat" bar with rectangle
            if isempty(iStr.SliderStep) || iStr.SliderStep(1) ~= 1
                iStr.ProgType = 'flat';

                % make bar from filling rectangle and surrounding line
                if xfigsngl.mlversion < 804
                    xchildren = rectangle('Position', [0, 0, eps, eps], ...
                        'FaceColor', iCFG, 'EdgeColor', iCFG, ...
                        'EraseMode', 'background', 'UserData', o, 'Visible', iVis);
                else
                    xchildren = rectangle('Position', [0, 0, eps, eps], ...
                        'FaceColor', iCFG, 'EdgeColor', iCFG, ...
                        'UserData', o, 'Visible', iVis);
                end

            % make "round" progress bars with graded colour image
            else
                iStr.ProgType = 'round';
                try
                    iStr.Value(end + 1) = 0;
                catch xfigerror
                    neuroelf_lasterr(xfigerror);
                end

                % how many colors
                numlines = max(xfigsngl.progbarface, iStr.SliderStep(2));
                numcols = fix((numlines + 1) / 2);

                % put together bar image
                imgpdata(1, numlines, 3) = uint8(0);
                for colc = 1:numcols
                    colb = (numcols - colc) / numcols;
                    colf = colc / numcols;
                    ridx = [colc (1 + numlines - colc)];
                    imgpdata(1, ridx, 1) = ...
                        uint8(fix(255 * (colb * iCBG(1) + colf * iCFG(1)) + 0.5));
                    imgpdata(1, ridx, 2) = ...
                        uint8(fix(255 * (colb * iCBG(2) + colf * iCFG(2)) + 0.5));
                    imgpdata(1, ridx, 3) = ...
                        uint8(fix(255 * (colb * iCBG(3) + colf * iCFG(3)) + 0.5));
                end

                % if direction is Y-axes, reshape image to Yx1
                if numel(iCap) < 1 || lower(iCap(1)) ~= 'y'
                    imgpdata = reshape(imgpdata, [numlines, 1, 3]);
                end

                % add image to "first" axes
                xchildren = image(imgpdata, 'Parent', hOut);
                set(xchildren, 'Visible', iVis, 'UserData', o);
                set(hOut, 'Tag', ['PBX_' oStr.Tag], 'UserData', o, 'Visible', 'off');

                % create "second" axes object for outline and caption
                hOut = axes(oStr);
                set(hOut, 'Units', 'normalized');
            end
            if xfigsngl.mlversion < 804
                xchildren = [xchildren, line([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], ...
                    'Parent', hOut, 'EraseMode', 'none', 'Color', iCBG * 0.75, ...
                    'UserData', o, 'Visible', iVis)];
            else
                xchildren = [xchildren, line([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], ...
                    'Parent', hOut, 'Color', iCBG * 0.75, 'UserData', o, ...
                    'Visible', iVis)];
            end

            % prepare Caption
            if numel(iCap) > 1 && iCap(2) == ':'
                ixCap = {iCap(1), iCap(3:end)};
            else
                ixCap = {'x', iCap};
            end
            myct = struct('Parent', hOut);
            myct.Units = 'normalized';
            myct.Position = [0.5 0.5];
            myct.HorizontalAlign = 'center';
            myct.VerticalAlign = 'middle';
            myct.FontAngle = iStr.FontItalic;
            myct.FontName = iStr.FontName;
            myct.FontSize = iStr.FontSize;
            myct.FontUnits = 'points';
            myct.FontWeight = iStr.FontWeight;
            myct.UserData = o;
            myct.Visible = iVis;
            if prod(iCFG) < 0.125 && all(iCFG < 0.25)
                itColor = [0.9325, 0.9325, 0.9325];
            else
                itColor = iCBG * 0.125;
            end
            myct.Color = itColor;

            % if tex interpreter is requested (prefix: "tex:") use it!
            if numel(ixCap{2}) > 3 && strcmpi(ixCap{2}(1:4), 'tex:')
                myct.Interpreter = 'tex';
                ixCap{2} = ixCap{2}(5:end);
            else
                myct.Interpreter = 'none';
            end

            % generate text object
            hTxt = text(myct.Position(1), myct.Position(2), ixCap{2}, myct);
            if ~isempty(ixCap{1}) && lower(ixCap{1}) =='y'
                set(hTxt, 'Rotation', 90);
                iStr.ProgDir = 'y';
            else
                iStr.ProgDir = 'x';
            end
            iStr.Caption = ixCap{2};

            % get size, set Units back to original units and init progress
            set(hOut, 'Units', iStr.Units, 'Visible', 'off');
            try
                oStr.progress = iStr.Value(end);
            catch xfigerror
                neuroelf_lasterr(xfigerror);
                oStr.progress = 0;
            end
            xchildren(end + 1) = hTxt;

    end

    % keep a note on what type was used
    iStr.xtype = iCTyp;

% "only" a MATLAB uicontrol type
else

    % deal with radiobutton / radiogroup click events
    iStr.Callbacks{1} = iStr.Callback;
    if strcmp(iRTyp, 'radiobutton') && ~isempty(iStr.RGroup)
        iStr.Callback = 'rgroupclick(xfigure(gcbo));';

    % otherwise standard behaviour
    else
        iStr.Callbacks = {iStr.Callback, '', iStr.CallbackDblClick};
        iStr.Callback = 'doubleclick(xfigure(gcbo));';
    end

    % check MinMaxTop according to edit/multiedit type
    if strcmp(iCTyp, 'edit')
        iStr.MinMaxTop = [0, 0, 0];
    elseif strcmp(iCTyp, 'multiedit')
        iStr.MinMaxTop = [0, 2, 0];
    elseif numel(iStr.MinMaxTop) > 2
        iStr.MinMaxTop = iStr.MinMaxTop(1:3);
    else
        iStr.MinMaxTop = [0, 1, 1];
    end

    % start building struct for MATLAB uicontrol() call
    oStr = struct( ...
        'Parent',             iPar, ...
        'Units',              iStr.Units, ...
        'Position',           iPos, ...
        'ButtonDownFcn',      cbclick, ...
        'DeleteFcn',          '', ...
        'Enable',             iStr.Enabled, ...
        'FontAngle',          iStr.FontItalic, ...
        'FontName',           iStr.FontName, ...
        'FontSize',           iStr.FontSize, ...
        'FontWeight',         iStr.FontWeight, ...
        'ForegroundColor',    iCFG, ...
        'HorizontalAlign',    iStr.HAlign, ...
        'Interruptible',      iStr.Interrupts, ...
        'ListboxTop',         iStr.MinMaxTop(3), ...
        'Max',                iStr.MinMaxTop(2), ...
        'Min',                iStr.MinMaxTop(1), ...
        'SelectionHighlight', iStr.Selectable, ...
        'Style',              iRTyp, ...
        'Tag',                utag, ...
        'TooltipString',      iStr.ToolTip, ...
        'UserData',           iStr.UserData, ...
        'Visible',            iVis ...
    );
    oStr.Callback = iStr.Callback;
    if ~isempty(iStr.Value)
        oStr.Value = iStr.Value;
    end

    % background color specified ?
    if hasBGColor
        oStr.BackgroundColor = iCBG;
    end

    % make images loadable for button and toggle UIControls
    if any(strcmp(iCTyp, {'button', 'toggle'})) && ~isempty(iStr.Caption)

        % caption = filename?
        if ischar(iStr.Caption) && any(iStr.Caption(:)' == '.') && ...
            exist(iStr.Caption(:)', 'file') == 2
            try
                oStr.CData = double(imread(iStr.Caption(:)')) / 256;
                iStr.ImageFile = iStr.Caption(:)';
                iStr.Caption = '';
            catch xfigerror
                warning('neuroelf:xfigure:badImageFile', ...
                    'Not a valid image file: %s (%s).', iStr.Caption(:)', xfigerror.message);
            end

        % binary caption = image
        elseif isnumeric(iStr.Caption) && ~isempty(iStr.Caption) && ndims(iStr.Caption) == 3
            oStr.CData = double(iStr.Caption);
            if any(oStr.CData(:) > 1)
                oStr.CData = oStr.CData ./ 256;
            end
            iStr.ImageData = iStr.Caption;
            iStr.Caption = '';

        % otherwise copy caption
        else
            oStr.String = iStr.Caption;
        end

    % otherwise simply copy caption
    else
        oStr.String = iStr.Caption;
    end

    % set style and do your job
    hOut = uicontrol(oStr);
    iStr.xtype = '';
end

% update UIControl unique counter and make sure we're deleted correctly
o.H = hOut;
o.T = 2;
o.X = makeostruct(xfigsngl.objtypes.uicontrol);
oStr.MLHandle = hOut;

% uicontext menu
if ~isempty(iStr.ContextMenu)
    uicm = findobj('Type', 'uicontextmenu', 'Tag', iStr.ContextMenu);
    if ~isempty(uicm)
        try
            set(hOut, 'UIContextMenu', uicm(1));
        catch xfigerror
            neuroelf_lasterr(xfigerror);
        end
    end
end

% any Enable-groups ?
iGroups = ne_methods.splittocellc(iStr.EGroups, ',; ', true, true);
for iGroupC = numel(iGroups):-1:1
    iGroup = ne_methods.ddeblank(iGroups{iGroupC});
    if isrealvarname(iGroup)
        if isfield(xo.X.figprops.egroups, iGroup)
            xo.X.figprops.egroups.(iGroup)(end + 1) = o;
        else
            xo.X.figprops.egroups.(iGroup) = o;
        end
    else
        iGroups(iGroupC) = [];
    end
end
iStr.EGroups = iGroups;

% any RadioGroup ?
iGroup = iStr.RGroup;
if ~isempty(iGroup)
    if isfield(xo.X.figprops.rgroups, iGroup)
        xo.X.figprops.rgroups.(iGroup)(end + 1) = o;
    else
        xo.X.figprops.rgroups.(iGroup) = o;
    end
end

% any SlideXY-group ?
iGroup = iStr.SGroups;
if ~isempty(iGroup)
    if isfield(xo.X.figprops.sgroups, iGroup)
        xo.X.figprops.sgroups.(iGroup)(end + 1) = o;
    else
        xo.X.figprops.sgroups.(iGroup) = o;
    end
end

% a page named ?
iPg  = unique(iStr.Page);
myVG = deblank(iStr.VGroups);
if ~isempty(iPg)

    % just add the Pages to the VGroups field

    % positive page number -> groups: pageno and anypage
    if iPg(1) > 0
        myVG = [myVG sprintf(',UICPage%d', iPg(:)')];
        myVG = [myVG ',UICPage_any'];

    % number = -1 -> groups: allpages
    elseif iPg(1) == -1
        myVG = [myVG ',UICPage_all'];
        iStr.Page = [];

    % number < -1 or number == 0 , same as no pages -> no groups
    else
        iStr.Page = [];
    end

    % set back iStr.VGroups for correct handling
    if ~isempty(myVG) && myVG(1) == ','
        myVG(1) = [];
    end

    % then add to figure's array
    xo.X.figprops.pages = union(xo.X.figprops.pages, iStr.Page);
end

% any Visible-groups (or pages)
iGroups = ne_methods.splittocellc(myVG, ',; ', true, true);
for iGroupC = numel(iGroups):-1:1
    iGroup = deblank(iGroups{iGroupC});
    if isrealvarname(iGroup)
        if isfield(xo.X.figprops.vgroups, iGroup)
            xo.X.figprops.vgroups.(iGroup)(end+1) = o;
        else
            xo.X.figprops.vgroups.(iGroup) = o;
        end
    else
        iGroups(iGroupC) = [];
    end
end
iStr.VGroups = iGroups;

% any valid resize spec
if numel(iStr.ResizeSpec) ~= 2 || ~isrealvarname(iStr.ResizeSpec{1}) || ...
   ~isa(iStr.ResizeSpec{2}, 'double') || numel(iStr.ResizeSpec{2}) ~= 8
    iStr.ResizeSpec = {};
else
    iStr.ResizeSpec = {iStr.ResizeSpec{1}(:)', iStr.ResizeSpec{2}(:)'};
end

% complete object represetation
o.X.callbacks = iStr.Callbacks;
o.X.loadprops = iStr;
o.X.timeclick = {0, 0};
if isfield(oStr, 'progress')
    o.X.uicprops(1).nextupdate = now;
    o.X.uicprops(1).progress = oStr.progress;
    o.X.uicprops(1).updrate = xfigsngl.progupdrate;
end
o.X.uicprops(1).xchildren = xchildren;
xfigmlup(end + 1) = hOut;
xfigures(end + 1) = o;

% set info to parent object
xo.X.figprops.lilookup = -1;
if ~isempty(iStr.ResizeSpec) && ~isempty(iStr.Tag)
    xo.X.figprops.rszuics(end+1) = {o};
end

% if progress bar has initial value issue command
if iRSpec && strcmp(iCTyp, 'xprogress') && ~isempty(iStr.Value)
    progress(o, iStr.Value(1));

% if listbox and no line is to be selected
elseif strcmp(iCTyp, 'listbox') && ~isempty(iStr.Value) && iStr.Value(1) == 0
    set(hOut, 'Value', []);
end

% add to global tag lookup struct
xfigsngl.tags.(['UIC_' iStr.Tag]) = o;

% set handle visibility
set(hOut, 'HandleVisibility', xfigsngl.hvisible);

% get new OUID
dfcn = sprintf('delete(xfigure(hxdouble(''%s'')));', hxdouble(double(hOut)));
set(hOut, 'DeleteFcn', dfcn);
if ~isempty(iStr.ContextMenu)
    set(hOut, 'Callback', ['setcontent(xfigure(gcbo));' cbclick]);
end
