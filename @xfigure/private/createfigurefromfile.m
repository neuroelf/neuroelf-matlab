function o = createfigurefromfile(xo, tfgfile, varargin)
%XFIGURE::CREATEFIGUREFROMFILE  Create a figure based on a TFG file.
%   F = CREATEFIGUREFROMFILE(xfigure, FILENAME) creates the figure from
%   TFG file FILENAME.

% Version:  v1.1
% Build:    16042110
% Date:     Apr-21 2016, 10:15 AM EST
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

% TFG parser is a NeuroElf function
global ne_methods xfigures;

% only valid for the ROOT object
if numel(xo) ~= 1 || xo.T ~= 0
    error('neuroelf:xfigure:badObjectType', ...
        'CreateFigureFromFile can only be used for the ROOT object.');
end

% check filename
if nargin < 2 || ~ischar(tfgfile) || isempty(tfgfile) || exist(tfgfile(:)', 'file') ~= 2
    error('neuroelf:xfigure:fileNotFound', ...
        'The specified figure file was not found.');
end
tfgfile = tfgfile(:)';

% input options
if nargin < 3 || ~isstruct(varargin{1})
    figopts = struct;
    for vac = 3:nargin
        if ischar(varargin{vac-2}) && isrealvarname(varargin{vac-2})
            figopts.(varargin{vac-2}(:)') = true;
        end
    end
else
    figopts = varargin{1};
end

% mat/figstruct file?
if ~isfield(figopts, 'IsStruct') || isempty(figopts.IsStruct) || ...
   ~figopts.IsStruct(1)

    % read file and prepare options
    try
        fgstr = ne_methods.tfgparse(tfgfile);
    catch xfigerror
        error('neuroelf:xfigure:invalidTFGFile', ...
            'The given file (%s) is no valid TFG file (%s).', tfgfile, xfigerror.message);
    end
else
    try
        fgstr = load(iStr);
        if ~isfield(fgstr, 'TFG')
            error('INVALID_TFGMAT');
        end
        fgstr = fgstr.TFG;
        try
            if isfield(figopts, 'Evaluate') && numel(figopts.Evaluate) == 1 && ...
                islogical(figopts.Evaluate) && figopts.Evaluate
                fgstr = ne_methods.tfgparse(fgstr, struct('evaluate', true));
            end
        catch xfigerror
            rethrow(xfigerror);
        end
    catch xfigerror
        error('neuroelf:xfigure:invalidTFGMAT', ...
            'Error loading TFG struct file %s (%s).', tfgfile, xfigerror.message);
    end
end

% check loaded file
if ~isstruct(fgstr) || numel(fgstr) ~= 1 || ~isfield(fgstr, 'FIGURE') || ...
   ~isfield(fgstr, 'UICONTROLS') || ~isfield(fgstr, 'MENU') || ~isfield(fgstr, 'CONTEXTMENUS') || ...
   ~isstruct(fgstr.FIGURE) || ~isstruct(fgstr.UICONTROLS) || ~isstruct(fgstr.UIRESIZE) || ...
   ~isstruct(fgstr.MENU) || ~isstruct(fgstr.CONTEXTMENUS) || isempty(fgstr.FIGURE) || ...
   ~isfield(fgstr.FIGURE, 'Position') || ~isnumeric(fgstr.FIGURE.Position) || ...
    any(isnan(fgstr.FIGURE.Position(:)) | isinf(fgstr.FIGURE.Position(:))) || ...
    numel(fgstr.FIGURE.Position) ~= 4 || ~isfield(fgstr.UICONTROLS, 'Position') || ...
   ~isfield(fgstr.UICONTROLS, 'Type') || ~isfield(fgstr.UIRESIZE, 'Tag') || ...
   ~isfield(fgstr.UIRESIZE, 'Reference') || ~isfield(fgstr.UIRESIZE, 'RelPosition') || ...
   ~isfield(fgstr.MENU, 'Callback') || ~isfield(fgstr.MENU, 'Caption') || ...
   ~isfield(fgstr.MENU, 'Level') || ~isfield(fgstr.CONTEXTMENUS, 'IsCM') || ...
   ~isfield(fgstr.CONTEXTMENUS, 'Callback') || ~isfield(fgstr.CONTEXTMENUS, 'Caption') || ...
   ~isfield(fgstr.CONTEXTMENUS, 'Level') || ~isfield(fgstr.CONTEXTMENUS, 'Tag')
    error('neuroelf:xfigure:badTFGStruct', ...
        'Missing crucial information in TFG file (%s).', tfgfile);
end

% memorize visible flag, set to off during creation
if isfield(fgstr.FIGURE,'Visible')
    nVisible = fgstr.FIGURE.Visible;
else
    nVisible = 'on';
end
fgstr.FIGURE.Visible = 'off';

% tell CreateFigure(...) what filename we've got
fgstr.FIGURE.CFilename = tfgfile;

% override FieldLink / FieldLinkIni / OnLoad from command line
if isfield(figopts, 'FieldLinkCont')
    fgstr.FIGURE.FieldLinkCont = figopts.FieldLinkCont;
end
if isfield(figopts, 'FieldLinkSpec')
    fgstr.FIGURE.FieldLinkSpec = figopts.FieldLinkSpec;
end
if isfield(figopts, 'OnLoad')
    fgstr.FIGURE.OnLoad = figopts.OnLoad;
end

% create object and bark out on failure
try
    o = createfigure(xfigure, fgstr.FIGURE);
    figmpos = numel(xfigures);
    if ~isxfigure(o, 1)
        error('INVALID_FIGURE_OBJECT');
    end
catch xfigerror
    error('neuroelf:xfigure:createFigureFailed', 'Figure creation failed: %s.', xfigerror.message);
end

% add all tags of figure to a struct if no UserData in figure
if ~isfield(fgstr.FIGURE, 'UserData')
    UDTags = true;
    if isfield(fgstr.FIGURE, 'Tag') && isrealvarname(deblank(fgstr.FIGURE.Tag))
        UDTagStruct.(deblank(fgstr.FIGURE.Tag)) = o;
    else
        UDTagStruct = struct;
    end
else
    UDTags = false;
end

% process any uicontextmenus
if size(fgstr.CONTEXTMENUS, 1)
    numCMenus = numel(fgstr.CONTEXTMENUS);
    lastCMenu = xfigure;

    % iterate over context menus
    for ncount = 1:numCMenus
        thisMenu = fgstr.CONTEXTMENUS(ncount);

        % uicontextmenu object
        if ~isempty(thisMenu.IsCM) && thisMenu.IsCM(1)
            try
                lastCMenu(1:end) = [];
                hCMenu = adduicontextmenu(o, thisMenu);
                if ~isxfigure(hCMenu, 1)
                    error('INVALID_CONTEXTMENU');
                end
            catch xfigerror
                neuroelf_lasterr(xfigerror);
                warning('neuroelf:xfigure:badPropertyStruct', ...
                    'Couldn''t add context menu to figure.');
                continue;
            end
            lastCMenu(1) = hCMenu;

            % add tag to taglist (if OK)
            try
                if UDTags && isrealvarname(deblank(thisMenu.Tag))
                    UDTagStruct.(deblank(thisMenu.Tag)) = hCMenu;
                end
            catch xfigerror
                neuroelf_lasterr(xfigerror);
            end

        elseif ~isempty(lastCMenu)
            % must have a Caption and no Comment
            if (isfield(thisMenu,'Comment') && ~isempty(thisMenu.Comment)) || ...
                isempty(thisMenu.Caption)
                continue;
            end

            % handle UIMenu level
            hMCLevel = thisMenu.Level;
            if ischar(hMCLevel)
                try
                    hMCLevel = str2double(hMCLevel);
                catch xfigerror
                    neuroelf_lasterr(xfigerror);
                    continue;
                end
            elseif ~isnumeric(hMCLevel)
                continue;
            end
            if isempty(hMCLevel) || isnan(hMCLevel(1)) || isinf(hMCLevel(1)) || hMCLevel(1) < 1
                continue;
            end
            hMCLevel = fix(hMCLevel(1));
            if hMCLevel > numel(lastCMenu)
                continue;
            end

            % add UIMenu
            try
                hMenuItem = adduimenu(lastCMenu(hMCLevel), 'AddUIMenu', thisMenu);
                if ~isxfigure(hMenuItem, 1)
                    error('INVALID_MENUITEM');
                end
            catch xfigerror
                warning('neuroelf:xfigure:badPropertyStruct', ...
                    'Couldn''t add uimenu to context menu (%s).', xfigerror.message);
                continue;
            end

            % set new parent menu level object for lower levels
            lastCMenu(hMCLevel + 1) = hMenuItem;

            % add tag to taglist (if OK)
            try
                if UDTags && isrealvarname(deblank(thisMenu.Tag))
                    UDTagStruct.(deblank(thisMenu.Tag)) = hMenuItem;
                end
            catch xfigerror
                neuroelf_lasterr(xfigerror);
            end
        end
    end
end

% create resize info struct
rszinfo = struct;
rszused = false;
for nrsz = 1:numel(fgstr.UIRESIZE)
    thisResz = fgstr.UIRESIZE(nrsz);
    if isvarname(thisResz.Tag) && (~isfield(thisResz, 'Comment') || isempty(thisResz.Comment))
        rszinfo.(thisResz.Tag) = thisResz;
    end
end

% add UIControls if any
if size(fgstr.UICONTROLS, 1)
    numControls = numel(fgstr.UICONTROLS);

    % remember last button bar
    lastButtonBar = [];

    % set relative Position to 0,0 -> 0,0
    lastCtrlPos = [0, 0, 0, 0];
    for ncount = 1:numControls
        thisCtrl = fgstr.UICONTROLS(ncount);

        % must have a Type and no Comment
        if numel(thisCtrl.Type) > 0 && (~isfield(thisCtrl,'Comment') || ...
            isempty(thisCtrl.Comment))
            try
                % handle relative positioning
                if all(thisCtrl.Position(3:4) == 0)
                    thisCtrl.Position = thisCtrl.Position + lastCtrlPos;
                elseif any(thisCtrl.Position(3:4) < 0)
                    thisCtrl.Position(3:4) = abs(thisCtrl.Position(3:4));
                    if thisCtrl.Position(3) == 0
                        thisCtrl.Position(3) = lastCtrlPos(3);
                    elseif thisCtrl.Position(4) == 0
                        thisCtrl.Position(4) = lastCtrlPos(4);
                    end
                    thisCtrl.Position = [thisCtrl.Position(1:2) + lastCtrlPos(1:2), thisCtrl.Position(3:4)];
                end

                % add resize info if available
                try
                    if isfield(thisCtrl, 'Tag') && ~isempty(thisCtrl.Tag) && ...
                        isfield(rszinfo, thisCtrl.Tag)
                        thisCtrl.ResizeSpec = {rszinfo.(thisCtrl.Tag).Reference,  ...
                            rszinfo.(thisCtrl.Tag).RelPosition};
                        rszused = true;
                    end
                catch xfigerror
                    neuroelf_lasterr(xfigerror);
                end

                % add UIControl to figure (or button bar)
                if ~strcmpi(thisCtrl.Type, 'xbarbutton')
                    hControl = adduicontrol(o, thisCtrl);
                else
                    hControl = addbarbutton(lastButtonBar, thisCtrl);
                end
                if ~isxfigure(hControl, 1)
                    error('INVALID_UICONTROL');
                end

            % if couldn't add control say so
            catch xfigerror
                warning('neuroelf:xfigure:badPropertyStruct', ...
                    'Couldn''t add uicontrol to figure (%s).', xfigerror.message);
                continue;
            end

            % remember button bar
            if strcmpi(thisCtrl.Type, 'xbuttonbar')
                lastButtonBar = hControl;
            end

            % remember last position
            lastCtrlPos = thisCtrl.Position;

            % add tag to taglist (if OK)
            try
                if UDTags && isrealvarname(deblank(thisCtrl.Tag))
                    UDTagStruct.(deblank(thisCtrl.Tag)) = hControl;
                end
            catch xfigerror
                neuroelf_lasterr(xfigerror);
            end
        end
    end
end

% add UIMenus if any
if ~isempty(fgstr.MENU)
    numMenuLines = numel(fgstr.MENU);

    % setup parent objects array, parent for level 1 is figure !
    hMControl = o;
    for mcount = 1:numMenuLines
        thisMenu = fgstr.MENU(mcount);

        % must have a Caption and no Comment
        if (isfield(thisMenu,'Comment') && ~isempty(thisMenu.Comment)) || isempty(thisMenu.Caption)
            continue;
        end

        % handle UIMenu level
        hMCLevel = thisMenu.Level;
        if ischar(hMCLevel)
            try
                hMCLevel = str2double(hMCLevel);
            catch xfigerror
                neuroelf_lasterr(xfigerror);
                continue;
            end
        elseif ~isnumeric(hMCLevel)
            continue;
        end
        if isempty(hMCLevel) || isnan(hMCLevel(1)) || isinf(hMCLevel(1)) || hMCLevel(1) < 1
            continue;
        end
        hMCLevel = fix(hMCLevel(1));
        if hMCLevel > numel(hMControl)
            continue;
        end

        % add UIMenu
        try
            hMenuItem = adduimenu(hMControl(hMCLevel), thisMenu);
            if ~isxfigure(hMenuItem, 1)
                error('INVALID_MENUITEM');
            end
        catch xfigerror
            warning('neuroelf:xfigure:badPropertyStruct', ...
                'Couldn''t add uimenu to figure (%s).', xfigerror.message);
            continue;
        end

        % set new parent menu level object for lower levels
        hMControl(hMCLevel + 1) = hMenuItem;

        % add tag to taglist (if OK)
        try
            if UDTags && isrealvarname(deblank(thisMenu.Tag))
                UDTagStruct.(deblank(thisMenu.Tag)) = hMenuItem;
            end
        catch xfigerror
            neuroelf_lasterr(xfigerror);
        end
    end
end

% if no UserData given, fill with TagList
if UDTags
    set(o.H, 'UserData', struct('xfigure_UDTagStruct', UDTagStruct));
end

% context menu requested in figure
if isfield(fgstr.FIGURE, 'ContextMenu') && ~isempty(fgstr.FIGURE.ContextMenu) && ...
    ischar(fgstr.FIGURE.ContextMenu)
    if isfield(UDTagStruct, fgstr.FIGURE.ContextMenu(:)')
        try
            set(o.H, 'UIContextMenu', UDTagStruct.(fgstr.FIGURE.ContextMenu(:)').H);
        catch xfigerror
            neuroelf_lasterr(xfigerror);
        end
    end
end

% is FieldLink feature enabled ?
if isstruct(xfigures(figmpos).X.figprops.lgroups)

    % should we initially load any field groups named in file
    if isfield(fgstr.FIGURE, 'LoadFields') && ischar(fgstr.FIGURE.LoadFields)
        lfg = fgstr.FIGURE.LoadFields(:)';
    else
        lfg = '';
    end

    % override with selection made on command line
    if isfield(figopts, 'LoadFields') && ischar(figopts.LoadFields)
        lfg = figopts.LoadFields(:)';
    end

    % lfg = 'on' ?
    if strcmpi(lfg, 'on')
        lfg = 'all_groups';
    end

    % are we supposed to load any field groups ?
    if ~isempty(lfg)
        try
            loadfields(o, lfg, false);
        catch xfigerror
            neuroelf_lasterr(xfigerror);
        end
    end
end

% show any initial page from file ?
if isfield(fgstr.FIGURE, 'Page')
    pgnum = fgstr.FIGURE.Page;
else
    pgnum = 0;
end

% override with command line options
if isfield(figopts, 'Page')
    pgnum = figopts.Page;
end

% show a page?
if ~isempty(pgnum) && isnumeric(pgnum) && ~isnan(pgnum(1)) && ~isinf(pgnum(1))
    pgnum = fix(pgnum(1));
else
    pgnum = 0;
end

% if a page is requested then show it
if pgnum > 0
    try
        showpage(o, pgnum);
    catch xfigerror
        neuroelf_lasterr(xfigerror);
    end
end

% enable resize if requested
if rszused
    set(o.H, 'Resize', 'on', 'ResizeFcn', 'resize(xfigure(gcbf));');
end

% OnLoad method in either figopts or file
if isfield(fgstr.FIGURE,'OnLoad') && ischar(fgstr.FIGURE.OnLoad) && ...
   ~isempty(fgstr.FIGURE.OnLoad) && isempty(ne_methods.checksyntax(fgstr.FIGURE.OnLoad))
    try
        xfigurecallback(fgstr.FIGURE.OnLoad, o.H, o.H);
    catch xfigerror
        warning('neuroelf:xfigure:onLoadError', ...
            'Error executing OnLoad figure property: %s.', xfigerror.message);
    end
end

% set memorized Visible flag
if ~strcmpi(nVisible, 'off')
    set(o.H, 'Visible', 'on');
    figure(o.H);
    redrawfig(o.H);
end
