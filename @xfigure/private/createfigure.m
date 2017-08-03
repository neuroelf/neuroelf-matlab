function o = createfigure(xo, iStr)
%XFIGURE::CREATEFIGURE  Create figure from struct-based configuration.
%   FIG = CREATEFIGURE(xfigure, FIGSTRUCT) uses the information in
%   FIGSTRUCT to create the figure with handle FIG.

% Version:  v1.1
% Build:    16042110
% Date:     Apr-21 2016, 10:16 AM EST
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

% NeuroElf library global storage and references
global ne_methods xfigmlup xfigsngl xfigures;

% only valid if parent is ROOT object
if numel(xo) ~= 1 || xo.T ~= 0
    error('neuroelf:xfigure:invalidObjectType', ...
        'CreateFigure can only be used for the ROOT object.');
end

% test iStr
if ~isstruct(iStr) || numel(iStr) ~= 1
    error('neuroelf:xfigure:badPropertyStruct', ...
        'Figure properties must be of type struct.');
end

% perform some checks on struct but allow function handles
icallbacks = repmat({''}, 1, 9);
if isfield(iStr, 'CallbackKey') && isa(iStr.CallbackKey, 'function_handle')
    icallbacks{5} = iStr.CallbackKey;
end
if isfield(iStr, 'CallbackMDown') && isa(iStr.CallbackMDown, 'function_handle')
    icallbacks{6} = iStr.CallbackMDown;
end
if isfield(iStr, 'CallbackMMove') && isa(iStr.CallbackMMove, 'function_handle')
    icallbacks{7} = iStr.CallbackMMove;
end
if isfield(iStr, 'CallbackMUp') && isa(iStr.CallbackMUp, 'function_handle')
    icallbacks{8} = iStr.CallbackMUp;
end
if isfield(iStr, 'CallbackResize') && isa(iStr.CallbackResize, 'function_handle')
    icallbacks{9} = iStr.CallbackResize;
end
iStr = ne_methods.checkstruct(iStr, xfigsngl.optfig);
if ~isempty(iStr.CallbackKey)
    icallbacks{5} = iStr.CallbackKey;
end
if ~isempty(iStr.CallbackMDown)
    icallbacks{6} = iStr.CallbackMDown;
end
if ~isempty(iStr.CallbackMMove)
    icallbacks{7} = iStr.CallbackMMove;
end
if ~isempty(iStr.CallbackMUp)
    icallbacks{8} = iStr.CallbackMUp;
end
if ~isempty(iStr.CallbackResize)
    icallbacks{9} = iStr.CallbackResize;
end
if numel(iStr.Position) ~= 4 || any(isnan(iStr.Position(:)) | isinf(iStr.Position(:)))
    error('neuroelf:xfigure:badPropertyStruct', 'Bad figure property struct supplied.');
end

% use tag from iStr?
if ~isempty(iStr.Tag) && numel(iStr.Tag) < 28 && isrealvarname(iStr.Tag(:)')
    utag = iStr.Tag(:)';
else
    iStr.Tag = sprintf('FIG_%010.0f', floor(1e10 * rand(1, 1)));
    utag = ['xfigure_' iStr.Tag];
end

% background color
if numel(iStr.Color) == 1
    iCBG(1:3) = max(0, min(1, iStr.Color));
elseif numel(iStr.Color) == 3
    iCBG = max(0, min(1, iStr.Color));
else
    iCBG = [];
end

% close requestor
if isempty(iStr.CallbackClReq)
    iStr.CallbackClReq = 'closereq;';
end

% context menu
if ~isempty(iStr.ContextMenu)
    try
        uicm = findobj('Type', 'uicontextmenu', 'Tag', iStr.ContextMenu);
        if numel(uicm) ~= 1
            uicm = [];
        end
    catch xfigerror
        neuroelf_lasterr(xfigerror);
        uicm = [];
    end
else
    uicm = [];
end

% set required options
iPar = 0;
iPos = iStr.Position;
if all(iPos(1:2) == -1) && strcmpi(iStr.Units, 'pixels')
    iPos(1) = fix((xfigsngl.units.pixels(3) - iPos(3)) / 2 + 0.5001);
    iPos(2) = fix((xfigsngl.units.pixels(4) - iPos(4)) / 2 + 0.5001);
end
iPos(3:4) = max(0, iPos(3:4));
if numel(iStr.MinSize) == 2 && all(iStr.MinSize > 0) && ...
   ~any(isnan(iStr.MinSize) | isinf(iStr.MinSize))
    iStr.MinSize = iStr.MinSize(:)';
else
    iStr.MinSize = [];
end

% fill with those and other, fixed options
oStr = struct( ...
    'Parent',                iPar, ...
    'Units',                 iStr.Units, ...
    'Position',              iPos, ...
    'BackingStore',          iStr.BackingStore, ...
    'ButtonDownFcn',         '', ...
    'CloseRequestFcn',       iStr.CallbackClReq, ...
    'DeleteFcn',             '', ...
    'DoubleBuffer',          iStr.DoubleBuffer, ...
    'IntegerHandle',         iStr.IntegerHandle, ...
    'Interruptible',         iStr.Interrupts, ...
    'InvertHardCopy',        iStr.PrintBW, ...
    'KeyPressFcn',           icallbacks{5}, ...
    'Menu',                  'none', ...
    'MenuBar',               iStr.MenuBar, ...
    'Name',                  iStr.Title, ...
    'NumberTitle',           'off', ...
    'PaperUnits',            iStr.PaperUnits, ...
    'PaperOrientation',      iStr.PaperOrientation, ...
    'PaperPosition',         iStr.PaperPosition, ...
    'PaperSize',             iStr.PaperSize, ...
    'PaperType',             iStr.PaperType, ...
    'Resize',                iStr.Resizeable, ...
    'ResizeFcn',             icallbacks{9}, ...
    'Tag',                   utag, ...
    'Toolbar',               'none', ...
    'UIContextMenu',         uicm, ...
    'UserData',              iStr.UserData, ...
    'Visible',               'off', ...
    'WindowButtonDownFcn',   icallbacks{6}, ...
    'WindowButtonMotionFcn', icallbacks{7}, ...
    'WindowButtonUpFcn',     icallbacks{8});
if ~isempty(iCBG)
    oStr.Color = iCBG;
end
if ~strcmpi(iStr.Modal, 'on')
    oStr.WindowStyle = 'normal';
else
    oStr.WindowStyle = 'modal';
end

% create figure and fill fields
hOut = figure(oStr);
set(hOut, 'Position', iPos);
o = xfigure('new');
o.H = hOut;
o.T = 1;

% FieldLink requested
if ~isempty(iStr.FieldLinkCont) && ~isempty(iStr.FieldLinkSpec)

    flcont = iStr.FieldLinkCont(:)';
    flspec = iStr.FieldLinkSpec(:)';

    % support for single content
    if ischar(flcont) && (exist(flcont, 'file') == 2 || strcmpi(flcont, '_auto'))
        if strcmpi(flcont, '_auto') && isfield(iStr, 'CFilename')
            [tgfpath, tgfname] = fileparts(iStr.CFilename);
            flcont = [tgfpath filesep tgfname '.ini'];
            if exist(flcont, 'file') ~= 2
                warning('neuroelf:xfigure:fieldLinkFailed', ...
                    'Auto modus for FieldLinkCont failed.');
                iStr.FieldLinkCont = [];
                flcont = [];
            end
        end

        if ~isempty(flcont)
            flcont = xini(flcont, 'convert');
            flname = Filename(flcont);
            [trfpath{1:2}] = fileparts(flname);
            iStr.FieldLinkCont = struct(ne_methods.makelabel(trfpath{2}), flcont);
        end

    % input is already in struct format
    elseif isstruct(flcont)

        % make sure content is either a valid file or an xini object
        inis = fieldnames(flcont);
        for cc = numel(inis):-1:1
            if ischar(flcont.(inis{cc})) && exist(flcont.(inis{cc}), 'file') == 2
                iStr.FieldLinkCont.(inis{cc}) = xini(flcont.(inis{cc}), 'convert');
            end
            if ~isxini(flcont.(inis{cc}), 1)
                iStr.FieldLinkCont = rmfield(iStr.FieldLinkCont, inis{cc});
            end
        end

    % is it an xini handle ?
    elseif isxini(flcont, 1)
        [trfpath{1:2}] = fileparts(Filename(iStr.FieldLinkCont));
        iStr.FieldLinkCont = struct(ne_methods.makelabel(trfpath{2}), iStr.FieldLinkCont);

    % else set to empty array
    else
        iStr.FieldLinkCont = [];
    end

    % support for single specification file
    if ischar(flspec) && (exist(flspec, 'file') == 2 || strcmpi(flspec, '_auto'))
        if strcmpi(flspec, '_auto') && isfield(iStr, 'CFilename')
            [tgfpath, tgfname] = fileparts(iStr.CFilename);
            iStr.FieldLinkSpec = [tgfpath filesep tgfname '.fln'];
            if exist(iStr.FieldLinkSpec, 'file') ~= 2
                warning('neuroelf:xfigure:fieldLinkFailed', ...
                    'Auto modus for FieldLinkSpec failed.');
                iStr.FieldLinkSpec = {};
            end
        else
            iStr.FieldLinkSpec = {};
        end
        if ~isempty(iStr.FieldLinkSpec)
            iStr.FieldLinkSpec = {xini(iStr.FieldLinkSpec, 'convert')};
        end

    % input is already in cell format
    elseif iscell(flspec)

        % make sure cell content is either a valid file or an xini object
        for cc = numel(flspec):-1:1
            if ischar(flspec{cc}) && exist(flspec{cc},'file') == 2
                iStr.FieldLinkSpec{cc} = xini(iStr.FieldLinkSpec{cc}, 'convert');
            end
            if ~isxini(iStr.FieldLinkSpec{cc}, 1)
                iStr.FieldLinkSpec(cc) = [];
            end
        end

    % even not ? then make an empty cell array
    else
        iStr.FieldLinkCont = [];
    end

    % check whether we still have content ?
    if isempty(iStr.FieldLinkCont) || isempty(iStr.FieldLinkSpec)
        iStr.FieldLinkCont = struct;
        iStr.FieldLinkSpec = {};
    end

    % find groups
    myLGroups = [];
    for cc = 1:numel(iStr.FieldLinkSpec)
        gnames = IniSections(iStr.FieldLinkSpec{cc});
        for gc=1:numel(gnames)
            fnames = IniSectionSettings(iStr.FieldLinkSpec{cc}, gnames{gc});
            myLGroups.(gnames{gc}) = {cc, fnames};
        end
    end
    iStr.FieldLinkGroups = myLGroups;

% make sure fields are used well
else
    iStr.FieldLinkCont   = struct;
    iStr.FieldLinkSpec   = {};
    iStr.FieldLinkGroups = [];
end

% fill additional internal object representation
o.X = makeostruct(1);
o.X.callbacks = icallbacks;
o.X.figprops(1).cpage = -2;
o.X.figprops.egroups  = struct;
o.X.figprops.lgroups  = iStr.FieldLinkGroups;
o.X.figprops.lilookup = -1;
o.X.figprops.linkcont = iStr.FieldLinkCont;
o.X.figprops.linkspec = iStr.FieldLinkSpec;
o.X.figprops.llookup  = 0;
o.X.figprops.pages    = [];
o.X.figprops.rgroups  = struct;
o.X.figprops.rszuics  = {};
o.X.figprops.sgroups  = struct;
o.X.figprops.vgroups  = struct;
o.X.loadprops = iStr;
xfigmlup(end + 1) = hOut;
xfigures(end + 1) = o;

% add to global tag lookup struct
xfigsngl.tags.(['FIG_' iStr.Tag]) = o;

% set correct visible state
set(hOut, 'Visible', iStr.Visible);

% set some Fcn settings
if ~isempty(uicm)
    btdfcn = sprintf('setcontent(xfigure(hxdouble(''%s'')));', hxdouble(double(hOut)));
else
    btdfcn = '';
end
set(hOut, ...
    'ButtonDownFcn', [btdfcn, iStr.CallbackClick(:)'], ...
    'DeleteFcn', sprintf('delete(xfigure(hxdouble(''%s'')));', hxdouble(double(hOut))));
