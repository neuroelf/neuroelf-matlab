classdef xfigure < handle
%XFIGURE  eXtended Figure (UI handle) class with additional properties.
%
%   XFIG = XFIGURE creates a factory object, which can be used to access
%   general methods and settings.
%
%   Other objects returned by XFIG are of valid graphics handle subtypes.
%
%   -> ROOT            (the factory itself, corresponds to MATLAB handle 0)
%   -> Figure          (MATLAB figure objects)
%   -> UIControl       (MATLAB uicontrols objects)
%   -> UIMenu          (MATLAB uimenu objects)
%   -> UIContextMenu   (MATLAB uicontextmenu objects)
%
% The following constructors are available:
%
% xfigure
% xfigure(matlabhandle [, 'delete'])
% xfigure(filename [, options])
% xfigure(objecttag)
%
% Using: checkstruct, checksyntax, gluetostring, mainver, makelabel,
%        splittocellc, subget, tfgparse.

% xfigures          object array with .X properties
%       .callbacks      = cell array for objects callbacks
%          {1}            Callback / CallbackClick
%          {2}            CallbackClReq
%          {3}            CallbackDblClick
%          {4}            CallbackDelete
%          {5}            CallbackKey
%          {6}            CallbackMDown
%          {7}            CallbackMMove
%          {8}            CallbackMUp
%          {9}            CallbackResize
%       .deletefcn      = original CallbackDelete value (replaced!)
%       .figprops       = struct with sub fields
%          .cpage       = currently shown page for figures [set to -2]
%          .egroups     = SetGroupEnabled groups (with lookup pointers
%                         UIControls: positive/UIMenus: negative doubles)
%          .lgroups     = groups of linked fields for load and update
%          .lilookup    = struct for lookups (max 31 chars!)
%          .linkcont    = FieldLinkCont xini object handles struct
%          .linkspec    = FieldLinkSpec xini object handles array
%          .llookup     = 1x1 struct array for lookups (max 31 chars!)
%          .pages       = list of available pages
%          .rgroups     = RadioGroupSetOne groups (with lookup pointers)
%          .rszuics     = cell array with affected UICs for CallbackResize
%          .sgroups     = SlideGroupXY groups (with lookup pointers)
%          .vgroups     = SetGroupVisible groups (with lookup pointers
%                         UIControls: positive/UIMenus: negative doubles)
%       .loadprops      = struct with input fields at creation time
%       .timeclick      = time index of last OnClick event (for DblClick)
%       .uicprops       = struct with extended properties

% Version:  v1.1
% Build:    16060812
% Date:     Jun-08 2016, 12:57 PM EST
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

% hidden properties
properties (Access = 'private', Hidden = true)
    H = [];     % MATLAB internal handle
    T = -1;     % numeric type (for internal branching)
    X = struct; % extended properties (type dependent)
end


% methods (no indentation for readability)
methods (Access = 'public', Hidden = true)


% constructor
function xo = xfigure(varargin)

    % requires neuroelf and global factory
    global xfigsngl xfigures xfigmlup xfiguremeth;

    % check for class initialization -> do if necessary
    if isempty(xfigsngl) || ~isstruct(xfigsngl) || ~xfigsngl.is_init

        % initialize
        xfiginit('v1.0;16040912;20160409123808');

        % create root object
        xo.H = xfigmlup(1);
        xo.T = 0;
        xo.X = struct('loadprops', xfigsngl.rp_init);
        xfigures = xo;

        % also call methods
        methods(xo);
    end

    % no or empty first argument
    if nargin < 1 || isempty(varargin{1})
        xo = xfigures(1);
        return;

    % new and empty (unregistered) object
    elseif ischar(varargin{1}) && strcmpi(varargin{1}(:)', 'new')
        xo.H = -1;
        xo.T = -1;
        xo.X = struct;
        return;
    end

    % for objects
    if isa(varargin{1}, 'xfigure') && numel(varargin{1}) == 1
        xo = varargin{1};

    % for numeric input
    elseif isa(varargin{1}, 'double') && numel(varargin{1}) == 1 && ...
        ~isinf(varargin{1}) && ~isnan(varargin{1})

        % internal handle lookup first
        myhnd = varargin{1};
        xo = [];
        ilup = find(double(xfigmlup) == myhnd);
        if ~isempty(ilup)

            % handle no longer exists
            if ~ishandle(xfigmlup(ilup(1)))

                % remove from lookups and internal field
                xfigmlup(ilup) = [];
                xfigures(ilup)  = [];
                xo(:) = [];
                return;
            end

            % otherwise use found object
            xo = xfigures(ilup(1));
        end

        % special case -> MATLAB handle, child of xfigure object
        if isempty(ilup) && ishandle(myhnd)

            % check type first
            myhtype = get(myhnd, 'Type');
            if isfield(xfigsngl.xobjtypes, myhtype)

                % get UserData
                pv = get(myhnd, 'UserData');

                % if that is a xfigure
                if numel(pv) == 1 && isxfigure(pv, 1)
                    xo = pv;

                % otherwise
                else
                    xo(:) = [];
                    return;
                end
            end
        end

        % error if lookup failed
        if isempty(xo)
            if nargin < 2 || ~ischar(varargin{2}) || ...
               ~any(strcmpi(varargin{2}(:)', {'delete', 'isvalid'}))
                error('neuroelf:xfigure:lookupFailed', 'Lookup of type double failed.');
            end
            xo(:) = [];
            return;
        end

        % second argument type char and == 'delete'
        if nargin > 1 && ischar(varargin{2}) && strcmpi(varargin{2}(:)', 'delete')

            % delete all objects found and return
            try
                delete(xo);
            catch xfigerror
                neuroelf_lasterr(xfigerror);
            end
            return;
        end

    % new MATLAB R2014b classes
    elseif numel(varargin{1}) == 1 && ...
       (isa(varargin{1}, 'matlab.ui.control.UIControl') || ...
        isa(varargin{1}, 'matlab.ui.container.Menu') || ...
        isa(varargin{1}, 'matlab.graphics.axis.Axes') || ...
        isa(varargin{1}, 'matlab.ui.Figure'))

        % look-up
        myhnd = varargin{1};
        ilup = find(xfigmlup == myhnd);
        if ~isempty(ilup)

            % handle no longer exists
            if ~ishandle(myhnd)

                % remove from lookups and internal field
                xfiggc();
                xo(:) = [];
                return;
            end

            % otherwise make object from it
            xo = xfigures(ilup(1));

        % error out
        else
            error('neuroelf:xfigure:lookupFailed', 'Error looking up MATLAB UI handle.');
        end

    % char constructor
    elseif ischar(varargin{1}) && ~isempty(varargin{1})

        % is it a file
        tchar = varargin{1}(:)';
        if exist(tchar, 'file') == 2

            % create figure from file then
            try
                xo = createfigurefromfile(xfigures(1), tchar, varargin{2:nargin});
                if ~isxfigure(xo, 1) || numel(xo) ~= 1
                    error('FIGURE_NOT_CREATED');
                end
                return;
            catch xfigerror
                error('neuroelf:xfigure:figureCreationFailed', ...
                    'Couldn''t create figure from TFG file (%s): %s.', tchar, ...
                    xfigerror.message);
            end
        end

        % not a valid TAG
        if any(tchar == filesep | tchar == '.') && ~isrealvarname(tchar)
            error('neuroelf:xfigure:invalidTag', 'Bad TAG name for lookup (%s).', tchar);
        end

        % try object lookups in Tag table, iterate over lookup types
        xo = [];
        for olc = 1:numel(xfigsngl.objtlup)
            tTag = [xfigsngl.objtlup{olc} tchar];

            % is field in tag list
            if isfield(xfigsngl.tags, tTag)

                % get reference
                xo = xfigsngl.tags.(tTag);
                break;
            end
        end

        % error if lookup failed
        if numel(xo) ~= 1 || ~isxfigure(xo, 1)
            error('neuroelf:xfigure:lookupFailed', 'Error looking up xfigure object (%s).', tchar)
        end

    % other input type
    else
        error('neuroelf:xfigure:badArgument', 'Invalid input argument class: %s.', class(varargin{1}));
    end

    % method call
    if nargin > 1 && numel(xo) == 1 && ischar(varargin{2}) && ...
       ~isempty(varargin{2}) && isfield(xfiguremeth.m, lower(varargin{2}(:)'))

        % try to pass on
        try
            xmeth = xfiguremeth.m.(lower(varargin{2}(:)'));
            feval(xmeth{1}, xo, varargin{3:end});
        catch xfigerror
            rethrow(xfigerror);
        end
        return;
    end

    % if nothing else to do return
    if nargin > 1
        error('neuroelf:xfigure:invalidConstructorSyntax', 'Invalid constructor syntax.');
    end
end

% single-object delete
function delete(xo)
    if numel(xo.H) == 1 && ishandle(xo.H) && (isa(xo.H, 'double') || isvalid(xo.H)) && xo.T ~= 0
        if xo.T > 0 && numel(xo.X.callbacks) > 3 && ~isempty(xo.X.callbacks{4})
            set(xo.H, 'DeleteFcn', '');
            try
                evalin('base', xo.X.callbacks{4});
            catch xfigerror
                warning(xfigerror.message);
            end
            xo.X.callbacks{4} = '';
        end
        if ~strcmpi(get(xo.H, 'BeingDeleted'), 'on')
            delete(xo.H);
        end
    end
end

% delete
function deletefigure(xo)

    % global reference storage
    global xfigmlup xfigures;

    % remove invalid handles
    xo(~isvalid(xo)) = [];
    if isempty(xo)
        return;
    end

    % get handles
    h = [xo.H];

    % iterate over objects
    for oc = 1:numel(xo)

        % don't do anything if no longer valid
        if ~ishandle(h(oc))
            continue;
        elseif ~isnumeric(h) && ~isvalid(h(oc))
            continue;
        end

        % anything but the root object
        if xo(oc).T > 0

            % call delete
            delete(xo(oc));
        end
    end

    % remove
    [h, ih] = intersect(xfigmlup(:), h(:));
    xfigmlup(ih) = [];
    xfigures(ih) = [];
end

% properties (for tab-completion)
function p = properties(xo)
    global xfiguremeth;
    m = xfiguremeth.m;
    if numel(xo) == 1
        try
            p = fieldnames(get(xo.H));
            pm = methods(xo);
            for pmc = 1:numel(pm)
                pm(pmc) = m.(pm{pmc})(3);
            end
            p = [p(:); pm(:)];
        catch xfigerror
            neuroelf_lasterr(xfigerror);
            p = {};
        end
    elseif isempty(xo)
        p = {};
    else
        p = properties(xo(1));
        for c = 1:numel(xo)
            p = intersect(p, properties(xo(c)));
            if isempty(p)
                return;
            end
        end
    end
end


% end of methods
end

% end of classdef
end
