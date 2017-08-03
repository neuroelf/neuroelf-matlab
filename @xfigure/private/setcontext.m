% handling of context menu events
case {'setcontext'}

if nargin < 3 || ...
    isempty(varargin{3}) || ...
   ~any(strcmpi(class(varargin{3}), {'xfigure', 'double'}))
    error( ...
        'xfigure:BadArgument', ...
        'Setting the context menu object requires more input.' ...
    );
end
hCMObj = varargin{3}(1);

% possible GUI handle
if isa(hCMObj, 'double')

    % is it a handle at all
    if ~ishandle(hCMObj)
        error( ...
            'xfigure:BadArgument', ...
            'Illegal GUI handle passed as reference.' ...
        );
    end

    % context menu property OK?
    try
        ttype = lower(get(hCMObj, 'Type'));
        if ~isfield(xfigsngl.uixtypes, ttype)
            error('BAD_OBJECT_TYPE');
        end
        testcm = get(hCMObj, 'UIContextMenu');

        % only set if there is a context menu!
        if isempty(testcm)
            return;
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error( ...
            'xfigure:BadArgument', ...
            'Illegal GUI handle passed as reference.' ...
        );
    end

    % try looking up object
    hlup = find(xfigmlup == double(hCMObj));
    if isempty(hlup)
        error( ...
            'xfigure:LookupFailed', ...
            'Requested object is not under xfigure control.' ...
        );
    end

    % set object
    xfigsngl.contextobject = xfigure(0, 'makeobj', ...
        xfig_ilup(hlup(1)), ...
        xfigmlup(hlup(1)), ...
        xfig_type(hlup(1)));

% xfigure object
else

    % is this a valid handle
    if ~isxfigure(hCMObj, 1) || ...
       ~any([xfigsngl.objtypes.figure, ...
             xfigsngl.objtypes.uicontrol] == hCMObj.type)
        error( ...
            'xfigure:BadArgument', ...
            'Illegal xfigure handle passed as reference.' ...
        );
    end

    % return on empty context menu
    try
        testcm = get(hCMObj.mhnd, 'UIContextMenu');
        if isempty(testcm)
            return;
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error( ...
            'xfigure:BadArgument', ...
            'Illegal xfigure handle passed as reference.' ...
        );
    end

    % set object
    xfigsngl.contextobject = hCMObj;
end
