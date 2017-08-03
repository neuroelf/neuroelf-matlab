function i_loadfield(hxfigure, flinkini, flinktype, withcbs)

% global library
global ne_methods;

% get handle and callback
try
    tmphnm = hxfigure.H;
    tmpcbk = get(tmphnm, 'Callback');
    docbk  = false;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    return;
end

% and current value to set
try
    tmpval = [];
    eval(['tmpval=flinkini{1}.' flinkini{2} '.' flinkini{3} ';']);
catch ne_eo;
    disp(['xfigure::LoadField:*i_loadfield(...) -> couldn''t read ' ...
          flinkini{2} '.' flinkini{3} ' from ' ...
          strrep(Filename(flinkini{1}),'\','\\') ': ' ne_eo.message '.']);
    return;
end

% on/off switch
switch lower(flinktype{1}), case {'b', 'bool', 'oo', 'onoff'}

    try
        % ini setting is true
        if all(tmpval)
            set(tmphnm, 'Value', get(tmphnm, 'Max'));
            docbk = true;

        % not true
        else
            set(tmphnm, 'Value', get(tmphnm, 'Min'));
            if ~strcmpi(get(tmphnm, 'Style'), 'radiobutton')
                docbk = true;
            end
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

% single line strings
case {'c', 'char', 'chararray', 'string'}

    % reject bad input
    if ~ischar(tmpval) || ...
        numel(tmpval) ~= length(tmpval)
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Invalid datatype for <string> link.']);
        return;
    end

    % set value
    set(tmphnm, 'String', tmpval(:)');

% indices
case {'i', 'index'}

    % accept only numeric input
    if ~isnumeric(tmpval)
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Invalid datatype for <index> link.']);
        return;
    end

    % set value
    try
        set(tmphnm, 'Value', tmpval(:)');
        docbk = true;
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

% stringlist
case {'l', 'list', 'stringlist'}

    % make character lists
    if ischar(tmpval)
        if isempty(tmpval)
            tmpval = {''};
        elseif strcmpi(get(tmphnm, 'Style'), 'popupmenu') && ...
            length(tmpval) == numel(tmpval) && ...
            any(tmpval == '|')
            tmpval = ne_methods.splittocell(tmpval, '|');
        else
            tmpval = cellstr(tempval);
        end
    elseif iscell(tmpval) && ...
        isempty(tmpval) && ...
        strcmpi(get(tmphnm, 'Style'), 'popupmenu') && ...
        isempty(tmpval)
        tmpval = {''};
    end

    % reject non-cell fields
    if ~iscell(tmpval)
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Invalid datatype for <stringlist> link.']);
        return;
    end
    tmpval = tmpval(:);

    % default value ?
    if numel(flinktype) > 1 && ...
        isempty(tmpval) || ...
       (length(tmpval) == 1 && ...
        isempty(tmpval{1}))
        try
            tmpval = cellstr(flinktype{2});
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    end

    % set value
    try
        set(tmphnm, 'String', tmpval);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

% single number
case {'n', 'num', 'numstr', 'numeric'}

    % reject bad input
    if ~isnumeric(tmpval) || ...
        numel(tmpval) ~= 1
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Invalid datatype or content for <numeric> link.']);
        return;
    end

    % set value
    try
        set(tmphnm, 'String', num2str(tmpval));
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

% numeric array (MxN matrices only)
case {'na', 'numarray'}

    % reject bad input
    if ~isnumeric(tmpval) || ...
        ndims(tmpval) > 2 || ...
       (any(size(tmpval) == 0) && ...
        ~all(size(tmpval) == 0))
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Invalid datatype or content for <numarray> link.']);
        return;
    end

    % format value
    tempout = strrep(strrep(any2ascii(tmpval, 8), ',', ', '), ';', '; ');

    % set value
    try
        set(tmphnm, 'String', tempout(2:(end-1)));
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

% negated on/off
case {'negb', 'negbool', 'negoo', 'negonoff'}

    try
        % ini setting is false
        if ~any(tmpval)
            set(tmphnm, 'Value', get(tmphnm, 'Max'));
            docbk = true;

        % not true
        else
            set(tmphnm, 'Value', get(tmphnm, 'Min'));
            if ~strcmpi('radiobutton', get(tmphnm, 'Style'))
                docbk = true;
            end
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

% set (on specific value to true)
case {'s', 'set'}

    % reject bad link spec
    if length(flinktype) < 2
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Bad spec for <set> link.']);
        return;
    end

    try
        % ini setting matches
        if (numel(flinktype{2}) == numel(tmpval) || ...
            numel(flinktype{2}) == 1 || ...
            numel(tmpval) == 1) && ...
           all(flinktype{2} == tmpval(1))
            set(tmphnm, 'Value', get(tmphnm, 'Max'));
            docbk = true;

        % doesn't match
        else
            set(tmphnm, 'Value', get(tmphnm, 'Min'));
            if ~strcmpi(get(tmphnm, 'Style'), 'radiobutton')
                docbk = true;
            end
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

% sub-indexed on/off
case {'sb', 'subbool', 'soo', 'subonoff'}

    % reject bad link spec
    if length(flinktype) < 2 || ...
       ~isnumeric(flinktype{2}) || ...
        isempty(flinktype{2}) || ...
        isnan(flinktype{2}(1)) || ...
        isinf(flinktype{2}(1)) || ...
        flinktype{2}(1) < 1 || ...
        flinktype{2}(1) > numel(tmpval)
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Bad spec for <subonoff> link.']);
        return;
    end

    try
        % ini setting true
        if tmpval(fix(flinktype{2}(1)))
            set(tmphnm, 'Value', get(tmphnm, 'Max'));
            docbk = true;

        % false
        else
            set(tmphnm, 'Value', get(tmphnm, 'Min'));
            if ~strcmpi(get(tmphnm, 'Style'), 'radiobutton')
                docbk = true;
            end
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

% sub-indexed indices
case {'si', 'subindex'}

    % reject bad link spec
    if length(flinktype) < 2 || ...
       ~isnumeric(flinktype{2}) || ...
        isempty(flinktype{2}) || ...
        isnan(flinktype{2}(1)) || ...
        isinf(flinktype{2}(1)) || ...
        flinktype{2}(1) < 1 || ...
        flinktype{2}(1) > numel(tmpval)
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Bad spec for <subindex> link.']);
        return;
    end

    % reject bad input value
    if ~isnumeric(tmpval)
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Invalid datatype for <subindex> link.']);
        return;
    end

    % set value
    try
        set(tmphnm, 'Value', tmpval(fix(flinktype{2}(1))));
        docbk = true;
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

% sub-indexed single number
case {'sn', 'subnum', 'subnumstr', 'subnumeric'}

    % reject bad link spec
    if length(flinktype) < 2 || ...
       ~isnumeric(flinktype{2}) || ...
        isempty(flinktype{2}) || ...
        isnan(flinktype{2}(1)) || ...
        isinf(flinktype{2}(1)) || ...
        flinktype{2}(1) < 1 || ...
        flinktype{2}(1) > numel(tmpval)
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Bad spec for <subnumeric> link.']);
        return;
    end

    % reject bad input value
    if ~isnumeric(tmpval)
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Invalid datatype for <subnumeric> link.']);
        return;
    end

    % set value
    try
        set(tmphnm, 'String', num2str(tmpval(fix(flinktype{2}(1)))));
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

% sub-indexed set
case {'ss', 'subset'}

    % reject bad link spec
    if length(flinktype) < 3 || ...
       ~isnumeric(flinktype{2}) || ...
        isempty(flinktype{2}) || ...
        isnan(flinktype{2}(1)) || ...
        isinf(flinktype{2}(1)) || ...
        flinktype{2}(1) < 1 || ...
        flinktype{2}(1) > numel(tmpval)
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Bad spec for <subnumeric> link.']);
        return;
    end

    % get correct value
    tmpval = tmpval(fix(flinktype{2}));
    if iscell(tmpval)
        try
            tmpval = [tmpval{:}];
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            return;
        end
    end

    try
        if flinktype{3} == tmpval(flinktype{2})
            set(tmphnm, 'Value', get(tmphnm, 'Max'));
            docbk = true;
        else
            set(tmphnm, 'Value', get(tmphnm, 'Min'));
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        try
            set(tmphnm, 'Value', get(tmphnm, 'Min'));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    end

% unknown link type -> display warning
otherwise
    disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
          'Unknown link type <' flinktype{1}(:)' '>.']);
end

% do callback ?
if withcbs && ...
    docbk && ...
   ~isempty(tmpcbk)
    try
        xfigurecallback(tmpcbk, ...
            findmlparent(tmphnm, 'figure', struct( ...
            'r', 0,           'root',            0, ...
            'f', 1, 'fig', 1, 'figure',          1, ...
            'c', 2, 'uic', 2, 'uicontrol',       2, ...
            'm', 3, 'uim', 3, 'uimenu',          3, ...
            'x', 4, 'uix', 4, 'uicontextmenu',   4  ...
            )), tmphnm, true);
    catch ne_eo;
        disp(['xfigure::LoadField:*i_loadfield(...) -> ' ...
              'Error executing Callback:' char(10) ...
              'Callback: ' tmpcbk char(10) ...
              'Lasterr:  ' ne_eo.message]);
    end
end
