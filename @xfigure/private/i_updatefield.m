function usuccess = i_updatefield(hxfigure, fieldlink, flinktype)

% neuroelf library
global ne_methods;

% suppose we didn't succeed
usuccess = false;

% get handles
try
    fieldhnm = hxfigure.mhnd;
    fieldini = fieldlink{1};
    fieldfld = [fieldlink{2} '.' fieldlink{3}];
    fieldevl = ['fieldini.' fieldfld];
catch ne_eo;
    disp(['xfigure::UpdateField:*i_updatefield(...) -> ' ...
          'Bad arguments, calling convention error (' ne_eo.message ').']);
    return;
end

% get current value (assume empty matrix for non-existant)
try
    eval(['fieldval=' fieldevl ';']);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    fieldval = [];
    try
        eval([fieldevl '=[];']);
    catch ne_eo;
        disp(['xfigure::UpdateField:*i_updatefield(...) -> ' ...
              'Illegal fieldname: ' fieldevl ' (' ne_eo.message ').']);
        return;
    end
end

% on/off switch
switch lower(flinktype{1}), case {'b', 'bool', 'oo', 'onoff'}
    try
        bval = get(fieldhnm, 'Value');
        if all(bval == get(fieldhnm, 'Max'))
            bval = 1;
        else
            bval = 0;
        end
        eval([fieldevl '=bval;']);
        usuccess = true;
    catch ne_eo;
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Couldn''t set <onoff> field ' fieldfld ' (' ne_eo.message ').']);
    end

% single line strings
case {'c', 'char', 'chararray', 'string'}
    try
        eval([fieldevl '=get(fieldhnm,''String'');']);
        usuccess = true;
    catch ne_eo;
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Couldn''t set <string> field ' fieldfld ' (' ne_eo.message ').']);
    end

% indices
case {'i', 'index'}
    try
        eval([fieldevl '=get(fieldhnm,''Value'');']);
        usuccess = 1;
    catch ne_eo;
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Couldn''t set <index> field ' fieldfld ' (' ne_eo.message ').']);
    end

% string list
case {'l', 'list', 'stringlist'}
    try
        cstr = get(fieldhnm, 'String');
        if ischar(cstr)
            if isempty(cstr)
                cstr = {};
            elseif strcmpi(get(fieldhnm, 'Style'), 'popupmenu') && ...
                length(cstr) == numel(cstr) && ...
                any(cstr == '|')
                cstr = ne_methods.splittocell(cstr, '|');
            else
                cstr = cellstr(cstr);
            end
        end
        eval([fieldevl '=cstr(:);']);
        usuccess = 1;
    catch ne_eo;
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Couldn''t set <stringlist> field ' fieldfld ' (' ne_eo.message ').']);
    end

% single number (default: 0)
case {'n', 'num', 'numstr', 'numeric'}
    try
        cstr = ne_methods.splittocellc(get(fieldhnm, 'String'), ' ,;', true, true);
        if ~isempty(cstr) && ...
            isempty(cstr{1})
            cstr(1) = [];
        end
        if isempty(cstr)
            tval = 0;
        else
            try
                tval = str2double(cstr{1});
                tval = tval(1);
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                tval = 0;
            end
        end
        eval([fieldevl '=tval;']);
        usuccess = 1;
    catch ne_eo;
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Couldn''t set <numeric> field ' fieldfld ' (' ne_eo.message ').']);
    end

% numeric array (MxN matrices only)
case {'na', 'numarray'}
    try
        cstr = get(fieldhnm, 'String');
        if ~ischar(cstr) || ...
            length(cstr) ~= numel(cstr) || ...
            any(cstr == '|')
            error('INVALID_NUMARRAY_STRING');
        end
        eval([fieldevl '=[' cstr '];']);
        usuccess = 1;
    catch ne_eo;
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Couldn''t set <numarray> field ' fieldfld ' (' ne_eo.message ').']);
    end

% negated on/off switch
case {'negb', 'negbool', 'negoo', 'negonoff'}
    try
        bval = get(fieldhnm, 'Value');
        if ~all(bval == get(fieldhnm, 'Max'))
            bval = 0;
        else
            bval = 1;
        end
        eval([fieldevl '=bval;']);
        usuccess = 1;
    catch ne_eo;
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Couldn''t set <negonoff> field ' fieldfld ' (' ne_eo.message ').']);
    end

% set
case {'s', 'set'}
    if length(flinktype) < 2
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Bad spec for <set> link.']);
        return;
    end
    try
        cval = get(fieldhnm, 'Value');
        if all(cval == get(fieldhnm, 'Max'))
            eval([fieldevl '=flinktype{2};']);
        end
        usuccess = 1;
    catch ne_eo;
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Couldn''t set <set> field ' fieldfld ' (' ne_eo.message ').']);
    end

% sub-indexed on/off switch
case {'sb', 'subbool', 'soo', 'subonoff'}

    % reject bad link spec
    if length(flinktype) < 2 || ...
       ~isnumeric(flinktype{2}) || ...
        isempty(flinktype{2}) || ...
       ~isnumeric(fieldval) || ...
        isnan(flinktype{2}(1)) || ...
        isinf(flinktype{2}(1)) || ...
        flinktype{2}(1) < 1
        disp(['xfigure::UpdateField:*i_updatefield(...) -> ' ...
              'Bad spec for <subonoff> link.']);
        return;
    end

    try
        cval = get(fieldhnm, 'Value');
        if all(cval == get(fieldhnm, 'Max'))
            fieldval(fix(flinktype{2}(1))) = 1;
        else
            fieldval(fix(flinktype{2}(1))) = 0;
        end
        eval([fieldevl '=fieldval;']);
        usuccess = 1;
    catch ne_eo;
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Couldn''t set <subonoff> field ' fieldfld ' (' ne_eo.message ').']);
    end

% sub-indexed indices
case {'si', 'subindex'}

    % reject bad link spec
    if length(flinktype) < 2 || ...
       ~isnumeric(flinktype{2}) || ...
        isempty(flinktype{2}) || ...
       ~isnumeric(fieldval) || ...
        isnan(flinktype{2}(1)) || ...
        isinf(flinktype{2}(1)) || ...
        flinktype{2}(1) < 1
        disp(['xfigure::UpdateField:*i_updatefield(...) -> ' ...
              'Bad spec for <subindex> link.']);
        return;
    end

    try
        cval = get(fieldhnm, 'Value');
        if numel(cval) ~= 1
            disp(['xfigure::UpdateField:*i_updatefield(...) -> ' ...
                  'Bad Value for <subindex> link field ' fieldfld '.']);
            return;
        end
        fieldval(fix(flinktype{2}(1))) = cval;
        eval([fieldevl '=fieldval;']);
        usuccess = 1;
    catch ne_eo;
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Couldn''t set <subindex> field ' fieldfld ' (' ne_eo.message ').']);
    end

% sub-indexed number (default: 0)
case {'sn','subnum','subnumstr','subnumeric'}

    % reject bad link spec
    if length(flinktype) < 2 || ...
       ~isnumeric(flinktype{2}) || ...
       isempty(flinktype{2}) || ...
      ~isnumeric(fieldval) || ...
       isnan(flinktype{2}(1)) || ...
       isinf(flinktype{2}(1)) || ...
       flinktype{2}(1) < 1
        disp(['xfigure::UpdateField:*i_updatefield(...) -> ' ...
              'Bad spec for <subnumeric> link.']);
        return;
    end

    try
        cstr = ne_methods.splittocellc(get(fieldhnm, 'String'), ' ,;', true, true);
        if ~isempty(cstr) && ...
            isempty(cstr{1})
            cstr(1) = [];
        end
        if isempty(cstr)
            tval = 0;
        else
            try
                tval = str2double(cstr{1});
                tval = tval(1);
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                tval = 0;
            end
        end
        fieldval(fix(flinktype{2}(1))) = tval;
        eval([fieldevl '=fieldval;']);
        usuccess = 1;
    catch ne_eo;
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Couldn''t set <subnumeric> field ' fieldfld ' (' ne_eo.message ').']);
    end

% sub-indexed set
case {'ss','subset'}

    % reject bad link spec
    if length(flinktype) < 3 || ...
       ~isnumeric(flinktype{2}) || ...
        isempty(flinktype{2}) || ...
       ~isnumeric(fieldval) || ...
        isnan(flinktype{2}(1)) || ...
        isinf(flinktype{2}(1)) || ...
        flinktype{2}(1) < 1
        disp(['xfigure::UpdateField:*i_updatefield(...) -> ' ...
              'Bad spec for <subset> link.']);
        return;
    end

    try
        cval = get(fieldhnm, 'Value');
        if all(cval == get(fieldhnm, 'Max'));
            fieldval(fix(flinktype{2}(1))) = flinktype{3};
            eval([fieldevl '=fieldval;']);
        end
        usuccess = 1;
    catch ne_eo;
        disp(['xfigure::UpdateField*:i_updatefield(...) -> ' ...
              'Couldn''t set <subset> field ' fieldfld ' (' ne_eo.message ').']);
    end
end
