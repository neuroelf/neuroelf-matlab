function s = mstring(xo, iStr, varargin)

% only valid for dropdown or listbox uicontrols
if numel(xo) ~= 1 || xo.T ~= 2 || ~any(strcmpi(get(xo.H, 'Style'), {'popupmenu', 'listbox'}))
    error('neuroelf:xfigure:invalidObjectType', ...
        'MString is only valid for DropDown or ListBox UIControls.');
end

% which strings to insert where
if nargin < 2 || ~isnumeric(iStr) || isempty(iStr) || any(isnan(iStr))
    if nargin < 2
        positions = get(xo.H, 'Value');
    elseif ischar(iStr)
        positions = get(xo.H, 'Value');
        if isempty(positions)
            positions = 1;
        end
    elseif iscell(iStr)
        positions = ones(1, numel(iStr)) * Inf;
    else
        error('neuroelf:xfigure:badArgument', ...
            'First argument to MString must either be pos or setstr.');
    end
    if iscell(iStr)
        istring = iStr(:)';
    elseif ischar(iStr)
        istring = cellstr(iStr);
    else
        istring = {};
    end
else
    positions = fix(iStr(:)');
    if nargin > 2 && (iscell(varargin{1}) || ischar(varargin{1}))
        istring = varargin{1};
        if iscell(istring)
            istring = istring(:)';
        else
            istring = cellstr(istring);
        end
    else
        istring = {};
    end
end

% get current string
rstring = get(xo.H, 'String');
if ~iscell(rstring)
    if ~isempty(rstring)
        rstring = cellstr(rstring);
    else
        rstring = {};
    end
    waschar = true;
else
    waschar = false;
end
rsize = numel(rstring);

% no insertion
if isempty(istring)
    positions(positions < 1 | positions > rsize) = [];
    s = rstring(positions);
    if length(positions) == 1
        s = s{1};
    end

% check for insertion
else
    if (islogical(varargin{end}) || isnumeric(varargin{end})) && ...
       ~isempty(varargin{end}) && varargin{end}(1)
        doinsert = true;
        positions(positions > rsize) = rsize + 1;
    else
        doinsert = false;
        positions(positions < 1 | positions > rsize) = [];
    end
    isize = numel(istring);
    psize = numel(positions);

    if ~doinsert
        icount = 1;
        while icount <= isize && icount <= psize
            rstring(positions(icount)) = istring(icount);
            icount = icount + 1;
        end

    else
        icount = min(isize, psize);
        for icount = icount:-1:1
            if isempty(rstring)
                rstring = istring(icount);
            else
                tpos = positions(icount);
                if tpos < 2
                    rstring = [istring(icount); rstring];
                elseif tpos > length(rstring)
                    rstring = [rstring; istring(icount)];
                else
                    rstring = [rstring(1:(tpos-1)); istring(icount); rstring(tpos:end)];
                end
            end
        end
        tposition = max(min([positions, length(rstring)]), 1);
    end

    for cc = 1:length(rstring)
        if isempty(rstring{cc})
            rstring{cc} = '';
        end
    end

    % conversion
    if waschar
        rstring = char(rstring);
    end

    % set string and value
    if doinsert
        set(xo.H, 'String', rstring, 'Value', tposition);
    else
        set(xo.H, 'String', rstring);
    end
end
