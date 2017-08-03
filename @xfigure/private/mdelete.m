function mdelete(xo, iStr, varargin)

% only valid for dropdown or listbox uicontrols
if numel(xo) ~= 1 || xo.T ~= 2 || ~any(strcmpi(get(xo.H, 'Style'), {'popupmenu', 'listbox'}))
    error('neuroelf:xfigure:invalidObjectType', ...
        'MDelete is only valid for DropDown or ListBox UIControls.');
end

% require numeric requests
if nargin < 2 || isempty(iStr) || ~isnumeric(iStr)
    error('neuroelf:xfigure:badArgument', 'MDelete requires a numeric argument.');
end

% get current content and make sure it's a cell array
tstring = get(xo.H, 'String');
if ~iscell(tstring)
    tstring = cellstr(tstring);
    waschar = true;
else
    waschar = false;
end
tsize = numel(tstring);
iStr  = fix(iStr(:)');

% no outofbounds, please
iStr(iStr < 1 | iStr > tsize | isnan(iStr) | isinf(iStr)) = [];
if isempty(iStr)
    return;
end

% delete requested items and set selection accordingly
tstring(iStr) = [];
islistbox = strcmpi(get(xo.H, 'Style'), 'listbox');
if islistbox && get(xo.H, 'Max') > (get(xo.H, 'Min') + 1)
    tvalue = intersect(get(xo.H, 'Value'), (1:numel(tstring)));
    if isempty(tvalue) && numel(tstring) > 0
        tvalue = numel(tstring);
    end
else
    if ~isempty(tstring)
        tvalue = max(1, min(numel(tstring), iStr(1)));
    else
        tvalue = [];
    end
end

% convert and set string and value
if islistbox
    set(xo.H, 'ListboxTop', max(1, min(get(xo.H, 'ListboxTop', numel(tstring)))));
end
if waschar
    tstring = char(tstring);
end
set(xo.H, 'String', tstring, 'Value', tvalue);
