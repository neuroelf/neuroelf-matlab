function mmove(xo, iStr, varargin)

% only valid for dropdown or listbox uicontrols
if numel(xo) ~= 1 || xo.T ~= 2 || ~any(strcmpi(get(xo.H, 'Style'), {'popupmenu', 'listbox'}))
    error('neuroelf:xfigure:invalidObjectType', ...
        'MMove is only valid for DropDown or ListBox UIControls.');
end

% check input arguments
if nargin < 2 || isempty(iStr) || ~isnumeric(iStr) || fix(iStr(1)) == 0 || ...
    isinf(iStr(1)) || isnan(iStr(1))
    error('neuroelf:xfigure:badArgument', ...
        'MMove requires a numeric argument.');
end

% what to do
if nargin > 2 && isnumeric(varargin{1}) && ~isempty(varargin{1})
    mindex  = unique(varargin{1}(:));
    msource = false;
else
    mindex  = sort(get(xo.H, 'Value'));
    msource = true;
end

% get current content
tstring = get(xo.H, 'String');
if ~iscell(tstring)
    tstring = cellstr(tstring);
    waschar = true;
else
    waschar = false;
end
tlength = length(tstring);
iStr    = fix(iStr(:)');

% make sure selection is valid
mindex(mindex < 1 | mindex > tlength | isnan(mindex) | isinf(mindex)) = [];
if isempty(mindex)
    return;
end
mlength = numel(mindex);

% where to move selection
moveto = fix(iStr(1));

% prepare target array and selection
tgstring = {};
tindex   = mindex + moveto;
while any(tindex < 1)
    tindex = tindex + 1;
end
while any(tindex > tlength)
    tindex = tindex - 1;
end
if all(tindex == mindex)
    return;
end

% what remains without moved strings
rindex = setdiff(1:tlength, mindex(:)');
oindex = setdiff(1:tlength, tindex(:)');
rlength = numel(rindex);

% compile and convert new array, set value
tgstring(oindex(1:rlength)) = tstring(rindex(1:rlength));
tgstring(tindex(1:mlength)) = tstring(mindex(1:mlength));
if waschar
    tgstring = char(tgstring);
end
set(xo.H, 'String', tgstring);
if msource
    set(xo.H, 'Value', tindex);
end
