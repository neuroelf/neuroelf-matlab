function addstring(xo, iStr, varargin)

% only valid for dropdown or listbox uicontrols
if numel(xo) ~= 1 || xo.T ~= 2 || ~any(strcmpi(get(xo.H, 'Style'), {'popupmenu', 'listbox'}))
    error('neuroelf:xfigure:invalidObjectType', ...
        'AddString is only valid for DropDown or ListBox UIControls.');
end

% only accept valid insertions
if ~ischar(iStr) && ~iscell(iStr)
    error('neuroelf:xfigure:badArgument', ...
        'AddString requires a CHAR or CELL argument to add.');
end

% positional argument given
if nargin < 3
    if ischar(iStr)
        pos = ones(1, size(iStr, 1)) * Inf;
    else
        pos = ones(1, numel(iStr)) * Inf;
    end
else
    pos = varargin{1}(:)';
end

% route through MString method
try
    mstring(xo, pos, iStr, 1);
catch xfigerror
    rethrow(xfigerror);
end
