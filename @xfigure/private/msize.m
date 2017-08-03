function s = msize(xo)

% only valid for dropdown or listbox uicontrols
if numel(xo) ~= 1 || xo.T ~= 2 || ~any(strcmpi(get(xo.H, 'style'), {'popupmenu', 'listbox'}))
    error('neuroelf:xfigure:invalidObjectType', ...
        'MSize is only valid for DropDown or ListBox UIControls.');
end

% get correct length
tstring = get(xo.H, 'String');
if iscell(tstring)
    s = numel(tstring);
else
    s = size(tstring, 1);
end
