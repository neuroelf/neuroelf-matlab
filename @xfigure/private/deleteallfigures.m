function xo = deleteallfigures(xo, varargin)

% global storage
global xfigmlup xfigures;

% we need all handles...
shht = get(0, 'ShowHiddenHandles');
set(0, 'ShowHiddenHandles', 'on');

% find figures and delete them
figs = findobj;
figtype = get(figs, 'Type');
dfigs = false(size(figtype));
for ft = {'figure', 'axes', 'uicontrol', 'uimenu', 'image'}
    dfigs = dfigs | strcmpi(figtype, ft{1});
end

% delete all but root
set(figs(dfigs), 'DeleteFcn', '');
delete(findobj('type', 'figure'));

% re-init arrays
xfigmlup(2:end) = [];
xfigures(2:end)  = [];

% set
set(0, 'ShowHiddenHandles', shht);
