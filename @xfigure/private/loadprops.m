function p = loadprops(xo)

% only valid for uicontrol objects
if numel(xo) ~= 1 || xo.T == 0
    warning('neuroelf:xfigure:invalidObjectType', 'LoadProps not valid for ROOT object.');
    p = struct;
else
    p = xo.X.loadprops;
end
