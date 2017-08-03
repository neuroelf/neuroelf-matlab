function redraw(xo)

% check
if numel(xo) ~= 1
    error('neuroelf:xfigure:badObject', 'Requires 1x1 object.');

% for figures
elseif xo.T == 1
    redrawfig(xo.H);

% for non root objects
elseif xo.T > 0
    redrawfig(ancestor(xo.H, 'figure'));

% give an error for the root object
else
    error('neuroelf:xfigure:invalidObjectType', ...
        'Redraw is only valid for non-root objects.');
end
