function setbarydata(xo, ydata)
%XFIGURE::SETBARYDATA  Update Y data of bars

% requires bar
if numel(xo) ~= 1 || xo.T ~= -1 || ~strcmpi(xo.X.loadprops.Type, 'barplot')
    error('neuroelf:xfigure:invalidObject', 'SetBarYData requires barplot object.');
end
p = xo.X.loadprops;
if nargin < 2 || ~isa(ydata, 'double') || numel(ydata) ~= (p.NumBars * p.NumGroups)
    error('neuroelf:xfigure:invalidYData', 'Invalid YData supplied.');
end

% prepare YData
ydata = ydata';
yp = ydata(:)';
yp = [zeros(1, numel(yp)); yp; yp; zeros(1, numel(yp))];
yp = yp(:);
set(xo.H, 'Vertices', [p.VerticesX, yp]);
xo.X.loadprops.VerticesY = yp;
