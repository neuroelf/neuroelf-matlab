function o = adduicontextmenu(xo, iStr)

% global settings and registered storage
global xfigmlup xfigsngl xfigures;

% only valid for figure and uipanel objects
if numel(xo) ~= 1 || (xo.T ~= 1 && xo.T ~= 2)
    error('neuroelf:xfigure:invalidObjectType', ...
        'Only figures can be parents of uicontextmenus.');
end

% test iStr
if ~isstruct(iStr) || numel(iStr) ~= 1
    error('neuroelf:xfigure:badPropertyStruct', ...
        'Uicontextmenu properties must be of type struct.');
end

% perform some checks on struct
iStr = checkstruct(iStr, xfigsngl.optuix);

% use tag from iStr?
if ~isempty(iStr.Tag) && numel(iStr.Tag) < 28 && isrealvarname(iStr.Tag(:)')
    utag = iStr.Tag(:)';
else
    iStr.Tag = sprintf('UIX_%010.0f', uuid);
    utag = ['xfigure_' iStr.Tag];
end

% get new object
o = xfigure('new');

% prepare struct for call (still needs DeleteFcn!)
oStr = struct('Parent', xo.H, 'Callback', iStr.Callback, 'Interruptible', iStr.Interrupts, ...
    'Tag', utag, 'UserData', iStr.UserData);

% create object and fill/update fields as necessary
o.H = uicontextmenu(oStr);
o.T = 4;

% complete object representation
oObj = makeostruct(xfigsngl.objtypes.uicontextmenu);
oObj.callbacks = {'', '', '', iStr.CallbackDelete};
oObj.loadprops = iStr;
o.X = oObj;
xfigmlup(end + 1) = o.H;
xfigures(end + 1)  = o;

% finally, add to global tag lookup struct
xfigsngl.tags.(['UIX_' iStr.Tag]) = cout;
