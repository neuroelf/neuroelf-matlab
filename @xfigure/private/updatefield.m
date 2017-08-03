% write a value to a linked xini
case {'updatefield'}

% only valid for uicontrols
if hFigType ~= xfigsngl.objtypes.uicontrol
    error( ...
        'xfigure:InvalidObjectType', ...
        'UpdateField is only valid for UIControls.' ...
    );
end

% no link activated
if ~isstruct(iFObj.figprops.lilookup) || ...
    isempty(fieldnames(iFObj.figprops.lilookup))
    error( ...
        'xfigure:NoFieldLinkSupport', ...
        'Parent figure has no FieldLink support enabled.' ...
    );
end

% no Tag -> return
if isempty(iObj.figprops.loadprops.Tag)
    return;
end

% uicontrol is not linked
if ~isfield(iFObj.figprops.lilookup, iObj.loadprops.Tag)
    warning( ...
        'xfigure:NoFieldLink', ...
        'UIControl with tag %s is not linked.', ...
        iObj.loadprops.Tag ...
    );
    return;
end

% update field
varargout{1} = i_updatefield(xo, iFObj.figprops.lilookup.(iObj.loadprops.Tag){1:2});
