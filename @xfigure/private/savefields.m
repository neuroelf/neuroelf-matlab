% storing state of fields into xini object
case {'savefields'}

% only valid for figures and uicontrols
if ~any([xfigsngl.objtypes.figure, ...
         xfigsngl.objtypes.uicontrol] == hFigType)
    error( ...
        'xfigure:InvalidObjectType', ...
        'SaveFields is only valid for figures or UIControls.' ...
    );
end

% get correct object
if hFigType == xfigsngl.objtypes.figure
    rObj = iObj;
else
    rObj = iFObj;
end

% any linked fields
if ~isstruct(rObj.lilookup) || ...
    isempty(fieldnames(rObj.lilookup))
    return;
end

% update fields to xini first
usuccess = true;
if ischar(iStr) || ...
    iscell(iStr)
    usuccess = xfigure(rObj.handle, 'updatefields', iStr);
end

% then write the files back to disk
contfiles = fieldnames(rObj.linkcont);
for fc = 1:length(contfiles)
    ifile = rObj.linkcont.(contfiles{fc});
    try
        if exist(Filename(ifile), 'file') == 2
            WriteIniFile(ifile);
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        usuccess = false;
    end
end

% report result
varargout{1} = usuccess;
