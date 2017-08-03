% write all linked values to linked xini object
case {'updatefields'}

% only valid for figures
if hFigType ~= xfigsngl.objtypes.figure
    error( ...
        'xfigure:InvalidObjectType', ...
        'UpdateFields is only valid for figures.' ...
    );
end

% check groups
if ~ischar(iStr) || ...
    numel(deblank(iStr)) < 3
    error( ...
        'xfigure:BadArgument', ...
        'No valid group selection given.' ...
    );
end
iStr = iStr(:)';

% group selection
if strcmpi(iStr, 'all_groups')
    ugroups = fieldnames(iObj.figprops.lgroups);
else
    ugroups = splittocellc(deblank(iStr), ',; ', true, true);
end
usuccess = true;

% iterate over groups
for gc = 1:length(ugroups)

    % only go on if group exists
    try
        gcontrols = iObj.figprops.lgroups.(ugroups{gc});
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        continue;
    end

    % iterate over controls
    gcontrols = gcontrols{2};
    for cc = 1:length(gcontrols)
        try
            glookup   = iObj.figprops.lilookup.(gcontrols{cc});
            cusuccess = xfigure(glookup{end}, 'UpdateField');
            if ~cusuccess
                usuccess = false;
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            usuccess = false;
        end
    end
end

% successful ?
varargout{1} = usuccess;
