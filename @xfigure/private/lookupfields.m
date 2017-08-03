% perform one-time field lookup
case {'lookupfields'}

% only for figures...
if hFigType ~= xfigsngl.objtypes.figure
    error( ...
        'xfigure:InvalidObjectType', ...
        'LookupFields is only valid for figures.' ...
    );
end

% don't lookup twice !
if isstruct(iObj.figprops.lilookup)
    return;
end

% make lilookup field an empty struct...
xfigures(ihPos).figprops.lilookup = struct;

% if the FieldLink feature isn't properly enabled return
if ~isstruct(iObj.figprops.lgroups)
    return;
end

% go through all groups registered from FieldLinkSpec files
groupname = fieldnames(iObj.figprops.lgroups);
for gc = 1:numel(groupname)
    inispec = iObj.figprops.lgroups.(groupname{gc});
    specobj = iObj.figprops.linkspec{inispec{1}};

    % go through all UIC tag names in group
    specfield = inispec{2};
    for fc = 1:numel(specfield)
        speclab = deblank(specfield{fc});

        % unknown UIC Tag
        if ~isvarname(speclab) || ...
           ~isfield(xfigsngl.tags,['UIC_' speclab])
            continue;
        end

        % get object and check whether it belongs to this figure
        uicobject = xfigsngl.tags.(['UIC_' speclab]);
        if ~any(find(get(hFigMHnd, 'Children') == uicobject.mhnd))
            continue;
        end

        % get ini and test for availability
        try
            speccont = specobj.(groupname{gc}).(speclab);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            speccont = {};
        end
        if ~iscell(speccont) || ...
            numel(speccont) < 2 || ...
           ~iscell(speccont{1}) || ...
           ~iscell(speccont{2})
            warning( ...
                'xfigure:BadLinkSpec', ...
                'Link spec too short for UIControl (%s).', ...
                speclab ...
            );
            continue;
        end

        % try to resolve link
        speclink = speccont{1};
        if ~isfield(iObj.figprops.linkcont, speclink{1})
            continue;
        end
        inicont = iObj.figprops.linkcont.(speclink{1});
        try
            eval('inicont.(speclink{2}).(speclink{3});');
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            continue;
        end

        % combine link target and store it
        linktarget = {{inicont, speclink{2:3}}, speccont{2}, uicobject};
        xfigures(ihPos).figprops.lilookup.(speclab) = linktarget;
    end
end
