% loading field contents from link
case {'loadfields'}

% only valid for figures...
if hFigType ~= xfigsngl.objtypes.figure
    error( ...
        'xfigure:InvalidObjectType', ...
        'LoadFields is only valid for figures.' ...
    );
end

% group name(s)
if nargin < 3 || ...
   ~ischar(iStr) || ...
    isempty(iStr)
    error( ...
        'xfigure:BadArgument', ...
        'LoadFields requires the fieldgroup(s) to be named.' ...
    );
end

% with callback
if nargin < 4 || ...
   ~islogical(varargin{4}) || ...
    isempty(varargin{4}) || ...
    any(varargin{4})
    withcbs = true;
else
    withcbs = false;
end

% link disabled
if ~isstruct(iObj.lgroups)
    return;
end

% links not yet looked up
if ~isstruct(iObj.lilookup)
    xfigure(xo, 'LookupFields');
    iObj = xfigures(ihPos);
end

% get shortcuts and group names from input
lgroups  = iObj.figprops.lgroups;
lilookup = iObj.figprops.lilookup;
if strcmpi(iStr(:)', 'all_groups')
    fgnames = fieldnames(lgroups);
else
    fgnames = ne_methods.splittocellc(deblank(iStr(:)'), ',; ', true, true);
end

% iterate over given group names
for fgc = 1:numel(fgnames)

    % discard non-existing groups
    if ~isfield(lgroups, fgnames{fgc})
        continue;
    end

    % get group spec
    groupspec   = lgroups.(fgnames{fgc});
    uicfields   = groupspec{2};

    % iterate over uicontrols
    for fc = 1:numel(uicfields)

        % discard non-existing uicontrols
        if ~isfield(lilookup, uicfields{fc})
            continue;
        end

        % retrieve link information
        flink    = lilookup.(uicfields{fc});
        flinkuic = flink{3};

        % no longer available
        if ~isxfigure(flinkuic, 1)
            continue;
        end

        % perform load
        i_loadfield(flinkuic, flink{1:2}, withcbs);
    end
end
