function delete(xo)

% global reference storage
global xfigmlup xfigures;

% remove invalid handles
xo(~isvalid(xo)) = [];
if isempty(xo)
    return;
end

% get handles
h = [xo.H];

% iterate over objects
for oc = 1:numel(xo)
    
    % don't do anything if no longer valid
    if ~ishandle(h(oc))
        continue;
    elseif ~isnumeric(h) && ~isvalid(h(oc))
        continue;
    end

    % anything but the root object
    if xo(c).T > 0

        % remove DeleteFcn
        set(xo(c).H, 'DeleteFcn', '');

        % execute original DeleteFcn
        if numel(xo(c).X.callbacks) > 3 && ~isempty(xo(c).X.callbacks{4})
            try
                evalin('base', xo(c).X.callbacks{4});
            catch xfigerror
                warning(xfigerror.message);
            end
        end

        % delete this (if not already in process)
        try
            if ~strcmpi(get(h(oc), 'BeingDeleted'), 'on')
                delete(h(oc));
            end
        catch xfigerror
            neuroelf_lasterr(xfigerror);
        end
    end
end

% remove
[h, ih] = intersect(xfigmlup(:), h(:));
xfigmlup(ih) = [];
xfigures(ih) = [];
