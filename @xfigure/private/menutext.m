function menutxt = menutext(h, ind)

% global library
global ne_methods;

menutxt = '';
if nargin < 2 || isempty(ind)
    ind = '';
end
if ~ishandle(h) || ~isvalid(h)
    return;
end
t = lower(get(h, 'Type'));
if ~any(strcmp(t, {'uicontextmenu', 'uimenu'}))
    return;
end
try
    c = get(h, 'Children');
catch xfigerror
    neuroelf_lasterr(xfigerror);
    c = [];
end
if strcmp(t, 'uimenu')
    try
        mp = ne_methods.subget(h, {'Label', 'Separator', 'Callback', 'Enable'});
    catch xfigerror
        neuroelf_lasterr(xfigerror);
        mp = struct('Label', '', 'Separator', 'off', 'Callback', '');
    end
    if strcmpi(mp.Separator, 'on')
        menutxt = [menutxt ind '------------' char(10)];
    end
    if ~isempty(mp.Callback)
        if ischar(mp.Callback)
            mp.Callback = [' -> [' mp.Callback ']'];
        elseif iscell(mp.Callback)
            if numel(mp.Callback) > 1
                try
                    cbargs = any2ascii(mp.Callback(2:end));
                    cbargs = cbargs(2:end-1);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    cbargs = '';
                end
            else
                cbargs = '';
            end

            if ischar(mp.Callback{1})
                mp.Callback = [' -> [' mp.Callback{1} '(' cbargs ')]'];
            elseif isa(mp.Callback{1}, 'function_handle')
                mp.Callback = [' -> [@' func2str(mp.Callback{1}) '(' cbargs ')]'];
            else
                mp.Callback = ' (#! bad Callback cell array !)';
            end
        elseif isa(mp.Callback, 'function_handle')
            mp.Callback = [' -> [@' func2str(mp.Callback) ']'];
        else
            mp.Callback = ' (#! bad Callback property !)';
        end
    else
        mp.Callback = '<no callback>';
    end
    if ~isempty(mp.Label)
        if strcmpi(mp.Enable, 'on')
            menutxt = [menutxt ind mp.Label mp.Callback char(10)];
        else
            menutxt = [menutxt ind mp.Label '(disabled/gray) ' mp.Callback char(10)];
        end
    end
    ind = [ind '  '];
end
for cc = length(c):-1:1
    menutxt = [menutxt menutext(c(cc), ind)];
end
