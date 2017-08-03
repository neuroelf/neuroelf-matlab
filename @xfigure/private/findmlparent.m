function [mlh, ht] = findmlparent(mlh,t,xfp_type)
if ~isnumeric(t)
    try
        t = xfp_type.(lower(t));
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        mlh = [];
        ht  = 0;
        return;
    end
end
try
    mlh = get(mlh, 'Parent');
    ht  = xfp_type.(lower(get(mlh, 'Type')));
    while ht ~= 1 && ...
       ~any(t == ht)
        mlh = get(mlh, 'Parent');
        try
            ht = xfp_type.(lower(get(mlh, 'Type')));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            mlh = [];
            ht  = 0;
            return;
        end
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    mlh = [];
    ht  = 0;
    return;
end
if ~any(t == 1) && ...
    ht == 1
    mlh = [];
    ht  = 0;
    return;
end
