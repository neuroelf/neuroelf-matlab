function fo = findobject(xo, tag, type)

% global config
global xfigsngl;

% we need valid input (Tag name)
if numel(xo) ~= 1 || xo.T ~= 0 || ~isrealvarname(tag)
    error('neuroelf:xfigure:badArgument', ...
        'FindObject requires a valid Tag argument.');
end

% only tag
if nargin < 3 || ~ischar(type)
    findtag = ['UIC_' tag];
else
    findtag = deblank(type(:)');
    if isvarname(findtag)
        switch lower(tag)
        case {'c', 'uic', 'uicontrol'}
            findtag = ['UIC_' findtag];
        case {'f', 'fig', 'figure'}
            findtag = ['FIG_' findtag];
        case {'m', 'uim', 'uimenu'}
            findtag = ['UIM_' findtag];
        case {'x', 'uix', 'uicontextmenu'}
            findtag = ['UIX_' findtag];
        otherwise
            error('neuroelf:xfigure:badArgument', ...
                'Invalid search type: %s.', tag);
        end
    else
        error('neuroelf:xfigure:badArgument', ...
            'Invalid tag: %s.', type(:)');
    end
end

% try lookup
try
    fo = xfigsngl.tags.(findtag);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    warning('neuroelf:xfigure:lookupFailed', ...
        'Lookup of Tag failed: %s.', tag);
end
