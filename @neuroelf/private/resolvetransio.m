function d = resolvetransio(d)
% resolvetransio  - resolve transio from data
%
% FORMAT:       d = resolvetransio(d)
%
% Input fields:
%
%       d           struct, cell or transio
%
% Output fields:
%
%       d           input field with resolved data

% argument check
if nargin < 1
    error( ...
        'neuroelf:TooFewArguments', ...
        'Missing argument.' ...
    );
end

% get class
dc = builtin('class', d);

% return if not necessary
if isempty(d) || ...
   ~any('cstx' == lower(dc(1))) || ...
    any(strcmpi(dc, {'char', 'single', 'xfigure'}))
    return;
end

% transio itself
if istransio(d)
    d = resolve(d);

% cell array
elseif iscell(d)

    % iterate over cells
    try
        for cc = 1:numel(d)
            if ~isempty(d{cc})
                d{cc} = resolvetransio(d{cc});
            end
        end
    catch ne_eo;
        rethrow(ne_eo);
    end

% struct array
elseif isstruct(d) || ...
   (numel(d) == 1 && ...
    isxff(d, true))

    % iterate over structs with fields
    if isstruct(d)
        f = fieldnames(d);
    else
        f = fieldnames(getcont(d));
    end
    try
        for sc = 1:numel(d)
            for fc = 1:numel(f)
                if ~isempty(d(sc).(f{fc}))
                    dc = builtin('class', d(sc).(f{fc}));
                    if any('cstx' == lower(dc(1))) && ...
                       ~any(strcmpi(dc, {'char', 'single', 'xfigure'}))
                        d(sc).(f{fc}) = resolvetransio(d(sc).(f{fc}));
                    end
                end
            end
        end
    catch ne_eo;
        rethrow(ne_eo);
    end

% otherwise, leave unchanged
end
