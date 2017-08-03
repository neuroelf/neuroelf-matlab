function o = setprop(xo, varargin)
%XFIGURE::SETPROP  Set object properties (passed on to underlying object)

% global references
global xfigsngl;

% try block
try
    % without arguments
    if nargin == 1
        if numel(xo) == 1
            o = set(xo.H);
        else
            o = set([xo.H]);
        end

    % single argument (list of values)
    elseif nargin == 2 && ~isstruct(varargin{1})
        if numel(xo) == 1
            o = set(xo.H, varargin{1});
        else
            o = set([xo.H], varargin{1});
        end

    % regular cases (at least a setting, value pair)
    else
        iStr = varargin{1};
        o = [];

        % struct
        if isstruct(iStr) && numel(iStr) == 1
            cfields = fieldnames(iStr);
            for fc = 1:numel(cfields)
                setprop(xo, cfields{fc}, iStr.(cfields{fc}));
            end
            return;
        end

        % regular mode
        if ischar(iStr)
            if numel(xo) == 1
                otype = lower(get(xo.H, 'Type'));
                if any(strcmp(otype, {'axes', 'image'}))
                    rtype = lower(xo.X.loadprops.xtype);
                else
                    rtype = otype;
                end
                try
                    if isfield(xfigsngl.aliases{xo.T + 2}, lower(iStr))
                        iStr = lower(xfigsngl.aliases{xo.T + 2}.(lower(iStr)));
                    else
                        iStr = lower(iStr);
                    end
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    iStr = lower(iStr);
                end

                % handle different objects
                switch (otype)
                    case 'axes'

                        % successful
                        setsf = false;

                        % which property
                        switch iStr

                            % callback
                            case 'callback'

                                % try to set ButtonDownFcn
                                try
                                    if isempty(varargin{2})
                                        set(xo.X.uicprops.xchildren(1), 'ButtonDownFcn', '');
                                    else
                                        set(xo.X.uicprops.xchildren(1), 'ButtonDownFcn', ...
                                            {@xclick, xo});
                                        xo.X.loadprops.XClick = varargin{2};
                                    end
                                    setsf = true;
                                catch xfigerror
                                    neuroelf_lasterr(xfigerror);
                                    warning('neuroelf:xfigure:setPropError', ...
                                        'Error setting extended Callback: %s', ...
                                        xfigerror.message);
                                end

                            % image CData
                            case 'cdata'
                                if ~strcmp(rtype, 'ximage')
                                    error('neuroelf:xfigure:invalidObjectType', ...
                                        'Bad %s property: CData.', rtype);
                                end
                                try
                                    set(xo.X.uicprops.xchildren(1), 'CData', varargin{2});
                                    setsf = true;
                                    try
                                        set(xo.H, 'XLim', [0.5, 0.5 + size(varargin{2}, 2)]);
                                        set(xo.H, 'YLim', [0.5, 0.5 + size(varargin{2}, 1)]);
                                    catch xfigerror
                                        neuroelf_lasterr(xfigerror);
                                    end
                                catch xfigerror
                                    rethrow(xfigerror);
                                end

                            % enabled
                            case 'enable'

                                % for ximage
                                if strcmp(rtype, 'ximage')

                                    % get data
                                    imgpdata = xo.X.loadprops.ImageData;

                                    % disabled
                                    if nargin > 2 && ...
                                            ((ischar(varargin{2}) && strcmpi(varargin{2}(:)', 'off')) || ...
                                             (islogical(varargin{2}) && isequal(varargin{2}, false)) || ...
                                             (isa(varargin{2}, 'double') && isequal(varargin{2}, 0)))
                                         imgpdata = uint8(round(0.3 .* double(imgpdata) + ...
                                             155 .* xfigsngl.figbgcolor(1) + 25));
                                         imgbtd = '';
                                    else
                                         imgbtd = {@xclick, xo};
                                    end

                                    % set
                                    try
                                        set(xo.X.uicprops.xchildren(1), 'CData', imgpdata, ...
                                            'ButtonDownFcn', imgbtd);
                                        setsf = true;
                                    catch xfigerror
                                        neuroelf_lasterr(xfigerror);
                                    end

                                % for xlabels
                                elseif strcmp(rtype, 'xlabel')
                                    ncolor = xo.X.loadprops.xcolor;
                                    try
                                        if ~strcmpi(varargin{2}, 'on')
                                            bcolor = get(xo.H, 'Color');
                                            if ~isnumeric(bcolor) || isempty(bcolor)
                                                bcolor = xfigsngl.figbgcolor;
                                            end
                                            ncolor = (1.5 * bcolor + ncolor) * 0.4;
                                        end
                                        set(xo.X.uicprops.xchildren(1), 'Color', ncolor);
                                        setsf = true;
                                    catch xfigerror
                                        neuroelf_lasterr(xfigerror);
                                    end

                                % for xbarbuttons
                                elseif strcmp(rtype, 'xbarbutton')
                                end

                            % position
                            case {'position', 'units'}
                                try
                                    set(xo.H, iStr, varargin{2});
                                    setsf = true;
                                    if strcmp(rtype, 'xprogress')
                                        xfigure(xo, 'ProgressBar', NaN);
                                    end
                                catch xfigerror
                                    neuroelf_lasterr(xfigerror);
                                end

                            % otherwise
                            otherwise

                                % first try on children
                                for cobj = xo.X.uicprops.xchildren(:)'
                                    try
                                        set(cobj, iStr, varargin{2});
                                        setsf = true;
                                        if setsf
                                            break;
                                        end
                                    catch xfigerror
                                        neuroelf_lasterr(xfigerror);
                                    end
                                end

                                % yet on axis
                                if ~setsf
                                    try
                                        if isfield(xo.X.loadprops, 'ButtonList') && ...
                                           ~isempty(xo.X.loadprops.ButtonList)
                                            setprop([xo.X.loadprops.ButtonList{:}], iStr, varargin{2});
                                        else
                                            set(xo.H, iStr, varargin{2});
                                        end
                                        setsf = true;
                                    catch xfigerror
                                        neuroelf_lasterr(xfigerror);
                                    end
                                end
                        end

                        % error if not successful
                        if ~setsf
                           error('neuroelf:xfigure:badPropertyOrValue', ...
                               'Couldn''t set property %s on %s object type.', iStr, rtype);
                        end

                    % images
                    case 'image'
                        disp('image');

                    % other objects
                    otherwise

                        try
                            set(xo.H, iStr, varargin{2});
                        catch xfigerror
                            neuroelf_lasterr(xfigerror);
                            warning('neuroelf:xfigure:badPropertyOrValue', ...
                                'Couldn''t set property %s on %s object type.', iStr, rtype);
                        end
                end

            % more than one object, default
            else
                set([xo.H], varargin{:});
            end
            return;

        % invalid input
        else
            error('xfigure:BadArgument', 'Bad/missing argument provided for Set.');
        end        
    end

% deal with errors
catch xfigerror
    rethrow(xfigerror);
end
