function o = getprop(xo, varargin)
%XFIGURE::GETPROP  Get object properties (passed on to underlying object)

% global references
global xfigsngl;

% try block
try
    % without arguments
    if nargin == 1
        if numel(xo) == 1
            o = get(xo.H);
        else
            o = get([xo.H]);
        end
    else
        iStr = varargin{1};
        if numel(xo) == 1
            otype = lower(get(xo.H, 'Type'));
            if strcmp(otype, 'axes')
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
            catch xfigerror
                neuroelf_lasterr(xfigerror);
                iStr = lower(iStr);
            end
            
            % handle different objects
            switch (otype)

                % for special objects
                case 'axes'

                    % progess
                    switch (iStr)
                        case 'progress'
                            if ~strcmp(rtype, 'xprogress')
                                error('neuroelf:xfigure:invalidObjectType', ...
                                    'Invalid uicontrol/%s object property: Progress.', rtype);
                            end
                            o = iObj.uicprops.progress;

                        case 'visible'
                            try
                                o = get(xo.X.uicprops.xchildren(1), 'Visible');
                            catch xfigerror
                                error('neuroelf:xfigure:badChild', ...
                                    'Couldn''t get Visible property for axes child (%s).', ...
                                    xfigerror.message);
                            end

                        % otherwise
                        otherwise
                            try
                                o = get(xo.H, iStr);
                            catch xfigerror
                                neuroelf_lasterr(xfigerror);
                                getsf = false;
                                for cc = xo.X.uicprops.xchildren
                                    try
                                        o = get(cc, iStr);
                                        getsf = true;
                                        break;
                                    catch xfigerror
                                        neuroelf_lasterr(xfigerror);
                                    end
                                end
                                if ~getsf
                                    error('neuroelf:xfigure:invalidProperty', ...
                                        'Couldn''t get %s property from axes child.', iStr);
                                end
                            end
                    end

                % otherwise
                otherwise
                    try
                        o = get(xo.H, iStr);

                         % invalid values ? (Matlab R2009b patch)
                         if strcmp(iStr, 'value')
                             switch (otype)
                                 case 'dropdown'
                                     o = unique(max(1, o));
                             end
                         end
                    catch xfigerror
                        neuroelf_lasterr(xfigerror);
                        error('neuroelf:xfigure:invalidProperty', ...
                            'Invalid %s property: %s.', rtype, iStr);
                    end
            end

        % try to retrieve multiple props
        elseif isstruct(iStr) && numel(iStr) == 1
            o = iStr;
            cfields = fieldnames(o);
            for fc = 1:numel(cfields)
                try
                    o.(cfields{fc}) = get(xo.H, cfields{fc});
                catch xfigerror
                    neuroelf_lasterr(xfigerror);
                    warning('neuroelf:xfigure:invalidProperty', ...
                        'Invalid property (%s) for object of type %s.', cfields{fc}, ...
                        get(xo.H, 'Type'));
                    o = rmfield(o, cfields{fc});
                end
            end

        % multiple objects (all bets are off)
        else
            o = get([xo.H], varargin{1:end});
        end
    end
catch xfigerror
    rethrow(xfigerror);
end
