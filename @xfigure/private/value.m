function v = value(xo, varargin)
%XFIGURE::VALUE  Retrieve default value for object.

% only valid for single object
if numel(xo) ~= 1
    error('neuroelf:xfigure:multipleObjectError', 'Cannot retrieve value for multiple objects.');
end

% switch over type
switch (xo.T)

    % root
    case 0
        v = 'Off';
        if ~all(get(0, 'ScreenSize') == 1)
            v = 'On';
        end

    % figure
    case 1
        v = get(xo.H, 'Name');

    % uicontrol
    case 2

        % string-based type
        if any(strcmpi(get(xo.H, 'Style'), {'edit', 'frame', 'pushbutton', 'text'}))
            v = get(xo.H, 'String');
        else
            v = get(xo.H, 'Value');
            if nargin > 1
                try
                    v = v(varargin{:});
                catch xfigerror
                    warning('neuroelf:xfigure:subValueError', ...
                        'Couldn''t pass on subsref: %s', xfigerror.message);
                end
            end
        end

    % uimenu
    case 3
        v = get(hFigMHnd, 'Label');

    case 4
        v = 'Off';
        if strcmpi(get(0, 'Visible'), 'on')
            v = 'On';
        end

    otherwise
        error('neuroelf:xfigure:invalidObjType', ...
            'Invalid object type. No default value specified.');
end
