function o = rgroupclick(xo, varargin)
%XFIGURE::RGROUPCLICK  Update radio group correctly

% only valid for uicontrols
if numel(xo) ~= 1 || xo.T ~= 2 || ~strcmpi(get(xo.H, 'Type'), 'uicontrol') || ...
   ~strcmpi(get(xo.H, 'Style'), 'radiobutton')
    error('neuroelf:xfigure:invalidObjectType', ...
        'The RadioButtonClick event is only valid for RadioButtons.');
end

% set all other buttons to off first !
radiogroupsetone(xo);

% no further callback named -> return
if isempty(xo.X.callbacks) || isempty(xo.X.callbacks{1})
    return;
end

% do callback
todocb = xo.X.callbacks{1};
if ~isempty(todocb)
    try
        if isa(todocb, 'function_handle')
            feval(todocb, varargin{:});
        elseif iscell(todocb) && isa(todocb{1}, 'function_handle');
            feval(todocb{1}, xo.H, 0, todocb{2:end}, varargin{:});
        elseif ischar(todocb)
            mygcbf = ancestor(xo.H, 'figure');
            assignin('base', 'gcbf', mygcbf);
            assignin('base', 'gcbo', xo.H);
            assignin('base', 'gcf',  mygcbf);
            assignin('base', 'this', xo);
            evalin('base', todocb);
            evalin('base', 'clear gcbc gcbf gcbo gcf this;', '');
        end
        try
            xo.X.prevprops = get(xo.H);
        catch xfigerror
            neuroelf_lasterr(xfigerror);
        end
    catch xfigerror
        evalin('base', 'clear gcbc gcbf gcbo gcf this;', '');
        warning('neuroelf:xfigure:rGroupClickFailed', ...
            'Error executing UIControl callback #%.0f: %s.\n%s%s', ...
            1, char(todocb), 'Error message: ', xfigerror.message);
    end
end
