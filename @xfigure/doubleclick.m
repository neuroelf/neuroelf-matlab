function doubleclick(xo, varargin)

% global config
global xfigsngl;

if numel(xo) ~= 1 || xo.T ~= 2
    warning('neuroelf:xfigure:invalidObjectType', 'DoubleClick is only valid for UIControls.');
    return;
end

% do nothing on empty callback list
if isempty(xo.X.callbacks)
    return;
end

% what callback should be processed
docb = 1;
if numel(xo.X.callbacks) > 2 && ~isempty(xo.X.callbacks{3})
    lastclick = xo.X.timeclick;
    mynow = now;
    drawnow;
    myvalue = get(xo.H, 'Value');
    xo.X.timeclick = {mynow, myvalue};
    if ((mynow - lastclick{1}) <= (xfigsngl.dblclickint) && isequal(myvalue, lastclick{2}))
        docb = 3;
    end
end
todocb = xo.X.callbacks{docb};

if ~isempty(todocb)
    try
        if isa(todocb, 'function_handle')
            feval(todocb, varargin{:});
        elseif iscell(todocb) && isa(todocb{1}, 'function_handle')
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
        if docb == 1
            try
                xo.X.prevprops = get(xo.H);
            catch xfigerror
                neuroelf_lasterr(xfigerror);
            end
        end
    catch xfigerror
        evalin('base', 'clear gcbc gcbf gcbo gcf this;', '');
        warning('neuroelf:xfigure:doubleClickFailed', ...
            'Error executing UIControl callback #%.0f: %s.\n%s%s', ...
            docb, char(todocb), 'Error message: ', xfigerror.message);
    end
end
