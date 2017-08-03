function docallback(xo, iStr)

% check
if numel(xo) ~= 1 || xo.T ~= 2
    warning('neuroelf:xfigure:invalidObjectType', ...
        'DoCallback is only valid for UIControls.');
    return;
end

% get callback string and evaluate in base workspace
cbstr = get(xo.H, 'Callback');
if ~isempty(cbstr)
    fp = ancestor(xo.H, 'figure');
    assignin('base', 'gcbf', fp);
    assignin('base', 'gcbo', xo.H);
    assignin('base', 'gcf',  fp);
    assignin('base', 'this', xo);
    try
        evalin('base', cbstr);
        xo.X.prevprops = get(xo.H);
    catch xfigerror
        warning('neuroelf:xfigure:callbackFailed', ...
            'Error executing UIControl callback: %s.', xfigerror.message);
    end
    evalin('base', 'clear gcbf gcbo gcf this;', '');
end
