function xclick(xo, varargin)

% override input
if nargin > 2 && isa(varargin{2}, 'xfigure') && numel(varargin{2}) == 1
    xo = varargin{2};
    varargin = varargin(3:end);
end

% check
if numel(xo) ~= 1 || xo.T ~= 2 || ~strcmpi(get(xo.H, 'Type'), 'axes')
    error('neuroelf:xfigure:badArgument', 'Invalid Object for XClick.');
end

% nothing to do
if ~isfield(xo.X.loadprops, 'XClick') || isempty(xo.X.loadprops.XClick)
    return;
end
xc = xo.X.loadprops.XClick;

% udpate image
imgcdata = uint8(max(0, round(0.7 .* double(xo.X.loadprops.ImageData) - 16)));

% what to do
try
    ph = ancestor(xo.H, 'figure');
    set(xo.X.uicprops.xchildren(1), 'CData', imgcdata);
    drawnow;
    if ischar(xc)
        assignin('base', 'gcbf', ph);
        assignin('base', 'gcbo', xo.X.uicprops.xchildren(1));
        assignin('base', 'gcf', ph);
        assignin('base', 'this', xo);
        evalin('base', xc(:)');
        evalin('base', 'clear gcbf gcbo gcf this;');
    elseif isa(xc, 'function_handle')
        feval(xc, varargin{:});
    elseif iscell(xc)
        if ischar(xc{1})
            assignin('base', 'gcbf', ph);
            assignin('base', 'gcbo', xo.X.uicprops.xchildren(1));
            assignin('base', 'gcf', ph);
            assignin('base', 'this', xo);
            evalin('base', xc{1}(:)');
            evalin('base', 'clear gcbf gcbo gcf this;');
        elseif isa(xc{1}, 'function_handle')
            feval(xc{1}, xo.X.uicprops.xchildren(1), [], xc{2:end});
        end
    end
catch xfigerror
    neuroelf_lasterr(xfigerror);
    warning('neuroelf:xfigure:xClickError', 'Error executing callback: %s.', ...
        xfigerror.message);
end
set(xo.X.uicprops.xchildren(1), 'CData', xo.X.loadprops.ImageData);
drawnow;
