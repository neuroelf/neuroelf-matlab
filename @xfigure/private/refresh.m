function varargout = refresh(xo, varargin)
%XFIGURE::REFRESH  Redraw figure content.

% first argument must be xfigure.
if ~isa(xo, 'xfigure')
    error('neuroelf:xfigure:badArgument', 'Invalid call to xfigure::refresh.');
end

% safely pass on yo redraw
try
    if nargout > 0
        [varargout{1:nargout}] = redraw(xo, varargin{:});
    else
        redraw(xo, varargin{:});
    end
catch xfigerror
    rethrow(xfigerror);
end
