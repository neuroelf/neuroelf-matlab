function tf = isenabled(xo)

% check
if numel(xo) ~= 1 || ~any((2:4) == xo.T)
    error('neuroelf:xfigure:invalidObject', 'Invalid object for IsEnabled.');
end

% default: true
tf = true;
try

    % override
    if ~strcmpi(get(xo.H, 'Enable'), 'on')
        tf = false;
    end
catch xfigerror
    neuroelf_lasterr(xfigerror);
end
