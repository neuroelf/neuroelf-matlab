function tf = isvisible(xo)

% check
if numel(xo) ~= 1 || xo.T == 0
    error('neuroelf:xfigure:invalidObject', 'Invalid object for IsVisible.');
end

tf = true;
try
    if ~strcmpi(get(xo.H, 'Visible'), 'on')
        tf = false;
    end
catch xfigerror
    neuroelf_lasterr(xfigerror);
end
