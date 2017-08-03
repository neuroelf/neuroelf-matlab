function tf = isactive(xo)

% check
if numel(xo) ~= 1 || ~any((2:4) == xo.T)
     error('neuroelf:xfigure:invalidObject', 'Call only valid for uicontrols, uimenus.');
end

% preset output
tf = false;
try
    if xo.T == 2
        tf = isequal(get(xo.H, 'Value'), get(xo.H, 'Max'));
    else
        tf = strcmpi(get(xo.H, 'Checked'), 'on');
    end
catch xfigerror
    neuroelf_lasterr(xfigerror);
end
