function po = parentobject(xo)

% global storage
global xfigures;

% check
if numel(xo) ~= 1
    error('neuroelf:xfigure:badInput', 'Requires 1x1 input.');
end

% for figures
if xo.T <= 1
    po = xfigures(1);

% for anything else
else
    try
        po = xfigure(get(xo.H, 'Parent'));
    catch xfigerror;
        neuroelf_lasterr(xfigerror);
        error('neuroelf:xfigure:objectDisappeared', ...
            'The xfigure/MATLAB UI object is no longer available.');
    end
end
