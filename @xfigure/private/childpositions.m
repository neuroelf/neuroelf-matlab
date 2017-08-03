function [ph, vh, p] = childpositions(xo)

% only for figures
if numel(xo) ~= 1 || xo.T ~= 1
    error('neuroelf:xfigure:invalidObject', 'ChildPositions requires a 1x1 figure object.');
end

% get children
shh = get(0, 'ShowHiddenHandles');
set(0, 'ShowHiddenHandles', 'on');
ch = get(xo.H, 'Children');
set(0, 'ShowHiddenHandles', shh);

% get types
ct = get(ch, 'Type');

% remove all but uicontrols and axes
keepa = strcmpi(ct, 'axes');
keeph = (keepa | strcmpi(ct, 'uicontrol'));
ph = ch(keeph);
if nargout < 2
    return;
end

% visibility handles
vh = ph;
keepa = find(keepa(keeph));
ch = get(vh(keepa), 'Children');
for hc = 1:numel(keepa)
    if ~isempty(ch{hc})
        vh(keepa(hc)) = ch{hc}(1);
    end
end
if nargout < 3
    return;
end

% get positions
p = get(ph, 'Position');
if iscell(p)
    p = cat(1, p{:});
end
