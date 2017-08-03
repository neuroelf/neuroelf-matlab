function t = isvalidxfigure(xo)
%XFIGURE::ISVALIDXFIGURE  Return valid status for xfigure object

% first pass
t = isvalid(xo);

% single object
if numel(t) == 1 && t
    t = (ishandle(xo.H) && (isnumeric(xo.H) || isvalid(xo.H)));
    return;
end

% then pass on (to also test underlying MATLAB UI handle)
if any(t(:))
    t(t(:)) = cellfun(@ishandle, {xo(t(:)).H});
    tn = t;
    t(t(:)) = cellfun(@isnumeric, {xo(t(:)).H});
    t(tn(:) & ~t(:)) = cellfun(@isvalid, {xo(tn(:) & ~t(:)).H});
end
