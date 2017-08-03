function t = tagstruct(xo)
%XFIGURE::TAGSTRUCT  Return TagStruct from UserData struct.

% only valid for figures
if numel(xo) ~= 1 || xo.T ~= 1
    error('neuroelf:xfigure:invalidObjectType', 'TagStruct is only valid for figures.');
end

% get (and set) struct
udt = get(xo.H, 'UserData');
if isstruct(udt) && ~isempty(udt) && isfield(udt, 'xfigure_UDTagStruct')
    t = udt.xfigure_UDTagStruct;
else
    t = struct;
end
