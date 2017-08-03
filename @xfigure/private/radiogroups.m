function g = radiogroups(xo)
%XFIGURE::RADIOGROUPS  Return the radiogroups (names) for a figure.
%   GROUPS = RADIOGROUPS(FIG) returns the list of radio (button) group
%   names for figure FIG.

% only for figures
if numel(xo) ~= 1 || xo.T ~= 1
    error('neuroelf:xfigure:invalidObjectType', 'RadioGroups is only valid for figures.');
end

% return radio groups
g = xo.X.figprops.rgroups;
