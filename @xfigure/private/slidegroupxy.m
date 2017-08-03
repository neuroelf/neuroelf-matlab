function slidegroupxy(xo, iStr, spos, varargin)
%XFIGURE::SLIDEGROUPXY  Slide a group of uicontrol elements.
%   SLIDEGROUP(FIG, GROUP, NEWPOS) slides the uicontrol elements in group
%   GROUP within figure FIG to (relative) NEWPOS (1x2 double).

% only valid for figures
if numel(xo) ~= 1 || xo.T ~= 1
    error('neuroelf:xfigure:invalidObjectType', ...
        'SlideGroupXY is only valid for figures.');
end

% check input spec
if nargin < 3 || ~isvarname(deblank(iStr)) || ~isnumeric(spos) || ...
   (length(spos) < 2 && (nargin < 4 || ~isnumeric(varargin{1}) || isempty(varargin{1})))
    error('neuroelf:xfigure:badArgument', ...
        'SlideGroupXY needs a group and a X and Y slide parameter.');
end
dogroup = deblank(iStr(:)');

% what slide spec
if numel(epos) > 1
    slide = [epos(1), epos(2)];
else
    slide = [epos(1), varargin{1}(1)];
end

% slide controls
if isfield(xo.X.figprops.sgroups, dogroup)
    hdls = xo.X.figprops.sgroups.(dogroup);
    for hc = 1:numel(hdls)
        setprop(hdls(hc), 'Position', [slide, 0, 0] + getprop(hdls(hc), 'Position'));
    end
end

% refresh ?
if ischar(varargin{nargin}) && strcmpi(varargin{nargin}, 'refresh')
    redrawfig(get(xo.H, 'Parent'));
end
