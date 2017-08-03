function o = addbarbutton(xo, iStr)
%XFIGURE::ADDBARBUTTON  Add a bar button element to a xbuttonbar uicontrol.
%   B = ADDBARBUTTON(XBB, BSTRUCT) adds a bar button to xbuttonbar XBB with
%   the settings in BSTRUCT.

% Version:  v1.1
% Build:    16042110
% Date:     Apr-21 2016, 10:18 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% global references and storage
global xfigmlup xfigsngl xfigures;

% only valid for figure and uipanel objects
if numel(xo) ~= 1 || xo.T ~= 2 || ~strcmpi(xo.X.loadprops.Type, 'xbuttonbar')
    error('neuroelf:xfigure:invalidObjectType', ...
        'Only xbuttonbar objects can be parents of bar buttons.');
end

% test iStr
if ~isstruct(iStr) || numel(iStr) ~= 1 || ~isfield(iStr, 'Caption')
    error('neuroelf:xfigure:badPropertyStruct', ...
        'BarButton properties must be of type struct with appropriate fields.');
end
if ~isfield(iStr, 'Position') || ~isnumeric(iStr.Position) || numel(iStr.Position) ~= 4 || ...
    any(isinf(iStr.Position(:)) | isnan(iStr.Position(:)) | iStr.Position(:) < 0)
    error('neuroelf:xfigure:badPropertyStruct', ...
        'BarButton properties requires valid Position argument.');
end
iStr.Position = iStr.Position(:)';
iPos = iStr.Position;

% use tag from iStr?
if ~isempty(iStr.Tag) && numel(iStr.Tag) < 28 && isrealvarname(iStr.Tag(:)')
    utag = iStr.Tag(:)';
else
    iStr.Tag = sprintf('BBT_%010.0f', floor(1e10 * rand(1, 1)));
    utag = ['xfigure_' iStr.Tag];
end

% get new object
o = xfigure('new');
o.T = -1; % keep as placeholder

% get self position
sPos = get(xo.H, 'Position');

% get background color
bgcol = get(get(xo.H, 'Parent'), 'Color');
if ~isnumeric(bgcol)
    bgcol = xfigsngl.figbgcolor;
end
bgimg = uint8(repmat(reshape(round(127.5 * (1 + bgcol)), [1, 1, 3]), iPos([4, 3])));

% border
if ~isfield(iStr, 'EdgeColor') || ~isa(iStr.EdgeColor, 'double') || numel(iStr.EdgeColor) ~= 3 || ...
    any(isinf(iStr.EdgeColor) | isnan(iStr.EdgeColor) | iStr.EdgeColor < 0 | iStr.EdgeColor > 1)
    iStr.EdgeColor = [0, 0, 0];
else
    iStr.EdgeColor = iStr.EdgeColor(:)';
end

% load image if necessary
if ischar(iStr.Caption)
    try
        iStr.Caption = imread(iStr.Caption);
        if size(iStr.Caption, 3) ~= 3
            iStr.Caption = repmat(iStr.Caption(:, :, 1), [1, 1, 3]);
        end
    catch xfigerror;
        neuroelf_lasterr(xfigerror);
        warning('neuroelf:xfigure:badImageFilename', 'Invalid image filename.');
        iStr.Caption = bgimg;
    end
end

% add frame?
capsz = size(iStr.Caption);
fr = floor(0.5 * (iPos([4, 3]) - capsz(1:2)));
if capsz(1) < iPos(4) && capsz(2) < iPos(3)
    cap = bgimg;
    cap(fr(1)+1:fr(1)+capsz(1), fr(2)+1:fr(2)+capsz(2), :) = iStr.Caption;
    cap(1, :, :) = 160;
    cap(end, :, :) = 160;
    cap(:, 1, :) = 160;
    cap(:, end, :) = 160;
    iStr.Caption = cap;
elseif capsz(1) > iPos(4) || capsz(2) > iPos(3)
    fr = min(0, fr);
    iStr.Caption = iStr.Caption(-fr(1)+1:-fr(1)+iPos(4), -fr(2)+1:-fr(2)+iPos(3), :);
    capsz = size(iStr.Caption);
    fr = floor(0.5 * (iPos([4, 3]) - capsz(1:2)));
    if capsz(1) < iPos(4)
        cap = bgimg;
        cap(fr(1)+1:fr(1)+capsz(1), fr(2)+1:fr(2)+capsz(2), :) = iStr.Caption;
        cap(1, :, :) = 160;
        cap(end, :, :) = 160;
        iStr.Caption = cap;
    elseif capsz(2) < iPos(3)
        cap = bgimg;
        cap(fr(1)+1:fr(1)+capsz(1), fr(2)+1:fr(2)+capsz(2), :) = iStr.Caption;
        cap(:, 1, :) = 160;
        cap(:, end, :) = 160;
        iStr.Caption = cap;
    end
end
iStr.ImageEnabled = iStr.Caption;

% generate different "styles" (disabled, clicked)
iStr.ImageDisabled = uint8(round(0.5 .* double(bgimg) + 0.3 .* repmat(mean(iStr.Caption, 3), [1, 1, 3])));
iStr.ImageClicked = uint8(max(0, round(0.75 .* double(iStr.Caption) - 16)));

% enabled
if ~isfield(iStr, 'Enabled') || ~ischar(iStr.Enabled) || ~strcmpi(iStr.Enabled(:)', 'off')
    cap = iStr.ImageEnabled;
else
    cap = iStr.ImageDisabled;
end

% create plot
o.H = surface('XData', (0.5 + iPos(1) - sPos(1)) + [0, iPos(3)], ...
    'YData', (0.5 + iPos(2) - sPos(2)) + [iPos(4); 0], 'CData', cap, ...
    'ZData', zeros(2, 2), 'EdgeColor', 'none', 'FaceColor', 'texturemap', 'Parent', xo.H);
dfcn = sprintf('delete(xfigure(hxdouble(''%s'')));', hxdouble(double(o.H)));
set(o.H, 'DeleteFcn', dfcn, 'HandleVisibility', xfigsngl.hvisible, 'Tag', utag, ...
    'UserData', iStr);

% output properties
iStr.Type = 'barbutton';
o.X.loadprops = iStr;
xfigmlup(end + 1) = o.H;
xfigures(end + 1) = o;

% add to parent list
xo.X.loadprops.ButtonList{end+1} = o;
o.X.loadprops.xtype = 'xbarbutton';

% add to global tag lookup struct
xfigsngl.tags.(['BBT_' iStr.Tag]) = o;
