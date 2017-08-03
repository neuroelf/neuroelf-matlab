function o = addbarplot(xo, iStr)
%XFIGURE::ADDBARPLOT  Add a bar plot element to a xaxes uicontrol.
%   BP = ADDBARPLOT(XAXES, BPSTRUCT) adds a bar plot to xaxes XAXES with
%   the settings in BPSTRUCT.

% Version:  v1.1
% Build:    16042110
% Date:     Apr-21 2016, 10:18 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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
if numel(xo) ~= 1 || xo.T ~= 2 || ~strcmpi(xo.X.loadprops.Type, 'xaxes')
    error('neuroelf:xfigure:invalidObjectType', ...
        'Only xaxes can be parents of bar plots.');
end

% test iStr
if ~isstruct(iStr) || numel(iStr) ~= 1 || ~isfield(iStr, 'NumBars') || ~isfield(iStr, 'NumGroups')
    error('neuroelf:xfigure:badPropertyStruct', ...
        'Barplot properties must be of type struct with appropriate fields.');
end

% check props further
if ~isa(iStr.NumBars, 'double') || numel(iStr.NumBars) ~= 1 || ...
    isinf(iStr.NumBars) || isnan(iStr.NumBars) || iStr.NumBars < 1
    error('neuroelf:xfigure:badProperty', 'Invalid NumBars property');
end
iStr.NumBars = round(iStr.NumBars);
if ~isa(iStr.NumGroups, 'double') || numel(iStr.NumGroups) ~= 1 || ...
    isinf(iStr.NumGroups) || isnan(iStr.NumGroups) || iStr.NumGroups < 1
    error('neuroelf:xfigure:badProperty', 'Invalid NumGroups property');
end
iStr.NumGroups = round(iStr.NumGroups);

% get new object
o = xfigure('new');
o.T = -1; % keep as placeholder

% total width per group
if ~isfield(iStr, 'GroupWidth') || ~isa(iStr.GroupWidth, 'double') || numel(iStr.GroupWidth) ~= 1 || ...
    isinf(iStr.GroupWidth) || isnan(iStr.GroupWidth) || iStr.GroupWidth <= 0 || iStr.GroupWidth >= 1
    iStr.GroupWidth = 0.8;
end
if ~isfield(iStr, 'ExtraBins') || ~isa(iStr.ExtraBins, 'double') || numel(iStr.ExtraBins) ~= 1 || ...
    isinf(iStr.ExtraBins) || isnan(iStr.ExtraBins) || iStr.ExtraBins < 0
    iStr.ExtraBins = 0;
end

% compute bar width
bw = iStr.GroupWidth / (iStr.NumBars + iStr.ExtraBins);
iStr.BarWidth = bw;

% spacing
if ~isfield(iStr, 'BarSpacing') || ~isa(iStr.BarSpacing, 'double') || numel(iStr.BarSpacing) ~= 1 || ...
    isinf(iStr.BarSpacing) || isnan(iStr.BarSpacing) || iStr.BarSpacing < 0
    iStr.BarSpacing = 0.1;
else
    iStr.BarSpacing = min(0.5, iStr.BarSpacing);
end
bs = bw * iStr.BarSpacing;
iStr.BarSpacingRel = bs;

% compute X positions
xp = 0.5 * (2 + bs + bw * iStr.ExtraBins - iStr.GroupWidth);
xp = repmat(xp + bw .* (0:(iStr.NumBars-1)), 2, 1);
xp = [xp; (bw - bs) + xp];
xp = xp(:);
if iStr.NumGroups > 1
    xp = xp(:, ones(1, iStr.NumGroups)) + ones(numel(xp), 1) * (0:(iStr.NumGroups-1));
    xp = xp(:);
end
iStr.VerticesX = xp;
iStr.BarPos = squeeze(mean(reshape(xp, [4, iStr.NumBars, iStr.NumGroups]), 1))';

% no Y positions given
if ~isfield(iStr, 'YData') || ~isa(iStr.YData, 'double') || numel(iStr.YData) ~= (iStr.NumBars * iStr.NumGroups)
    iStr.YData = zeros(iStr.NumGroups, iStr.NumBars);
end
yp = iStr.YData';
yp = yp(:)';
yp = [zeros(1, numel(yp)); yp; yp; zeros(1, numel(yp))];
yp = yp(:);
iStr.VerticesY = yp;

% use tag from iStr?
if isfield(iStr, 'Tag') && ischar(iStr.Tag) && ~isempty(iStr.Tag) && numel(iStr.Tag) < 28 && isrealvarname(iStr.Tag(:)')
    utag = iStr.Tag(:)';
else
    iStr.Tag = sprintf('BAR_%010.0f', floor(1e10 * rand(1, 1)));
    utag = ['xfigure_' iStr.Tag];
end

% foreground color
if ~isfield(iStr, 'EdgeColor') || ~isa(iStr.EdgeColor, 'double') || numel(iStr.EdgeColor) ~= 3 || ...
    any(isinf(iStr.EdgeColor) | isnan(iStr.EdgeColor) | iStr.EdgeColor < 0 | iStr.EdgeColor > 1)
    iStr.EdgeColor = [0, 0, 0];
else
    iStr.EdgeColor = iStr.EdgeColor(:)';
end

% create plot
o.H = patch('Faces', reshape(1:numel(xp), 4, iStr.NumBars * iStr.NumGroups)', ...
    'Vertices', [xp, yp], 'Parent', xo.H);
dfcn = sprintf('delete(xfigure(hxdouble(''%s'')));', hxdouble(double(o.H)));
set(o.H, 'DeleteFcn', dfcn, 'EdgeColor', iStr.EdgeColor, ...
    'HandleVisibility', xfigsngl.hvisible, 'Tag', utag);

% re-set axes X limit
if strcmpi(get(xo.H, 'XLimMode'), 'auto')
    set(xo.H, 'XLim', [0.5, 0.5 + iStr.NumGroups]);
end

% output properties
iStr.Type = 'barplot';
o.X.loadprops = iStr;
xfigmlup(end + 1) = o.H;
xfigures(end + 1) = o;

% add to global tag lookup struct
xfigsngl.tags.(['BAR_' iStr.Tag]) = o;
