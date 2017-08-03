function o = adduimenu(xo, iStr)
%XFIGURE::ADDUIMENU  Add a uimenu to a figure, uimenu, or uicontextmenu.
%   UIM = ADDUIMENU(XFIG, UIMSTRUCT) adds the uimenu with settings in
%   UIMSTRUCT to the xfigure object in XFIG, which must be of type figure,
%   uimenu, or uicontextmenu.

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
global ne_methods xfigmlup xfigsngl xfigures;

% only valid for figures or UIMenus
if numel(xo) ~= 1 || ~any([1, 3, 4] == xo.T)
    error('neuroelf:xfigure:invalidObjectType', ...
        'Bad object type (%s) for AddUIMenu(...).', get(xo.H, 'Type'));
end

% test iStr
if ~isstruct(iStr) || numel(iStr) ~= 1
    error('neuroelf:xfigure:badPropertyStruct', ...
        'Uimenu properties must be of type struct.');
end

% perform some checks on struct
iStr = ne_methods.checkstruct(iStr, xfigsngl.optuim);
if isempty(iStr.Caption)
    error('neuroelf:xfigure:badPropertyStruct', ...
        'Bad uimenu property struct supplied.');
end

% get new object
o = xfigure('new');
o.T = 3;
o.X = makeostruct(xfigsngl.objtypes.uimenu);

% get parent object
if xo.T == 1
    po = xo;
else
    po = xfigure(ancestor(xo.H, 'figure'));
end

% use tag from iStr?
if ~isempty(iStr.Tag) && numel(iStr.Tag) < 28 && isrealvarname(iStr.Tag(:)')
    utag = iStr.Tag(:)';
else
    iStr.Tag = sprintf('UIM_%010.0f', floor(1e10 * rand(1, 1)));
    utag = ['xfigure_' iStr.Tag];
end

% accelerator
if ~isempty(iStr.Accelerator)
    iStr.Accelerator = iStr.Accelerator(1);
end

% text color
if numel(iStr.Color) == 1
    iStr.Color(1:3) = iStr.Color;
elseif numel(iStr.Color) ~= 3
    iStr.Color = [0, 0, 0];
end
iCFG = max(0, min(1, iStr.Color));
iPos = iStr.Position;

% any Enable-groups ?
iGroups = ne_methods.splittocellc(iStr.EGroups, ',; ', true, true);
for iGroupC = numel(iGroups):-1:1
    iGroup = deblank(iGroups{iGroupC});
    if isrealvarname(iGroup)
        if isfield(po.X.figprops.egroups, iGroup)
            po.X.figprops.egroups.(iGroup)(end + 1) = o;
        else
            po.X.figprops.egroups.(iGroup) = o;
        end
    else
        iGroups(iGroupC) = [];
    end
end
iStr.EGroups = iGroups;

% any Visible-groups (or pages)
iGroups = ne_methods.splittocellc(iStr.VGroups, ',; ', true, true);
for iGroupC = numel(iGroups):-1:1
    iGroup = deblank(iGroups{iGroupC});
    if isrealvarname(iGroup)
        if isfield(po.X.figprops.vgroups, iGroup)
            po.X.figprops.vgroups.(iGroup)(end + 1) = o;
        else
            po.X.figprops.vgroups.(iGroup) = o;
        end
    else
        iGroups(iGroupC) = [];
    end
end
iStr.VGroups = iGroups;

% prepare output structure
oStr = struct( ...
    'Parent',          xo.H, ...
    'Label',           iStr.Caption, ...
    'Accelerator',     iStr.Accelerator, ...
    'Checked',         iStr.Checked, ...
    'DeleteFcn',       '', ...
    'Enable',          iStr.Enabled, ...
    'ForegroundColor', iCFG, ...
    'Separator',       iStr.Separator, ...
    'Tag',             utag, ...
    'UserData',        iStr.UserData, ...
    'Visible',         iStr.Visible);
oStr.Callback = iStr.Callback;
if ~isempty(iPos) && ~isnan(iPos(1)) && ~isinf(iPos(1))
    oStr.Position = fix(iPos(1));
end

% create object and fill/update fields as necessary
hOut = uimenu(oStr);
o.H = hOut;

% complete object representation
o.X.callbacks = {iStr.Callback, '', '', iStr.CallbackDelete};
o.X.loadprops = iStr;
xfigmlup(end + 1) = hOut;
xfigures(end + 1) = o;

% finally, add to global tag lookup struct
xfigsngl.tags.(['UIM_' iStr.Tag]) = o;

% set handle visibility
set(hOut, 'HandleVisibility', xfigsngl.hvisible);

% and return appropriate object
set(hOut, 'DeleteFcn', sprintf('delete(xfigure(hxdouble(''%s'')));', ...
    hxdouble(double(hOut))));
