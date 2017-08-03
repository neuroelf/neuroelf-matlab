function varargout = setgroupvisible(xo, iStr, vmode, varargin)
%XFIGURE::SETGROUPVISIBLE  Sets a group of elements (VGroup) visible.
%   SETGROUPVISIBLE(FIG, GROUP) sets the group GROUP in figure FIG
%   visible.
%
%   SETGROUPVISIBLE(FIG, GROUP, VMODE) sets the visibiltiy of the group
%   GROUP in figure FIG to VMODE (which can also be 1x1 numeric).

% Version:  v1.1
% Build:    16041811
% Date:     Apr-18 2016, 11:51 AM EST
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

% global NeuroElf methods
global ne_methods;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% only valid for figures
if numel(xo) ~= 1 || xo.T ~= 1
    error('neuroelf:xfigure:invalidObjectType', ...
        'SetGroupVisible is only valid for figures.');
end

% check group spec
if ~ischar(iStr) || isempty(iStr(:))
    error('neuroelf:xfigure:invalidObjectType', ...
        'SetGroupVisible requires a valid group name.');
end
dogroup = deblank(iStr(:)');
vismode = 'on';

% determine target vis mode
if nargin > 2 && ((isnumeric(vmode) && ~vmode(1)) || ...
   (ischar(vmode) && strcmpi(vmode, 'off')))
    vismode = 'off';
end

% all groups ?
if strcmpi(dogroup, 'all_groups')
    gnames = fieldnames(xo.X.figprops.vgroups);
    for gc = 1:length(gnames)
        setgroupvisible(xo, gnames{gc}, vismode);
    end
    return;
elseif any(dogroup == ';' | dogroup == ',' | dogroup == ' ')
    gnames = ne_methods.splittocellc(dogroup, ';, ', true, true);
    for gc = 1:length(gnames)
        if isvarname(gnames{gc})
            setgroupvisible(xo, gnames{gc}, vismode);
        end
    end
    return;
elseif ~isvarname(dogroup)
    error('neuroelf:xfigure:invalidObjectType', ...
        'SetGroupVisible requires a valid group name.');
end

% single group
if isfield(xo.X.figprops.vgroups, dogroup)
    hdls = xo.X.figprops.vgroups.(dogroup);
    for hc = 1:numel(hdls)
        setprop(hdls(hc), 'Visible', vismode);
    end
end

% refresh graphics
if (nargin > 2 && ischar(vmode) && strcmpi(vmode(:)', 'refresh')) || ...
   (nargin > 3 && ischar(varargin{end}) && strcmpi(varargin{end}, 'refresh'))
    redrawfig(xo.H);
end
