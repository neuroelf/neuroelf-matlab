% FUNCTION ne_showsatpage: show a page on a satellite window
function varargout = ne_showsatpage(varargin)

% Version:  v1.1
% Build:    16052921
% Date:     May-29 2016, 9:48 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% global variable
global ne_gcfg;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% check input arguments (either page number or '+1', '-1')
if nargin < 4 || ~ischar(varargin{3}) || ~isfield(ne_gcfg.cc, varargin{3}) || ...
   ((~isa(varargin{4}, 'double') || numel(varargin{4}) ~= 1) && ...
    (~ischar(varargin{4}) || ~any(strcmp(varargin{4}, {'+1', '-1'}))))
    return;
end
tsat = varargin{3}(:)';

% get config and handles
ch = ne_gcfg.cc.(tsat);
cc = ch.Config;
cini = ne_gcfg.c.ini;

% what about crosshair display
if cc.chair
    sl = 'on';
else
    sl = 'off';
end
cbv = 'off';
zbv = 'off';

% show requested page
ch.Satellite.ShowPage(varargin{4}, 'norefresh');

% then get current page (allowing for string call)
ne_gcfg.cc.(tsat).Config.page = ch.Satellite.ShowPage('cur');
cc = ne_gcfg.cc.(tsat).Config;

% update window
switch (cc.sattype)
    case {'slice'}
        ne_setsatslicepos(0, 0, tsat);
    case {'surf'}
        ne_setsurfpos(0, 0, tsat);
end

% for first page
if cc.page == 1

    % take show flag for 3-slice crosshairs
    sl1 = sl;

    % and set zoomed crosshairs invisible
    sl2 = 'off';

    % and surface view axes and content
    sra = 'off';
    if cini.Statistics.ShowThreshBars && cini.Statistics.ShowThreshText
        cbv = 'on';
    end

% for second page
elseif cc.page == 2

    % vice versa
    sl1 = 'off';
    sl2 = sl;
    sra = 'off';
    if cini.Statistics.ShowThreshBars && cini.Statistics.ShowThreshText
        zbv = 'on';
    end

% for third page (surface view)
elseif cc.page == 3

    % both slicing views are off
    sl1 = 'off';
    sl2 = 'off';

    % surface axes and children on
    sra = 'on';

% for any other page
else

    % everything is off
    sl1 = 'off';
    sl2 = 'off';
    sra = 'off';
end

% make the settings
set([ch.CorLineX, ch.CorLineY, ch.SagLineX, ch.SagLineY, ...
     ch.TraLineX, ch.TraLineY], 'Visible', sl1);
set([ch.ZoomLineX, ch.ZoomLineY], 'Visible', sl2);
set([get(ch.CorLineX, 'Parent'), get(ch.SagLineX, 'Parent'), ...
     get(ch.TraLineX, 'Parent'), get(ch.ZoomLineX, 'Parent')], 'Visible', 'off');
if strcmpi(cc.sattype, 'surf')
    set(ch.Surface, 'Visible', sra);
    src = get(ch.Surface, 'Children');
    if ~isempty(src)
        set(src, 'Visible', sra);
    end
end
set(ch.CorStatsText, 'Visible', cbv);
set(ch.ZoomStatsText, 'Visible', zbv);
