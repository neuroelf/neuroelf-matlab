% FUNCTION ne_setstatthrcolor: toggle between LUT and RGB schemes (for VMPs)
function varargout = ne_setstatthrcolor(varargin)

% Version:  v1.1
% Build:    16031522
% Date:     Mar-15 2016, 10:38 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% get configuration and handles
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% and StatsVar and map number
if nargin < 4 || ...
   ~ischar(varargin{4}) || ...
   ~strcmpi(varargin{4}(:)', 'smp')
    stvar = cc.StatsVar;
    stvix = cc.StatsVarIdx;
    ch = ch.Stats;
    typt = {'cmp', 'glm', 'hdr', 'head', 'vmp', 'vtc'};
else
    stvar = cc.SurfStatsVar;
    stvix = cc.SurfStatsVarIdx;
    ch = ch.SrfStats;
    typt = {'fsmf', 'smp'};
end

% which button
button = varargin{3};
if strcmpi(button, 'lut')
    button = 'UseLUT';
    otherb = 'UseRGB';
else
    button = 'UseRGB';
    otherb = 'UseLUT';
end

% set button states
ch.(button).Value = 1;
ch.(otherb).Value = 0;

% only valid for VMPs with single map selection
if isxff(stvar, typt) && ...
    numel(stvix) == 1

    % determine whether or not to use RGBs
    usergb = double(ch.UseRGB.Value > 0);

    % then re-color the VMP map
    try
        stvar.Map(stvix).UseRGBColor = usergb;
        stvar.Map(stvix).OverlayColors = [];
        stvar.SetColors(stvix, 'xauto');
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end

% update screen
if any(strcmp(typt, 'vmp'))
    ne_setslicepos;
else
    ne_setcsrfstatmap;
end
