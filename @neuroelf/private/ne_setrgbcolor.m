% FUNCTION ne_setrgbcolor: get one of the RGB colors in a VMP map
function varargout = ne_setrgbcolor(varargin)

% Version:  v1.1
% Build:    16031523
% Date:     Mar-15 2016, 11:55 PM EST
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

% only valid with one of the four button settings
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
   ~any(strcmp(varargin{3}, {'+', '++', '-', '--'}))
    return;
end

% get current configuration
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% get correct stvar/stvix
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

% check stvar/stvix
if numel(stvar) ~= 1 || ...
   ~isxff(stvar, typt) || ...
    numel(stvix) ~= 1
    return;
end

% field and button depends on source (string) of call
switch (varargin{3})

    % lower positive threshold
    case {'+'}
        field = 'RGBLowerThreshPos';
        button = 'RGBLPos';

    % upper positive threshold
    case {'++'}
        field = 'RGBUpperThreshPos';
        button = 'RGBUPos';

    % lower negative threshold
    case {'-'}
        field = 'RGBLowerThreshNeg';
        button = 'RGBLNeg';

    % upper negative threshold
    case {'--'}
        field = 'RGBUpperThreshNeg';
        button = 'RGBUNeg';
end

% allow for errors
try

    % call colorpicker with current color
    bcolor = min(255, max(0, colorpicker(stvar.Map(stvix).(field), field)));

    % set back to field
    stvar.Map(stvix).(field) = bcolor(:)';

    % and also set button background color
    ch.(button).BackgroundColor = (1 / 255) .* bcolor(:)';

    % if the map uses RGBcolor scheme
    if stvar.Map(stvix).UseRGBColor > 0

        % delete current colors
        stvar.Map(stvix).OverlayColors = [];

        % and use @VMP::SetColors to re-create LUT from RGB colors
        stvar.SetColors(stvix, 'xauto');
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% update
if ~any(strcmp(typt, 'smp'))
    ne_setslicepos;
else
    ne_setcsrfstatmap;
end
