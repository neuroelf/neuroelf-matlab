% FUNCTION ne_tcbtdown: react on click event (down) for TCPlot
function ne_tcbtdown(varargin)

% Version:  v1.1
% Build:    16052621
% Date:     May-26 2016, 9:21 PM EST
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

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% invalid call
if nargin < 3 || ~ischar(varargin{3}) || ~isfield(ne_gcfg.cc, varargin{3}(:)')
    return;
end
iSat = varargin{3}(:)';
cch = ne_gcfg.cc.(iSat);

% only if nothing is waiting
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

% take note that the button is pressed!
ne_gcfg.c.btdown = gcbf;
ne_gcfg.c.btdoup = false;
if nargin < 1 || numel(varargin{1}) ~= 1 || ~ishandle(varargin{1})
    downtest = ne_gcfg.c.btdown;
else
    downtest = varargin{1};    
end

% for main window
if downtest == cch.SatelliteMLH

    % record where and what modifiers at the time
    ne_gcfg.cc.(iSat).Config.mpos.down = cch.Satellite.CurrentPoint;
    ne_gcfg.cc.(iSat).Config.mpos.mods = ne_gcfg.cc.(iSat).Config.mods;
    ne_gcfg.cc.(iSat).Config.tcpos = cch.TCPlot.Position;
    tcpsz = ne_gcfg.cc.(iSat).Config.tcpos(3:4);
    ne_gcfg.cc.(iSat).Config.tcpos(3:4) = ne_gcfg.cc.(iSat).Config.tcpos(1:2) + tcpsz;

    % update?
    fPos = ne_gcfg.cc.(iSat).Config.tcpos;
    nPos = ne_gcfg.cc.(iSat).Config.mpos.down;
    cp = nPos - fPos(1, 1:2);
    tcvar = ne_gcfg.cc.(iSat).Config.tcvar;
    if isempty(cc.mods) && all(nPos >= (fPos(1, 1:2) - 1)) && all(nPos <= (fPos(1, 3:4) + 1))
        if isxff(tcvar, {'fmr', 'hdr', 'head', 'vtc'})
            tcrtv = tcvar.RunTimeVars;

            % then get the number of volumes and adapt the position
            if ~isfield(tcrtv, 'AvgVTC') || ~islogical(tcrtv.AvgVTC) || ~tcrtv.AvgVTC
                nvol = tcvar.NrOfVolumes;
                tsvalue = min(nvol, max(1, round((nvol + 2) * cp(1) / tcpsz(1)) - 1));
                ch.Coord.TempSlider.Value = tsvalue;
            else
                nvol = tcrtv.NrOfVolumesPerTC;
                tsvalue = min(nvol, max(1, ((nvol + 2) * cp(1) / (1.06 * tcpsz(1))) + 1));
                tcvar.RunTimeVars.SubMapVol = tsvalue;
            end
        end
        ne_setslicepos;
    else
        ne_gcfg.c.btdown = -1;
    end
end

% unblock
ne_gcfg.c.incb = false;
