% FUNCTION ne_tcbtmove: update window with pressed button for TCPlot
function ne_tcbtmove(varargin)

% Version:  v1.1
% Build:    16052621
% Date:     May-26 2016, 9:22 PM EST
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

% invalid call
if nargin < 3 || ~ischar(varargin{3}) || ~isfield(ne_gcfg.cc, varargin{3}(:)')
    return;
end
iSat = varargin{3}(:)';
cch = ne_gcfg.cc.(iSat);

% only once
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

% hide ToolTip (menu)
ne_gcfg.h.ToolTip.Visible = 'off';

% tcplot figure
if ne_gcfg.c.btdown == cch.SatelliteMLH

    % do nothing if same position
    last = get(cch.SatelliteMLH, 'CurrentPoint');
    if isequal(last, cch.Config.mpos.last)
        ne_gcfg.c.incb = false;
        return;
    end
    ne_gcfg.cc.(iSat).Config.mpos.last = last;
    fPos = cch.Config.tcpos;
    oPos = cch.Config.mpos.down;
    nPos = last;
    cp = nPos - fPos(1, 1:2);
    tcvar = cch.Config.tcvar;
    tcpsz = cch.TCPlot.Position(3:4);

    if all(oPos >= (fPos(1, 1:2) - 1)) && all(oPos <= (fPos(1, 3:4) + 1))
        if isxff(tcvar, {'fmr', 'hdr', 'head', 'vtc'})
            tcrtv = tcvar.RunTimeVars;

            % then get the number of volumes and adapt the position
            if ~isfield(tcrtv, 'AvgVTC') || ~islogical(tcrtv.AvgVTC) || ~tcrtv.AvgVTC
                nvol = tcvar.NrOfVolumes;
                tsvalue = min(nvol, max(1, round((nvol + 2) * cp(1) / tcpsz(1)) - 1));
                ne_gcfg.h.Coord.TempSlider.Value = tsvalue;
            else
                nvol = tcrtv.NrOfVolumesPerTC;
                tsvalue = min(nvol, max(1, ((nvol + 2) * cp(1) / (1.06 * tcpsz(1))) + 1));
                tcvar.RunTimeVars.SubMapVol = tsvalue;
            end
        end
        ne_setslicepos;
    end
end

% lift restriction
ne_gcfg.c.incb = false;

% check if it was up inbetween
if ne_gcfg.c.btdoup
    ne_gcfg.c.btdown = [];
    ne_gcfg.c.btdoup = false;
    ne_gcfg.cc.(iSat).Config.mpos.ddat = {};
end
