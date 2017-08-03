% FUNCTION ne_vmp_applyfdr: apply FDR to selection of maps
function varargout = ne_vmp_applyfdr(varargin)

% Version:  v0.9c
% Build:    12011212
% Date:     Jan-05 2012, 12:51 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, Jochen Weber
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

% only valid if single VMP
stvar = ne_gcfg.fcfg.StatsVar;
if numel(stvar) ~= 1 || ...
   ~isxff(stvar, 'vmp')
    return;
end

% get maps
maps = stvar.Map;
mapnames = {maps.Name};
mapnames = mapnames(:);
mapsel = true([numel(maps), 1]);
for mc = 1:numel(maps)
    if ~isequal(maps(mc).Type, 1) && ...
       ~isequal(maps(mc).Type, 2) && ...
       ~isequal(maps(mc).Type, 4)
        mapsel(mc) = false;
    else
        mapnames{mc} = sprintf('Map %02d: %s', mc, mapnames{mc});
    end
end
if ~any(mapsel)
    return;
end
mapsel = find(mapsel);

% show a selector
[csel, csok] = listdlg( ...
    'PromptString', 'Please select maps to apply FDR to...', ...
    'ListString',   mapnames(mapsel), ...
    'InitialValue', (1:numel(mapsel))', ...
    'ListSize',     [min(600, max(300, 10 .* size(char(mapnames(mapsel)), 2))), 200]);
if ~isequal(csok, 1) || ...
    isempty(csel)
    return;
end
mapsel = mapsel(csel);

% apply FDR
mf = ne_gcfg.h.MainFig;
mfp = mf.Pointer;
mf.Pointer = 'watch';
drawnow;
try
    stvar.PerformFDR(struct('mapsel', mapsel));
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    mf.Pointer = mfp;
    drawnow;
    return;
end

% update
ne_openfile(0, 0, stvar);
mf.Pointer = mfp;
drawnow;
