% FUNCTION ne_vmp_combine: combine current with other VMPs
function varargout = ne_vmp_combine(varargin)

% Version:  v1.1
% Build:    16033111
% Date:     Mar-31 2016, 11:36 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, 2016, Jochen Weber
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
if numel(stvar) ~= 1 || ~isxff(stvar, 'vmp')
    return;
end

% get handle and layout of vmp
vh = objectid(stvar);
vl = stvar.Layout;

% get list of all VMPs that are available for a match
svar = ne_gcfg.h.StatsVar.UserData;
svari = false(size(svar, 1), 1);
for vc = 1:numel(svari)
    if strcmpi(svar{vc, 2}, 'vmp') && ~isequal(objectid(svar{vc, 4}), vh) && ...
        isequal(svar{vc, 4}.Layout, vl)
        svari(vc) = true;
    end
end

% no match
if ~any(svari)
    return;
end

% create list of names and objects
svarn = svar(svari, 1);
svaro = svar(svari, 4);

% make sure all objects have a good name
for vc = 1:numel(svarn)
    if isempty(svaro{vc}.FilenameOnDisk)
        svarn{vc} = sprintf('<unsaved.vmp> (VMP %d)', vc);
    else
        svarn{vc} = svaro{vc}.FilenameOnDisk;
    end
end

% show a selector
[csel, csok] = listdlg( ...
    'PromptString', 'Please select VMPs to combine with...', ...
    'ListString',   svarn, ...
    'InitialValue', [], ...
    'ListSize',     [min(600, max(300, 10 .* size(char(svarn), 2))), 200]);
if ~isequal(csok, 1) || isempty(csel)
    return;
end

% get maps of objects
svaro = svaro(csel);
for vc = 1:numel(svaro)
    svaro{vc} = svaro{vc}.Map(:);
end

% add to current list of maps
stvar.Map = catstruct(stvar.Map(:), svaro{:})';
stvar.NrOfMaps = numel(stvar.Map);

% update
ne_openfile(0, 0, stvar);
