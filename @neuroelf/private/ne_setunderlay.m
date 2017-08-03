% FUNCTION ne_setunderlay: set underlay object
function ne_setunderlay(varargin)

% Version:  v1.1
% Build:    16031616
% Date:     Mar-16 2016, 4:22 PM EST
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

% only valid if valid xff object
svar = ne_gcfg.fcfg.SliceVar;
if numel(svar) ~= 1 || ...
   ~isxff(svar, true)
    return;
end
% unset current underlay
svar.SetHandle('Underlay', []);

% get list of all objects
svars = ne_gcfg.h.SliceVar.UserData;
svari = false(size(svars, 1), 1);
for vc = 1:numel(svari)
    if numel(svars{vc, 4}) == 1 && ...
        isxff(svars{vc, 4}, {'hdr', 'head', 'mgh', 'vmr', 'vtc'}) && svars{vc, 4} ~= svar
        svari(vc) = true;
    end
end

% no match
if ~any(svari)
    ne_setcvar;
    return;
end

% create list of names and objects
svarn = svars(svari, 1);
svaro = svars(svari, 4);

% make sure all objects have a good name
for vc = 1:numel(svarn)
    svarn{vc} = svaro{vc}.FilenameOnDisk(true);
    if isempty(svarn{vc})
        svarn{vc} = sprintf('<unsaved> (object #%d)', vc);
    end
end

% show a selector
[svarpath, svarname] = fileparts(svar.FilenameOnDisk(true));
if isempty(svarname)
    svarname = 'current slicing object';
end
[csel, csok] = listdlg( ...
    'PromptString',  sprintf('Please select underlay object for %s...', svarname), ...
    'ListString',    svarn, ...
    'InitialValue',  [], ...
    'ListSize',      [min(600, max(300, 10 .* size(char(svarn), 2))), 200], ...
    'SelectionMode', 'single');
if ~isequal(csok, 1) || ...
    isempty(csel)
    ne_setcvar;
    return;
end

% update Underlay handle
svar.SetHandle('Underlay', svaro{csel});

% update
ne_setcvar;
