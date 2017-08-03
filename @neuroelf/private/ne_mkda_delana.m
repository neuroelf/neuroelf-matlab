% PUBLIC FUNCTION ne_mkda_delana: remove an analysis from the list
function varargout = ne_mkda_delana(varargin)

% Version:  v0.9c
% Build:    11120211
% Date:     Nov-08 2011, 1:08 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
hFig = ne_gcfg.h.MKDA.MKDAFig;
ch = ne_gcfg.h.MKDA.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get content of PLP dropdown
plps = ch.PLPs;
plpud = plps.UserData;
plpid = plps.Value;
try
    plp = plpud{plpid, 3};
    if numel(plp) ~= 1 || ...
       ~isxff(plp, 'plp')
        error( ...
            'neuroelf:GUI:BadPLP', ...
            'Bad PLP object.' ...
        );
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end
rtv = plp.RunTimeVars;
anas = rtv.MKDAAnalyses;
anaidx = ch.Analyses.Value;

% invalid call
if size(anas, 1) == 1 && ...
    strcmp(anas{1, 1}, '<no analyses defined>')
    hFig.SetGroupEnabled('HasAnas', 'off');
    return;
end

% ask for "are you sure?"
vans = questdlg('Are you sure you wish to delete the current analysis?', ...
    'NeuroElf - user request', 'No', 'Yes', 'No');
if ~isequal(vans, 'Yes')
    return;
end

% last analysis
if size(anas, 1) == 1

    % simply replace name and disable group
    plp.RunTimeVars.MKDAAnalyses{1, 1} = '<no analyses defined>';
    hFig.SetGroupEnabled('HasAnas', 'off');

% otherwise
else

    % remove selection
    plp.RunTimeVars.MKDAAnalyses(anaidx, :) = [];

    % update value, string and then set ana
    if anaidx >= size(anas, 1)
        anaidx = anaidx - 1;
    end
    ch.Analyses.Value = anaidx;
    ch.Analyses.String = plp.RunTimeVars.MKDAAnalyses(:, 1);
    ne_mkda_setana;
end

% save configuration
if ~isempty(plp.FilenameOnDisk)
    try
        plp.SaveRunTimeVars;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
