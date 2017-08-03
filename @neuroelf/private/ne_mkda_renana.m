% PUBLIC FUNCTION ne_mkda_renana: rename an analysis in the list
function varargout = ne_mkda_renana(varargin)

% Version:  v0.9c
% Build:    13012611
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

% ask for new name
newname = inputdlg({'Rename analysis into:     '}, 'NeuroElf - user request', ...
    1, anas(anaidx, 1));
if ~iscell(newname) || ...
    numel(newname) ~= 1 || ...
   ~ischar(newname{1}) || ...
    isempty(ddeblank(newname{1}))
    return;
end

% update name in RunTimeVars and String
plp.RunTimeVars.MKDAAnalyses{anaidx, 1} = newname{1};
ch.Analyses.String = plp.RunTimeVars.MKDAAnalyses(:, 1);

% save configuration
if ~isempty(plp.FilenameOnDisk)
    try
        plp.SaveRunTimeVars;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
