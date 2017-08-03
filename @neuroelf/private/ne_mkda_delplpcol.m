% PUBLIC FUNCTION ne_mkda_delplpcol: delete PLP column
function varargout = ne_mkda_delplpcol(varargin)

% Version:  v0.9c
% Build:    11120710
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

% get column names
clnames = plp.ColumnNames(:);

% ask for which to remove (other than essentials!)
clnames(multimatch(lower(clnames), {'study'; 'x'; 'y'; 'z'}) > 0) = [];
[todel, lok] = listdlg( ...
    'ListString',     clnames, ...
    'SelectionMode',  'multiple', ...
    'ListSize',       [10 * size(char(clnames), 2), min(400, max(200, 10 * numel(clnames)))], ...
    'InitialValue',   [], ...
    'Name',           'NeuroElf - user input', ...
    'PromptString',   'Please select the columns you wish to delete:');
if ~isequal(lok, 1) || ...
    isempty(todel)
    return;
end

% remove those columns
clnames = clnames(todel);
remcols = find(multimatch(plp.ColumnNames(:), clnames) > 0);
plp.ColumnNames(remcols) = [];
plp.Points(:, remcols) = [];

% check that no column remains selected
for ac = 1:size(rtv.MKDAAnalyses, 1)
    ana = rtv.MKDAAnalyses{ac, 2};
    if any(strcmpi(ana.ContColumn, clnames))
        rtv.MKDAAnalyses{ac, 2}.ContColumn = '';
    end
    if any(strcmpi(ana.StudyColumn, clnames))
        rtv.MKDAAnalyses{ac, 2}.StudyColumn = '';
    end
    cps = lower(ana.CndParts);
    clnames = lower(clnames);
    ff = false;
    for cc = 1:numel(cps)
        for nc = 1:numel(clnames)
            if ~isempty(strfind(cps{cc}, ['$' clnames{nc}]))
                ff = true;
                break;
            end
        end
        if ff
            break;
        end
    end
    if ff
        rtv.MKDAAnalyses{ac, 2}.CndParts = {};
    end
end
plp.RunTimeVars = rtv;

% update PLP
ne_mkda_setplp;

% store saved flag
plp.RunTimeVars.Saved = false;
