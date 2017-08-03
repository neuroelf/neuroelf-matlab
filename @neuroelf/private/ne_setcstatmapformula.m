% FUNCTION ne_setcstatmapformula: compute formula and create new map
function varargout = ne_setcstatmapformula(varargin)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-10 2011, 4:52 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% make sure the current statsvar is a VMP
stvar = cc.StatsVar;
if numel(stvar) ~= 1 || ...
   ~isxff(stvar, 'vmp') || ...
    isempty(stvar.Map)
    return;
end

% get current selection
stvix = cc.StatsVarIdx;
formdesc = '#i := i-th map';
if ~isempty(stvix)
    formdesc = [formdesc ', $i := i-th selected map'];
end

% get formula and other values
formset = {'  #1', sprintf('  %d', numel(stvar.Map) + 1), ...
    '  ', '  ', '  ', '  ', '  n'};
try
    formset = ddeblank(inputdlg({ ...
        ['Formula: ' formdesc], ...
        'Target map number:', ...
        'Target map name (leave blank to use formula):', ...
        'Map type (either of t, r, CC, or F; leave blank to copy from 1st):', ...
        'DF1 (leave blank to copy from first):', ...
        'DF2 (leave blank to copy from first):', ...
        'Use p-values (conversion required) instead of stats, (y)es/(n)o:'}, ...
        'NeuroElf GUI - Map settings', 1, formset));
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end

% check formula fields (as much as possible)
if ~iscell(formset) || ...
    numel(formset) ~= 7
    return;
end
maptypes = {'t', 'r', 'CC', 'F', '', '', '', '', '', '', ...
    '%', 'z_ica', 'TH', '', 'beta', 'prob', '', '', '', ...
    'MD', 'FA'};
if isempty(formset{1}) || ...
   ~any(formset{1} == '$' | formset{1} == '#') || ...
    isempty(formset{2}) || ...
    any(formset{2} < '0' | formset{2} > '9') || ...
    numel(formset{4}) > 1 || ...
   (~isempty(formset{4}) && ...
    ~any(strcmpi(formset{4}, maptypes))) || ...
    isempty(formset{7}) || ...
    ~any(lower(formset{7}(1)) == 'ny')
    return;
end
formset{2} = str2double(formset{2});
if isempty(formset{3})
    formset{3} = formset{1};
end
if formset{2} < 1 || ...
    formset{2} > (numel(stvar.Map) + 1)
    return;
end

% perform command
opts = struct( ...
    'mapsel',  stvix(:)', ...
    'name',    formset{3}, ...
    'pvalues', (lower(formset{7}(1)) == 'y'), ...
    'target',  formset{2});
stvar.ComputeFormula(formset{1}, opts);

% echo
if ne_gcfg.c.echo
    ne_echo('vmp', 'ComputeFormula', opts);
end

% adapt type and DF if requested
if ~isempty(formset{4})
    stvar.Map(formset{2}).Type = findfirst(strcmpi(maptypes, formset{4}));
end
if ~isempty(formset{5}) && ...
   ~isempty(regexpi(formset{5}, '^\d+$')) && ...
    formset{5}(1) > '0'
    stvar.Map(formset{2}).DF1 = str2double(formset{5});
end
if ~isempty(formset{6}) && ...
   ~isempty(regexpi(formset{6}, '^\d+$')) && ...
    formset{6}(1) > '0'
    stvar.Map(formset{2}).DF2 = str2double(formset{6});
end

% bring up new map
stvar.SetColors(formset{2}, 'xauto');
ne_openfile(0, 0, stvar);
ne_gcfg.fcfg.StatsVarIdx = formset{2};
ne_gcfg.h.StatsVarMaps.Value = formset{2};
ne_setcstatmap;
