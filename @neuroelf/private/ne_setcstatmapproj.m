% FUNCTION ne_setcstatmapproj: project current selection to other subjects
function varargout = ne_setcstatmapproj(varargin)

% Version:  v1.1
% Build:    16031016
% Date:     Mar-10 2016, 4:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% check stats var
if nargin < 3 || ~ischar(varargin{3}) || isempty(varargin{3}) || ~strcmpi(varargin{3}(:)', 'surfstatsvar')
    stv = true;
    stvar = cc.StatsVar;
    csel = ch.StatsVarMaps.Value;
    ssel = ch.StatsVarMaps.String;
else
    stv = false;
    stvar = cc.SurfStatsVar;
    csel = ch.SurfStatsVarMaps.Value;
    ssel = ch.SurfStatsVarMaps.String;
end
if numel(stvar) ~= 1 || ~isxff(stvar, 'glm')
    return;
end

% get current selection
if ~iscell(ssel)
    ssel = cellstr(ssel);
end
ssel = unique(regexprep(ssel(csel), 'Subject .*\:\s', ': '));

% get available map names
avl = stvar.MapNames;

% create logical array for index
avb = false(size(avl));

% set in index
for sc = 1:numel(ssel)
    avb = avb | ~cellfun('isempty', regexpi(avl, [ssel{sc} '$']));
end

% set in value and update
if stv
    ch.StatsVarMaps.Value = find(avb);
    ne_setcstatmap;
else
    ch.SurfStatsVarMaps.Value = find(avb);
    ne_setcsrfstatmap;
end
