% FUNCTION ne_movestatmap: move selected stats maps
function varargout = ne_movestatmap(varargin)

% Version:  v0.9d
% Build:    14052912
% Date:     May-29 2014, 12:59 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

% allow special treatment of RFX-GLMs
if nargin > 3 && ...
    ischar(varargin{4}) && ...
   ~isempty(varargin{4}) && ...
    isfield(cc, varargin{4}(:)')
    stvar = cc.(varargin{4}(:)');
else
    stvar = cc.StatsVar;
end
if numel(stvar) ~= 1 || ...
   ~isxff(stvar, true) || ...
    nargin < 3 || ...
   ~isa(varargin{3}, 'double') || ...
    numel(varargin{3}) ~= 1
       return;
end
sttyp = stvar.Filetype;
if strcmpi(sttyp, 'glm') && ...
    stvar.ProjectTypeRFX > 0 && ...
    varargin{3} == 0

    % get selected subject IDs
    if stvar.ProjectType > 1
        stvix = cc.SurfStatsVarIdx;
        sids = ch.SurfStatsVarMaps.String(stvix);
    else
        stvix = cc.StatsVarIdx;
        sids = ch.StatsVarMaps.String(stvix);
    end
    if isempty(stvix)
        return;
    end
    if ~iscell(sids)
        sids = cellstr(sids);
    end
    sids = strrep(unique(regexprep(sids, '\:.*$', '')), 'Subject ', '');

    % then ask for removal of subjects!
    if numel(sids) == 1
        a = questdlg(sprintf('Do you want to remove subject ''%s'' from GLM %s?', ...
            sids{1}, stvar.FilenameOnDisk), 'NeuroElf GUI - request', 'Yes', 'No', 'No');
    else
        a = questdlg(sprintf('Do you want to remove %d subjects from GLM %s?', ...
            numel(sids), stvar.FilenameOnDisk), 'NeuroElf GUI - request', 'Yes', 'No', 'No');
    end

    % remove?
    resaved = true;
    if strcmpi(a, 'yes')
        stvar.RemoveSubject(sids);

        % offer save as
        resaved = false;
        oldname = stvar.FilenameOnDisk;
        stvar.SaveAs;
        if ~strcmp(oldname, stvar.FilenameOnDisk)
            resaved = true;

        % if not resaved
        else

            % clear from lastglm to prevent saving of RunTimeVars
            if isstruct(ne_gcfg.h.CM) && ...
                isfield(ne_gcfg.h.CM, 'CMFig') && ...
                isxfigure(ne_gcfg.h.CM.CMFig, true) && ...
                isfield(ne_gcfg.h.CM.CMFig.UserData, 'lastglm') && ...
                stvar == ne_gcfg.h.CM.CMFig.UserData.lastglm
                ne_gcfg.h.CM.CMFig.UserData.lastglm = [];
            end
            if isstruct(ne_gcfg.h.RM) && ...
                isfield(ne_gcfg.h.RM, 'RMFig') && ...
                isxfigure(ne_gcfg.h.RM.RMFig, true) && ...
                isfield(ne_gcfg.h.RM.RMFig.UserData, 'lastglm') && ...
                stvar == ne_gcfg.h.RM.RMFig.UserData.lastglm
                ne_gcfg.h.RM.RMFig.UserData.lastglm = [];
            end
        end

        % re-initialize list
        if stvar.ProjectType > 1
            ne_setcsrfstats;
        else
            ne_setcstats;
        end
    end

    % not saved?
    if ~resaved
        uiwait(warndlg(['Please re-save under the GLM under a different name ' ...
            'to keep the original covariates and groups!'], 'NeuroElf - info', 'modal'));
    end

    % return!
    return;
end

% rest of code only valid if currently selected file is either VMP or HEAD
if ~any(strcmpi(sttyp, {'head', 'cmp', 'smp', 'vmp'}))
    return;
end

% also only valid if maps are selected
if strcmpi(sttyp, 'smp')
    surfs = true;
    stvix = cc.SurfStatsVarIdx;
else
    surfs = false;
    stvix = cc.StatsVarIdx;
end
if isempty(stvix)
    return;
end

% get requested objects
switch (lower(sttyp))
    case {'head'}
        field = 'Brick';
    case {'cmp', 'smp', 'vmp'}
        field = 'Map';
end
try
    maps = stvar.(field);
catch ne_eo;
    uiwait(warndlg(['An error occurred: ' ne_eo.message], 'NeuroElf GUI - warning', 'modal'));
    return;
end

% remove maps
if varargin{3} == 0

    % remove from container
    maps(stvix) = [];

    % and re-set
    stvar.(field) = maps;

    % then update list
    if surfs

        % clear index
        ch.SurfStatsVarMaps.Value = [];
        ch.SurfStatsVarMaps.String = stvar.MapNames(ne_gcfg.c.extmapnames);

        % update map selection
        ne_setcsrfstatmap;
    else

        % clear index
        cc.StatsVarIdx = [];
        ch.StatsVarMaps.Value = [];
        ch.StatsVarMaps.String = stvar.MapNames(ne_gcfg.c.extmapnames);

        % update screen
        ne_setcstatmap;
    end

    % and return early
    return;
end

% compute new indices
nmaps = numel(maps);
idxto = varargin{3};
sttix = stvix + idxto;
while any(sttix < 1)
    sttix = sttix + 1;
end
while any(sttix > nmaps)
    sttix = sttix - 1;
end

% no movement
if all(sttix == stvix)
    return;
end

% where to move to
ridx = setdiff(1:nmaps, stvix(:)');
oidx = setdiff(1:nmaps, sttix(:)');
tmaps = maps;
tmaps(oidx) = maps(ridx);
tmaps(sttix) = maps(stvix);

% re-set field
stvar.(field) = tmaps;

% and update List (names)
if surfs
    ch.SurfStatsVarMaps.String = stvar.MapNames(ne_gcfg.c.extmapnames);
    ch.SurfStatsVarMaps.Value = sttix;
    ne_setcsrfstatmap;
else
    ch.StatsVarMaps.String = stvar.MapNames(ne_gcfg.c.extmapnames);
    ch.StatsVarMaps.Value = sttix;
    ne_gcfg.fcfg.StatsVarIdx = sttix;
end
