function varargout = ne_clustertable(varargin)
% ne_clustertable  - create or add to cluster table
%
% FORMAT:       ne_clustertable(SRC, EVT [, type [, addflag]])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       type        if set to 'smp' use for SMP, otherwise StatsVar
%       addflag     if set to 'add' or 'new' override config setting
%
% Example:
%
%   VMP_OBJECT.Browse(1);
%   ne_clustertable(0, 0, 'vmp', 'new');
%   VMP_OBJECT.Browse(2);
%   ne_clustertable(0, 0, 'vmp', 'add');

% Version:  v1.1
% Build:    16031522
% Date:     Mar-15 2016, 10:22 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% do not call multiple times
if any(strcmp('clustertable', ne_gcfg.c.blockcb))
    return;
end
ne_gcfg.c.blockcb{end+1} = 'clustertable';

% get config and handles, as well as StatsVar and index
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;
ci = ne_gcfg.c.ini;
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
   ~strcmpi(varargin{3}(:)', 'smp')
    sfvar = {};
    stvar = cc.StatsVar;
    stvix = cc.StatsVarIdx;
    typt = {'cmp', 'hdr', 'head', 'vmp'};
else
    sfvar = {cc.SurfVar};
    stvar = cc.SurfStatsVar;
    stvix = cc.SurfStatsVarIdx;
    typt = {'fsmf', 'smp'};
end
if nargin < 4 || ...
   ~ischar(varargin{4}) || ...
   ~any(strcmpi(varargin{4}(:)', {'add', 'new'}))
    addtolist = ci.Statistics.ClusterTableAdd;
elseif strcmpi(varargin{4}(:)', 'add')
    addtolist = true;
else
    addtolist = false;
end

% for anything but VMPs with a single map selection
if ~isxff(stvar, typt) || ...
    numel(stvix) ~= 1 || ...
   ((any(strcmp(typt, 'fsmf')) || any(strcmp(typt, 'smp'))) && ...
    (~isxff(sfvar{1}, 'srf') || ...
     sfvar{1}.NrOfVertices ~= stvar.NrOfVertices))

    % unblock further calls and return
    ne_gcfg.c.blockcb(strcmp('clustertable', ne_gcfg.c.blockcb)) = [];
    return;
end

% enable hourglass
chg = ch.MainFig.Pointer;
ch.MainFig.Pointer = 'watch';

% build options struct for @VMP::ClusterTable
% .localmax   - guess a useful threshold to split clusters
% .localmin   - minimum size of split sub-clusters (local maxima)
% .localmsz   - print sizes of sub-clusters
% .mni2tal    - switch from UI, convert coordinates for text table
% .tdclient   - run tdclient(NGM) on coordinates
% .pbar       - progress bar handle
ctopt = struct( ...
    'clconn',   cc.clconn, ...
    'icbm2tal', (ch.Stats.ICBM2TAL.Value > 0), ...
    'localmsz', cc.localmaxsz, ...
    'lupcrd',   ci.Statistics.LookupCoord, ...
    'sorting',  cc.clsort, ...
    'tdclient', (ch.Stats.TDClient.Value > 0));
if cc.localmax
    if any(strcmp(stvar.Filetype, {'cmp', 'vmp'}))
        ctopt.localmax = floor(100 ./ (stvar.Resolution .^ 2));
        ctopt.localmin = round(125 ./ (stvar.Resolution .^ 3));
    elseif ~any(strcmp(typt, 'fsmf')) && ~any(strcmp(typt, 'smp'))
        cfr = stvar.CoordinateFrame;
        ctopt.localmax = floor(100 / prod(cfr.Resolution));
        ctopt.localmin = round(125 / prod(cfr.Resolution));
    else
        ctopt.localmax = cc.localmaxsrf;
    end
end

% set progress bar visible
cprog = ne_progress(0, 0, {true, 0, 'Cluster table'});

% echo
if ne_gcfg.c.echo
    ne_echo(lower(stvar.Filetype), 'ClusterTable', sfvar{:}, stvix, [], ctopt);
end
ctopt.pbar = ch.Progress;

% create clustertable (ttab), and get VOI with clusters
[clusters, ttab, clmap, clvoi] = ...
    stvar.ClusterTable(sfvar{:}, stvix, [], ctopt);

% set progress bar invisible again
ne_progress(0, 0, cprog);

% rest depends on type
if ~any(strcmp(typt, 'fsmf')) && ~any(strcmp(typt, 'smp'))

    % if currently a VOI object is loaded
    if isxff(ne_gcfg.voi, true)

        % keep VOIs in list
        oldvois = ne_gcfg.voi.VOI;

        % then clear it
        voih = handles(ne_gcfg.voi);
        if isfield(voih, 'RGBImage') && ...
            numel(voih.RGBImage) == 1 && ...
            isxff(voih.RGBImage, 'hdr')
            voih.RGBImage.ClearObject;
        end
        ne_gcfg.voi.ClearObject;
    else
        oldvois = [];
    end

    % and keep the new one
    ne_gcfg.voi = clvoi;

    % add old VOIs?
    if addtolist && ...
       ~isempty(oldvois)
        clvoi.VOI = joinstructs(oldvois(:)', clvoi.VOI(:)');
        fcc = numel(oldvois) + 1;
    else
        fcc = 1;
    end
    clvoi.NrOfVOIs = numel(clvoi.VOI);

    % iterate over VOIs
    vnames = clvoi.VOINames;
    for cc = fcc:numel(vnames)

        % add "local max." text if appropriate
        if strcmpi(clusters(cc+1-fcc).localmax, 'l')
            lmaxs = ' (local max.)';
        else
            lmaxs = '';
        end

        % and also give number of voxels
        vnames{cc} = sprintf('%s - %d voxel(s)%s', vnames{cc}, ...
            clvoi.VOI(cc).NrOfVoxels, lmaxs);
    end

    % set cluster selection
    ch.Clusters.ListboxTop = 1;
    ch.Clusters.Value = [];
    ch.Clusters.String = vnames;

    % and if clusters were found
    if ~isempty(vnames)

        % also enable table control
        ch.Clusters.Enable = 'on';

    % or if not
    else

        % then disable
        ch.Clusters.Enable = 'off';
    end

% for SMP/POIs
else

    % if currently a POI object is loaded
    if isxff(ne_gcfg.poi, true)

        % keep VOIs in list
        oldpois = ne_gcfg.poi.POI;
        oldpoin = ne_gcfg.poi.NrOfMeshVertices;

        % then clear it
        ne_gcfg.poi.ClearObject;
    else
        oldpois = [];
    end

    % and keep the new one
    ne_gcfg.poi = clvoi;

    % add old VOIs?
    if addtolist && ...
       ~isempty(oldpois) && ...
        oldpoin == clvoi.NrOfVertices
        clvoi.POI = joinstructs(oldpois(:)', clvoi.POI(:)');
        fcc = numel(oldpois) + 1;
    else
        fcc = 1;
    end
    clvoi.NrOfPOIs = numel(clvoi.POI);

    % get names of VOI regions (clusters), which contain coordinates
    pnames = clvoi.POINames;

    % iterate over VOIs
    for cc = fcc:numel(pnames)

        % add "local max." text if appropriate
        if strcmpi(clusters(cc).localmax, 'l')
            lmaxs = ' (local max.)';
        else
            lmaxs = '';
        end

        % and also give number of voxels
        pnames{cc} = sprintf('%s - %d vertices%s', pnames{cc}, ...
            clvoi.POI(cc).NrOfVertices, lmaxs);
    end

    % set cluster selection
    ch.ClustersSrf.ListboxTop = 1;
    ch.ClustersSrf.Value = [];
    ch.ClustersSrf.String = pnames;

    % and if clusters were found
    if ~isempty(pnames)

        % also enable table control
        ch.ClustersSrf.Enable = 'on';

    % or if not
    else

        % then disable
        ch.ClustersSrf.Enable = 'off';
    end
end

% set text to edit control
ch.ClusterTable.String = ttab;

% in minimized state also on console
if ~ne_gcfg.fcfg.fullsized
    disp(' ');
    disp(ttab);
end

% un-hour glass
ch.MainFig.Pointer = chg;

% unblock callback
ne_gcfg.c.blockcb(strcmp('clustertable', ne_gcfg.c.blockcb)) = [];

% and update window
if ~any(strcmp(typt, 'fsmf')) && ~any(strcmp(typt, 'smp'))
    ne_setslicepos;
end
