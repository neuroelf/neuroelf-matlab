function varargout = ne_addatlascluster(varargin)
% ne_addatlascluster: add (Talairach) atlas cluster to list of VOIs
%
% FORMAT:       [vois, voio] = ne_addatlascluster(SRC, EVT, cluster)
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       cluster     either single cluster name (TDDB label) or list
%
% Output fields:
%
%       vois        VOI structures added to VOI object
%       voio        VOI object (stored also in global ne_gcfg.voi)
%
% Example:
%
%   [vois, voi] = ne_addatlascluster(0, 0, 'amygdala')
%   amygdala_tc = VTC_OBJECT.VOITimeCourse(voi);

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:31 AM EST
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
c = ne_gcfg.c;
atl = c.atlas;
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% no region requested
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~any(strcmpi(makelabel(varargin{3}(:)'), atl.tal.shorts))

    % request region(s)
    [regsel, selok] = listdlg( ...
        'ListString', atl.tal.labels, 'SelectionMode', 'multiple', ...
    'ListSize', [480, 320], 'InitialValue', (1:numel(atl.tal.labels))', ...
    'Name', 'Please select the (Talairach) labels to add as VOI clusters...');
    if isequal(selok, 0) || ...
        isempty(regsel)
        return;
    end
    region = atl.tal.labels(regsel);

% specific region requested
else
    region = atl.tal.labels(strcmpi(makelabel(varargin{3}(:)'), atl.tal.shorts));
end

% prepare for processing
dotal2icbm = (ch.Stats.ICBM2TAL.Value > 0);
voi = ne_gcfg.voi;
if numel(voi) ~= 1 || ...
   ~isxff(voi, 'voi')
    voi = xff('new:voi');
end
vois = voi.VOI;
voistr = ch.Clusters.String;
if isempty(voistr)
    voistr = {};
else
    if ~iscell(voistr)
        voistr = cellstr(voistr);
    end
end
ti = min(numel(vois), numel(voistr)) + 1;
vois(ti+numel(region)-1).Name = '';
voistr(end+1:ti+numel(region)-1) = {''};
mfp = ch.MainFig.Pointer;
ch.MainFig.Pointer = 'watch';

% process
for rc = 1:numel(region)

    % get coordinates
    rcoords = tdlocal2(6, region{rc});

    % determine color
    voicolor = 1 ./ (std(rcoords, [], 1) + sqrt(eps));
    voicolor = round(255 .* voicolor ./ sqrt(sum(voicolor .* voicolor)));

    % split into right and left coords
    if any(diff(unique(rcoords(:, 1))) > 2)
        lcoords = rcoords(rcoords(:, 1) <= 0, :);
        rcoords(rcoords(:, 1) <= 0, :) = [];
    else
        lcoords = zeros(0, 3);
    end

    % convert to ICBM
    if dotal2icbm
        if ~isempty(lcoords)
            lcoords = unique(round(tal2icbm(lcoords)), 'rows');
            lcenter = round(mean(lcoords));
            ldist = sqrt(sum((lcoords - ones(size(lcoords, 1), 1) * lcenter) .^ 2, 2));
            [ldist, lsort] = sort(ldist);
            lcoords = lcoords(lsort, :);
        end
        if ~isempty(rcoords)
            rcoords = unique(round(tal2icbm(rcoords)), 'rows');
            rcenter = round(mean(rcoords));
            rdist = sqrt(sum((rcoords - ones(size(rcoords, 1), 1) * rcenter) .^ 2, 2));
            [rdist, rsort] = sort(rdist);
            rcoords = rcoords(rsort, :);
        end
    end

    % add to VOI object
    if ~isempty(lcoords)
        vois(ti).Name = sprintf('left %s', region{rc});
        vois(ti).Voxels = lcoords;
        vois(ti).NrOfVoxels = size(lcoords, 1);
        vois(ti).Color = voicolor;
        voistr{ti} = sprintf('LH_%s_%d_%d_%d - %d voxel(s)', ...
            makelabel(region{rc}), lcenter(1), lcenter(2), lcenter(3), ...
            size(lcoords, 1));
        ti = ti + 1;
        if ~isempty(rcoords)
            vois(ti).Name = sprintf('right %s', region{rc});
            vois(ti).Voxels = rcoords;
            vois(ti).NrOfVoxels = size(rcoords, 1);
            vois(ti).Color = voicolor;
            voistr{ti} = sprintf('RH_%s_%d_%d_%d - %d voxel(s)', ...
                makelabel(region{rc}), rcenter(1), rcenter(2), rcenter(3), ...
                size(rcoords, 1));
            ti = ti + 1;
        end
    else
        vois(ti).Name = sprintf('%s', region{rc});
        vois(ti).Voxels = rcoords;
        vois(ti).NrOfVoxels = size(rcoords, 1);
        vois(ti).Color = voicolor;
        voistr{ti} = sprintf('BH_%s_%d_%d_%d - %d voxel(s)', ...
            makelabel(region{rc}), rcenter(1), rcenter(2), rcenter(3), ...
            size(rcoords, 1));
        ti = ti + 1;
    end
end

% back to arrays
voi.VOI = vois;
ne_gcfg.voi = voi;
ch.Clusters.String = voistr;
ch.MainFig.Pointer = mfp;
drawnow;

% return vois and VOI object if requested
if nargout > 0
    varargout{1} = vois;
    if nargout > 1
        varargout{2} = voi;
    end
end
