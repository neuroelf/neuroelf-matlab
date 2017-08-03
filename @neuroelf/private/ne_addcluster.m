function ne_addcluster(varargin)
% ne_addcluster  - add cluster from marking to VOI list
%
% FORMAT:       ne_addcluster(SRC, EVT, source [, srcopts])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       source      source to get cluster from, currently MUST be 'ana'
%       srcopts     optional settings in 1x3 cell array OF STRINGS
%        {1}        1x1 code to get data from (e.g. '240')
%        {2}        1xN name for new ROI (e.g. 'New ROI')
%        {3}        1x3 RGB color (e.g. '255 0 0' for red)
%
% Example:
%
%   neuroelf_gui('setdrawmode', 3, {'8', '240', 's', '0 255', '0'});
%   neuroelf_gui('draw', [-20, -4, -22]);
%   neuroelf_gui('setdrawmode', 1);
%   ne_addcluster(0, 0, 'ana', {'240', 'left-Amygdala-sph8mm', '255 0 0'});

% Version:  v1.0
% Build:    14091314
% Date:     Sep-13 2014, 2:54 PM EST
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

% global variable and most used fields
global ne_gcfg;
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% only valid if current selection is a VMR, and voi is valid also
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
   ~any(strcmpi(varargin{3}, {'ana'})) || ...
    numel(cc.SliceVar) ~= 1 || ...
   ~isxff(cc.SliceVar, 'vmr') || ...
    numel(ne_gcfg.voi) ~= 1 || ...
   ~isxff(ne_gcfg.voi, 'voi')
    return;
end

% which code to select from, name and color of ROI
if nargin < 4 || ...
   ~iscell(varargin{4}) || ...
    numel(varargin{4}) ~= 3
    udata = {sprintf('  %d', cc.paint.code), 'New ROI', ...
        sprintf('  %d', floor(255.999 * rand(1, 3)))};
    udata = inputdlg( ...
        {'Color code to define ROI from:', 'Name of ROI:', 'ROI color (RGB):'}, ...
        'NeuroElf - Anatomical ROI definition', 1, udata);
    if ~iscell(udata) || ...
        numel(udata) ~= 3
        return;
    end
else
    udata = varargin{4}(:);
end
try
    udcode = str2double(udata{1});
    udname = udata{2}(:)';
    udrgbc = eval(['[' udata{3}(:)' ']']);
    if numel(udcode) ~= 1 || ...
        isinf(udcode) || ...
        isnan(udcode) || ...
        udcode < 0 || ...
        udcode > 255 || ...
        isempty(udname) || ...
        numel(udrgbc) ~= 3 || ...
        any(isinf(udrgbc) | isnan(udrgbc) | udrgbc < 0 | udrgbc > 255)
        return;
    end
    udcode = floor(udcode);
    udrgbc = udrgbc(:)';
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end

% no voxels, return
if ~any(cc.SliceVar.VMRData(:) == udcode)
    errordlg('ROI is empty, not adding to VOI.', 'NeuroElf GUI - error', 'modal');
    return;
end

% cluster
[cs, cv, l, c] = clustercoordsc(cc.SliceVar.VMRData == udcode);
if numel(cs) > 1
    split = questdlg('Split multiple clusters into separate VOIs?', ...
        'NeuroElf - user request');
    if ~ischar(split) || ...
       ~any(strcmpi(split, {'no', 'yes'}))
        return;
    end

    % split
    if strcmpi(split, 'yes')
        return;
    end
end

% get the coordinates of voxels in the VMR that match the current code
roivox = find(cc.SliceVar.VMRData(:) == udcode);

% convert to TAL coords for ROI
roivox = bvcoordconv(roivox, 'bvx2tal', cc.SliceVar.BoundingBox);

% sort by distance to center
meanvox = round(mean(roivox, 1));
voxdist = sum((roivox - meanvox(ones(size(roivox, 1), 1), :)) .^ 2, 2);
[voxdist, sortidx] = sort(voxdist);
roivox = roivox(sortidx, :);

% add to VOI
newvoi = emptystruct(fieldnames(ne_gcfg.voi.VOI), [1, 1]);
newvoi.Name = udname;
newvoi.Color = udrgbc;
newvoi.Voxels = roivox;
newvoi.NrOfVoxels = size(roivox, 1);
ne_gcfg.voi.VOI(end+1) = newvoi;
ne_gcfg.voi.NrOfVOIs = numel(ne_gcfg.voi.VOI);

% echo
if ne_gcfg.c.echo
    ne_echo({['voi.VOI(end + 1) = ' ...
        'struct(''Name'', ''%s'', ''Color'', %s, ''Voxels'', ' ...
        'roivox, ''NrOfVoxels'', size(roivox, 1));'], ...
        udname, any2ascii(udrgbc)});
end

% add name to list of VOI names
clnames = ch.Clusters.String;
if ~iscell(clnames)
    clnames = cellstr(clnames);
end
if numel(clnames) >= ne_gcfg.voi.NrOfVOIs
    clnames(ne_gcfg.voi.NrOfVOIs:end) = [];
    ch.Clusters.Value = [];
end
clnames{end+1} = sprintf('%s (%d voxels around [%d, %d, %d])', ...
    udname, size(roivox, 1), round(meanvox));
ch.Clusters.String = clnames;
ch.Clusters.Enable = 'on';
