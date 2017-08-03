function plp = plp_ReduceToMKDAMask(xo, vmp, opts)
% PLP::ReduceToMKDAMask  - reduce coordinates to those within a mask
%
% FORMAT:       rplp = plp.MKDA(vmp [, opts])
%
% Input fields:
%
%       vmp         VMP object containing the meta map (for clusters)
%       opts        1x1 struct with options
%        .cond      conditional statement, default from VMP.Map.RunTimeVars
%        .fixcoords fix coordinates to nearest integer (VMP space, false)
%        .map       map number for mask (default: 1)
%        .mapthresh thresholds to apply (default: [0.05, 0.01])
%        .unique    only use unique points (within study, default: true)
%        .unitcols  column names for study and/or statistical unit
%                   default from plp.RunTimeVars.StudyColumns
%        .unitsel   selection of units (default: from VMP)
%
% Output fields:
%
%       rplp        reduced PLP object (new object!)
%
% Using: findfirst, lsqueeze, multimatch, psetdists.

% Version:  v1.1
% Build:    16021210
% Date:     Feb-12 2016, 10:22 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2013, 2014, 2016, Jochen Weber
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

% neuroelf library
global ne_methods;
findfirst = ne_methods.findfirst;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'plp') || ...
    numel(vmp) ~= 1 || ~xffisobject(vmp, true, 'vmp')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
cn = bc.ColumnNames(:);
vmpc = vmp.C;
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'cond') || ~ischar(opts.cond) || isempty(opts.cond)
    opts.cond = 0;
else
    opts.cond = opts.cond(:)';
end
if ~isfield(opts, 'fixcoords') || ~islogical(opts.fixcoords) || numel(opts.fixcoords) ~= 1
    opts.fixcoords = false;
end
if ~isfield(opts, 'map') || ~isa(opts.map, 'double') || numel(opts.map) ~= 1 || ...
    isinf(opts.map) || isnan(opts.map) || opts.map < 1 || opts.map > numel(vmpc.Map)
    opts.map = 1;
else
    opts.map = fix(opts.map);
end
map = vmpc.Map(opts.map);
if ~isfield(map, 'RunTimeVars') || ~isstruct(map.RunTimeVars) || numel(map.RunTimeVars) ~= 1 || ...
   ~isfield(map.RunTimeVars, 'AlphaSimMKDA') || ~isfield(map.RunTimeVars, 'Condition') || ...
   ~isfield(map.RunTimeVars, 'UnitColumn') || ~isfield(map.RunTimeVars, 'UnitID')
    error('neuroelf:xff:badArgument', 'Invalid VMP object (not a valid META map).');
end
rtv = map.RunTimeVars;
asthresh = rtv.AlphaSimMKDA{2};
if isa(opts.cond, 'double')
    opts.cond = rtv.Condition;
end
uc = rtv.UnitColumn;
uids = unique(rtv.UnitID(:))';
if ~isfield(opts, 'mapthresh') || ~isa(opts.mapthresh, 'double') || isempty(opts.mapthresh)
    opts.mapthresh = [0.05, 0.01];
end
opts.mapthresh = opts.mapthresh(:);
opts.mapthresh(isinf(opts.mapthresh) | isnan(opts.mapthresh) | opts.mapthresh <= 0) = [];
mts = asthresh(:, 1)';
mtd = opts.mapthresh * ones(size(mts)) - ones(numel(opts.mapthresh), 1) * mts;
opts.mapthresh = mts(min(abs(mtd), [], 1) < 1e-7);
if ~isfield(opts, 'unique') || ~islogical(opts.unique) || numel(opts.unique) ~= 1
    opts.unique = true;
end
if ~isfield(opts, 'unitcols') || ~iscell(opts.unitcols) || ~any(numel(opts.unitcols) == [1, 2]) || ...
    any(ne_methods.multimatch(opts.unitcols(:), bc.RunTimeVars.StudyColumns(:)) == 0)
    opts.unitcols = bc.RunTimeVars.StudyColumns(:);
else
    opts.unitcols = opts.unitcols(:);
    if ~any(strcmpi(opts.unitcols, cn{uc}))
        opts.unitcols(end+1) = cn(uc);
    end
end
if ~isfield(opts, 'unitsel') || ~isa(opts.unitsel, 'double') || isempty(opts.unitsel) || ...
    any(isinf(opts.unitsel(:)) | isnan(opts.unitsel(:)))
    opts.unitsel = uids;
else
    opts.unitsel = unique(opts.unitsel(:))';
end

% cluster meta-VMP (map 1) at two thresholds
cl = cell(numel(opts.mapthresh), 1);
for cc = 1:numel(cl)
    astm = (asthresh(:, 1) == opts.mapthresh(cc));
    [cl{cc}, clt] = aft_ClusterTable(vmp, opts.map, ...
        asthresh(astm, 2), struct('minsize', asthresh(astm, 3)));
    cl{cc} = cat(1, cl{cc}.rwcoords);
end

% get combined list of coordinates (VMP-resolution mask)
coords = cat(1, cl{:});
coords = unique(coords, 'rows');

% get list of points (including study and contrast number)
lcn = lower(cn);
ucn = cn{uc};
xc = findfirst(strcmp(lcn, 'x'));
yc = findfirst(strcmp(lcn, 'y'));
zc = findfirst(strcmp(lcn, 'z'));
sc = zeros(1, numel(opts.unitcols));
for cc = 1:numel(sc)
    sc(cc) = findfirst(strcmpi(lcn, opts.unitcols{cc}));
end
csel = [xc, yc, zc, sc, setdiff(1:numel(cn), [xc, yc, zc, sc])];
cn = bc.ColumnNames(csel);
cn = cn(:);
uc = findfirst(strcmp(cn, ucn));
xyz = bc.Points(:, csel);

% generate new PLP
plp = aft_CopyObject(xo);
plpc = plp.C;
plpc.ColumnNames = cn;
plpc.Points = xyz;

% minimal distance column
dc = findfirst(strcmp(cn, 'MaxMaskDist'));
if isempty(dc)
    plpc.ColumnNames{end+1} = 'MaxMaskDist';
    dc = numel(plpc.ColumnNames);
    plpc.RunTimeVars.ColumnIsText.MaxMaskDist = false;
end

% extend list with minimal distance to any (closest) voxel in mask
plpc.Points(:, dc) = min(ne_methods.psetdists(xyz(:, 1:3), coords), [], 2);
plp.C = plpc;

% select points (with max distance of 5 mm)
try
    sel = plp_Select(plp, ['(' opts.cond ') & $MaxMaskDist <= 5']);
catch xfferror
    rethrow(xfferror);
end

% apply selection
plpc.Points = plpc.Points(sel, :);

% apply unit selection
usel = plpc.Points(:, uc);
plpc.Points(~any(usel * ones(size(uids)) == ones(size(usel)) * uids, 2), :) = [];

% prune list
[xyzs, xyzi] = sortrows(plpc.Points(:,1:3));
sd = sum(abs(diff(plpc.Points(xyzi, 1:3))), 2);
ti = xyzi(find(sd == 0) + 0);
qi = xyzi(find(sd == 0) + 1);
txyz = plpc.Points(ne_methods.lsqueeze(([ti(:), qi(:)])'), 4);
plpc.Points(sort(qi(txyz(1:2:end) == txyz(2:2:end))), :) = [];
plpc.NrOfPoints = size(plpc.Points, 1);
plp.C = plpc;

% extend list once more and store 1 / sqrt(number of studies per contrast)
fc = findfirst(strcmp(cn, 'SCFactor'));
if isempty(fc)
    plpc.ColumnNames{end+1} = 'SCFactor';
    fc = numel(plpc.ColumnNames);
    plpc.RunTimeVars.ColumnIsText.SCFactor = false;
end
plpc.Points(:, fc) = 1;
study = plpc.Points(:, 4);
ustudy = unique(study);
for sc = 1:numel(ustudy)
    plpc.Points(study == ustudy(sc), fc) = ...
        1 ./ sqrt(numel(unique(plpc.Points(study == ustudy(sc), 5))));
end
plp.C = plpc;
