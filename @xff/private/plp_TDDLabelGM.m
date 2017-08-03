function [labels, lf] = plp_TDDLabelGM(xo, opts)
% PLP::TDDLabelGM  - label points according to ICBM-TalairachDD
%
% FORMAT:       [labels, lf] = plp.TDDLabelGM([opts])
%
% Input fields:
%
%       opts        1x1 struct with optional settings
%        .radius    radius to scan for gray-matter labels
%
% Output fields:
%
%       labels      Px1 cell array with chosen labels
%       lf          Px1 double array with fractional value of frequency
%
% Using: findfirst, histcount, multimatch.

% Version:  v1.1
% Build:    16021210
% Date:     Feb-12 2016, 10:17 AM EST
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
histcount = ne_methods.histcount;

% requires TD image in ICBM format
tpath = [neuroelf_path('tal'), filesep];
if exist([tpath, 'talairach_ICBMnorm.nii'], 'file') ~= 2
    error('neuroelf:xff:missingFile', ...
        'This method requires talairach_ICBMnorm.nii. Please run neuroelf_makefiles(''all'').');
end
tobj = {[]};

% load TD image (in ICBM space)
try
    tobj{1} = xff([tpath, 'talairach_ICBMnorm.nii']);
    tc = tobj{1}.C;
    if istransio(tc.VoxelData)
        tc.VoxelData = resolve(tc.VoxelData);
        tobj{1}.C = tc;
    end

    % get labels
    tallabels = tc.RunTimeVars.TALLabels(:);
    usetlabel = false(size(tallabels));
    maxlabel = numel(tallabels);

    % get graymatter status
    graymatter = ~cellfun('isempty', regexpi(tallabels, 'gr.y\s+m.tter'));
catch xfferror
    clearxffobjects(tobj);
    rethrow(xfferror);
end

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'plp')
    clearxffobjects(tobj);
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
cn = lower(bc.ColumnNames(:));
xc = findfirst(strcmp(cn, 'x'));
yc = findfirst(strcmp(cn, 'y'));
zc = findfirst(strcmp(cn, 'z'));
if numel([xc, yc, zc]) ~= 3 || numel(cn) ~= size(bc.Points, 2)
    clearxffobjects(tobj);
    error('neuroelf:xff:badObject', ...
        'PLP object must have X, Y, and Z columns and correct number of column labels.');
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'radius') || ~isa(opts.radius, 'double') || ...
    numel(opts.radius) ~= 1 || isinf(opts.radius) || isnan(opts.radius)
    opts.radius = 8;
else
    opts.radius = max(1, min(opts.radius, 15));
end

% pre-generate output
labels = cell(size(bc.Points, 1), 1);
lf = zeros(size(labels));

% generate spherical coordinates
r = ceil(opts.radius);
[x, y, z] = ndgrid(-r:r, -r:r, -r:r);
s = [x(:), y(:), z(:)];
s(sum(s .* s, 2) > (opts.radius * opts.radius), :) = [];
sones = ones(size(s, 1), 1);

% get points
p = bc.Points(:, [xc, yc, zc]);
np = size(p, 1);

% null label
nl = findfirst(strcmp(bc.Labels(:), '*,*,*,*,*'));
if isempty(nl)
    bc.Labels{end+1} = '*,*,*,*,*';
    nl = numel(bc.Labels);
    bc.NrOfLabels = nl;
end
labeln = zeros(np, 1);

% label column
lc = findfirst(strcmp(cn, 'tddngm'));
if isempty(lc)
    bc.ColumnNames{end+1} = 'TDDNGM';
    lc = numel(bc.ColumnNames);
    bc.Points(:, end+1) = nl;
    bc.NrOfColumns = lc;
end

% iterate over points
for pc = 1:np

    % sample TDD with sphere
    sl = aft_SampleData3D(tobj{1}, sones * p(pc, :) + s);
    sl(isnan(sl) | sl == 0) = [];

    % nothing left?
    if isempty(sl)
        continue;
    end

    % histogram over labels
    slh = histcount(sl, 1, maxlabel, 1);

    % remove entries that are not gray matter
    slhi = find(slh);
    slh(slhi(~graymatter(slhi))) = 0;

    % get max and max position
    [slhm, slhmi] = max(slh(:));

    % store max (position = label number)
    labeln(pc) = slhmi;

    % mark as used
    usetlabel(slhmi) = true;

    % and compute frequency (among gray matter labels)
    lf(pc) = slhm / sum(slh);
end

% lookup labels
labels(labeln == 0) = {'*,*,*,*,*'};
labels(labeln > 0) = tallabels(labeln(labeln > 0));

% compress label indices
usedlabels = find(usetlabel);

% match agains existing labels
lm = ne_methods.multimatch(tallabels(usedlabels), bc.Labels(:));

% generate target function
tfunc = zeros(1, maxlabel);

% assign existing labels
tfunc(usedlabels(lm > 0)) = lm(lm > 0);

% generate new labels
numlab = numel(bc.Labels);
tfunc(usedlabels(lm == 0)) = (numlab+1):(numlab+sum(lm == 0));
bc.Labels(end+1:end+sum(lm == 0)) = tallabels(usedlabels(lm == 0));
bc.NrOfLabels = numel(bc.Labels);
labeln(labeln > 0) = tfunc(labeln(labeln > 0));
labeln(labeln == 0) = nl;

% set in object
bc.Points(:, lc) = labeln;

% update PLP
xo.C = bc;

% clear TDD object
clearxffobjects(tobj);
