function [bcorrs, blist, conds, tvi] = glm_VOICorrelations(xo, vfile, opts)
% GLM::VOICorrelations  - return CxCxS cross-correlation matrices
%
% FORMAT:       [bcorr, betas, conds, tvi] = glm.VOICorrelations(voi, opts)
%
% Input fields:
%
%       voi         VOI object
%       opts        optional settings
%        .conds     condition selection (single regexp or cell array)
%        .maxvox    maximum number of voxels (default: Inf)
%        .res1mm    check that VOI has a 1mm resolution (default: false)
%        .voisel    VOI selection (default: 1, must be single number!)
%
% Output fields:
%
%       bcorr       Condition-by-Condition-by-Subject correlation matrices
%       betas       VxCxS betas extracted from voxels
%       conds       Cx1 selected conditions
%       tvi         Sx1 voxel indices the betas were extracted from
%
% Note: it is HIGHLY recommended to issue glm.LoadTransIOData; first!!
%
% Using: findfirst, histcount, multimatch.

% Version:  v1.1
% Build:    16020409
% Date:     Feb-04 2016, 9:41 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2014, 2016, Jochen Weber
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

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm') || ...
   ((numel(vfile) ~= 1 || ~xffisobject(vfile, true, 'voi')) && ...
    (~isa(vfile, 'double') || size(vfile, 2) ~= 3 || isempty(vfile) || ...
    any(isinf(vfile(:)) | isnan(vfile(:)) | vfile(:) < -128 | vfile(:) > 128)))
    error('neuroelf:xff:badArgument', 'Invalid object handle in call.');
end

% must be RFX
bc = xo.C;
if bc.ProjectTypeRFX ~= 1
    error('neuroelf:xff:badObject', 'Only valid for RFX-GLMs.');
end

% get VOI contents
vc = vfile.C;

% and get subject predictors and map size
conds = glm_SubjectPredictors(xo);
nsp = numel(conds);
mne = numel(bc.GLMData.RFXGlobalMap);

% check options
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'conds') || isempty(opts.conds) || (~ischar(opts.conds) && ~iscell(opts.conds))
    opts.conds = 1:numel(conds);
elseif ischar(opts.conds)
    opts.conds = find(~cellfun('isempty', regexpi(conds, opts.conds(:)')));
else
    opts.conds = find(ne_methods.multimatch(conds, opts.conds(:), true) > 0);
end
opts.conds = opts.conds(:)';
conds = conds(opts.conds, 1);
ngp = numel(conds);
if ~isfield(opts, 'maxvox') || ~isa(opts.maxvox, 'double') || numel(opts.maxvox) ~= 1 || ...
    isnan(opts.maxvox) || opts.maxvox <= 1
    opts.maxvox = Inf;
else
    opts.maxvox = max(8, opts.maxvox);
end
if ~isfield(opts, 'res1mm') || ~islogical(opts.res1mm) || numel(opts.res1mm) ~= 1
    opts.res1mm = false;
end
if ~isfield(opts, 'voisel') || ~isa(opts.voisel, 'double') || numel(opts.voisel) ~= 1 || ...
    isinf(opts.voisel) || isnan(opts.voisel) || opts.voisel < 1 || opts.voisel > numel(vc.VOI)
    opts.voisel = 1;
else
    opts.voisel = round(opts.voisel);
end

% check resolution
numvox = size(vc.VOI(opts.voisel).Voxels, 1);
if opts.res1mm
    detres = min(min( ...
        min(diff(unique(vc.VOI(opts.voisel).Voxels(:, 1)))), ...
        min(diff(unique(vc.VOI(opts.voisel).Voxels(:, 2))))), ...
        min(diff(unique(vc.VOI(opts.voisel).Voxels(:, 3)))));
    if detres > 1
        voi_SortCoords(vfile, opts.voisel, 'a', 1);
        vc = vfile.C;
        vc.VOI(opts.voisel).Voxels(ceil(numvox * (detres ^ 3) + 1):end, :) = [];
        vc.VOI(opts.voisel).NrOfVoxels = size(vc.VOI(opts.voisel).Voxels, 1);
        vfile.C = vc;
    end
end

% hard-limit number of voxels
numvox = min(opts.maxvox, vc.VOI(opts.voisel).NrOfVoxels);

% get VOI Betas
[blist, vui, vbi] = glm_VOIBetas(xo, vfile, struct('vl', opts.voisel));

% for each subject
vui = vbi;
for sc = 1:numel(vbi)

    % ensure that the same number of voxels is indeed available
    numvox = min(numvox, numel(unique(vbi{sc}(vbi{sc} > 0))));

    % get unique voxel indices (in same order)
    vui{sc} = unique(vbi{sc}(vbi{sc} > 0), 'stable');
end

% create target voxel indices
tvi = repmat({zeros(numvox, 1)}, numel(vbi), 1);

% for each subject
for sc = 1:numel(tvi)

    % compute histogram over voxels
    voxh = histcount(vbi{sc}, 1, mne, 1);

    % target index
    ti = 1;

    % repeat until all indices are filled
    while ti <= numvox

        % find most-implicated voxel(s)
        [mvv, mvx] = max(voxh);

        % only one
        matches = find(voxh == mvv);
        if numel(matches) == 1

            % add
            tvi{sc}(ti) = mvx;
            voxh(mvx) = 0;

        % multiple voxels
        else

            % find closest one
            svi = findfirst(any(repmat(vui{sc}, 1, numel(matches)) == repmat(matches(:)', numel(vui{sc}), 1), 2));

            % add
            tvi{sc}(ti) = vui{sc}(svi);
            voxh(vui{sc}(svi)) = 0;
        end

        % next voxel
        ti = ti + 1;
    end
end

% create betas array
blist = zeros(numvox, ngp, numel(tvi));
bcorrs = zeros(ngp, ngp, numel(tvi));

% access and store
for sc = 1:numel(tvi)

    % reshape maps
    maps = reshape(bc.GLMData.Subject(sc).BetaMaps, mne, nsp);

    % then access
    blist(:, :, sc) = maps(tvi{sc}, opts.conds);
    bcorrs(:, :, sc) = corrcoef(blist(:, :, sc));
end
