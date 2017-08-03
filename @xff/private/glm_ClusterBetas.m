function [betas, cl, tables, vois] = glm_ClusterBetas(xo, vmp, opts)
% GLM::ClusterBetas  - create tables and extracts from GLM/VMPs
%
% FORMAT:       [betas, cl, tables, vois] = glm.ClusterBetas(vmp, [opts])
%
% Input fields:
%
%       vmp         VMP object to use for clustering (empty: GUI request)
%       opts        optional settings
%        .clconn    cluster connectivity ('face', {'edge'}, 'vertex)
%        .contrasts contrast definition, by default, each unique beta
%        .localmax  break down larger clusters threshold (default: Inf)
%        .maps      which maps of the VMP to cluster from (default: all)
%        .maxrad    maximum radius around peak (Inf, use 0 for peak voxel)
%        .minsize   size threshold override (don't use map settings)
%        .mni2tal   flag indicating that VMP is in MNI space (for TAL)
%        .savevois  save VOI objects to disk (automatic filenames)
%        .showneg   show negative tail override (don't use map settings)
%        .showpos   show positive tail override (don't use map settings)
%        .tdclient  add tdclient to tables (only in effect for nargout > 2)
%        .thresh    stats threshold override (don't use map specific)
%
% Output fields:
%
%       betas       cell array of SxCxV arrays with beta/contrast extracts
%       cl          cell array with ClusterTable outputs
%       tables      cell array with text tables
%       vois        cell array with VOI objects
%
% Note: only works for RFX GLMs!

% Version:  v1.1
% Build:    16020314
% Date:     Feb-03 2016, 2:27 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm')
    error('neuroelf:xff:badArgument', 'Bad object or argument in call.');
end
bc = xo.C;
if bc.ProjectTypeRFX ~= 1
    error('neuroelf:xff:badArgument', 'Only valid for RFX GLMs.');
end
reqvmp = false;
if nargin < 2 || numel(vmp) ~= 1 || ~xffisobject(vmp, true, 'vmp')
    reqvmp = true;
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'contrasts')
    opts.contrasts = [];
end
if ~isfield(opts, 'maps') || ~isa(opts.maps, 'double')
    opts.maps = [];
else
    opts.maps = unique(round(opts.maps(:)'));
    opts.maps(isinf(opts.maps) | isnan(opts.maps)) = [];
end
if ~isfield(opts, 'maxrad') || ~isa(opts.maxrad, 'double') || numel(opts.maxrad) ~= 1 || ...
    isinf(opts.maxrad) || isnan(opts.maxrad) || opts.maxrad < 0 || opts.maxrad > 256
    opts.maxrad = Inf;
end
if ~isfield(opts, 'savevois') || ~islogical(opts.savevois) || numel(opts.savevois) ~= 1
    opts.savevois = false;
end
if ~isfield(opts, 'thresh') || ~isa(opts.thresh, 'double') || numel(opts.thresh) ~= 1
    opts.thresh = 0;
end

% request VMP?
if reqvmp
    vmp = xff('*.vmp', 'Please select a VMP files to create tables from...');
    if ~xffisobject(vmp, true, 'vmp')
        if xffisobject(vmp, true)
            delete(vmp);
        end
        warning('neuroelf:xff:badObject', 'No VMP was selected.');
        betas = {};
        cl = {};
        tables = {};
        vois = {};
        return;
    end
end

% perform map selection
vmpc = vmp.C;
vmpm = vmpc.Map;
if isempty(opts.maps)
    opts.maps = 1:numel(vmpm);
else
    opts.maps = intersect(1:numel(vmpm), opts.maps);
end
nm = numel(opts.maps);

% init output variables
betas = cell(1, nm);
cl = cell(1, nm);
tables = cell(1, nm);
vois = cell(1, nm);

% iterate over maps
for mc = 1:nm

    % get table, clusters, etc.
    [cl{mc}, tables{mc}, sthreshmap, vois{mc}] = ...
        vmp_ClusterTable(vmp, opts.maps(mc), opts.thresh, opts);

    % patch VOIs with radius
    if ~isinf(opts.maxrad)

        % get voi content
        voic = vois{mc}.C;

        % iterate over VOIs in object
        for vc = 1:numel(voic.VOI)

            % compute distance from peak voxel
            dist = ones(size(voic.VOI(vc).Voxels, 1), 1) * ...
                voic.VOI(vc).Voxels(1, :) - voic.VOI(vc).Voxels;
            dist = sqrt(sum(dist .* dist, 2));

            % remove "over-far" voxels from VOI
            voic.VOI(vc).Voxels(dist > VOI_radius, :) = [];
            voic.VOI(vc).NrOfVoxels = size(voic.VOI(vc).Voxels, 1);
        end

        % store back
        vois{mc}.C = voic;
    end

    % save VOI?
    if opts.savevois
        aft_SaveAs(vois{mc}, sprintf('%s_map%02d.voi', ...
            strrep(vmp.FilenameOnDisk, '.vmp', ''), opts.maps(mc)));
    end
end

% extract contrast values from GLM
for ec = 1:nm
    betas{ec} = glm_VOIBetas(xo, vois{ec}, struct('c', opts.contrasts));
end

% clear VOIs if not used
if nargout < 4
    clearxffobjects(vois);
end

% clear VMP if requested
if reqvmp
    delete(vmp);
end
