function tablesextract
% tablesextract  - create tables and extracts from GLM/VMPs
%
% FORMAT:       tablesextract
%
% for now, this is more a script than a function. no input/output arguments
% are used.
%
% variables assigned (into BASE workspace!) during the script:
%
% - vmp_tables: cell array with text tables (char)
% - vmp_vois:   cell array with VOI objects (one per table)
% - glm_betas:  cell array with extracted values
%
% NOTE: you *MUST* adapt the following part of this script for usefulness!!

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:30 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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

% This script-like function was written for Ethan Kross, UMich.

% options as struct
options = struct( ...
    'localmax', split_to_sub_clusters_size, ...
    'minsize',  k_threshold, ...
    'mni2tal',  mni_to_tal, ...
    'showneg',  show_neg, ...
    'showpos',  show_pos, ...
    'tdclient', TAL_lookup);

% init stored variables
vmp_tables = {};
vmp_vois   = {};
glm_betas  = {};
mpc = 0;

% request GLM file
glm = xff('*.glm', 'Please select the GLM file for data extraction...');
vmp = 1;

% continue until no more selected
while ~isempty(vmp)

    % load vmp
    vmp = xff('*.vmp', 'Please select a VMP files to create tables from...');
    if ~isxff(vmp, 'vmp')
        if isxff(vmp)
            vmp.ClearObject;
            continue;
        end
        break;
    end

    % iterate over maps
    if isempty(map_sel)
        my_map_sel = 1:numel(vmp.Map);
    else
        my_map_sel = intersect(1:numel(vmp.Map), map_sel(:)');
    end
    for mc = my_map_sel

        % patch options for defaults
        if isempty(k_threshold)
            if vmp.Map(mc).EnableClusterCheck
                opts.minsize = vmp.Map(mc).ClusterSize;
            else
                opts.minsize = 1;
            end
        end
        if isempty(show_neg)
            opts.showneg = (vmp.Map(mc).ShowPositiveNegativeFlag > 1);
        end
        if isempty(show_pos)
            opts.showneg = (mod(vmp.Map(mc).ShowPositiveNegativeFlag, 2) > 0);
        end

        % get table, clusters, etc.
        [clusters_struct, text_table, threshmap, clusters_voi] = ...
            vmp.ClusterTable(mc, threshold, options);

        % store table
        mpc = mpc + 1;
        vmp_tables{mpc} = text_table;

        % patch VOIs with radius
        if ~isinf(VOI_radius)

            % iterate over VOIs in object
            for vc = 1:numel(clusters_voi.VOI)

                % compute distance from peak voxel
                dist = ones(size(clusters_voi.VOI(vc).Voxels, 1), 1) * ...
                    clusters_voi.VOI(vc).Voxels(1, :) - ...
                    clusters_voi.VOI(vc).Voxels;
                dist = sqrt(sum(dist .* dist, 2));

                % remove "over-far" voxels from VOI
                clusters_voi.VOI(vc).Voxels(dist > VOI_radius, :) = [];
                clusters_voi.VOI(vc).NrOfVoxels = ...
                    size(clusters_voi.VOI(vc).Voxels, 1);
            end
        end

        % store VOIs in array
        vmp_vois{mpc} = clusters_voi;

        % save VOI?
        if save_VOIs
            clusters_voi.SaveAs(sprintf('%s_map%02d.voi', ...
                strrep(vmp.FilenameOnDisk, '.vmp', ''), mc));
        end
    end

    % clear vmp in memory
    vmp.ClearObject;
end

% if GLM ok
if numel(glm) == 1 && ...
    isxff(glm, 'glm')

    % extract contrast values from GLM
    glm_betas = cell(1, mpc);
    for ec = 1:mpc
        glm_betas{ec} = glm.VOIBetas(vmp_vois{ec}, struct('c', glm_contrasts));
    end
end

% clear GLM
if numel(glm) == 1 && ...
    isxff(glm)
    glm.ClearObject;
end

% assign into BASE WS
assignin('base', 'vmp_tables', vmp_tables);
assignin('base', 'vmp_vois', vmp_vois);
assignin('base', 'glm_betas', glm_betas);
