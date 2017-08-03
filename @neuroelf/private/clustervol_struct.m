function cl = clustervol_struct(img, threshold, k)
% clustervol_struct  - create cluster struct (cell array) for TW tools
%
% FORMAT:       cl = clustervol_struct(img, thresholds, ks)
%
% Input fields:
%
%       img         filename of image (e.g. t-map) to cluster
%       thresholds  list of thresholds (1xT double array)
%                   values in the range of [-0.5 .. 0] are swapped in sign
%                   and interpreted as a probability threshold
%       ks          list of extent thresholds (1xT double array)
%
% Output fields:
%
%       cl          1xT cell array holding 1xC structs each with fields
%        .title     customary text
%        .threshold numeric threshold value
%        .M         matrix to compute real-world coordinates from voxels
%        .dim       image dimensions used to compute indices from voxels
%        .voxSize   used to compute volume of cluster
%        .name      customary text
%        .Z         extracted stats
%        .XYZmm     real-world coordinates of voxel coordinates
%        .XYZ       voxel coordinates
%        .numVox    number of voxels
%        .mm_center center of cluster (in mm, rounded)
%
% Note: These structures are used for the tools by Tor Wager! This
%       function uses spm_vol and spm_read_vols (for compatibility!)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 1:55 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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
if nargin < 3 || ...
   ~ischar(img) || ...
   ~isa(threshold, 'double') || ...
    isempty(threshold) || ...
   ~isa(k, 'double') || ...
    numel(threshold) ~= numel(k) || ...
    any(k(:) < 0)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument in call.' ...
    );
end

% try reading the image
try
    spm_defaults;
    v = spm_vol(img);
    y = spm_read_vols(v);
catch ne_eo;
    error( ...
        'neuroelf:SPMError', ...
        'Error using SPM function: ''%s''.', ...
        ne_eo.message ...
    );
end

% create output argument as cell array
clt = {cell2struct(cell(1,1,11), ...
    {'title', 'threshold', 'M', 'dim', 'voxSize', 'name', ...
     'Z', 'XYZmm', 'XYZ', 'numVox', 'mm_center'}, 3)};
cl = clt(ones(1, numel(threshold)));

if any(threshold < 0 & threshold > -0.5)
    yt = -y;
    yt(y == 0) = -1;
end

% iterate over thresholds
for tc = 1:numel(threshold)

    % cluster volume
    if threshold(tc) < 0 && ...
        threshold(tc) > -0.5
        [clt, cs] = clustervol(yt, -threshold(tc), k(tc), struct('mat', v.mat));
    else
        [clt, cs] = clustervol(y, threshold(tc), k(tc), struct('mat', v.mat));
    end

    % set initial fields of struct
    cl{tc}.threshold = threshold(tc);
    cl{tc}.M = v.mat;
    cl{tc}.dim = v.dim;
    cl{tc}.voxSize = diag(v.mat(1:3, 1:3))';
    cl{tc} = cl{tc}(ones(1, numel(cs)));

    % iterate over found clusters
    for cc = 1:numel(cs)

        % set further fields
        nv = numel(cs(cc).values);
        cl{tc}(cc).title = sprintf('Cluster of %d voxels from %s', nv, img);
        cl{tc}(cc).name = '';
        cl{tc}(cc).Z = cs(cc).values';
        cl{tc}(cc).XYZmm = cs(cc).rwcoords';
        cl{tc}(cc).XYZ = cs(cc).coords';
        cl{tc}(cc).numVox = nv;
        cl{tc}(cc).mm_center = round(mean(cs(cc).rwcoords));
    end
end
