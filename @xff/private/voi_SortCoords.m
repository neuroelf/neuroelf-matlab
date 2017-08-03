function xo = voi_SortCoords(xo, vl, so, res)
% VOI::SortCoords  - sort coordinates by distance from the center
%
% FORMAT:       [voi = ] voi.SortCoords([vl, [so [, res]]])
%
% Input fields:
%
%       vl         list of VOIs to sort (default: all)
%       so         sort order: {'ascend'}, 'descend'
%       res        force to resolution (default: not)
%
% Output fields:
%
%       voi        altered VOI object

% Version:  v1.1
% Build:    16020409
% Date:     MaY-05 2014, 3:59 PM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'voi')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
numvois = numel(bc.VOI);
if numvois == 0
    return;
end
if nargin < 2 || ~isa(vl, 'double') || isempty(vl)
    vl = 1:numvois;
elseif any(isinf(vl(:)) | isnan(vl(:)) | vl(:) < 1 | vl(:) > numvois)
    error('neuroelf:xff:badArgument', 'Invalid sub-VOI list.');
else
    vl = unique(round(vl(:)))';
end
if nargin < 3 || ~ischar(so) || isempty(so) || ~any(lower(so(1)) == 'ad')
    so = 'a';
else
    so = lower(so(1));
end
if nargin < 4 || ~isa(res, 'double') || numel(res) ~= 1 || ...
    isinf(res) || isnan(res) || res <= 0
    res = [];
else
    res = min(6, res);
end

% resolution given, detect given resolution
if ~isempty(res)
    detres = Inf;
    for vc = 1:numel(vl)
        for xc = 1:3
            detres = min(detres, min(diff(unique(bc.VOI(vl(vc)).Voxels(:, xc)))));
        end
    end
end

% sorting algorithm
for vc = 1:numel(vl)

    % get coordinates
    coords = bc.VOI(vl(vc)).Voxels;

    % remove bad entries (?)
    coords(any(isinf(coords) | isnan(coords), 2), :) = [];
    coordn = size(coords, 1);

    % find center
    center = mean(coords, 1);

    % resolution change
    if ~isempty(res) && res < detres

        % make room for more coordinates
        resfac = ceil(detres / res);
        residx = (1-resfac):(resfac-1);
        resk = numel(residx);
        rest = resk * resk * resk;

        % place into target
        ncoords = zeros(rest * coordn, 3);
        nci = 1;
        for xc = 1:resk
            for yc = 1:resk
                for zc = 1:resk
                    ncoords(nci:nci+coordn-1, :) = ...
                        [coords(:, 1) + residx(xc), coords(:, 2) + residx(yc), coords(:, 3) + residx(zc)];
                    nci = nci + coordn;
                end
            end
        end

        % replace
        coords = unique(ncoords, 'rows');
    end

    % compute distance
    distance = sqrt(sum((coords - center(ones(size(coords, 1), 1), :)) .^ 2, 2));

    % sort distance according to algorithm
    if so == 'a'
        [distance, dso] = sort(distance, 1);
    else
        [distance, dso] = sort(distance, 1, 'descend');
    end

    % sort coordinates
    bc.VOI(vl(vc)).Voxels = coords(dso, :);
    bc.VOI(vl(vc)).NrOfVoxels = size(coords, 1);

    % and if VoxelValues are present
    if isfield(bc.VOI, 'VoxelValues')

        % either discard them
        if ~isempty(res) && res < detres
            bc.VOI(vl(vc)).VoxelValues = [];

        % or sort them as well
        elseif numel(bc.VOI(vl(vc)).VoxelValues) == numel(dso)
            bc.VOI(vl(vc)).VoxelValues = bc.VOI(vl(vc)).VoxelValues(dso);
        end
    end
end

% store
xo.C = bc;
