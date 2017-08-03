function v = createmarkervmp(coords, opts)
% createmarkervmp  - create a VMP map with markers (electrode indicators)
%
% FORMAT:       mvmp = createmarkervmp(coords [, opts])
%
% Input fields:
%
%       coords      Cx3 coordinate list
%       opts        optional settings
%        .bbox      bounding box (default: full 256x256x256 matrix) in TAL
%        .negative  use negative scale flag, default: false
%        .radmax    maximal radius of spherical markers (default: 10)
%        .radvis    initially visible radius of spherical markers (4)
%        .saveas    if given, a .SaveAs(...) attempt will be made
%        .saveclr   flag, VMP cleared on successful saving (default: true)
%        .shape     either of 'box', {'sphere'}
%        .trans     4x4 transformation matrix required for BV's TAL space
%
% Output fields:
%
%       mvmp        VMP object with map (unless cleared after saving)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:15 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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
if nargin < 1 || ...
   ~isa(coords, 'double') || ...
    size(coords, 2) ~= 3 || ...
    any(isinf(coords(:)) | isnan(coords(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument supplied.' ...
    );
end
if nargin < 2 || ...
    numel(opts) ~= 1 || ...
   ~isstruct(opts)
    opts = struct;
end
if ~isfield(opts, 'bbox') || ...
   ~isa(opts.bbox, 'double') || ...
   ~isequal(size(opts.bbox), [2, 3]) || ...
    any(isinf(opts.bbox(:)) | isnan(opts.bbox(:))) || ...
    any(diff(opts.bbox) < 1) || ...
    any(opts.bbox(:) < 0 | opts.bbox(:) > 256)
    opts.bbox = [0, 0, 0; 256, 256, 256];
end
if ~isfield(opts, 'negative') || ...
   ~islogical(opts.negative) || ...
    numel(opts.negative) ~= 1
    opts.negative = false;
end
if ~isfield(opts, 'radmax') || ...
   ~isa(opts.radmax, 'double') || ...
    numel(opts.radmax) ~= 1 || ...
    isinf(opts.radmax) || ...
    isnan(opts.radmax) || ...
    opts.radmax < 1 || ...
    opts.radmax > 20
    opts.radmax = 10;
end
if ~isfield(opts, 'radvis') || ...
   ~isa(opts.radvis, 'double') || ...
    numel(opts.radvis) ~= 1 || ...
    isinf(opts.radvis) || ...
    isnan(opts.radvis) || ...
    opts.radvis < 1 || ...
    opts.radvis > opts.radmax
    opts.radvis = max(1, 0.5 * opts.radmax);
end
if ~isfield(opts, 'saveas') || ...
   ~ischar(opts.saveas)
    opts.saveas = '';
else
    opts.saveas = opts.saveas(:)';
end
if ~isfield(opts, 'saveclr') || ...
   ~islogical(opts.saveclr) || ...
    numel(opts.saveclr) ~= 1
    opts.saveclr = true;
end
if ~isfield(opts, 'shape') || ...
   ~ischar(opts.shape) || ...
   ~any(strcmpi(opts.shape(:)', {'box', 'sphere'}))
    opts.shape = 'sphere';
else
    opts.shape = lower(opts.shape(:)');
end
if ~isfield(opts, 'trans') || ...
   ~isa(opts.trans, 'double') || ...
   ~isequal(size(opts.trans), [4, 4]) || ...
    any(isinf(opts.trans(:)) | isnan(opts.trans(:))) || ...
    any(opts.trans(4, :) ~= [0, 0, 0, 1])
    opts.trans = [];
end

% create new VMP
v = xff('new:vmp');

% make global settings
v.XStart = opts.bbox(1, 2);
v.YStart = opts.bbox(1, 3);
v.ZStart = opts.bbox(1, 1);
v.XEnd = opts.bbox(2, 2);
v.YEnd = opts.bbox(2, 3);
v.ZEnd = opts.bbox(2, 1);
v.Resolution = 1;

% make map defaults
v.Map.Type = 1;
v.Map.LowerThreshold = opts.radmax - opts.radvis;
v.Map.UpperThreshold = opts.radmax;
v.Map.EnableClusterCheck = 0;

% set new map data to correct size
ms = diff(opts.bbox(:, [2, 3, 1]));
v.Map.VMPData = single(zeros(ms));

% transformation needed
if ~isempty(opts.trans)
    coords(:, 4) = 1;
    coords = (opts.trans * coords')';
    coords(:, 4) = [];
end

% convert coordinates to vmp-Matlab indices
cv = bvcoordconv(coords, 'tal2bvc', v.BoundingBox);

% iterate over coordinates
for cc = 1:size(cv, 1)

    % create coordinate grid
    [x,y,z] = ndgrid( ...
        cv(cc, 1)-10:cv(cc, 1)+10, ...
        cv(cc, 2)-10:cv(cc, 2)+10, ...
        cv(cc, 3)-10:cv(cc, 3)+10);

    % combine into one variable
    xyz = round([x(:), y(:), z(:)]);

    % remove "invalid" coordinates (< 1, > max)
    xyz(xyz(:, 1) < 1 | xyz(:, 1) > ms(1) | ...
        xyz(:, 2) < 1 | xyz(:, 2) > ms(2) | ...
        xyz(:, 3) < 1 | xyz(:, 3) > ms(3), :) = [];

    % compute distance to center for spheres
    if strcmp(opts.shape, 'sphere')
        d = sqrt(sum((xyz - ones(size(xyz, 1), 1) * cv(cc, :)) .^2, 2));

    % for boxes
    else
        d = max(abs(xyz - ones(size(xyz, 1), 1) * cv(cc, :)), [], 2);
    end

    % set voxels to value - radius
    ixyz = sub2ind(ms, xyz(:, 1), xyz(:, 2), xyz(:, 3));
    if opts.negative
        v.Map.VMPData(ixyz) = min(v.Map.VMPData(ixyz), d - opts.radmax);
    else
        v.Map.VMPData(ixyz) = max(v.Map.VMPData(ixyz), opts.radmax - d);
    end
end

% map name
if opts.negative
    v.Map.Name = sprintf('Marker VMP (%d points, negative)', size(coords, 1));
else
    v.Map.Name = sprintf('Marker VMP (%d points, positive)', size(coords, 1));
end

% SaveAs ?
if ~isempty(opts.saveas)
    try
        % save VMP
        v.SaveAs(opts.saveas);
        if opts.saveclr
            v.ClearObject;
            v = [];
        end
    catch ne_eo;
        % give warning
        warning( ...
            'neuroelf:SaveAsFailed', ...
            'Error using SaveAs method on VMP: %s.', ...
            ne_eo.message ...
        );
    end
end
