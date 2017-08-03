function cfr = nlf_CoordinateFrame(xo, vol)
% NLF::CoordinateFrame  - get coordinate frame from NeuroeLF file
%
% FORMAT:       cframe = nlf.CoordinateFrame([vol])
%
% Input fields:
%
%       vol         1x2 volume number (default: first)
%
% Output fields:
%
%       cframe      struct with fields
%        .DimX/Y/Z  spatial dimension of image
%        .DimT      number of volumes for multi-volume image, default 1
%        .Dimensions combined list
%        .ResX/Y/Z  spatial resolution (in TAL/MNI/real world mm)
%        .Resolution combined list
%        .Slice1Center  center coordinate of first slice
%        .SliceNCenter  center coordinate of last slice
%        .SystemOrigin  coordinate system origin
%        .RowDir    equivalent to AffineTransX(1:3)
%        .ColDir    equivalent to AffineTransY(1:3)
%        .SlcDir    equivalent to AffineTransZ(1:3)
%        .IsRadiological  flag whether coordinate system is left-handed
%        .Trf       4x4 quaternion matrix, so that
%                   cframe.Trf * voxel := mm
%          - and -  inv(cframe.Trf) * mm := voxel

% Version:  v1.1
% Build:    16020314
% Date:     Feb-03 2016, 2:34 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'nlf')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isa(vol, 'double') || numel(vol) ~= 2 || ...
    any(isinf(vol) | isnan(vol) | vol < 1)
    vol = [1, 1];
else
    vol = floor(vol(:)');
end
bc = xo.C;
cfr = struct;

% make sure this is some volume file
if ~strcmp(bc.Intent, 'tri') && isempty(regexpi(bc.DimMeaning, 'xyz'))
    error('neuroelf:xff:badSubType', ...
        'Invalid subtype for call to NLF::CoordinateFrame.');
end

% for all but surfaces
if ~strcmp(bc.Intent, 'tri')

    % get dimensions
    cfr.DimX = bc.Size(bc.DimMeaning == 'x');
    cfr.DimY = bc.Size(bc.DimMeaning == 'y');
    cfr.DimZ = bc.Size(bc.DimMeaning == 'z');
    itd = [];
    if any(bc.DimMeaning == 't')
        cfr.DimT = bc.Size(bc.DimMeaning == 't');
        if any(bc.IndivTransformDims == find(bc.DimMeaning == 't'))
            itd = find(bc.IndivTransformDims == find(bc.DimMeaning == 't'));
        end
    else
        cfr.DimT = 1;
    end
    vol(1) = min(vol(1), cfr.DimT);
    cfr.Dimensions = [cfr.DimX, cfr.DimY, cfr.DimZ, cfr.DimT];

    % get resolution
    cfr.ResX = bc.Resolution(bc.DimMeaning == 'x');
    cfr.ResY = bc.Resolution(bc.DimMeaning == 'y');
    cfr.ResZ = bc.Resolution(bc.DimMeaning == 'z');

    % set full resolution element
    cfr.Resolution = [cfr.ResX, cfr.ResY, cfr.ResZ];

    % use global transformation and, possibly additional transformation
    if ~isempty(itd)
        itu = [1, 1];
        itu(itd) = vol;
        mfilec = bc.IndivTransforms(:, :, itu(1), itu(2));
    else
        mfilec = bc.GlobalTransform;
    end

% for surfaces
else
end

% compute slice 1 and N center coordinates (Trf is one-based!)
sl1c = mfilec * [0.5 + cfr.DimX / 2; 0.5 + cfr.DimY / 2; 0.5; 1];
slnc = mfilec * [0.5 + cfr.DimX / 2; 0.5 + cfr.DimY / 2; 0.5 + cfr.DimZ; 1];

% get direction components
dcp = -mfilec(1:3, 1:3);
dcp = dcp ./ repmat(sqrt(sum(dcp .* dcp, 2)), [1, 3]);
dcp(4, :) = -cross(dcp(1, :), dcp(2, :));

% make settings
cfr.Slice1Center = sl1c(1:3)';
cfr.SliceNCenter = slnc(1:3)';
cfr.SystemOrigin = mfilec(1:3, 4)';
cfr.RowDir = dcp(1, :);
cfr.ColDir = dcp(2, :);
cfr.SlcDir = dcp(3, :);
cfr.IsRadiological = (sum(dcp(3, :) .* dcp(4, :)) > 0);
cfr.Trf = mfilec;
