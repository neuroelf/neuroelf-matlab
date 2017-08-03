function cfr = head_CoordinateFrame(xo, vol)
% HEAD::CoordinateFrame  - get coordinate frame from HEAD/BRIK image
%
% FORMAT:       cframe = head.CoordinateFrame([sbrick])
%
% Input fields:
%
%       sbrick      1x1 sub-brick number (default: last for 4-D files)
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
% Build:    16021412
% Date:     Feb-14 2016, 12:56 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'head')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isa(vol, 'double') || numel(vol) ~= 1 || isinf(vol) || isnan(vol) || vol < 1
    vol = Inf;
else
    vol = floor(vol);
end
bc = xo.C;
cfr = struct;

% get dimensions
cfr.DimX = bc.DataDimensions(1);
cfr.DimY = bc.DataDimensions(2);
cfr.DimZ = bc.DataDimensions(3);
cfr.DimT = numel(bc.Brick);
cfr.Dimensions = [cfr.DimX, cfr.DimY, cfr.DimZ, cfr.DimT];

% get resolution
bco = bc.Origin;
bcr = bc.Resolution;
cfr.ResX = abs(bcr(1));
cfr.ResY = abs(bcr(2));
cfr.ResZ = abs(bcr(3));

% set full resolution element
cfr.Resolution = [cfr.ResX, cfr.ResY, cfr.ResZ];

% get axesorder
axo = floor(0.5 .* bc.AxesOrientation);
if ~any(axo == 0) || ~any(axo == 1) || ~any(axo == 2)
    error('neuroelf:xff:badArgument', 'Invalid axes orientation combination.');
end

% combine for switching, possible values: 21, 15, 19, 7, 11, 5
axo = axo(1) + 3 * axo(2) + 9 * axo(3);

% compute 1-based transformation matrix (as used in SPM)
switch (axo)
    case 5
        mfilec = [ ...
              0   ,   0   , bcr(3), bco(3) - bcr(3); ...
              0   , bcr(2),   0   , bco(2) - bcr(2); ...
            bcr(1),   0   ,   0   , bco(1) - bcr(1); ...
              0   ,   0   ,   0   ,        1];
    case 7
        mfilec = [ ...
              0   ,   0   , bcr(3), bco(3) - bcr(3); ...
            bcr(1),   0   ,   0   , bco(1) - bcr(1); ...
              0   , bcr(2),   0   , bco(2) - bcr(2); ...
              0   ,   0   ,   0   ,        1];
    case 11
        mfilec = [ ...
              0   , bcr(2),   0   , bco(2) - bcr(2); ...
              0   ,   0   , bcr(3), bco(3) - bcr(3); ...
            bcr(1),   0   ,   0   , bco(1) - bcr(1); ...
              0   ,   0   ,   0   ,        1];
    case 15
        mfilec = [ ...
            bcr(1),   0   ,   0   , bco(1) - bcr(1); ...
              0   ,   0   , bcr(3), bco(3) - bcr(3); ...
              0   , bcr(2),   0   , bco(2) - bcr(2); ...
              0   ,   0   ,   0   ,        1];
    case 19
        mfilec = [ ...
              0   , bcr(2),   0   , bco(2) - bcr(2); ...
            bcr(1),   0   ,   0   , bco(1) - bcr(1); ...
              0   ,   0   , bcr(3), bco(3) - bcr(3); ...
              0   ,   0   ,   0   ,        1];
    case 21
        mfilec = [ ...
            bcr(1),    0  ,    0  , bco(1) - bcr(1); ...
              0   , bcr(2),    0  , bco(2) - bcr(2); ...
              0   ,    0  , bcr(3), bco(3) - bcr(3); ...
              0   ,    0  ,    0  ,        1];
    otherwise
        warning('neuroelf:xff:invalidObject', ...
            'Invalid axes orientation: %d-%d-%d.', bc.AxesOrientation);
end

% perform x/y-flip to get to TAL
mfilec(1:2, :) = -mfilec(1:2, :);

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
