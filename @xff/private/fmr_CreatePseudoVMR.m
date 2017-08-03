function vmr = fmr_CreatePseudoVMR(xo, vnum, flipped)
% FMR::CreatePseudoVMR  - create a PseudoVMR from one time point
%
% FORMAT:       vmr = fmr.CreateVMR([vnum, flipped])
%
% Input fields:
%
%       vnum        volume number, default: 1
%       flipped     boolean, flip along BV's Y axis, default: true
%
% Output fields:
%
%       vmr         VMR object

% Version:  v1.1
% Build:    16020310
% Date:     Feb-03 2016, 10:51 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get reference
bc = xo.C;
nvol = bc.NrOfVolumes;
if nargin < 2 || ~isa(vnum, 'double') || numel(vnum) ~= 1 || ...
    isinf(vnum) || isnan(vnum) || vnum < 1 || vnum > nvol
    vnum = 1;
else
    vnum = fix(real(vnum));
end
if nargin < 3 || numel(flipped) ~= 1 || (~isnumeric(flipped) && ~islogical(flipped)) || flipped
    flipped = true;
else
    flipped = false;
end

% slices not loaded?
if isempty(bc.Slice) || ~isstruct(bc.Slice) || ~isfield(bc.Slice, 'STCData')
    try
        fmr_LoadSTC(xo);
        bc = xo.C;
    catch xfferror
        error('neuroelf:xff:internalError', 'Error loading slices: ''%s''.', xfferror.message);
    end
end
try
    slarr = aft_GetVolume(xo, vnum);
catch xfferror
    error('neuroelf:xff:internalError', 'Error getting volume %d: %s.', volnum, xfferror.message);
end

% restrict to 0-225 range
mxs = double(max(slarr(:)));
ms = 225 / mxs;
for sc = 1:vdz
    slarr = uint8(round(ms * double(sl)));
end

% flip ?
if flipped
    slarr = slarr(:, :, vdz:-1:1);
end

% change type
slarr = shiftdim(slarr, 1);

% create VMR
vmr = xff('new:vmr');
vmrc = vmr.C;
vmrc.DimZ = vdn(1);
vmrc.DimX = vdn(2);
vmrc.DimY = vdz;
vmrc.VMRData = slarr;
vmrc.OffsetZ = fix(128 - vdn(1) / 2);
vmrc.OffsetX = fix(128 - vdn(2) / 2);
vmrc.OffsetY = fix(128 - vdz / 2);
vmrc.FramingCube = 256;
vmrc.Slice1CenterZ = -vdn(1) / 2;
vmrc.SliceNCenterZ = vdn(1) / 2;
vmrc.NRows = vdn(2);
vmrc.NCols = vdz;
vmrc.FoVRows = bc.FoVCols;
vmrc.FoVCols = vdz * (bc.SliceThickness + bc.SliceGap);
vmrc.SliceThickness = bc.InplaneResolutionX;
vmrc.GapThickness = 0;
vmrc.VoxResX = bc.InplaneResolutionY;
vmrc.VoxResY = bc.SliceThickness + bc.SliceGap;
vmrc.VoxResZ = bc.InplaneResolutionX;
vmrc.MaxOriginalValue = mxs;

% put back in global array
vmr.C = vmrc;
