function xo = vmr_Reframe(xo, bbox)
% VMR::Reframe  - reframe the VMR
%
% FORMAT:       vmr.Reframe(bbox)
%
% Input fields:
%
%       bbox        2x3 double bounding box (end included!)
%
% No output fields.

% Version:  v1.1
% Build:    16021216
% Date:     Feb-12 2016, 4:23 PM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr') || ...
   ~isa(bbox, 'double') || length(size(bbox)) ~= 2 || numel(bbox) ~= 6 || ...
    any(isinf(bbox(:)) | isnan(bbox(:)))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get contents
bc = xo.C;

% check bbox size
if size(bbox, 2) == 3
    bbox = bbox';
end

% get current array size
oOffsetX = bc.OffsetX;
oOffsetY = bc.OffsetY;
oOffsetZ = bc.OffsetZ;
oDimX = size(bc.VMRData, 1);
oDimY = size(bc.VMRData, 2);
oDimZ = size(bc.VMRData, 3);

% special case, autofind (min, max, fringe)
if all(bbox(:, 1) == bbox(:, 2))
    vd = bc.VMRData(:, :, :);
    vs = size(vd);
    nxc = [vs(2) * vs(3), 1];
    mnt = bbox(1);
    mxt = bbox(2);
    for x1 = 1:vs(1)
        if any(reshape(vd(x1, :, :), nxc) >= mnt) && any(reshape(vd(x1, :, :), nxc) <= mxt)
            break;
        end
    end
    for xe = vs(1):-1:x1
        if any(reshape(vd(xe, :, :), nxc) >= mnt) && any(reshape(vd(xe, :, :), nxc) <= mxt)
            break;
        end
    end
    sxc = x1:xe;
    nyc = [numel(sxc) * vs(3), 1];
    for y1 = 1:vs(2)
        if any(reshape(vd(sxc, y1, :), nyc) >= mnt) && any(reshape(vd(sxc, y1, :), nyc) <= mxt)
            break;
        end
    end
    for ye = vs(2):-1:y1
        if any(reshape(vd(sxc, ye, :), nyc) >= mnt) && any(reshape(vd(sxc, ye, :), nyc) <= mxt)
            break;
        end
    end
    syc = y1:ye;
    nzc = [numel(sxc) * numel(syc), 1];
    for z1 = 1:vs(3)
        if any(reshape(vd(sxc, syc, z1), nzc) >= mnt) && any(reshape(vd(sxc, syc, z1), nzc) <= mxt)
            break;
        end
    end
    for ze = vs(3):-1:z1
        if any(reshape(vd(sxc, syc, ze), nzc) >= mnt) && any(reshape(vd(sxc, syc, ze), nzc) <= mxt)
            break;
        end
    end
    x1 = max(1, x1 - bbox(3));
    xe = min(vs(1), xe + bbox(3));
    y1 = max(1, y1 - bbox(3));
    ye = min(vs(2), ye + bbox(3));
    z1 = max(1, z1 - bbox(3));
    ze = min(vs(3), ze + bbox(3));

    % reframe data
    vd = vd(x1:xe, y1:ye, z1:ze);
    if ~isempty(bc.VMRData16) && ndims(bc.VMRData16) == ndims(bc.VMRData) && ...
        all(size(bc.VMRData16) == size(bc.VMRData))
        bc.VMRData16 = bc.VMRData16(x1:xe, y1:ye, z1:ze);
    end
    bc.VMRData = vd;
    bc.OffsetX = oOffsetX + x1 - 1;
    bc.OffsetY = oOffsetY + y1 - 1;
    bc.OffsetZ = oOffsetZ + z1 - 1;
    bc.DimX = size(vd, 1);
    bc.DimY = size(vd, 2);
    bc.DimZ = size(vd, 3);
    xo.C = bc;
    return;
end

% try to get new array size
nOffsetX = bbox(1);
nOffsetY = bbox(2);
nOffsetZ = bbox(3);
nDimX = bbox(4) + 1 - nOffsetX;
nDimY = bbox(5) + 1 - nOffsetY;
nDimZ = bbox(6) + 1 - nOffsetZ;

% check new dims
if any([nDimX, nDimY, nDimZ] < 0)
    error('neuroelf:xff:badArgument', 'Invalid bounding box.');
end

% create target array
zVMR = bc.VMRData(1);
zVMR = zVMR - zVMR;
try
    nVMR = zVMR([]);
    nVMR(nDimX, nDimY, nDimZ) = zVMR;
catch xfferror
    error('neuroelf:xff:internalError', ...
        'Invalid dimensions for reframing: %s.', xfferror.message);
end

% decide on source/target indexing...
sfx = max(1, 1 + nOffsetX - oOffsetX);
tfx = max(1, 1 + oOffsetX - nOffsetX);
stx = min(oDimX, nDimX + nOffsetX - oOffsetX);
ttx = tfx + stx - sfx;
sfy = max(1, 1 + nOffsetY - oOffsetY);
tfy = max(1, 1 + oOffsetY - nOffsetY);
sty = min(oDimY, nDimY + nOffsetY - oOffsetY);
tty = tfy + sty - sfy;
sfz = max(1, 1 + nOffsetZ - oOffsetZ);
tfz = max(1, 1 + oOffsetZ - nOffsetZ);
stz = min(oDimZ, nDimZ + nOffsetZ - oOffsetZ);
ttz = tfz + stz - sfz;

% try data copy
try
    nVMR(tfx:ttx, tfy:tty, tfz:ttz) = bc.VMRData(sfx:stx, sfy:sty, sfz:stz);
catch xfferror
    error('neuroelf:xff:internalError', ...
        'Error calculating copy from/to indices: %s.', xfferror.message);
end

% also work on V16 ?
if bc.VMR8bit && ~isempty(bc.VMRData16) && ndims(bc.VMRData16) == ndims(bc.VMRData) && ...
    all(size(bc.VMRData16) == size(bc.VMRData))
    nVMR16 = uint16(0);
    nVMR16(1:nDimX, 1:nDimY, 1:nDimZ) = nVMR16(1);
    nVMR16(tfx:ttx, tfy:tty, tfz:ttz) = bc.VMRData16(sfx:stx, sfy:sty, sfz:stz);
    bc.VMRData16 = nVMR16;
end

% overwrite VMRData in object
bc.VMRData = nVMR;

% set fields
bc.DimX = size(nVMR, 1);
bc.DimY = size(nVMR, 2);
bc.DimZ = size(nVMR, 3);
bc.OffsetX = nOffsetX;
bc.OffsetY = nOffsetY;
bc.OffsetZ = nOffsetZ;

% set content
xo.C = bc;
