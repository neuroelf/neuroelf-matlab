function xo = vtc_Reframe(xo, bbox)
% VTC::Reframe  - reframe a VTC
%
% FORMAT:       [vtc] = vtc.Reframe(bbox)
%
% Input fields:
%
%       bbox        2x3 bounding box

% Version:  v1.1
% Build:    16021413
% Date:     Feb-14 2016, 1:05 PM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc') || ...
   ~isa(bbox, 'double') || length(size(bbox)) ~= 2 || numel(bbox) ~= 6 || ...
    any(isinf(bbox(:)) | isnan(bbox(:)))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get object reference
bc = xo.C;

% check bbox size
if size(bbox, 2) == 3
    bbox = bbox';
end

% get current array size
oOffsetX = bc.XStart;
oOffsetY = bc.YStart;
oOffsetZ = bc.ZStart;
oDimX = size(bc.VTCData, 2);
oDimY = size(bc.VTCData, 3);
oDimZ = size(bc.VTCData, 4);
oRes = bc.Resolution;

% try to get new array size
nOffsetX = bbox(1);
nOffsetY = bbox(2);
nOffsetZ = bbox(3);
nDimX = round((bbox(4) - nOffsetX) / oRes);
nDimY = round((bbox(5) - nOffsetY) / oRes);
nDimZ = round((bbox(6) - nOffsetZ) / oRes);

% check new dims
if any([nDimX, nDimY, nDimZ] < 0)
    error('neuroelf:xff:badArgument', 'Invalid bounding box.');
end

nVTC = uint16(0);
try
    nVTC(size(bc.VTCData, 1), nDimX, nDimY, nDimZ) = nVTC(1);
catch xfferror
    error('neuroelf:xff:internalError', 'Invalid dimensions for reframing: %s.', xfferror.message);
end

% decide on source/target indexing...
sfx = max(1, round(1 + (nOffsetX - oOffsetX) / oRes));
tfx = max(1, round(1 + (oOffsetX - nOffsetX) / oRes));
stx = min(oDimX, nDimX + round((nOffsetX - oOffsetX) / oRes));
ttx = tfx + stx - sfx;
sfy = max(1, round(1 + (nOffsetY - oOffsetY) / oRes));
tfy = max(1, round(1 + (oOffsetY - nOffsetY) / oRes));
sty = min(oDimY, nDimY + round((nOffsetY - oOffsetY) / oRes));
tty = tfy + sty - sfy;
sfz = max(1, round(1 + (nOffsetZ - oOffsetZ) / oRes));
tfz = max(1, round(1 + (oOffsetZ - nOffsetZ) / oRes));
stz = min(oDimZ, nDimZ + round((nOffsetZ - oOffsetZ) / oRes));
ttz = tfz + stz - sfz;

% try data copy
try
    nVTC(1:end, tfx:ttx, tfy:tty, tfz:ttz) = bc.VTCData(1:end, sfx:stx, sfy:sty, sfz:stz);
catch xfferror
    error('neuroelf:xff:internalError', ...
        'Error calculating copy from/to indices: %s.', xfferror.message);
end

% overwrite VMRData in object
bc.VTCData = nVTC;

% set fields
bc.XStart = nOffsetX;
bc.YStart = nOffsetY;
bc.ZStart = nOffsetZ;
bc.XEnd = nOffsetX + oRes * size(nVTC, 2);
bc.YEnd = nOffsetY + oRes * size(nVTC, 3);
bc.ZEnd = nOffsetZ + oRes * size(nVTC, 4);

% put back into storage
xo.C = bc;
