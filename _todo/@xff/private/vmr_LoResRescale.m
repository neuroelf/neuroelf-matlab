function xo2 = vmr_LoResRescale(xo, cutoff)
% VMR::LoResRescale  - bring 0.5mm hires back to VMR space
%
% FORMAT:       lores = vmr.LoResRescale([cutoff])
%
% Input fields:
%
%       cutoff      optional flag, minimize VMR box (default false)
%
% Output fields:
%
%       lores       1x1x1 mm VMR

% Version:  v1.1
% Build:    16021315
% Date:     Feb-13 2016, 3:49 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
if ~isfield(bc, 'VoxResX') || any([bc.VoxResX, bc.VoxResY, bc.VoxResZ] ~= 0.5)
    error('neuroelf:xff:invalidObject', 'Method only valid for 0.5mm ISOvoxel VMRs.');
end
if nargin < 2 || ~islogical(cutoff) || isempty(cutoff)
    cutoff = false;
else
    cutoff = cutoff(1);
end

% what outbox
if cutoff
    obox = [bc.OffsetX, bc.OffsetY, bc.OffsetZ];
    obox = [floor(obox / 2); floor((obox + size(bc.VMRData) - 1) / 2)];
else
    obox = [0, 0, 0; 255, 255, 255];
end

% make temporary copy
to = [];
try
    to = aft_CopyObject(xo);
    vmr_Reframe(to, [obox(1, :) .* 2; obox(2, :) .* 2 + 1]);
    toc = to.C;
    delete(to);
catch xfferror
    neuroelf_lasterr(xfferror);
    if ~isempty(to)
        delete(to);
    end
    error('neuroelf:xff:outOfMemory', 'Error creating temporary sampling object.');
end

% create output VMR
xo2 = xff('new:vmr');
vmr_Reframe(xo2, obox);
bc2 = xo2.C;
nd = uint16([]);
nd(1:(obox(2, 1) + 1), 1:(obox(2, 2) + 1), 1:(obox(2, 3) + 1)) = 0;

% sum elements
for x = 1:2
    for y = 1:2
        for z = 1:2
            nd = nd + uint16(toc.VMRData(x:2:end, y:2:end, z:2:end));
        end
    end
end
for zc = 1:size(nd, 3)
    nd(:, :, zc) = uint16(round(double(nd(:, :, zc)) ./ 8));
end
if strcmpi(class(toc.VMRData), 'uint8')
    nd = uint8(nd);
end
bc2.VMRData = nd;

% do the same for V16 data if present
if ~isempty(toc.VMRData16) && numel(size(toc.VMRData16)) == numel(size(toc.VMRData)) && ...
    all(size(toc.VMRData16) == size(toc.VMRData))
    nd = uint16([]);
    nd(1:(obox(2, 1) + 1), 1:(obox(2, 2) + 1), 1:(obox(2, 3) + 1)) = 0;
    for x = 1:2
        for y = 1:2
            for z = 1:2
                nd = nd + toc.VMRData16(x:2:end, y:2:end, z:2:end);
            end
        end
    end
    for zc = 1:size(nd, 3)
        nd(:, :, zc) = uint16(double(double(nd(:, :, zc)) ./ 8));
    end
    bc2.VMRData16 = nd;
end

% set offsets / dims
bc2.OffSetX = obox(1, 1);
bc2.OffSetY = obox(1, 2);
bc2.OffSetZ = obox(1, 3);
bc2.DimX = size(nd, 1);
bc2.DimY = size(nd, 2);
bc2.DimZ = size(nd, 3);
xo2.C = bc2;
