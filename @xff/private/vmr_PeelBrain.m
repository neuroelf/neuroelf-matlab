function xo2 = vmr_PeelBrain(xo, thresh, nre, nrd)
% VMR::PeelBrain  - cheap try to peel the brain
%
% FORMAT:       peeled = vmr.PeelBrain([thresh, nre, nrd])
%
% Input fields:
%
%       thresh      matter inclusion threshold (default 0.005)
%       nre         number of erodes (default 5)
%       nrd         number of post-dilations (default nre + 1)
%
% Output fields:
%
%       peeled      peeled brain
%
% Note: this function is experimental !
%
% Using: dilate3d, erode3d, floodfill3, smoothdata3.

% Version:  v1.1
% Build:    16021320
% Date:     Feb-13 2016, 8:25 PM EST
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

% neuroelf library
global ne_methods;
dilate3d = ne_methods.dilate3d;
erode3d  = ne_methods.erode3d;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
if nargin < 3 || numel(thresh) ~= 1 || ~isa(thresh, 'double') || ...
    isnan(thresh) || thresh < 0 || thresh > 0.5
    thresh = 0.005;
end
if nargin < 3 || numel(nre) ~= 1 || ~isa(nre, 'double') || ...
    isnan(nre) || nre < 1 || nre > 11 || nre ~= fix(nre)
    nre = 5;
end
if nargin < 4 || numel(nrd) ~= 1 || ~isa(nrd, 'double') || ...
    isnan(nrd) || nrd < 1 || nrd > 12 || nrd ~= fix(nrd)
    nrd = nre + 1;
end

% which data to use
if ~bc.VMR8bit || isempty(bc.VMRData16) || ...
    numel(size(bc.VMRData16)) ~= numel(size(bc.VMRData)) || ...
    any(size(bc.VMRData16) ~= size(bc.VMRData))
    vmrd = bc.VMRData(:, :, :);
else
    vmrd = bc.VMRData16(:, :, :);
end

% find adaptive threshold
mnd = min(vmrd(:));
mxd = max(vmrd(:));
vxs = (vmrd > mnd);
nvox = sum(vxs(:));
hn = cumsum(hist(single(vmrd(vxs(:))), single((mnd + 1):mxd)));
thresh = find(hn >= (thresh * nvox));
if isempty(thresh)
    thresh = mnd;
end
if strcmpi(class(vmrd), 'uint8')
    thresh = mnd + uint8(thresh(1) - 1);
else
    thresh = mnd + uint16(thresh(1) - 1);
end
vdm = (vmrd > thresh);

% build erode3d/dilate3d arguments
r = sqrt(5);
try
    res = [bc.VoxResX, bc.VoxResY, bc.VoxResZ];
catch xfferror
    neuroelf_lasterr(xfferror);
    res = ones(1, 3);
end
opsz = 1 + 2 * round(sqrt(5) ./ res);

% erode
for ec = 1:nre
    vdm = erode3d(vdm, true(opsz), r, res);
end

% search in center for chunk coordinate
dochunk = false;
hsz = round(size(vdm) / 2);
for c = 0:ceil(min(hsz) / 3)
    tdm = vdm((hsz(1) - c):(hsz(1) + c), (hsz(2) - c):(hsz(2) + c), (hsz(3) - c):(hsz(3) + c));
    if any(tdm(:))
        dochunk = true;
        tdmi = find(tdm);
        [cx, cy, cz] = ind2sub(size(tdm), tdmi(1));
        tc = c + 1;
        cx = cx + hsz(1) - tc;
        cy = cy + hsz(2) - tc;
        cz = cz + hsz(3) - tc;
        break;
    end
end

% if found, try to find biggest chunk
if dochunk
    vdmc = ne_methods.floodfill3(vdm, cx, cy, cz);
    if sum(vdmc(:)) > (0.5 * sum(vdm))
        vdm = vdmc;
    end
end

% dilate back
for ec = 1:nrd
    vdm = dilate3d(vdm, true(opsz), r, res);
end

% mask with original mask
vdm = ~(vdm & vxs);

% generate new file
xo2 = aft_CopyObject(xo);
xo2.F = '';

% set data of non-mask to 0
vmrd(vdm) = 0;

% smooth and again remove anything below a threshold
vdm = find(~vdm);
vmrs = ne_methods.smoothdata3(double(vmrd), [2, 2, 2]);
tvmr = (vmrd(vdm) < ((1 / 2) * vmrs(vdm)));
vmrd(vdm(tvmr)) = 0;
if ~bc.VMR8bit || isa(vmrd, 'uint8')
    xo2.C.VMRData = vmrd;
else
    xo2.C.VMRData16 = vmrd;
    if istransio(xo2.C.VMRData)
        xo2.C.VMRData = xo2.C.VMRData(:, :, :);
    end
    xo2.C.VMRData(vdm) = 0;
end
