function vmr = spherevmr(c, sl, ic)
% spherevmr  - create a 256x256x256 VMR for spherical region growing
%
% FORMAT:       vmr = spherevmr([center [, slope, intercept])
%
% Input fields:
%
%       center      optional center coordinate (default: 128, 128, 128)
%       slope       value slope (default: 1)
%       intercept   value for center coordinate (default: 0)
%
% Output fields:
%
%       vmr         VMR

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:27 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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
   ~isa(c, 'double') || ...
    numel(c) ~= 3 || ...
    any(isnan(c) | c < 0 | c > 255)
    c = [128, 128, 128];
else
    c = c(:)';
end
if nargin < 2 || ...
   ~isa(sl, 'double') || ...
    numel(sl) ~= 1 || ...
    isnan(sl) || ...
    sl < -100 || ...
    sl > 100 || ...
    sl == 0
    sl = 1;
end
if nargin < 3 || ...
   ~isa(ic, 'double') || ...
    numel(ic) ~= 1 || ...
    isnan(ic) || ...
    ic < -200 || ...
    ic > 200
    ic = 0;
end

% prepare coordinates for slice
[cx, cz] = ndgrid(single(0:255), single(0:255));
cx = [cx(:) - c(2), cz(:) - c(3)];
cz = single(1:256) - c(1);
co = ones(65536, 1);

% create v16 data
v16d = single([]);
v16d(256, 256, 256) = 0;

% iterate over slices
for zc = 1:256
    cdist = sl * sqrt(sum([cx, cz(zc) * co] .^ 2, 2)) + ic;
    cdist(cdist < 0) = 0;
    v16d(:, :, zc) = reshape(cdist, [256, 256]);
end

% put into output
vmr = xff('new:vmr');
vmr.VMRData16 = uint16(round(v16d));

% recompute VMR
mn = min(v16d(:));
if mn > 100
    v16d = v16d - mn;
end
mx = max(v16d(:));
if mx > 225
    v16d = v16d ./ (mx / 225);
end
vmr.VMRData = uint8(round(v16d));
