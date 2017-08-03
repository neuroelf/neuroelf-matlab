function vmr = mergehemivmrs(h1, h2)
% mergehemivmrs  - merge segmented VMRs of LH and RH
%
% FORMAT:       vmr = mergehemivmrs(h1, h2)
%
% Input fields:
%
%       h1, h2      segmented hemisphere VMR
%
% Output fields:
%
%       vmr         combined VMR

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:28 AM EST
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

% copy object
vmr = h1.CopyObject;

% get voxels to copy
copyb = conv3d(vmr.VMRData == 0, 1);

% copy data
vmr.VMRData(copyb) = h2.VMRData(copyb);

% find hemisphere cut
po = squeeze(sum(sum(vmr.VMRData(:, :, 113:144) >= 235, 1), 2));
[pm{1:2}] = min(po);
pm = pm{2} + 112;

% get adjacent slices
h1s = vmr.VMRData(:, :, pm - 1);
h2s = vmr.VMRData(:, :, pm + 1);

% get biggest "chunk"
[cc{1:2}] = clustercoords(cat(3, h1s, h2s) > 0, 'edge');
if isempty(cc{2})
    error( ...
        'neuroelf:InternalError', ...
        'Error generating cluster list in connected hemisphere.' ...
    );
end
cc = cc{2};
cn = 0;
ci = 0;
for c = 1:numel(cc)
    if cn < size(cc{c}, 1)
        ci = c;
        cn = size(cc{c}, 1);
    end
end

% reconnect by marking overlap
cc = cc{ci};
s1i = cc(:, 3) == 1;
mi = intersect( ...
    sub2ind(size(h1s), cc(s1i, 1), cc(s1i, 2)), ...
    sub2ind(size(h1s), cc(~s1i, 1), cc(~s1i, 2))) + ...
    (pm - 1) * numel(h1s);
vmr.VMRData(mi) = 240;

% prepare again
vmr.PrepareForReco;
