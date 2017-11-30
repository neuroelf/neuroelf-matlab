function vmr = voi_CreateVMR(xo, res, vspec, vmrval)
% VOI::CreateVMR  - create a VMR object of VOI voxels
%
% FORMAT:       vmr = voi.CreateVMR([res [, vspec [, vmrval]]])
%
% Input fields:
%
%       res         presumed VOI resolution (default: auto-detect)
%       vspec       voi index specification (list of indices, default: all)
%       vmrval      VMR voxel value (default: 240)
%
% Output fields:
%
%       vmr         VMR object with VOI voxels set to vmrval
%
% Using: bvcoordconv.

% Version:  v1.1
% Build:    17113010
% Date:     Nov-30 2017, 10:14 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2017, Jochen Weber
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
bvcoordconv = ne_methods.bvcoordconv;

% argument check
if nargin < 1 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'voi')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isa(res, 'double') || numel(res) ~= 1 || isinf(res) || isnan(res) || res < 1
    res = [];
else
    res = round(res);
end
if nargin < 3 || ~isa(vspec, 'double') || isempty(vspec) || ...
    any(isinf(vspec(:)) | isnan(vspec(:)) | vspec(:) < 1 | vspec(:) ~= fix(vspec(:)))
    vspec = 1:numel(bc.VOI);
else
    vspec = unique(vspec(:))';
end
if any(vspec > numel(bc.VOI))
    error('neuroelf:xff:badArgument', 'Selected VOI(s) out of bounds.');
end
if nargin < 4 || ~isa(vmrval, 'double') || numel(vmrval) ~= 1 || ...
    isinf(vmrval) || isnan(vmrval) || vmrval < 0 || vmrval > 255
    vmrval = 240;
else
    vmrval = round(vmrval);
end
vmrval = uint8(vmrval);

% coordinate space
natspc = ~strcmpi(bc.ReferenceSpace, 'tal');
if natspc
    voxres = [bc.OriginalVMRResolutionX, bc.OriginalVMRResolutionY, bc.OriginalVMRResolutionZ];
    voxoff = [bc.OriginalVMROffsetX, bc.OriginalVMROffsetY, bc.OriginalVMROffsetZ];
    convflag = 'bvi2bvx';
else
    convflag = 'tal2bvx';
end

% create VMR
vmr = xff('new:vmr');
vmrc = vmr.C;
bbo = aft_BoundingBox(vmr);

% iterate over selected VOIs
for vc = 1:numel(vspec)

    % get voxels
    vox = bc.VOI(vspec(vc)).Voxels;
    
    % resolution expansion
    if isempty(res)
        xd = min(diff(unique(vox(:, 1))));
        yd = min(diff(unique(vox(:, 2))));
        zd = min(diff(unique(vox(:, 3))));
        md = min([xd, yd, zd]);
    elseif res > 1
        md = res;
    else
        md = 1;
    end
    if md > 1
        cv = 0.25 + 0.5 * md;
        fv = round(cv - 0.5 * md);
        tv = round(cv + 0.5 * md) - 1;
        rx = md * md * md;
        nvox = size(vox, 1);
        xvox = zeros(rx * nvox, 3);
        xi = 1;
        for xd = fv:tv
            for yd = fv:tv
                for zd = fv:tv
                    xvox(xi:xi+nvox-1, :) = vox + ones(nvox, 1) * [xd, yd, zd];
                    xi = xi + nvox;
                end
            end
        end
        vox = xvox;
    end

    % native space
    if natspc

        % adapt voxels first
        nvox = ones(size(vox, 1), 1);
        vox = vox .* voxres(nvox, :) + voxoff(nvox, :);
    end

    % convert
    vox = bvcoordconv(vox, convflag, bbo);

    % remove out-of-box indices (now NaN!)
    vox(isnan(vox)) = [];

    % fill coordinates with 1's
    vmrc.VMRData(vox) = vmrval;
end

% set content
vmr.C = vmrc;
