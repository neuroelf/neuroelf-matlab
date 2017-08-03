function xo = vmr_OverlayVMP(xo, vmp, vmpnum, vmpi)
% VMR::OverlayVMP  - overlay a VMP into the VMR (data changed!)
%
% FORMAT:       [vmr] = vmr.OverlayVMP(vmp [, vmpnum, vmpi])
%
% Input fields:
%
%       vmp         VMP object with at least one Map
%       vmpnum      map selection (default: 1)
%       vmpi        interpolation: 'nearest', {'linear'}, 'cubic'
%
% Output fields
%
%       vmr         altered (!) object
%
% Note: since the overlay is done with the statistical LUT colors
%       the VMR object is changed !!! So, it does not work on V16 objects
%
% Using: flexinterpn_method.

% Version:  v1.1
% Build:    16021316
% Date:     Feb-13 2016, 4:43 PM EST
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

% neuroelf library
global ne_methods;

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr') || ...
    numel(vmp) ~= 1 || ~xffisobject(vmp, true, 'vmp')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
vmpc = vmp.C;
if ~bc.VMR8bit || isempty(vmpc.Map)
    error('neuroelf:xff:invalidObject', 'VMR must be 8-bit and VMP not empty.');
end

% map selection
if nargin < 3 || ~isa(vmpnum, 'double') || isempty(vmpnum) || ...
    any(isinf(vmpnum(:)) | isnan(vmpnum(:)) | vmpnum(:) < 1 | vmpnum(:) ~= fix(vmpnum(:)))
    vmpnum = 1;
else
    vmpnum = min(vmpnum(:), numel(vmpc.Map));
    if isempty(vmpnum)
        vmpnum = 1;
    end
    if numel(unique(vmpnum)) ~= numel(vmpnum)
        for mc = numel(vmpnum):-1:2
            if any(vmpnum(1:(mc - 1)) == vmpnum(mc))
                vmpnum(mc) = [];
            end
        end
    end
end
vmpnum = vmpnum(:)';

% interpolation
if nargin < 4 || ~ischar(vmpi) || isempty(vmpi) || ~any(strcmpi(vmpi(:)', {'cubic', 'linear', 'nearest'}))
    vmpi = 'linear';
else
    vmpi = lower(vmpi(:)');
end

% get selected maps (in reverse order)
map = vmpc.Map(vmpnum);

% get dimensions and coordinates right
if isfield(bc, 'OffsetX')
    vmro = [bc.OffsetX, bc.OffsetY, bc.OffsetZ];
else
    vmro = [0, 0, 0];
end
vmrs = size(bc.VMRData);
vmpo = [vmpc.XStart, vmpc.YStart, vmpc.ZStart] - vmro;
vmpr = vmpc.Resolution;
vmps = size(vmpc.Map(1).VMPData);

% interpolation needed at all?
vmpir = 1 / vmpr;
scxyz = [Inf, Inf, Inf; 1, 1, 1; vmpir, vmpir, vmpir; vmps + 0.95];
szx = numel(scxyz(2, 1):scxyz(3, 1):scxyz(4, 1));
szy = numel(scxyz(2, 2):scxyz(3, 2):scxyz(4, 2));
szz = numel(scxyz(2, 3):scxyz(3, 3):scxyz(4, 3));
tmap = uint8([]);
tmap(1:szx, 1:szy, 1:szz) = 0;

% iterate over maps
for mc = 1:numel(map)

    % get map values in required representation
    vmpd = ne_methods.flexinterpn_method(map(mc).VMPData, scxyz, 0, vmpi);

    % get indices of positive > thresh and set in tmap
    lt = map(mc).LowerThreshold;
    ut = map(mc).UpperThreshold;
    actv = (vmpd > lt);
    tmap(actv) = uint8(226 + floor(9.0001 * (min( vmpd(actv) - lt, ut - lt) / (ut - lt))));

    % get indices of negative < -thresh and set in tmap
    actv = (vmpd < -lt);
    tmap(actv) = uint8(236 + floor(9.0001 * (min(-vmpd(actv) - lt, ut - lt) / (ut - lt))));
end

% put tmap into the right position
clear vmpd;
xtmap = uint8([]);
xtmap(1:vmrs(1), 1:vmrs(2), 1:vmrs(3)) = 0;
scpo = vmro + max(-vmpo, 0) - min(vmpo, 0) + 1;
tcpo = max(vmpo, 1) + 1;
cpsz = min(size(tmap) + 1 - scpo , vmrs + 1 - tcpo) - 1;
xtmap(tcpo(1):(tcpo(1) + cpsz(1)), tcpo(2):(tcpo(2) + cpsz(2)), tcpo(3):(tcpo(3) + cpsz(3))) = ...
  tmap(scpo(1):(scpo(1) + cpsz(1)), scpo(2):(scpo(2) + cpsz(2)), scpo(3):(scpo(3) + cpsz(3)));

% set in VMRData
if ~isfield(bc, 'VMRDataBack')
    bc.VMRDataBack = bc.VMRData;
end
bc.VMRData = bc.VMRDataBack;
bc.VMRData(xtmap > 0) = xtmap(xtmap > 0);
xo.C = bc;
