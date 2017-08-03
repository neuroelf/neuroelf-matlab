function vmr = vtc_CreateFuncVMR(xo, vol, iptype)
% VTC::CreateFuncVMR  - create a (pseudo) VMR
%
% FORMAT:       vmr = vtc.CreateFuncVMR([vol, iptype])
%
% Input fields:
%
%       vol         volume number (default 1)
%       iptype      interpolation 'cubic', 'linear', {'nearest'}
%
% Output fields:
%
%       vmr         VMR object (in 256x256x256 frame)
%
% Using: flexinterpn_method.

% Version:  v1.1
% Build:    16021321
% Date:     Feb-13 2016, 9:00 PM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
szvtc = size(bc.VTCData);
if nargin < 2 || ~isa(vol, 'double') || numel(vol) ~= 1 || isnan(vol) || vol < 0 || vol > szvtc(1) || vol ~= fix(vol)
    vol = 1;
end
if nargin < 3 || ~ischar(iptype) || ~any(strcmpi(iptype(:)', {'cubic', 'linear', 'nearest'}))
    iptype = 'nearest';
else
    iptype = lower(iptype(:)');
end

% create VMR
vmr = xff('new:vmr');
vmrc = vmr.C;

% get data
vd = bc.VTCData(:, :, :, :);
if vol > 0
    vd = squeeze(vd(vol, :, :, :));
else
    vd = squeeze(mean(vd));
end
vr = bc.Resolution;
iv = 1 / vr;
is = 1 - iv;
ixyz = [Inf, Inf, Inf; is, is, is; iv, iv, iv; szvtc(2:4) + eps + iv];

% what interpolation
vd = ne_methods.flexinterpn_method(vd, ixyz, 0, iptype);
vds = size(vd) - 1;

% put into V16 data at offset
vmrc.VMRData16 = uint16(vmrc.VMRData);
vmrc.VMRData16(bc.XStart:(bc.XStart + vds(1)), bc.YStart:(bc.YStart + vds(2)), ...
    bc.ZStart:(bc.ZStart + vds(3))) = uint16(vd);

% put into object
vmr.C = vmrc;

% thresholding
vmr_LimitVMR(vmr, struct('recalc8b', true));
