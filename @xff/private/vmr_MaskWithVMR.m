function xo = vmr_MaskWithVMR(xo, maskvmr, thold)
% VMR::MaskWithVMR  - set voxels 0 where mask VMR below threshold
%
% FORMAT:       [vmr] = vmr.MaskWithVMR(maskvmr [, threshold])
%
% Input fields:
%
%       maskvmr     masking VMR object
%       threshold   threshold value, default: 11
%
% Output fields:
%
%       vmr         masked VMR

% Version:  v1.1
% Build:    16021315
% Date:     Feb-13 2016, 3:41 PM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr') || ...
    numel(maskvmr) ~= 1 || ~xffisobject(maskvmr, true, 'vmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc1 = xo.C;
bc2 = maskvmr.C;
if any(size(bc1.VMRData) ~= size(bc2.VMRData))
    error('neuroelf:xff:invalidObject', 'Dimension mismatch.');
end
if nargin < 3 || ~isnumeric(thold) || numel(thold) ~= 1 || isnan(thold)
    thold = 11;
else
    if thold < 0
        thold = 0;
    elseif bc2.VMR8bit && thold > 225
        thold = 226;
    elseif thold > 32767
        thold = 32768;
    end
end
thold = fix(thold);

% resolve transio first
if istransio(bc1.VMRData)
    bc1.VMRData = bc1.VMRData(:, :, :);
end
if istransio(bc1.VMRData16)
    bc1.VMRData16 = bc1.VMRData16(:, :, :);
end

% apply mask
msk = (bc2.VMRData < thold);
bc1.VMRData(msk) = 0;
if numel(bc1.VMRData16) == numel(msk)
    bc1.VMRData16(msk) = 0;
end
xo.C = bc1;
