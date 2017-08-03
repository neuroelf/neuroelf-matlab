function xo = vmr_OverlaySRF(xo, srf, ocolor)
% VMR::OverlaySRF  - overlay an SRF into the VMR (data changed!)
%
% FORMAT:       [vmr] = vmr.OverlaySRF(srf [, ocolor])
%
% Input fields:
%
%       srf         SRF object
%       ocolor      optional overlay color (default: 10, max positive LUT)
%
% Output fields
%
%       vmr         altered (!) object
%
% Note: since the overlay is done with the statistical LUT colors
%       the VMR object is changed !!! So, it does not work on V16 objects

% Version:  v1.1
% Build:    16021316
% Date:     Feb-13 2016, 4:50 PM EST
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

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr') || ...
    numel(srf) ~= 1 || ~xffisobject(srf, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
srfbc = srf.C;
if ~bc.VMR8bit
    error('neuroelf:xff:invalidObject', 'Only valid for 8-bit VMR objects.');
end
if nargin < 3 || ~isnumeric(ocolor) || numel(ocolor) ~= 1
    ocolor = 10;
else
    ocolor = round(real(double(ocolor)));
end
if isinf(ocolor) || isnan(ocolor) || ocolor < 0 || ocolor > 30
    ocolor = 10;
end
ocolor = uint8(ocolor + 225);

% get good coordinates
c = round(srfbc.VertexCoordinate);
vs = size(bc.VMRData);
c((c(:, 1) < 1 | c(:, 1) > vs(1) | c(:, 2) < 1 | c(:, 2) > vs(2) | c(:, 3) < 1 | c(:, 3) > vs(3)), :) = [];
c = sub2ind(vs, c(:, 1), c(:, 2), c(:, 3));

% set in VMRData
bc.VMRData(c) = ocolor;
xo.C = bc;
