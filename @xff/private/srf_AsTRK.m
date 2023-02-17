function trk = srf_AsTRK(xo)
% SRF::AsTRK  - create a TRK object from a surface
%
% FORMAT:       srf.AsTRK();
%
% No input fields.

% Version:  v1.1
% Build:    23021013
% Date:     Feb-10 2023, 1:47 PM EST
% Author:   Jochen Weber, NeuroElf.net
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2023, Jochen Weber
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

% check input arguments
if nargin < 1 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if ~isfield(bc.RunTimeVars, 'FiberStarts')
    error('neuroelf:xff:badArgument', 'Convertion to TRK requires Fiber-SRF.');
end
rtv = bc.RunTimeVars;
fs = rtv.FiberStarts;
vc = bc.VertexCoordinate(:, [3, 1, 2]);
vc(:, 3) = 256 - vc(:, 3);
vnn = ~isnan(vc(:, 1));
mn = floor(min(vc(vnn, :)));
vc = vc - repmat(mn - 1, size(vc, 1), 1);
mx = ceil(max(vc(vnn, :)));

% convert
trk = xff('new:trk');
rtv.xffID = trk.C.RunTimeVars.xffID;
rtv.SourceSRFFilename = xo.F;
trk.C.RunTimeVars = rtv;
trk.C.ImageVolumeDims = 2 .* ceil(0.5 .* (mx + 1));
trk.C.NrOfTracts = numel(fs) - 1;
t = repmat(trk.C.Tracts(1), trk.C.NrOfTracts, 1);
for tc = 1:numel(t)
    t(tc).Points = vc(fs(tc):fs(tc+1)-2, :);
    t(tc).NrOfPoints = size(t(tc).Points, 1);
    t(tc).Values = zeros(t(tc).NrOfPoints, 0);
end
trk.C.Tracts = t;
