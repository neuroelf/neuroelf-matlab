function vtc = hdr_Dyn3DToVTC(xo, vols)
% HDR::Dyn3DToVTC  - convert a dynamics 3D Analyze into VTC
%
% FORMAT:       [vtc] = hdr.Dyn3DToVTC([vols])
%
% Input fields:
%
%       vols        range of volumes in file (default: all)
%
% Output fields:
%
%       vtc         created VTC object (with slices loaded)

% Version:  v1.1
% Build:    16020515
% Date:     Feb-05 2016, 3:48 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2011, 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'hdr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
hdrc = xo.C;
hdrf = hdr_CoordinateFrame(xo);
hdrt = hdrf.Trf;
hdrfile = xo.F;
vsz = size(hdrc.VoxelData);
if numel(vsz) < 4
    vsz(end+1:4) = 1;
end
if nargin < 2 || ~isa(vols, 'double') || numel(vols) ~= max(size(vols)) || ...
    any(isinf(vols) | isnan(vols) | vols < 1 | vols > vsz(4) | vols ~= fix(vols))
    vols = 1:vsz(4);
else
    vols = unique(vols(:))';
end

% make some initial FMR project settings
vtc = xff('new:vtc');
vtcc = vtc.C;
vtcc.FileVersion = 3;
vtcc.NameOfSourceFMR = hdrfile;
vtcc.DataType = 2;
vtcc.NrOfVolumes = numel(vols);
vtcc.Resolution = 1;
vtcc.XStart = 0;
vtcc.XEnd = vsz(1);
vtcc.YStart = 0;
vtcc.YEnd = vsz(2);
vtcc.ZStart = 0;
vtcc.ZEnd = vsz(3);
vtcc.TR = 2000;
vtcc.VTCData = single(0);
vtcc.VTCData(numel(vols), vsz(1), vsz(2), vsz(3)) = 0;

% scaling
slope = hdrc.ImgDim.ScalingSlope;
offset = hdrc.ImgDim.ScalingIntercept;
if numel(slope) ~= 1 || isinf(slope) || isnan(slope) || slope == 0
    slope = 1;
end
if numel(offset) ~= 1 || isinf(offset) || isnan(offset)
    offset = 0;
end

% copy data
for vc = 1:numel(vols)
    if slope ~= 1 || offset ~= 0
        vtcc.VTCData(vc, :, :, :) = offset + slope .* double(reshape( ...
            hdrc.VoxelData(:, :, :, vols(vc)), [1, vsz(1:3)]));
    else
        vtcc.VTCData(vc, :, :, :) = reshape( ...
            hdrc.VoxelData(:, :, :, vols(vc)), [1, vsz(1:3)]);
    end
end

% add transformation information!
vtcc.RunTimeVars.AutoSave = true;
if isfield(hdrc.RunTimeVars, 'Discard')
    vtcc.RunTimeVars.Discard = hdrc.RunTimeVars.Discard;
end
vtcc.RunTimeVars.TrfPlus = hdrt * ...
    [ 0, -1,  0,  129; ...
      0,  0, -1,  129; ...
     -1,  0,  0,  129; ...
      0,  0,  0,    1];

% set content
vtc.C = vtcc;
