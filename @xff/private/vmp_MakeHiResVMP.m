function hvmp = vmp_MakeHiResVMP(xo, mapno)
% VMP::MakeHiResVMP  - convert one of the VMP maps to hi-res
%
% FORMAT:       hvmp = vmp.MakeHiResVMP([mapno])
%
% Input fields:
%
%       mapno       if given, only convert sub-selection
%
% Output fields:
%
%       hvmp        1x1x1 interpolated VMP
%
% Using: flexinterpn_method.

% Version:  v1.1
% Build:    16012316
% Date:     Jan-23 2016, 4:00 PM EST
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
flexinterpn_method = ne_methods.flexinterpn_method;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isa(mapno, 'double') || isempty(mapno) || ...
    any(isinf(mapno(:)) | isnan(mapno(:)) | mapno(:) < 1)
    mapno = 1:numel(bc.Map);
end
mapno = unique(round(mapno(:)));
mapno(mapno > numel(bc.Map)) = [];
if isempty(mapno)
    error('neuroelf:xff:badArgument', 'Invalid sub-selection of maps.');
end
numm = numel(mapno);

% get dims right
xs = bc.XStart;
xe = bc.XEnd;
ys = bc.YStart;
ye = bc.YEnd;
zs = bc.ZStart;
ze = bc.ZEnd;
br = bc.Resolution;
sr = 1 / br;
mapx = (xe - xs) / br;
mapy = (ye - ys) / br;
mapz = (ze - zs) / br;
dimx = xe - xs + 1;
dimy = ye - ys + 1;
dimz = ze - zs + 1;
dims = [dimx, dimy, dimz];

% special cases
if ~bc.NativeResolutionFile
    error('neuroelf:xff:invalidObject', 'This method is only available for NativeResolutionFiles.');
end
if br == 1
    hvmp = aft_CopyObject(xo);
    hbc = hvmp.C;
    for mc = 1:numel(hbc.Map)
        dl = (size(hbc.Map(mc).VMPData) > dims);
        ds = (size(hbc.Map(mc).VMPData) < dims);
        if dl(1)
            hbc.Map(mc).VMPData(dims(1)+1:end, :, :) = [];
        end
        if dl(2)
            hbc.Map(mc).VMPData(:, dims(1)+1:end, :) = [];
        end
        if dl(3)
            hbc.Map(mc).VMPData(:, :, dims(1)+1:end) = [];
        end
        if any(ds)
            hbc.Map(mc).VMPData(dims(1), dims(2), dims(3)) = 0;
        end
    end
    hvmp.C = hbc;
    return;
end

% resampling required, so create output object
hvmp = xff('new:vmp');
hbc = hvmp.C;
hbc.NativeResolutionFile = 0;
hbc.FileVersion = 4;
hbc.NrOfMaps = numm;
hbc.VMRDimX = 256;
hbc.VMRDimY = 256;
hbc.VMRDimZ = 256;
hbc.XStart = xs;
hbc.XEnd = xe;
hbc.YStart = ys;
hbc.YEnd = ye;
hbc.ZStart = zs;
hbc.ZEnd = ze;
hbc.Resolution = 1;

% create default map
hbc.Map.Type = 1;
hbc.Map.Name = '';
hbc.Map.VMPData = single([]);
hbc.Map.VMPData(1:dimx, 1:dimy, 1:dimz) = single(0);
hbc.Map(2:numm) = hbc.Map(1);

% get intepolation grid
cxyz = [Inf, Inf, Inf; 1, 1, 1; sr, sr, sr; mapx + 1, mapy + 1, mapz + 1];

% resample maps
hMap = hbc.Map;
oMap = bc.Map;
for cc = 1:numm

    % copy settings first
    mc = mapno(cc);
    hMap(cc).Type                 = oMap(mc).Type;
    hMap(cc).LowerThreshold       = oMap(mc).LowerThreshold;
    hMap(cc).UpperThreshold       = oMap(mc).UpperThreshold;
    hMap(cc).Name                 = oMap(mc).Name;
    hMap(cc).RGBLowerThreshPos    = oMap(mc).RGBLowerThreshPos;
    hMap(cc).RGBUpperThreshPos    = oMap(mc).RGBUpperThreshPos;
    hMap(cc).RGBLowerThreshNeg    = oMap(mc).RGBLowerThreshNeg;
    hMap(cc).RGBUpperThreshNeg    = oMap(mc).RGBUpperThreshNeg;
    hMap(cc).UseRGBColor          = oMap(mc).UseRGBColor;
    hMap(cc).TransColorFactor     = oMap(mc).TransColorFactor;
    hMap(cc).NrOfLags             = oMap(mc).NrOfLags;
    hMap(cc).MinLag               = oMap(mc).MinLag;
    hMap(cc).MaxLag               = oMap(mc).MaxLag;
    hMap(cc).CCOverlay            = oMap(mc).CCOverlay;
    hMap(cc).ClusterSize          = oMap(mc).ClusterSize * (br ^ 3);
    hMap(cc).EnableClusterCheck   = oMap(mc).EnableClusterCheck;
    hMap(cc).UseValuesAboveThresh = oMap(mc).UseValuesAboveThresh;
    hMap(cc).DF1                  = oMap(mc).DF1;
    hMap(cc).DF2                  = oMap(mc).DF2;
    hMap(cc).ShowPositiveNegativeFlag = oMap(mc).ShowPositiveNegativeFlag;
    hMap(cc).BonferroniValue      = oMap(mc).BonferroniValue;

    % interpolate content
    hMap(cc).VMPData(:, :, :) = ...
        flexinterpn_method(oMap(mc).VMPData(:, :, :), cxyz, 0, 'linear');
end

% store back
hbc.Map = hMap;
hvmp.C = hbc;
