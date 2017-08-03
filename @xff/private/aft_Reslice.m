function rso2 = aft_Reslice(xo, rso, vol, imeth)
% AFT::Reslice  - reslice one object into the space of another
%
% FORMAT:       robj = obj.Reslice(rsobj, [vol, [, imeth]])
%
% Input fields:
%
%       rsobj       reslicing space object
%       vol         volume number (default: 1);
%       imeth       method, e.g. {'linear'}, 'cubic', 'lanczos3', 'nearest'
%
% Output fields:
%
%       robj        resliced object
%
% TYPES: HDR, HEAD, MSK, VMR
%
% Note: This method requires the MEX file flexinterpn.
%
% Using: bvcoordconv, flexinterpn, flexinterpn_method.

% Version:  v1.1
% Build:    17040615
% Date:     Apr-06 2017, 3:34 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015, 2016, 2017, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'hdr', 'head', 'msk', 'vmr'}) || ...
    numel(rso) ~= 1 || ~xffisobject(rso, true, {'hdr', 'head', 'msk', 'vmr'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end

% copy object and get type information
rso2 = aft_CopyObject(rso);
st = lower(xo.S.Extensions{1});
tst = lower(rso2.S.Extensions{1});

% volume number
if nargin < 3 || isempty(vol)
    vol = 1;
elseif numel(vol) ~= 1 || ~isa(vol, 'double') || isinf(vol) || isnan(vol) || vol < 1 || vol ~= fix(vol)
    error('neuroelf:xff:badArgument', 'Invalid vol argument.');
end

% interpolation method
if nargin < 4 || ((~ischar(imeth) || ...
   ~any(strcmpi(imeth(:)', {'cubic', 'lanczos3', 'lanczos8', 'linear', 'nearest'}))) && ...
    (~isa(imeth, 'double') || numel(imeth) ~= size(imeth(1)) || mod(numel(imeth), 2) == 0 || ...
     any(isinf(imeth) | isnan(imeth))))
    imeth = 'linear';
elseif ischar(imeth)
    imeth = lower(imeth(:)');
else
    if numel(imeth) > 4096 && mod(numel(imeth), 4096) == 1
        ks = 4096;
    else
        ks = 1;
    end
end

% get data and construct sampling box
vdt = aft_GetVolume(xo, vol);
vds = size(vdt);
bb = [Inf, Inf, Inf; ones(2, 3); vds(1:3)];

% get required transformation
if strcmp(st, 'hdr')
    sh = hdr_CoordinateFrame(xo, vol);
    sh = sh.Trf;
elseif strcmp(st, 'head')
    sh = head_CoordinateFrame(xo, vol);
    sh = sh.Trf;
else
    sh = ne_methods.bvcoordconv(zeros(0, 3), 'bvc2tal', aft_BoundingBox(xo));
end
if strcmp(tst, 'hdr')
    th = hdr_CoordinateFrame(rso, vol);
    th = th.Trf;
elseif strcmp(tst, 'head')
    th = head_CoordinateFrame(rso, vol);
    th = th.Trf;
else
    th = ne_methods.bvcoordconv(zeros(0, 3), 'bvc2tal', aft_BoundingBox(rso));
end
tsh = sh \ th;

% interpolation
if ischar(imeth)
    y = ne_methods.flexinterpn_method(vdt(:, :, :, 1), bb, 0, tsh, imeth);
    for d5c = 2:size(vdt, 5)
        y(:, :, :, d5c) = ne_methods.flexinterpn_method(vdt(:, :, :, 1, d5c), bb, 0, tsh, imeth);
    end
else
    y = ne_methods.flexinterpn(vdt(:, :, :, 1), bb, imeth, ks, 0, tsh);
    for d5c = 2:size(vdt, 5)
        y(:, :, :, d5c) = ne_methods.flexinterpn(vdt(:, :, :, 1, d5c), bb, imeth, ks, 0, tsh);
    end
end

% store into output
switch (tst)
    case 'hdr'
        rso2.C.VoxelData = single(y);
        rso2.C.ImgDim.Dim(5) = size(y, 4);
        rso2.C.ImgDim.DataType = 16;
        rso2.C.ImgDim.BitsPerPixel = 32;
        rso2.C.ImgDim.ScalingSlope = 1;
        rso2.C.ImgDim.ScalingIntercept = 0;
    case 'head'
        rso2.C.Brick = sc.C.Brick(1);
        rso2.C.NrOfVolumes = 1;
        rso2.C.Brick.Data = single(y);
        rso2.C.Brick.ScalingFactor = 1;
    case 'msk'
        rso2.C.Mask = uint8(min(1, round(y(:, :, :, 1))));
    case 'vmr'
        rso2.C.VMRData = uint8(min(225, round(y(:, :, :, 1))));
end

% copy some RunTimeVars from initial object
rtv = xo.C.RunTimeVars;
if isfield(rtv, 'ScalingWindow') && isfield(rtv, 'ScalingWindowLim') && isfield(rtv, 'ScalingHist')
    rso2.C.RunTimeVars.ScalingWindow = rtv.ScalingWindow;
    rso2.C.RunTimeVars.ScalingWindowLim = rtv.ScalingWindowLim;
    rso2.C.RunTimeVars.ScalingHist = rtv.ScalingHist;
end
