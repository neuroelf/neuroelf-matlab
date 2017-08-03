function xo2 = vmr_Talairach(xo, imeth, tal, acpc, inverse)
% VMR::Talairach  - Un/Talairachize VMR
%
% FORMAT:       talvmr = vmr.Talairach(imeth, tal [, acpc, inverse]);
%
% Input fields:
%
%       imeth   interpolation method ('linear', 'cubic', 'lanczos3')
%       tal     TAL object
%       acpc    optional ACPC TRF file to go from/back to native space
%       inverse perform inverse (un-Tal, un-ACPC) operation
%
% Output fields:
%
%       tvmr    (un-) Talairachized VMR
%
% Note: only works on VMRs in 1mm or 0.5mm ISOvoxel resolution
%
% Using: acpc2tal, flexinterpn_method.

% TODO: piecewise interpolation !!

% Version:  v1.1
% Build:    16051017
% Date:     May-10 2016, 5:01 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
acpc2tal           = ne_methods.acpc2tal;
flexinterpn_method = ne_methods.flexinterpn_method;

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr') || ...
   ~ischar(imeth) || isempty(imeth) || ~any(strcmpi(imeth(:)', {'cubic', 'lanczos3', 'linear'})) || ...
    numel(tal) ~= 1 || ~xffisobject(tal, true, 'tal')
    error('neuroelf:xff:badArgument', 'Bad or missing argument.');
end
bc = xo.C;
if (bc.FileVersion > 1 && (any([bc.VoxResX, bc.VoxResY, bc.VoxResZ] ~= 1 & ...
    [bc.VoxResX, bc.VoxResY, bc.VoxResZ] ~= 0.5) || ...
   ~all([bc.VoxResY, bc.VoxResZ] == bc.VoxResX)))
    error('neuroelf:xff:invaldObject', 'Method only valid for 1mm or 0.5mm ISOvoxel VMRs.');
end
doovinv = false;
if nargin == 4 && islogical(acpc) && ~isempty(acpc)
    ovinv = acpc(1);
    doovinv = true;
end
if nargin < 4 || numel(acpc) ~= 1 || ~xffisobject(acpc, true, 'trf')
    acpc = [];
else
    acpcbc = acpc.C;
    if ~strcmpi(acpcbc.DataFormat, 'Matrix') || acpcbc.TransformationType ~= 2
        acpc = [];
    else
        acpc = acpcbc.TFMatrix;
    end
end
if nargin < 5 || ~islogical(inverse) || isempty(inverse)
    inverse = false;
else
    inverse = inverse(1);
end
if doovinv
    inverse = ovinv;
end
if inverse && ~isempty(acpc)
    try
        acpc = inv(acpc)';
    catch xfferror
        neuroelf_lasterr(xfferror);
        error('neuroelf:xff:badArgument', 'Couldn''t invert ACPC transformation matrix.');
    end
elseif ~isempty(acpc)
    acpc = acpc';
end
tvmrd = bc.VMRData(1);
tvmrd(1) = 0;

% create new VMR
xo2 = xff('new:vmr');
co2 = xo2.C;
co2.VMR8bit = bc.VMR8bit;
if ~inverse
    co2.VoxResInTalairach = true;
end
co2.VoxResVerified = bc.VoxResVerified;

% grids depend on
res = bc.VoxResX;
hires = false;
if res == 1
    cs = 256;
else
    cs = 512;
    co2.DimX = 512;
    co2.DimY = 512;
    co2.DimZ = 512;
    co2.FramingCube = 512;
    co2.VoxResX = 0.5;
    co2.VoxResY = 0.5;
    co2.VoxResZ = 0.5;
    if tal.C.FileVersion > 1 && tal.C.TALVMRFramingCube == 512
        hires = true;
    end
end
[xgrd, ygrd] = ndgrid(0:res:255.5, 0:res:255.5);
tc = [xgrd(:), ygrd(:), zeros(numel(xgrd), 1)];
zval = 0:res:255.5;
zc = [zeros(numel(zval), 2), zval(:)];

% the simple way only works for forward or single transformation !
if ~inverse || isempty(acpc)
    tc(:, [3, 1, 2]) = acpc2tal(tc(:, [3, 1, 2]), tal, ~inverse, hires);
    zc(:, [3, 1, 2]) = acpc2tal(zc(:, [3, 1, 2]), tal, ~inverse, hires);
end
if ~isempty(acpc)
    tc = tc - 127.5;
    zc = zc - 127.5;
end
zval = zc(:, 3)';
tc(:, 4) = 1;

% determine progress bar capabilities
try
    if ~inverse
        if isempty(acpc)
            step = 'Talairach transforming';
        else
            step = 'ACPC / Talairach transforming';
        end
    else
        if isempty(acpc)
            step = 'Un-Talairach transforming';
        else
            step = 'Un-Tal / Un-ACPC transforming';
        end
    end

    pbar = xprogress;
    xprogress(pbar, 'setposition', [80, 200, 640, 36]);
    xprogress(pbar, 'settitle', sprintf('%s VMR...', step));
    xprogress(pbar, 0, 'Setting up target VMR...', 'visible', 0, 1);
    nxp = now + 0.5 / 86400;
catch xfferror
    neuroelf_lasterr(xfferror);
    pbar = [];
end

% get sampling data
try
    if istransio(bc.VMRData)
        sdt = bc.VMRData(:, :, :);
    else
        sdt = bc.VMRData;
    end
catch xfferror
    delete(xo2);
    rethrow(xfferror);
end

% create output matrix
co2.VMRData = tvmrd;
tvmrd(1:(cs * cs), 1:cs) = tvmrd;

% get sampling offset and size
ssz = size(sdt) + 1;
try
    sof = [bc.OffsetX, bc.OffsetY, bc.OffsetZ];
catch xfferror
    neuroelf_lasterr(xfferror);
    sof = [0, 0, 0];
end

% scaling
if res == 0.5
    acpc(1:3, 1:3) = acpc(1:3, 1:3) .* 2;
end

% iterate over z slices
if ~isempty(pbar)
    xprogress(pbar, 0, 'Interpolating...', 'visible', 0, numel(zval));
end
for zc = 1:numel(zval)

    % apply transformation
    tc(:, 3) = zval(zc);
    if ~isempty(acpc)
        sc = (tc * acpc) + (127.5 / res);
    else
        sc = tc;
    end

    % apply UnTal after ACPC
    if inverse && ~isempty(acpc)
        sc(:, [3, 1, 2]) = acpc2tal(sc(:, [3, 1, 2]), tal, false, hires);
    end

    % MATLAB offset
    sc = sc + 1;

    % offset
    if any(sof)
        sc = [sc(:, 1) - sof(1), sc(:, 2) - sof(2), sc(:, 3) - sof(3)];
    end

    % good samples
    p = (sc(:, 1) >= 0 & sc(:, 1) <= ssz(1) & sc(:, 2) >= 0 & sc(:, 2) <= ssz(2) & sc(:, 3) >= 0 & sc(:, 3) <= ssz(3));

    % sample VMR
    if any(p)

        % sample whole volume
        tvmrd(p, zc) = round(flexinterpn_method(sdt, sc(p, 1:3), 0, imeth));
    end

    % progress bar
    if ~isempty(pbar) && now > nxp
        xprogress(pbar, zc);
        nxp = now + 0.5 / 86400;
    end
end

% make sure to limit 8bit VMR
if strcmpi(class(tvmrd), 'uint8')
    tvmrd = min(tvmrd, uint8(225));
end

% set some more VMR fields
tvmrd = reshape(tvmrd, [cs, cs, cs]);
co2.VMRData = tvmrd;
co2.VoxResInTalairach = double(~inverse);
xo2.C = co2;

% make sure that V16 files are FileVersion 1
if ~co2.VMR8bit || strcmpi(class(co2.VMRData), 'uint16')
    vmr_Update(xo2, 'FileVersion', [], 1);
    xo2.C.VMR8bit = false;
end

% progress bar
if ~isempty(pbar)
    closebar(pbar);
end
