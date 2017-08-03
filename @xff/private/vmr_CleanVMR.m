function xo2 = vmr_CleanVMR(xo, opts)
% VMR::CleanVMR  - try automatic cleaning (experimental !)
%
% FORMAT:       cleaned = vmr.CleanVMR([opts])
%
% Input fields:
%
%       opts        1x1 struct with optional settings
%        .framed    if true, assumes the VMR is already framed
%        .smkern    smoothing kernel (for border detection, 3mm)
%
% Output fields:
%
%       cleaned     new, cleaned VMR object
%
% Using: clustercoordsc, dilate3d, histcount, maxpos, minmaxmean,
%        smoothdata3.

% Version:  v1.1
% Build:    16021320
% Date:     Feb-13 2016, 8:25 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% argument check
if nargin < 1 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'framed') || ~islogical(opts.framed) || isempty(opts.framed)
    opts.framed = false;
else
    opts.framed = opts.framed(1);
end
if ~isfield(opts, 'smkern') || ~isa(opts.smkern, 'double') || numel(opts.smkern) ~= 1 || ...
    isinf(opts.smkern) || isnan(opts.smkern) || opts.smkern <= 0.5
    opts.smkern = 3;
else
    opts.smkern = min(8, opts.smkern);
end

% copy vmr and get double precision data
xo2 = aft_CopyObject(xo);
xo2.F = '';
if ~bc.VMR8bit || isempty(bc.VMRData16) || ...
    numel(size(bc.VMRData16)) ~= numel(size(bc.VMRData)) || ...
    any(size(bc.VMRData16) ~= size(bc.VMRData))
    cvmr = class(bc.VMRData);
    vmrd = double(bc.VMRData(:, :, :));
else
    cvmr = 'uint16';
    vmrd = double(bc.VMRData16(:, :, :));
end
isz = size(vmrd);

% if not framed, inspect exterior planes first
if ~opts.framed
    pxy1 = vmrd(:, :, 1);
    pxye = vmrd(:, :, end);
    pxz1 = squeeze(vmrd(:, 1, :));
    pxze = squeeze(vmrd(:, end, :));
    pyz1 = squeeze(vmrd(1, :, :));
    pyze = squeeze(vmrd(end, :, :));

    % build mean / std of each plane
    mxy1 = mean(pxy1(:));
    mxye = mean(pxye(:));
    mxz1 = mean(pxz1(:));
    mxze = mean(pxze(:));
    myz1 = mean(pyz1(:));
    myze = mean(pyze(:));
    sxy1 = std(pxy1(:));
    sxye = std(pxye(:));
    sxz1 = std(pxz1(:));
    sxze = std(pxze(:));
    syz1 = std(pyz1(:));
    syze = std(pyze(:));

    % find planes for which the mean is around the mean of means
    dmea = mean([mxy1, mxye, mxz1, mxze, myz1, myze]);
    mstd = 0.5 + mean([sxy1, sxye, sxz1, sxze, syz1, syze]) / 1.5;
    tv = [];
    if abs(mxy1 - dmea) < mstd
        tv = [tv; pxy1(:)];
    end
    if abs(mxye - dmea) < mstd
        tv = [tv; pxye(:)];
    end
    if abs(mxz1 - dmea) < mstd
        tv = [tv; pxz1(:)];
    end
    if abs(mxze - dmea) < mstd
        tv = [tv; pxze(:)];
    end
    if abs(myz1 - dmea) < mstd
        tv = [tv; pyz1(:)];
    end
    if abs(myze - dmea) < mstd
        tv = [tv; pyze(:)];
    end

    % any background found
    if ~isempty(tv)

        % set all voxels up to mean + 0.5 * std + 1 auf 0
        coff = mean(tv) + 0.5 * std(tv) + 1;
        vmrd(vmrd <= coff) = 0;
    else
        coff = 0.5;
    end

    % find remaining data
    nxc = [isz(2) * isz(3), 1];
    for x1 = 1:isz(1)
        if any(reshape(vmrd(x1, :, :), nxc) > 0)
            break;
        end
    end
    for xe = isz(1):-1:x1
        if any(reshape(vmrd(xe, :, :), nxc) > 0)
            break;
        end
    end
    sxc = x1:xe;
    nyc = [numel(sxc) * isz(3), 1];
    for y1 = 1:isz(2)
        if any(reshape(vmrd(sxc, y1, :), nyc) > 0)
            break;
        end
    end
    for ye = isz(2):-1:y1
        if any(reshape(vmrd(sxc, ye, :), nyc) > 0)
            break;
        end
    end
    syc = y1:ye;
    nzc = [numel(sxc) * numel(syc), 1];
    for z1 = 1:isz(3)
        if any(reshape(vmrd(sxc, syc, z1), nzc) > 0)
            break;
        end
    end
    for ze = isz(3):-1:z1
        if any(reshape(vmrd(sxc, syc, ze), nzc) > 0)
            break;
        end
    end
    ioff = [x1, y1, z1];
    isz = [xe, ye, ze];
    vmrd = vmrd(x1:xe, y1:ye, z1:ze);
else

    % set offsets to 1
    ioff = [1, 1, 1];
    coff = 0.5;
end

% smooth remainder
svmrd = ne_methods.smoothdata3(vmrd, opts.smkern([1, 1, 1]));

% compute minmaxmean
svmmm = ne_methods.minmaxmean(svmrd);

% build histogram (above coff)
mxi = 0;
histcount = ne_methods.histcount;
while mxi < 12
    hn = histcount(svmrd, coff, svmmm(2), 0.0101 * (svmmm(2) - coff));
    [mxn, mxi] = max(hn);
    if mxi < 12
        coff = coff + 0.125 * (svmmm(2) - coff);
    end
end

% if biggest cluster from center, only keep that
svmrd = (svmrd > coff);
[cs, cv] = ne_methods.clustercoordsc(svmrd);
svmrd = ne_methods.dilate3d(cv == ne_methods.maxpos(cs));
[cs, cv] = ne_methods.clustercoordsc(~svmrd);
svmrd = (cv == ne_methods.maxpos(cs));
vmrd(svmrd) = 0;

% back to old class
if strcmpi(cvmr, 'uint8')
    vmrd = uint8(vmrd);
else
    vmrd = uint16(vmrd);
end

% resolve transio first
if istransio(xo2.C.VMRData)
    xo2.C.VMRData = reshape(bc.VMRData(:), size(bc.VMRData));
end
if istransio(xo2.C.VMRData16)
    xo2.C.VMRData16 = reshape(bc.VMRData16(:), size(bc.VMRData16));
end

% where to put this
if ~xo2.C.VMR8bit || isa(vmrd, 'uint8')
    xo2.C.VMRData(:) = 0;
    xo2.C.VMRData(ioff(1):isz(1), ioff(2):isz(2), ioff(3):isz(3)) = vmrd;
else
    xo2.C.VMRData16(:) = 0;
    xo2.C.VMRData16(ioff(1):isz(1), ioff(2):isz(2), ioff(3):isz(3)) = vmrd;
    xo2.C.VMRData(xo2.C.VMRData16 == 0) = 0;
end
