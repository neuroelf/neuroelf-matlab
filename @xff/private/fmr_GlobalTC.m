function gtc = fmr_GlobalTC(xo)
% FMR::GlobalTC  - extract global time course
%
% FORMAT:       gtc = fmr.GlobalTC;
%
% No input fields.
%
% Output fields:
%
%       gtc         structure with global time courses
%        .allvoxslc one per slice (TxS)
%        .allvoxtot all voxels (Tx1)
%        .bgrndmsk  estimated background voxels mask (XxYxS)
%        .bgrndslc  only background voxels (TxS)
%        .bgrndtot  only background voxels (Tx1)
%        .fgrndmsk  estimated foreground voxels mask (XxYxS)
%        .fgrndslc  only foreground voxels (TxS)
%        .fgrndtot  only foreground voxels (Tx1)
%        .highimsk  estimated high intensity voxels mask (XxYxS)
%        .highislc  high intensity voxels, probably "gray matter" (TxS)
%        .highitot  high intensity voxels, probably "gray matter" (Tx1)
%
% Using: clustercoordsc, erode3d, floodfill3c, lsqueeze, maxpos, minmaxmean.

% Version:  v1.1
% Build:    16021413
% Date:     Feb-14 2016, 1:23 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get data
bc = xo.C;

% depending on file version
if bc.FileVersion > 4 && numel(bc.Slice) == 1

    % do work
    stcd = bc.Slice.STCData;
    if istransio(stcd)
        stcd = resolve(stcd);
    end

% otherwise pack and unpack
else

    stcd = bc.Slice(1).STCData;
    if istransio(stcd)
        stcd = resolve(stcd);
    end
    stcd(1, 1, 1, numel(bc.Slice)) = 0;
    for sc = 2:numel(bc.Slice)
        stcd(:, :, :, sc) = bc.Slice(sc).STCData(:, :, :);
    end
end

% permute for easier access
stcd = permute(stcd, [3, 1, 2, 4]);
nvol = size(stcd, 1);

% produce mean over time
mstcd = (1 / size(stcd, 1)) * squeeze(sum(stcd));
sstcd = size(mstcd);
mmm = ne_methods.minmaxmean(mstcd, 5);

% produce all voxel time courses
allvoxslc = (1 / prod(sstcd(1:2))) * squeeze(sum(sum(stcd, 2), 3));
allvox = (1 / sstcd(3)) * sum(allvoxslc, 2);

% try to get background estimate
bg = ne_methods.floodfill3c(mstcd < mmm(3), [1, 1, 1], 1) | ...
     ne_methods.floodfill3c(mstcd < mmm(3), sstcd, 1);
bg(mstcd == 0) = false;
bgm = ne_methods.minmaxmean(mstcd(bg), 5);
bgh = bgm(3) + 1.5 * sqrt(bgm(6));
bgl = bgm(3) - 1.5 * sqrt(bgm(6));
bg = (mstcd < bgh) & (mstcd > bgl) & bg;

% get background time course
bgrndslc = zeros(size(allvoxslc));
fndi = zeros(1, sstcd(3));
lsqueeze = ne_methods.lsqueeze;
for sc = 1:sstcd(3)
    vi = find(lsqueeze(bg(:, :, sc)));
    sltc = stcd(:, vi + (sc - 1) * prod(sstcd(1:2)));
    sltt = all(sltc > bgl & sltc < bgh, 1);
    bg(vi(~sltt) + (sc - 1) * prod(sstcd(1:2))) = false;
    bgrndslc(:, sc) = (1 / sum(sltt)) * sum(sltc(:, sltt), 2);
    fndi(sc) = sum(sltt);
end
bgrnd = (1 / sum(fndi)) * sum(bgrndslc .* (ones(nvol, 1) * fndi), 2);

% now try to get good forground estimate as well
bgh = bgh + 1.5 * sqrt(bgm(6));
fg = (mstcd > bgh);
[fgs, fgc] = ne_methods.clustercoordsc(fg, 1, floor(prod(sstcd) / 64));
if isempty(fgs)
    warning('neuroelf:xff:internalError', ...
        'No large enough brain found in dataset for computation.');
    gtc = struct('allvoxslc', allvoxslc, 'allvoxtot', allvox, 'bgrndmsk', bg, ...
        'bgrndslc', bgrndslc, 'bgrndtot', bgrnd, 'fgrndmsk', false(size(bg)), ...
        'fgrndslc', zeros(nvol, sstcd(3)), 'fgrndtot', zeros(nvol, 1), ...
        'highimsk', false(size(bg)), 'highislc', zeros(nvol, sstcd(3)), ...
        'highitot', zeros(nvol, 1));
    return;
end
eop = zeros(3, 3);
eop = logical(cat(3, eop, eop + 1, eop));
fg = ne_methods.erode3d(fgc == ne_methods.maxpos(fgs), eop);

% get foreground time course
fgrndslc = zeros(size(allvoxslc));
for sc = 1:sstcd(3)
    vi = find(lsqueeze(fg(:, :, sc)));
    sltc = stcd(:, vi + (sc - 1) * prod(sstcd(1:2)));
    sltt = all(sltc > bgh, 1);
    fg(vi(~sltt) + (sc - 1) * prod(sstcd(1:2))) = false;
    fgrndslc(:, sc) = (1 / sum(sltt)) * sum(sltc(:, sltt), 2);
    fndi(sc) = sum(sltt);
end
fgrnd = (1 / sum(fndi)) * sum(fgrndslc .* (ones(nvol, 1) * fndi), 2);

% high intensity
him = ne_methods.minmaxmean(mstcd(fg), 5);
hi = fg & (mstcd >= him(3)) & (mstcd <= (him(3) + sqrt(him(6))));
highislc = zeros(size(allvoxslc));
for sc = 1:sstcd(3)
    vi = find(lsqueeze(hi(:, :, sc)));
    sltc = stcd(:, vi + (sc - 1) * prod(sstcd(1:2)));
    sltt = all(sltc > him(3), 1);
    hi(vi(~sltt) + (sc - 1) * prod(sstcd(1:2))) = false;
    highislc(:, sc) = (1 / sum(sltt)) * sum(sltc(:, sltt), 2);
    fndi(sc) = sum(sltt);
end
highi = (1 / sum(fndi)) * sum(highislc .* (ones(nvol, 1) * fndi), 2);

% put into output
gtc = struct('allvoxslc', allvoxslc, 'allvoxtot', allvox, 'bgrndmsk', bg, ...
    'bgrndslc', bgrndslc, 'bgrndtot', bgrnd, 'fgrndmsk', fg, 'fgrndslc', fgrndslc, ...
    'fgrndtot', fgrnd, 'highimsk', hi, 'highislc', highislc, 'highitot', highi);
