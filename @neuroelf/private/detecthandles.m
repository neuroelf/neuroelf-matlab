function dhv = detecthandles(v)
% detecthandles  - detect handles in 3D segmentation
%
% FORMAT:       dhv = detecthandles(v)
%
% Input fields:
%
%       v           logical 3D volume; 0 = background, 1 = foreground
%
% Output fields:
%
%       dhv         uint8 volume with
%                   -  0 = background
%                   -  1 = foreground
%                   -  2 = holes in foreground (background clusters)
%                   -  3 = holes in background (foreground clusters)
%                   -  4 = foreground handle candidates (X-dir)
%                   -  5 = foreground handle candidates (Y-dir)
%                   -  6 = foreground handle candidates (Z-dir)

% todo
%                   -  7 = background handle candidates (X-dir)
%                   -  8 = background handle candidates (Y-dir)
%                   -  9 = background handle candidates (Z-dir)

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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
if nargin < 1 || ...
   ~islogical(v) || ...
    isempty(v) || ...
    size(v, 1) < 2 || ...
    size(v, 2) < 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
vs = size(v);
if numel(vs) < 3
    vs(3) = 1;
end
ms = vs + 2;
fgc = 'edge';
bgc = 'edge';

% temporarily put object into larger frame
v = flexmask(v, 0, 2, 0, ms, [1, 1, 1]);

% build sums over image
sx = squeeze(sum(sum(v, 3), 2));
px = find(sx(:)' > 0);
px = px(1):px(end);
ax = px(1) - 1;
bx = px(end) + 1;
sy = squeeze(sum(sum(v(px, :, :), 3), 1));
py = find(sy(:)' > 0);
py = py(1):py(end);
ay = py(1) - 1;
by = py(end) + 1;
sz = squeeze(sum(sum(v(px, py, :), 2), 1));
pz = find(sz(:)' > 0);
pz = pz(1):pz(end);
az = pz(1) - 1;
bz = pz(end) + 1;

% only consider part of object that contains useful information
v = v(ax:bx, ay:by, az:bz);
ms = size(v);
sx = sx(ax:bx);

% create output
dhv = uint8(v);

% find "center" voxel (sort of center of gravity)
mx = maxpos(sx);
my = maxpos(squeeze(sum(v(mx, :, :), 3)));
mz = floor(median(find(squeeze(v(mx, my, :)))));

% make sure foreground is a continuous object
fg = floodfill3(v, mx, my, mz, fgc);
dhv(v & ~fg) = uint8(3);

% remove holes in background (i.e. extra foreground clusters)
v(~fg) = false;

% make sure background is a continous object
mh = floor(size(v) / 2);
fg = floodfill3(~v, 1, mh(2), mh(3), bgc);
dhv(~v & ~fg) = uint8(2);

% remove holes in foreground (i.e. extra background clusters)
v(~fg) = true;

% find holes in foreground (i.e. fore- enclosing background)
csi = cell(1, 2);
v = ~v;
bf = ~floodfill3(v, 1, 1, 1, '5') & v;
pms = ms;
for sc = 2:(pms(3)-1)
    sd = bf(:, :, sc);
    [cl{1:2}] = clustercoords(sd, fgc);
    cl = cl{2};
    for cc = 1:numel(cl)
        tc = cl{cc};
        rg1 = floodfill3(v(:, :, 1:sc), tc(1, 1), tc(1, 2), sc, bgc);
        rg2 = floodfill3(v(:, :, sc:end), tc(1, 1), tc(1, 2), 1, bgc);
        if ~rg1(1, 1, end) || ...
            ~rg2(1, 1, 1)
            mi = sub2ind(pms, tc(:, 1), tc(:, 2), repmat(sc, [size(tc, 1), 1]));
            bf(mi) = false;
        end
    end
end
[cl{1:2}] = clustercoords(bf);
cl = cl{2};
for cc = 1:numel(cl)
    if ~all(cl{cc}(:, 3) == cl{cc}(1, 3))
        bb = false(size(bf));
        bb(sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3))) = true;
        smz = squeeze(sum(sum(bb, 2), 1));
        smz(smz == 0) = max(smz);
        msz = minpos(smz(:));
        cl{cc}(cl{cc}(:, 3) ~= msz, :) = [];
    end
end
bb = false(size(bf));
for cc = 1:numel(cl)
    bb(sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3))) = true;
end
dhv(bb) = uint8(9);
v = permute(v, [1, 3, 2]);
pms = ms([1, 3, 2]);
bf = ~floodfill3(v, 1, 1, 1, '5') & v;
for sc = 2:(pms(3)-1)
    sd = bf(:, :, sc);
    [cl{1:2}] = clustercoords(sd, fgc);
    cl = cl{2};
    for cc = 1:numel(cl)
        tc = cl{cc};
        rg1 = floodfill3(v(:, :, 1:sc), tc(1, 1), tc(1, 2), sc, bgc);
        rg2 = floodfill3(v(:, :, sc:end), tc(1, 1), tc(1, 2), 1, bgc);
        if ~rg1(1, 1, end) || ...
            ~rg2(1, 1, 1)
            mi = sub2ind(pms, tc(:, 1), tc(:, 2), repmat(sc, [size(tc, 1), 1]));
            bf(mi) = false;
        end
    end
end
[cl{1:2}] = clustercoords(bf);
cl = cl{2};
for cc = 1:numel(cl)
    if ~all(cl{cc}(:, 3) == cl{cc}(1, 3))
        bb = false(size(bf));
        bb(sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3))) = true;
        smz = squeeze(sum(sum(bb, 2), 1));
        smz(smz == 0) = max(smz);
        msz = minpos(smz(:));
        cl{cc}(cl{cc}(:, 3) ~= msz, :) = [];
    end
end
bb = false(size(bf));
for cc = 1:numel(cl)
    bb(sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3))) = true;
end
dhv(permute(bb, [1, 3, 2])) = uint8(8);
v = permute(v, [2, 3, 1]);
pms = pms([2, 3, 1]);
bf = ~floodfill3(v, 1, 1, 1, '5') & v;
for sc = 2:(pms(3)-1)
    sd = bf(:, :, sc);
    [cl{1:2}] = clustercoords(sd, fgc);
    cl = cl{2};
    for cc = 1:numel(cl)
        tc = cl{cc};
        rg1 = floodfill3(v(:, :, 1:sc), tc(1, 1), tc(1, 2), sc, bgc);
        rg2 = floodfill3(v(:, :, sc:end), tc(1, 1), tc(1, 2), 1, bgc);
        if ~rg1(1, 1, end) || ...
            ~rg2(1, 1, 1)
            mi = sub2ind(pms, tc(:, 1), tc(:, 2), repmat(sc, [size(tc, 1), 1]));
            bf(mi) = false;
        end
    end
end
[cl{1:2}] = clustercoords(bf);
cl = cl{2};
for cc = 1:numel(cl)
    if ~all(cl{cc}(:, 3) == cl{cc}(1, 3))
        bb = false(size(bf));
        bb(sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3))) = true;
        smz = squeeze(sum(sum(bb, 2), 1));
        smz(smz == 0) = max(smz);
        msz = minpos(smz(:));
        cl{cc}(cl{cc}(:, 3) ~= msz, :) = [];
    end
end
bb = false(size(bf));
for cc = 1:numel(cl)
    bb(sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3))) = true;
end
dhv(permute(bb, [3, 2, 1])) = uint8(7);
v = permute(~v, [3, 2, 1]);

% iterate over XY-slices (start with second slice)
bb = false(size(v));
pms = ms;
for sc = 2:(pms(3)-1)

    % create clusters of foreground
    sd = v(:, :, sc);
    [cl{1:2}] = clustercoords(sd, fgc);
    cl = cl{2};

    % continue if only one cluster (i.e. no handles possible !)
    if numel(cl) < 2
        continue;
    end

    % sort clusters according to their size
    cs = zeros(1, numel(cl));
    for cc = 1:numel(cl);
        cs(cc) = size(cl{cc}, 1);
    end
    [csi{1:2}] = sort(cs);
    cl = cl(csi{2});

    % use largest cluster as reference, check connectivity in both dirs
    clx = cl{end}(1, 1);
    cly = cl{end}(1, 2);
    for cc = 1:(numel(cl)-1)
        tc = cl{cc};
        rg1 = floodfill3(v(:, :, 1:sc), tc(1, 1), tc(1, 2), sc, fgc);
        rg2 = floodfill3(v(:, :, sc:end), tc(1, 1), tc(1, 2), 1, fgc);
        if rg1(clx, cly, end) && ...
            rg2(clx, cly, 1)
            mi = sub2ind(pms, tc(:, 1), tc(:, 2), repmat(sc, [size(tc, 1), 1]));
            bb(mi) = true;
        end
    end
end
[cl{1:2}] = clustercoords(bb);
bb(:) = false;
cl = cl{2};
for cc = 1:numel(cl)
    if ~all(cl{cc}(:, 3) == cl{cc}(1, 3))
        bb = false(size(bf));
        bb(sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3))) = true;
        smz = squeeze(sum(sum(bb, 2), 1));
        smz(smz == 0) = max(smz);
        msz = minpos(smz(:));
        cl{cc}(cl{cc}(:, 3) ~= msz, :) = [];
    end
end
for cc = 1:numel(cl)
    bb(sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3))) = true;
end
dhv(bb) = uint8(6);

% do the same for XZ- and YZ-planes
bf = permute(v, [1, 3, 2]);
pms = ms([1, 3, 2]);
bb = false(size(bf));
for sc = 2:(pms(2)-1)
    sd = bf(:, :, sc);
    [cl{1:2}] = clustercoords(sd, fgc);
    cl = cl{2};
    if numel(cl) < 2
        continue;
    end
    cs = zeros(1, numel(cl));
    for cc = 1:numel(cl);
        cs(cc) = size(cl{cc}, 1);
    end
    [csi{1:2}] = sort(cs);
    cl = cl(csi{2});
    clx = cl{end}(1, 1);
    cly = cl{end}(1, 2);
    for cc = 1:(numel(cl)-1)
        tc = cl{cc};
        rg1 = floodfill3(bf(:, :, 1:sc), tc(1, 1), tc(1, 2), sc, fgc);
        rg2 = floodfill3(bf(:, :, sc:end), tc(1, 1), tc(1, 2), 1, fgc);
        if rg1(clx, cly, end) && ...
            rg2(clx, cly, 1)
            mi = sub2ind(pms, tc(:, 1), tc(:, 2), repmat(sc, [size(tc, 1), 1]));
            bb(mi) = true;
        end
    end
end
[cl{1:2}] = clustercoords(bb);
bb(:) = false;
cl = cl{2};
for cc = 1:numel(cl)
    if ~all(cl{cc}(:, 3) == cl{cc}(1, 3))
        bb = false(size(bf));
        bb(sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3))) = true;
        smz = squeeze(sum(sum(bb, 2), 1));
        smz(smz == 0) = max(smz);
        msz = minpos(smz(:));
        cl{cc}(cl{cc}(:, 3) ~= msz, :) = [];
    end
end
for cc = 1:numel(cl)
    bb(sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3))) = true;
end
dhv(permute(bb, [1, 3, 2])) = uint8(5);

bf = permute(v, [3, 2, 1]);
pms = ms([3, 2, 1]);
bb = false(size(bf));
for sc = 2:(pms(1)-1)
    sd = bf(:, :, sc);
    [cl{1:2}] = clustercoords(sd, fgc);
    cl = cl{2};
    if numel(cl) < 2
        continue;
    end
    cs = zeros(1, numel(cl));
    for cc = 1:numel(cl);
        cs(cc) = size(cl{cc}, 1);
    end
    [csi{1:2}] = sort(cs);
    cl = cl(csi{2});
    clx = cl{end}(1, 1);
    cly = cl{end}(1, 2);
    for cc = 1:(numel(cl)-1)
        tc = cl{cc};
        rg1 = floodfill3(bf(:, :, 1:sc), tc(1, 1), tc(1, 2), sc, fgc);
        rg2 = floodfill3(bf(:, :, sc:end), tc(1, 1), tc(1, 2), 1, fgc);
        if rg1(clx, cly, end) && ...
            rg2(clx, cly, 1)
            mi = sub2ind(pms, tc(:, 1), tc(:, 2), repmat(sc, [size(tc, 1), 1]));
            bb(mi) = true;
        end
    end
end
[cl{1:2}] = clustercoords(bb);
bb(:) = false;
cl = cl{2};
for cc = 1:numel(cl)
    if ~all(cl{cc}(:, 3) == cl{cc}(1, 3))
        bb = false(size(bf));
        bb(sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3))) = true;
        smz = squeeze(sum(sum(bb, 2), 1));
        smz(smz == 0) = max(smz);
        msz = minpos(smz(:));
        cl{cc}(cl{cc}(:, 3) ~= msz, :) = [];
    end
end
for cc = 1:numel(cl)
    bb(sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3))) = true;
end
dhv(permute(bb, [3, 2, 1])) = uint8(4);

% finally try to cluster back and foreground holes where possible
v = false(size(dhv));
pms = size(v);
bb = (dhv == uint8(8)) | (dhv == uint8(9));
bf = (dhv == uint8(4));
[cl{1:2}] = clustercoords(bb | bf);
cl = cl{2};
for cc = 1:numel(cl)
    cli = sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3));
    bs = sum(bb(cli));
    fs = sum(bf(cli));
    if bs > 0 && ...
        fs > 0
        if bs > fs
            dhv(bs) = 0;
        elseif bs < fs
            dhv(fs) = 1;
        end
    end
end
bb = (dhv == uint8(7)) | (dhv == uint8(9));
bf = (dhv == uint8(5));
[cl{1:2}] = clustercoords(bb | bf);
cl = cl{2};
for cc = 1:numel(cl)
    cli = sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3));
    bs = sum(bb(cli));
    fs = sum(bf(cli));
    if bs > 0 && ...
        fs > 0
        if bs > fs
            dhv(bs) = 0;
        elseif bs < fs
            dhv(fs) = 1;
        end
    end
end
bb = (dhv == uint8(7)) | (dhv == uint8(8));
bf = (dhv == uint8(6));
[cl{1:2}] = clustercoords(bb | bf);
cl = cl{2};
for cc = 1:numel(cl)
    cli = sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3));
    bs = sum(bb(cli));
    fs = sum(bf(cli));
    if bs > 0 && ...
        fs > 0
        if bs > fs
            dhv(bs) = 0;
        elseif bs < fs
            dhv(fs) = 1;
        end
    end
end
bb = (dhv == uint8(7));
bf = (dhv == uint8(5)) | (dhv == uint8(6));
[cl{1:2}] = clustercoords(bb | bf);
cl = cl{2};
for cc = 1:numel(cl)
    cli = sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3));
    bs = sum(bb(cli));
    fs = sum(bf(cli));
    if bs > 0 && ...
        fs > 0
        if bs > fs
            dhv(bs) = 0;
        elseif bs < fs
            dhv(fs) = 1;
        end
    end
end
bb = (dhv == uint8(8));
bf = (dhv == uint8(4)) | (dhv == uint8(6));
[cl{1:2}] = clustercoords(bb | bf);
cl = cl{2};
for cc = 1:numel(cl)
    cli = sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3));
    bs = sum(bb(cli));
    fs = sum(bf(cli));
    if bs > 0 && ...
        fs > 0
        if bs > fs
            dhv(bs) = 0;
        elseif bs < fs
            dhv(fs) = 1;
        end
    end
end
bb = (dhv == uint8(9));
bf = (dhv == uint8(4)) | (dhv == uint8(5));
[cl{1:2}] = clustercoords(bb | bf);
cl = cl{2};
for cc = 1:numel(cl)
    cli = sub2ind(pms, cl{cc}(:, 1), cl{cc}(:, 2), cl{cc}(:, 3));
    bs = sum(bb(cli));
    fs = sum(bf(cli));
    if bs > 0 && ...
        fs > 0
        if bs > fs
            dhv(bs) = 0;
        elseif bs < fs
            dhv(fs) = 1;
        end
    end
end

% get only part we need
dhv = dhv(2:end-1, 2:end-1, 2:end-1);
dhv = flexmask(dhv, 0, 2, 0, vs, [ax - 1, ay - 1, az - 1]);
