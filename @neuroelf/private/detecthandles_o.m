function dhv = detecthandles(v)
% detecthandles_o  - detect handles in 3D segmentation
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
% Build:    10062205
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
mh = floor(ms / 2);
fgc = 'edge';
bgc = 'edge';

% temporarily put object into larger frame
v = flexmask(v, 0, 2, 0, ms, [1, 1, 1]);

% create output
dhv = uint8(v);

% build sums over image
sx = squeeze(sum(sum(v, 3), 2));
px = find(sx(:)' > 0);
px = px(1):px(end);
ax = px(1) - 1;
sy = squeeze(sum(sum(v(px, :, :), 3), 1));
py = find(sy(:)' > 0);
py = py(1):py(end);
ay = py(1) - 1;
sz = squeeze(sum(sum(v(px, py, :), 2), 1));
pz = find(sz(:)' > 0);
pz = pz(1):pz(end);
az = pz(1) - 1;

% find "center" voxel (sort of center of gravity)
mx = maxpos(sx);
my = maxpos(squeeze(sum(v(mx, py, pz), 3))) + ay;
mz = floor(median(find(squeeze(v(mx, my, pz)))) + az);

% make sure foreground is a continuous object
fg = floodfill3(v, mx, my, mz, fgc);
dhv(v & ~fg) = uint8(3);

% remove holes in background (i.e. extra foreground clusters)
v(~fg) = false;

% make sure background is a continous object
if ms(3) > max(ms(1:2))
    fg = floodfill3(~v, mh(1), mh(2), 1, bgc);
elseif ms(2) > ms(1)
    fg = floodfill3(~v, mh(1), 1, mh(3), bgc);
else
    fg = floodfill3(~v, 1, mh(2), mh(3), bgc);
end
dhv(~v & ~fg) = uint8(2);

% remove holes in foreground (i.e. extra background clusters)
fg = ~fg;
v(fg) = true;

% iterate over X-slices (start with second slice)
for sc = 2:(numel(px)-1)

    % create clusters of foreground
    scc = px(sc);
    sd = v(scc, py, pz);
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
    cly = cl{end}(1, 2);
    clz = cl{end}(1, 3);
    for cc = 1:(numel(cl)-1)
        tc = cl{cc};
        rg1 = floodfill3(v(px(1):scc, py, pz), sc, tc(1, 2), tc(1, 3), fgc);
        rg2 = floodfill3(v(scc:px(end), py, pz), 1, tc(1, 2), tc(1, 3), fgc);
        if rg1(end, cly, clz) && ...
            rg2(1, cly, clz)
            mi = sub2ind(ms, repmat(scc, [size(tc, 1), 1]), tc(:, 2) + ay, tc(:, 3) + az);
            dhv(mi) = uint8(4);
        end
    end
end

% do the same for Y and Z
for sc = 2:(numel(py)-1)
    scc = py(sc);
    sd = v(px, scc, pz);
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
    clz = cl{end}(1, 3);
    for cc = 1:(numel(cl)-1)
        tc = cl{cc};
        rg1 = floodfill3(v(px, py(1):scc, pz), tc(1, 1), sc, tc(1, 3), fgc);
        rg2 = floodfill3(v(px, scc:py(end), pz), tc(1, 1), 1, tc(1, 3), fgc);
        if rg1(clx, end, clz) && ...
            rg2(clx, 1, clz)
            mi = sub2ind(ms, tc(:, 1) + ax, repmat(scc, [size(tc, 1), 1]), tc(:, 3) + az);
            dhv(mi) = uint8(5);
        end
    end
end
for sc = 2:(numel(pz)-1)
    scc = pz(sc);
    sd = v(px, py, scc);
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
        rg1 = floodfill3(v(px, py, pz(1):scc), tc(1, 1), tc(1, 2), sc, fgc);
        rg2 = floodfill3(v(px, py, scc:pz(end)), tc(1, 1), tc(1, 2), 1, fgc);
        if rg1(clx, cly, end) && ...
            rg2(clx, cly, 1)
            mi = sub2ind(ms, tc(:, 1) + ax, tc(:, 2) + ay, repmat(scc, [size(tc, 1), 1]));
            dhv(mi) = uint8(6);
        end
    end
end

% use this information to find candidates for background handles (todo...)

% get only part we need
dhv = dhv(2:end-1, 2:end-1, 2:end-1);
