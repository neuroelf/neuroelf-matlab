function v = shrinkfg(v)
% shrinkfg  - shrink foreground where no changes to topology are made

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

% expand v
sv = size(v);
v = flexmask(v, 0, 2, 0, sv + 2);
if sum(sum(v(2, :, :), 3), 2) > (0.5 * sv(2) * sv(3))
    v(1, :, :) = true;
end
if sum(sum(v(end-1, :, :), 3), 2) > (0.5 * sv(2) * sv(3))
    v(end, :, :) = true;
end
if sum(sum(v(:, 2, :), 3), 1) > (0.5 * sv(1) * sv(3))
    v(:, 1, :) = true;
end
if sum(sum(v(:, end-1, :), 3), 1) > (0.5 * sv(1) * sv(3))
    v(:, end, :) = true;
end
if sum(sum(v(:, :, 2), 2), 1) > (0.5 * sv(1) * sv(2))
    v(:, :, 1) = true;
end
if sum(sum(v(:, :, end-1), 2), 1) > (0.5 * sv(1) * sv(2))
    v(:, :, end) = true;
end
vs = size(v);

% find border worth checking
e = v & ~conv3d(v, 1);
e([1, vs(1)], :, :) = false;
e(:, [1, vs(1)], :) = false;
e(:, :, [1, vs(1)]) = false;

% for each voxel
[x, y, z] = ind2sub(vs, find(e));
xf = x - 1;
xt = x + 1;
yf = y - 1;
yt = y + 1;
zf = z - 1;
zt = z + 1;

fsi = reshape(1:27, [3, 3, 3]);
fsx = fsi(1, :, :);
msx = fsi(2, :, :);
lsx = fsi(3, :, :);
fsy = fsi(:, 1, :);
msy = fsi(:, 2, :);
lsy = fsi(:, 3, :);
fsz = fsi(:, :, 1);
msz = fsi(:, :, 2);
lsz = fsi(:, :, 3);
fsx = fsx(:)';
msx = msx(:)';
lsx = lsx(:)';
fsy = fsy(:)';
msy = msy(:)';
lsy = lsy(:)';
fsz = fsz(:)';
msz = msz(:)';
lsz = lsz(:)';
fsxs = fsx;
msxs = msx;
lsxs = lsx;
fsys = fsy;
msys = msy;
lsys = lsy;
fszs = fsz;
mszs = msz;
lszs = lsz;
cfx = 13;
cfy = 11;
cfz = 5;
clx = 15;
cly = 17;
clz = 23;
fsxs(fsxs == 13) = [];
msxs(msxs == 14) = [];
lsxs(lsxs == 15) = [];
fsys(fsys == 11) = [];
msys(msys == 14) = [];
lsys(lsys == 17) = [];
fszs(fszs == 5) = [];
mszs(mszs == 14) = [];
lszs(lszs == 23) = [];

vval = [1, 2, 4, 8, 16, 32, 64, 128];
vconn = [ ...
    0, ...
    2, 8, 16, 64, ...
    3, 6, 9, 20, 40, 96, 144, 192, ...
    11, 22, 104, 208, ...
    7, 41, 148, 224, ...
    15, 23, 43, 105, 150, 212, 232, 240, ...
    47, 151, 233, 244, ...
    31, 107, 214, 248, ...
    63, 111, 159, 215, 235, 246, 249, 252, ...
    191, 239, 247, 253, ...
    255, ...
    24, 66, ...
    25, 28, 56, 67, 70, 98, 152, 194, ...
    29, 99, 184, 198, ...
    57, 71, 156, 226, ...
    60, 102, 153, 195, ...
    61, 103, 157, 185, 188, 199, 227, 230, ...
    189, 231];

% check minicubes
for c = 1:numel(x)

    mc = v(xf(c):xt(c), yf(c):yt(c), zf(c):zt(c));
    fdx = mc(fsx);
    ldx = mc(lsx);
    fdy = mc(fsy);
    ldy = mc(lsy);
    fdz = mc(fsz);
    ldz = mc(lsz);
    fdxs = mc(fsxs);
    mdxs = mc(msxs);
    ldxs = mc(lsxs);
    fdys = mc(fsys);
    mdys = mc(msys);
    ldys = mc(lsys);
    fdzs = mc(fszs);
    mdzs = mc(mszs);
    ldzs = mc(lszs);
    if (~any(fdx) && ...
        (all(ldx) || ...
         (mc(clx) && ...
          any(vconn == vval * mdxs(:)) && ...
          (all(mdxs == ldxs))))) || ...
       (~any(ldx) && ...
        (all(fdx) || ...
         (mc(cfx) && ...
          any(vconn == vval * mdxs(:)) && ...
          (all(mdxs == fdxs))))) || ...
       (~any(fdy) && ...
        (all(ldy) || ...
         (mc(cly) && ...
          any(vconn == vval * mdys(:)) && ...
          (all(mdys == ldys))))) || ...
       (~any(ldy) && ...
        (all(fdy) || ...
         (mc(cfy) && ...
          any(vconn == vval * mdys(:)) && ...
          (all(mdys == fdys))))) || ...
       (~any(fdz) && ...
        (all(ldz) || ...
         (mc(clz) && ...
          any(vconn == vval * mdzs(:)) && ...
          (all(mdzs == ldzs))))) || ...
       (~any(ldz) && ...
        (all(fdz) || ...
         (mc(cfz) && ...
          any(vconn == vval * mdzs(:)) && ...
          (all(mdzs == fdzs)))))
        v(x(c), y(c), z(c)) = 0;
        continue;
    end
end

% put back in original size
v = v(2:end-1, 2:end-1, 2:end-1);
