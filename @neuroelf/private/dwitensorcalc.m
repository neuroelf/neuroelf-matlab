function [dt, b0, fa, evals, evecs, faev1, b2] = dwitensorcalc(dwi, bvals, bvecs, opts)
% dwitensorcalc  - compute diffusion tensor and other stats from DWI data
%
% FORMAT:       [dt, b0, fa, evals, evecs, faev1, bm] = dwitensorcalc(dwi, bvals, bvecs, opts)
%
% Input fields:
%
%       dwi         XxYxZxV diffusion weighted volumes
%       bvals       Vx1 B-values
%       bvecs       Vx3 B-vectors
%       opts        optional settings
%        .pbar      progress bar object
%        .res       data size dimension (default: [2, 2, 2])
%        .thresh    threshold applied to masking (default: 0.75) of Otsu
%        .trf       4x4 transformation matrix (used to determine rotations)
%
% Output fields:
%
%       dt          XxYxZx9 diffusion tensor (full 3x3 matrix packed)
%       b0          B0 estimate (from all data)
%       fa          FA (fractional anisotropy) estimate from eigenvalues
%       evals       eigenvalues (XxYxZx3)
%       evecs       eigenvectors (XxYxZx3x3, fifth dimension is vectors)
%       faev1       FA-weighted first eigenvector (suitable for RGB display)
%       bm          B-value/vector-based regression matrix
%
% Note: this function has been developed based on DTISearch/Streamline
%       available at http://www.mathworks.com/matlabcentral/fileexchange/34008

% Version:  v1.0
% Build:    15100617
% Date:     Oct-06 2015, 5:40 PM EST
% Author:   Chang Chia-Hao <oh75420 (at) gmail.com>
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
% Source:   http://www.mathworks.com/matlabcentral/fileexchange/34008

% Copyright (c) 2011 - 2012, Chang Chia-Hao; (c) 2014, Jochen Weber
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

% main argument check
if nargin < 3 || ...
    isempty(dwi) || ...
   ~isnumeric(dwi) || ...
    ndims(dwi) ~= 4 || ...
    ~isa(bvals, 'double') || ...
    numel(bvals) ~= size(dwi, 4) || ...
    any(isinf(bvals(:)) | isnan(bvals(:)) | bvals(:) < 0) || ...
    ~isa(bvecs, 'double') || ...
    ~isequal(size(bvecs), [size(dwi, 4), 3]) || ...
    any(isinf(bvecs(:)) | isnan(bvecs(:)) | abs(bvecs(:)) > 1) || ...
    any(all(bvecs(bvals ~= 0, :) == 0, 2))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% get size of volume
vxs = size(dwi);
vxs3 = prod(vxs(1:3));

% force type
if ~isa(dwi, 'single') && ...
   ~isa(dwi, 'double')
    dwi = double(dwi);
end

% check options
if nargin < 4 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end

% progress bar
if ~isfield(opts, 'pbar') || ...
    numel(opts.pbar) ~= 1 || ...
   (~isa(opts.pbar, 'xprogress') && ...
    ~isxfigure(opts.pbar, true))
    pbar = [];
else
    pbar = opts.pbar;
    pbarvis = pbar.Visible;
    pbar.Progress(0, 'Calculating DTI tensor and eigenvalues/vectors');
    pbar.Visible = 'on';
    drawnow;
end

% resolution
if ~isfield(opts, 'res') || ...
   ~isa(opts.res, 'double') || ...
    numel(opts.res) ~= 3 || ...
    any(isinf(opts.res) | isnan(opts.res) | opts.res <= 0)

    % assume 2mm cubed
    opts.res = [2, 2, 2];

% or
else
    
    % make a 1x3 vector
    opts.res = opts.res(:)';
end

% threshold (relative to Otso level, between 0.5 and 1.5)
if ~isfield(opts, 'thresh') || ...
   ~isa(opts.thresh, 'double') || ...
    numel(opts.thresh) ~= 1 || ...
    isinf(opts.thresh) || ...
    isnan(opts.thresh) || ...
    opts.thresh < 0.5
    opts.thresh = 0.75;
else
    opts.thresh = min(1.5, opts.thresh);
end

% transformation matrix
if ~isfield(opts, 'trf') || ...
   ~isa(opts.trf, 'double') || ...
   ~isequal(size(opts.trf), [4, 4]) || ...
    any(isinf(opts.trf(:)) | isnan(opts.trf(:))) || ...
    any(opts.trf(4, :) ~= [0, 0, 0, 1])

    % default to 3 cardinal axes
    opts.trf = diag(opts.res(:));

    % centering the dataset such that [0, 0, 0] is in the middle
    opts.trf(1:3, 4) = -opts.res(:) .* (1 + 0.5 .* vxs(1:3)');
    opts.trf(4, :) = [0, 0, 0, 1];
end

% B0 (no gradient applied) volumes
B0 = (bvals == 0);
B1 = ~B0;

% normalize gradients
gg = sqrt(sum(bvecs .* bvecs, 2));
g = bvecs ./ repmat(gg, [1, 3]);

% b-matrix (calculated from gradients and b-value)
b = bvals(:, ones(1, 6)) .* [ ...
    g(:, 1) .^ 2, ...
    2 .* g(:, 1) .* g(:, 2), ...
    2 .* g(:, 1) .* g(:, 3), ...
    g(:, 2) .^ 2, ...
    2 .* g(:, 2) .* g(:, 3), ...
    g(:, 3) .^ 2];

% force B0 to 0
b(B0, :) = 0;

% apply some ad-hoc masking to reduce computational effort
maxval = dwi(:, :, :, B1);
dwimask = mean(maxval, 4);
minval = sort(dwimask(:));
maxval = minval(floor(0.975 * numel(minval)));
minval = minval(ceil(0.025 * numel(minval)));
dwimask = limitrangec(dwimask, minval, maxval, minval);
minval = min(dwimask(:));
maxval = max(dwimask(:));
dwimask = (dwimask - minval) .* (1 / ((maxval - minval) + sqrt(eps)));

% find best threshold to separate into two classes
level = otsulevel(dwimask);

% apply mask
dwimask = (dwimask >= (opts.thresh * level));

% cluster foreground and background (to avoid separate "salt and pepper")
% [L, NUM] = bwlabeln(dwimask, 6);
% SL = L(L ~= 0);
% sd = hist(SL, 1:NUM);
% [M, I] = max(sd);
% dwimask(L ~= I) = false;
% dwimask = double(dwimask);
% dwimask = (imfill(dwimask, 6, 'holes') > 0);
[cs, dwimask] = clustercoordsc(dwimask, 1, 5);
[maxsize, maxpos] = max(cs(:));
dwimask = (dwimask == maxpos);
[cs, dwimask] = clustercoordsc(~dwimask, 1, 5);
[maxsize, maxpos] = max(cs(:));
dwimask = (dwimask ~= maxpos);
dwimask = dwimask(:);

% number of elements
nummask = sum(dwimask);

% generate matrix holding computed tensor values and estimate of B0
dtidata = nan(6, nummask);
dwib0 = nan(nummask, 1);
% chi_sq = dwib0;
% chi_sq_iqr = chi_sq;

% linear regression matrix (requires all-1 column for mean value)
b2 = [ones(size(dwi, 4), 1), -b];

% and transpose (for faster multiplication)
b2t = b2';

% permute into 1-by-volumes-by-totalvoxelcount array
dwi = permute(reshape(dwi, [1, vxs3, vxs(4)]), [1, 3, 2]);

% then mask data
dwi = dwi(:, :, dwimask);

% and make sure no values < 1 (log of DWI is used!!)
dwi(dwi < 1) = 1;

% square values (for estimation)
dwisq = dwi .* dwi;

% replicate values in first dimension (again for faster multiplication)
dwisq = repmat(dwisq, 7, 1);

% compute how many voxels are computed during single iteration
vstep = max(1, floor(2e7 / (vxs(4) * vxs(4))));

% replicate regression matrix accordingly
b2 = repmat(b2, [1, 1, vstep]);
b2t = repmat(b2t, [1, 1, vstep]);

% iterate over "packets" of (masked) data
for vc = 1:vstep:nummask

    % "to" index for packet
    vt = min(vc + vstep - 1, nummask);

    % number of voxels being processed
    vnum = vt + 1 - vc;

    % if less then step (possible on last iteration only)
    if vnum < vstep

        % shorten regression matrix
        b2(:, :, vnum+1:end) = [];
        b2t(:, :, vnum+1:end) = [];
    end

    % compute (weighted) estimates (using log of values)
    X = mtimesnd(mtimesnd(invnd(mtimesnd(b2t .* dwisq(:, :, vc:vt), ...
        b2)), b2t .* dwisq(:, :, vc:vt)), ...
        reshape(log(dwi(:, :, vc:vt)), [vxs(4), 1, vnum]));

    % compute fit (exp() back to original scale)
    fit = exp(mtimesnd(b2, X));

    % compute residual
    % resid = abs(squeeze(fit) - squeeze(dwi(:, :, vc:vt)));

    % store B0 estimate
    dwib0(vc:vt) = squeeze(fit(1, 1, :));
    
    % process residual
    % sresid = sort(resid);
    % ppos = 1 + [0.25, 0.5, 0.75] .* (vxs(4) - 1);
    % ppos2 = ceil(ppos);
    % ppos = floor(ppos);
    % y = 0.5 .* (sresid(ppos, :) + sresid(ppos2, :));
    % chi_sq(vc:vt) = y(2, :);
    % chi_sq_iqr(vc:vt) = y(3, :) - y(1, :);
    
    % store remaining data as tensor estimates
    dtidata(:, vc:vt) = squeeze(X(2:end, 1, :));

    % progress information
    if ~isempty(pbar)
        pbar.Progress(vt / nummask);
    end
end

% assign to first output
dt = nan(vxs3, 9);
dt(dwimask, [1, 2, 3, 5, 6, 9]) = dtidata';

% copy into according 3x3 matrix indices
dt(:, 4) = dt(:, 2);
dt(:, 7) = dt(:, 3);
dt(:, 8) = dt(:, 6);
dt = reshape(dt, [vxs(1:3), 9]);

% reset progress bar
if ~isempty(pbar)
    pbar.Visible = pbarvis;
end

% no more to do?
if nargout < 2
    return;
end

% assign to next output
b0 = nan(vxs3, 1);
b0(dwimask) = dwib0;
b0 = reshape(b0, vxs(1:3));

% no more to do?
if nargout < 3
    return;
end

% eigenvalue decomposition, get individual elements (faster access!)
a1 = dtidata(1, :);
a2 = dtidata(2, :);
a3 = dtidata(3, :);
a4 = dtidata(4, :);
a5 = dtidata(5, :);
a6 = dtidata(6, :);

% compute i1, i2, and i3
i1 = (1 / 3) .* (a1 + a4 + a6);
i2 = a1 .* a4 + a1 .* a6 + a4 .* a6;
i2 = i2 - (a2 .* a2 + a3 .* a3 + a5 .* a5);
i3 = a1 .* a4 .* a6 + 2 .* (a2 .* a3 .* a5);
i3 = i3 - (a6 .* a2 .* a2 + a4 .* a3 .* a3 + a1 .* a5 .* a5);

% then compute the factors to compute values
v = (i1 .* i1) - (1 / 3) .* i2;
s = (i1 .* i1 .* i1) - 0.5 .* (i1 .* i2 - i3);
t = min(1, max(-1, s ./ (v .* sqrt(v))));
t(isnan(t)) = 0;

% compute angle
phi = (1 / 3) .* acos(t);

% then store eigenvalues
d = zeros(3, nummask);
d(1, :) = i1 + 2 .* (sqrt(v) .* cos(phi));
d(2, :) = i1 - 2 .* (sqrt(v) .* cos(pi / 3 + phi));
d(3, :) = i1 - 2 .* (sqrt(v) .* cos(pi / 3 - phi));

% prepare eigenvector array
v = zeros([3, 3, nummask]);

% run along vector dimension
for ic = 1:3

    % computation
    a = a1 - d(ic, :);
    b = a4 - d(ic, :);
    c = a6 - d(ic, :);
    v(1, ic, :) = (a2 .* a5 - b .* a3) .* (a3 .* a5 - c .* a2);
    v(2, ic, :) = (a3 .* a5 - c .* a2) .* (a3 .* a2 - a .* a5);
    v(3, ic, :) = (a2 .* a5 - b .* a3) .* (a3 .* a2 - a .* a5);
end

% normalize each vector to length 1
qv = repmat(sqrt(sum(v(:, 1, :) .* v(:, 1, :), 1)), [3, 1, 1]);
qv(qv == 0) = 1;
v(:, 1, :) = v(:, 1, :) ./ qv;
qv = repmat(sqrt(sum(v(:, 2, :) .* v(:, 2, :), 1)), [3, 1, 1]);
qv(qv == 0) = 1;
v(:, 2, :) = v(:, 2, :) ./ qv;
qv = repmat(sqrt(sum(v(:, 3, :) .* v(:, 3, :), 1)), [3, 1, 1]);
qv(qv == 0) = 1;
v(:, 3, :) = v(:, 3, :) ./ qv;

% apply (inverse) transformation matrix
itrf = inv(opts.trf(1:3, 1:3));
itrf = itrf ./ repmat(sqrt(sum(itrf .* itrf, 2)), 1, 3);
v = mtimesnd(repmat(itrf, [1, 1, nummask]), v);

% make sure highest value (in vector) points in same (position) direction
[maxv, maxpos] = max(abs(v), [], 1);
maxpos = maxpos(:)' + (0:3:(numel(v)-1));
v = reshape(v, 3, 3 * nummask);
maxpos = (v(maxpos) < 0);
v(:, maxpos) = -v(:, maxpos);
v = reshape(v, [3, 3, nummask]);

% compute fractional anisotropy
fa = nan(vxs3, 1);
fa(dwimask) = min(1, max(0, fracanisotrop(d)));
fa = reshape(fa, vxs(1:3));

% no more to do?
if nargout < 4
    return;
end

% store eigenvalues
evals = NaN(vxs3, 3);
evals(dwimask, :) = d';
evals = reshape(evals, [vxs(1:3), 3]);

% no more to do?
if nargout < 5
    return;
end

% store eigenvectors
evecs = NaN([vxs3, 3, 3]);
evecs(dwimask, :, :) = permute(v, [3, 1, 2]);
evecs = reshape(evecs, [vxs(1:3), 3, 3]);

% no more to do?
if nargout < 6
    return;
end

% FA-weighted first eigenvector
faev1 = repmat(fa, [1, 1, 1, 3]) .* evecs(:, :, :, 1:3, 1);

% return regression matrix
if nargout > 6
    b2 = b2(:, :, 1);
end



% sub fuction fracanisotrop - compute fractional anisotropy
function fa = fracanisotrop(ev)

% make sure computation is valid
cmask = all(~isnan(ev), 1);
l1 = ev(1, cmask);
l2 = ev(2, cmask);
l3 = ev(3, cmask);

% computation
d1 = l1 - l2;
d1 = d1 .* d1;
d2 = l2 - l3;
d2 = d2 .* d2;
d1 = d1 + d2;
d2 = l3 - l1;
d2 = d2 .* d2;
d1 = d1 + d2;
d1 = sqrt(d1);
l1 = l1 .* l1;
l2 = l2 .* l2;
l3 = l3 .* l3;
l1 = l1 + l2;
l1 = sqrt(2 .* (l1 + l3));
fa = NaN(numel(cmask), 1);
fa(cmask) = d1(:) ./ l1(:);



% sub function invnd - matrix inversion along 3rd+ dim
function mi = invnd(mm)

% get size
ms = size(mm);
ms3 = prod(ms(3:end));

% prepare output
mi = zeros(size(mm));

% reshape (single 3rd dim!)
mm = reshape(mm, [ms(1:2), ms3]);

% iterate
for c = 1:ms3
    mi(:, :, c) = inv(mm(:, :, c));
end



% sub function mtimesnd - matrix multiplication along 3rd+ dim
function m = mtimesnd(f1, f2)

% get sizes
s1 = size(f1);
s2 = size(f2);

% determine output size
os = [s1(1), s2(2), s1(3:end)];
osx = prod(os(3:end));
ost = [1, 1, osx];

% determine output class
if isa(f1, 'single') && ...
    isa(f2, 'single')
    m = single(0);
    m(prod(os)) = 0;
    m = reshape(m, os);
else
    m = zeros(os);
end

% reshape if necessary
if numel(s1) > 3
    f1 = reshape(f1, [s1(1:2), osx]);
    f2 = reshape(f2, [s2(1:2), osx]);
    m  = reshape(m,  [os(1:2), osx]);
end

% iterate either along first two dims
if (s1(1) * s2(2)) < osx
    for c1 = 1:s1(1)
        for c2 = 1:s2(2)

            % perform multiplication and sum
            m(c1, c2, :) = m(c1, c2, :) + reshape(sum( ...
                reshape(f1(c1, :, :), s1(2), osx) .* ...
                reshape(f2(:, c2, :), s1(2), osx), 1), ost);
        end
    end
    
% or along 3rd+ dim
else
    for c1 = 1:osx
        m(:, :, c1) = f1(:, :, c1) * f2(:, :, c1);
    end
end

% reshape if necessary
if numel(os) > 3
    m = reshape(m, os);
end



% sub function otsulevel - estimate separation value (Otsu level) for image
function l = otsulevel(d)

% sort data in single precision
d = sort(single(d(:)));

% build "running sum of squares" and "running mean" to compute variance
dnum = single(1 ./ ((1:numel(d))'));
dnumm = single(1 ./ ((0:numel(d)-1)'));

% variance of class one (for each threshold)
rmean = cumsum(d);
rmeansq = cumsum(d .* d);
rvarc1 = (rmeansq - rmean .* rmean .* dnum) .* dnumm;

% reverse order
d = d(end:-1:1, 1);

% to compute variance of class two (for each threshold)
rmean = cumsum(d);
rmeansq = cumsum(d .* d);
rvarc2 = (rmeansq - rmean .* rmean .* dnum) .* dnumm;

% find level at the minimum of the sum of variances
rvarc1 = rvarc1 + rvarc2(end:-1:1, 1);
rvarc1(isnan(rvarc1)) = Inf;
[minvar, minpos] = min(rvarc1);
l = 0.5 * sum(d([0, 1] + (numel(d) - minpos)));
