function [ttrf, trp, rsl] = rbalign(v1, v2, opts)
% rbalign  - rigid-body alignment of one volume to a volume stack
%
% FORMAT:       [tfm, trp, aligned] = rbalign(target, totarget, opts)
%
% Input fields:
%
%       target      numeric volume data (will be used as target)
%       totarget    3d/4d numeric data (will be matched to target)
%       opts        optional settings
%        .hold      hold value for unsamplable voxels (default: 0)
%        .interpe   estim interpolation, {'linear'}, 'cubic', 'lanczos3'
%        .interpr   reslice interpolation, 'linear', 'cubic', {'lanczos3'}
%        .mask      boolean 3D volume, where bad values are rejected ([])
%        .params    1x6 boolean vector, 3 trans, 3 rotation, default: all
%        .pbar      xfigure::ProgressBar or xprogress object handle ([])
%        .prange    progress bar range (default: [0, 1])
%        .preuse    re-use params of volume (N-1) for Nth volume (true)
%        .robust    use robust regression to reject outliers (false)
%        .rtplot    real-time plot of param estimates (false)
%        .smooth    smoothing kernel for v2 (before detection, 0)
%        .smpl      interpolation sampling (in spacing units, voxelsize)
%        .trfv1     v1 voxel2mm mapping quaternion (default: eye(4))
%        .trfv2     v2 voxel2mm mapping quaternion (default: eye(4))
%        .tsmooth   smoothing kernel for v1 (before match)
%
% Output fields:
%
%        tfm        transformation matrices
%        trp        transformation parameters
%        aligned    corrected data

% Version:  v1.0
% Build:    16010821
% Date:     Jan-08 2016, 9:29 PM EST
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

% argument check
if nargin < 2 || ...
   ~isnumeric(v1) || ...
    ndims(v1) ~= 3 || ...
   ~isnumeric(v2) || ...
   (ndims(v2) ~= 3 && ...
    ndims(v2) ~= 4)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
nvol = size(v2, 4);
if nargin < 3 || ...
    numel(opts) ~= 1 || ...
   ~isstruct(opts)
    opts = struct;
end
if ~isfield(opts, 'hold') || ...
   ~isa(opts.hold, 'double') || ...
    numel(opts.hold) ~= 1
    opts.hold = 0;
end
if ~isfield(opts, 'interpe') || ...
   ~ischar(opts.interpe) || ...
   ~any(strcmpi(opts.interpe(:)', {'linear', 'cubic', 'lanczos3', 'lanczos5'}))
    opts.interpe = 'linear';
end
if ~isfield(opts, 'interpr') || ...
   ~ischar(opts.interpr) || ...
   ~any(strcmpi(opts.interpr(:)', {'linear', 'cubic', 'lanczos3', 'lanczos5'}))
    opts.interpr = 'lanczos3';
end
if ~isfield(opts, 'mask') || ...
   ~islogical(opts.mask) || ...
   ~isequal(size(opts.mask), size(v1))
    opts.mask = false(0, 0);
end
if ~isfield(opts, 'params') || ...
   ~islogical(opts.params) || ...
    numel(opts.params) ~= 6
    opts.params = true(1, 6);
else
    opts.params = opts.params(:)';
end
if ~any(opts.params)
    error( ...
        'neuroelf:BadArgument', ...
        'At least one of the parameters must be free.' ...
    );
end
if ~isfield(opts, 'pbar') || ...
    numel(opts.pbar) ~= 1 || ...
   ~any(strcmpi(class(opts.pbar), {'xfigure', 'xprogress'}))
    opts.pbar = [];
end
if ~isfield(opts, 'prange') || ...
   ~isa(opts.prange, 'double') || ...
    numel(opts.prange) ~= 2 || ...
    any(isinf(opts.prange) | isnan(opts.prange) | opts.prange < 0 | opts.prange > 1)
    opts.prange = [0, 1];
else
    opts.prange = sort(opts.prange(:))';
    if opts.prange(2) == opts.prange(1)
        opts.prange(2) = opts.prange(1) + sqrt(eps);
    end
end
if ~isfield(opts, 'preuse') || ...
   ~islogical(opts.preuse) || ...
    numel(opts.preuse) ~= 1
    opts.preuse = true;
end
if ~isfield(opts, 'robust') || ...
   ~islogical(opts.robust) || ...
    numel(opts.robust) ~= 1
    opts.robust = false;
end
if ~isfield(opts, 'rtplot') || ...
   ~islogical(opts.rtplot) || ...
    numel(opts.rtplot) ~= 1
    opts.rtplot = false;
end
if ~isfield(opts, 'smooth') || ...
   ~isa(opts.smooth, 'double') || ...
    numel(opts.smooth) ~= 3 || ...
    any(isinf(opts.smooth) | isnan(opts.smooth) | opts.smooth < 0)
    opts.smooth = [-1, -1, -1];
else
    opts.smooth = min(20, opts.smooth(:)');
end
if ~isfield(opts, 'smpl') || ...
   ~isa(opts.smpl, 'double') || ...
    numel(opts.smpl) ~= 3 || ...
    any(isinf(opts.smpl) | isnan(opts.smpl) | opts.smpl < 0.5)
    opts.smpl = [-1, -1, -1];
end
if ~isfield(opts, 'trfv1') || ...
   ~isa(opts.trfv1, 'double') || ...
   ~all(size(opts.trfv1) == 4) || ...
    numel(opts.trfv1) ~= 16 || ...
    any(isinf(opts.trfv1(:)) | isnan(opts.trfv1(:))) || ...
    any(opts.trfv1(4, :) ~= [0, 0, 0, 1])
    opts.trfv1 = eye(4);
end
if ~isfield(opts, 'trfv2') || ...
   ~isa(opts.trfv2, 'double') || ...
    size(opts.trfv2, 1) ~= 4 || ...
    size(opts.trfv2, 2) ~= 4 || ...
   ~any(size(opts.trfv2, 3) == [1, nvol]) || ...
    any(isinf(opts.trfv2(:)) | isnan(opts.trfv2(:))) || ...
    any(any(opts.trfv2(4, :, :) ~= repmat([0, 0, 0, 1], [1, 1, size(opts.trfv2, 3)])))
    opts.trfv2 = eye(4);
end
if size(opts.trfv2, 3) ~= nvol
    opts.trfv2 = opts.trfv2(:, :, ones(1, nvol));
end
otrf = opts.trfv2(:, :, 1);
if ~isfield(opts, 'tsmooth') || ...
   ~isa(opts.tsmooth, 'double') || ...
    numel(opts.tsmooth) ~= 3 || ...
    any(isinf(opts.tsmooth) | isnan(opts.tsmooth) | opts.tsmooth < 0)
    opts.tsmooth = opts.smooth;
else
    opts.tsmooth = min(20, opts.tsmooth(:)');
end
if any(opts.smpl < 0)
    sm = sqrt(sum(opts.trfv1(1:3,1:3) .^ 2));
    opts.smpl(opts.smpl < 0) = sm(opts.smpl < 0);
end

% smooth v1
if any(opts.smooth < 0)
    sms = sqrt(sum(opts.trfv1(1:3, 1:3) .^ 2));
    opts.smooth(opts.smooth < 0) = 2 * sms(opts.smooth < 0);
end
if any(opts.tsmooth < 0)
    sms = sqrt(sum(opts.trfv2(1:3, 1:3, 1) .^ 2));
    opts.tsmooth(opts.tsmooth < 0) = 2 * sms(opts.tsmooth < 0);
end
if any(opts.smooth)
    opts.smooth = (1 ./ sqrt(sum(opts.trfv2(1:3, 1:3, 1) .^ 2))) .* opts.smooth;
end

% progress bar
v1s = size(v1) + 0.25;
v2s = [size(v2), 1];
ttrf = opts.trfv2;
if any(size(v1) ~= v2s(1:3)) || ...
   ~all(lsqueeze(v1 == v2(:, :, :, 1)))
    fi = 1;
else
    fi = 2;
end
pbar = opts.pbar;
pbmin = opts.prange(1);
pbdif = opts.prange(2) - opts.prange(1);
if opts.robust
    mcef = 15;
else
    mcef = 4;
end
maxn = 5 + mcef * nvol;
if nargout > 2
    maxn = maxn + nvol;
end
if nvol > 1 && ...
    isempty(pbar)
    try
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 200, 640, 36]);
        xprogress(pbar, 'settitle', 'Rigid-body alignment...');
        xprogress(pbar, pbmin, sprintf('%d volumes: preparation...', nvol), ...
            'visible', 0, opts.prange(2));
        pbstr = '';
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        pbar = [];
    end
elseif ~isempty(pbar)
    pbar.Progress(pbmin, 'rbalign: preparation...');
    pbstr = 'rbalign: ';
end

% plot ?
if opts.rtplot
    itrf = inv(opts.trfv2(:, :, 1));
    rtf = figure;
    drawnow;
    set(0, 'CurrentFigure', rtf);
    set(rtf, 'Visible', 'on');
    drawnow;
    axes('parent', rtf);
    rtp = zeros(nvol, 6);
    rtl = plot(rtp);
    drawnow;
end

% compute gradient components
in3 = [Inf, Inf, Inf];
on3 = ones(1, 3);
v1 = double(v1);
if any(opts.tsmooth > 0)
    v1 = smoothdata3(v1, (1 ./ sqrt(sum(opts.trfv1(1:3, 1:3) .^ 2))) .* opts.tsmooth);
end
dgx = flexinterpn_method(v1, [in3;  1.125,  1.125, 1.125; on3; v1s], 0, opts.interpr);
dgy = dgx;
dgz = dgx;
cxyz = flexinterpn_method(v1, [in3;  1.125,  1.125, 0.875; on3; v1s], 0, opts.interpr);
dgx = dgx + cxyz;
dgy = dgy + cxyz;
dgz = dgz - cxyz;
cxyz = flexinterpn_method(v1, [in3;  1.125,  0.875, 1.125; on3; v1s], 0, opts.interpr);
dgx = dgx + cxyz;
dgy = dgy - cxyz;
dgz = dgz + cxyz;
cxyz = flexinterpn_method(v1, [in3;  1.125,  0.875, 0.875; on3; v1s], 0, opts.interpr);
dgx = dgx + cxyz;
dgy = dgy - cxyz;
dgz = dgz - cxyz;
if ~isempty(pbar)
    pbar.Progress(pbmin + pbdif / maxn);
end
cxyz = flexinterpn_method(v1, [in3;  0.875,  1.125, 1.125; on3; v1s], 0, opts.interpr);
dgx = dgx - cxyz;
dgy = dgy + cxyz;
dgz = dgz + cxyz;
cxyz = flexinterpn_method(v1, [in3;  0.875,  1.125, 0.875; on3; v1s], 0, opts.interpr);
dgx = dgx - cxyz;
dgy = dgy + cxyz;
dgz = dgz - cxyz;
cxyz = flexinterpn_method(v1, [in3;  0.875,  0.875, 1.125; on3; v1s], 0, opts.interpr);
dgx = dgx - cxyz;
dgy = dgy - cxyz;
dgz = dgz + cxyz;
cxyz = flexinterpn_method(v1, [in3;  0.875,  0.875, 0.875; on3; v1s], 0, opts.interpr);
dgx = dgx - cxyz;
dgy = dgy - cxyz;
dgz = dgz - cxyz;
if ~isempty(pbar)
    pbar.Progress(pbmin + pbdif * 2 / maxn);
end

% get coordinate step size
skipc = (1 ./ sqrt(sum(opts.trfv1(1:3,1:3).^2))) .* opts.smpl;

% get sample coordinates
[cxyz, d, dlim] = ndgrid(0.75:skipc(1):v1s(1), 0.75:skipc(2):v1s(2), 0.75:skipc(3):v1s(3));
cxyz = [cxyz(:), d(:), dlim(:)];

% mask coordinates
if ~isempty(opts.mask)
    d = flexinterpn_method(double(opts.mask), cxyz, 0, 'linear');
    cxyz(d < 0.5, :) = [];
end

% get data and gradients at coordinates
d = flexinterpn_method(v1, cxyz, 0, opts.interpe);
if ~isempty(pbar)
    pbar.Progress(pbmin + pbdif * 3 / maxn);
end
dgx = flexinterpn_method(dgx, cxyz, 0, opts.interpe);
dgy = flexinterpn_method(dgy, cxyz, 0, opts.interpe);
dgz = flexinterpn_method(dgz, cxyz, 0, opts.interpe);
if ~isempty(pbar)
    pbar.Progress(pbmin + 4 * pbdif / maxn);
end

% exclude data that has no real value
% rvx = [];
% cxyz(rvx, :) = [];
% d(rvx) = [];
% dgx(rvx) = [];
% dgy(rvx) = [];
% dgz(rvx) = [];
dlim = 0.25 + v2s(1:3);

% put into matrix for covariance
v1c = zeros(size(cxyz, 1), sum(opts.params));
ps = find(opts.params);
np = zeros(1, 6);
for pc = 1:numel(ps)
    pl = np;
    pl(ps(pc)) = 1e-6;
    dxyz = coords(cxyz, pl, opts.trfv1, opts.trfv1);
    v1c(:, pc) = -1e6 .* sum((dxyz - cxyz) .* [dgx dgy dgz], 2);
end

% would normally require orthogonalization with the original data!
% v1c = v1c - d * (pinv(d' * d) * d' * v1c);

% progress bar
if ~isempty(pbar)
    pbar.Progress(pbmin + 5 * pbdif / maxn);
end

% iterate over remaining images
for ic = fi:nvol

    % progress bar
    if ~isempty(pbar)
        pbar.Progress(pbmin + pbdif * (5 + mcef * (ic - 1)) / maxn, ...
            sprintf('%sdetect motion in vol %d/%d...', pbstr, ic, nvol));
    end

    % get smoothed version of data
	v2sm = smoothdata3(double(v2(:, :, :, ic)), opts.smooth);
    v2smm = uint8(v2sm >= 0.5);
    maxiter = 101;
    ss = Inf * ones(1, maxiter);
	pss = Inf;
    stablec = 0;
    if opts.preuse
        trf = repmat(otrf, [1, 1, maxiter]);
    else
        trf = repmat(opts.trfv2(:, :, max(1, ic - 1)), [1, 1, maxiter]);
    end
	while maxiter > 1
		dxyz = coords(cxyz, zeros(1, 6), opts.trfv1, trf(:, :, maxiter));
		msk = all(dxyz >= 0.5, 2) & (dxyz(:, 1) <= dlim(1)) & ...
            (dxyz(:, 2) <= dlim(2)) & (dxyz(:, 3) <= dlim(3));
        if sum(msk) < 32
            error( ...
                'neuroelf:InternalError', ...
                'Overlap too small.' ...
            );
        end
        msk = msk & flexinterpn_method(v2sm, dxyz, 0, 'linear') >= 0.5;
        f = flexinterpn_method(v2sm, dxyz(msk, :), 0, opts.interpe);
        cm1 = v1c(msk, :);
        dm1 = d(msk);
        sc = sum(dm1) / sum(f);
        dm1 = dm1 - f * sc;

        % which kind of regression
        if opts.robust
            soln = fitrobustbisquare(cm1, dm1);
        else
            soln = transmul(cm1) \ (cm1' * dm1);
        end

		p = zeros(1, 6);
		p(opts.params) = p(opts.params) + soln';

        maxiter = maxiter - 1;
		trf(:, :, maxiter) = inv(spmtrf(p(1:3), p(4:6))) * trf(:, :, maxiter + 1);
		ss(maxiter) = sum(dm1 .^ 2) / numel(dm1);
        if (pss - ss(maxiter)) / pss < 1e-6
            stablec = stablec + 1;
            if stablec > 2
                break;
            end
        else
            stablec = 0;
        end
        pss = ss(maxiter);
	end
    ttrf(:, :, ic) = trf(:, :, minpos(ss));
    if opts.rtplot
        mot = ttrf(:, :, ic) * itrf;
        pmot = spmitrf(mot);
        rtp(ic, :) = [pmot{1}, (180 / pi) * pmot{2}];
        for rtlc = 1:6
            set(rtl(rtlc), 'YData', rtp(:, rtlc));
        end
        drawnow;
    end
    if opts.preuse
        otrf = ttrf(:, :, ic);
    end
end

% create second output
if nargout > 1
    trp = ttrf;
    itrf = inv(opts.trfv2(:, :, 1));
    for ic = 1:size(trp, 3)
        trp(:, :, ic) = trp(:, :, ic) * itrf;
    end

    % create third output
    if nargout > 2
        if v2s(4) == 1
            rsl = v2(:, :, :);
            if ~isa(rsl, 'single') && ...
               ~isa(rsl, 'double')
                rsl = single(rsl);
            end
        else
            if ~isa(v2, 'single') && ...
               ~isa(v2, 'double')
                rsl = single(0);
                rsl(numel(v2)) = 0;
                rsl = reshape(rsl, size(v2));
                for ic = fi:nvol
                    rsl(:, :, :, ic) = v2(:, :, :, ic);
                end
            else
                rsl = v2(:, :, :, :);
            end
        end
        for ic = fi:nvol
            rsl(:, :, :, ic) = flexinterpn_method(rsl(:, :, :, ic), ...
                [Inf, Inf, Inf; 1, 1, 1; 1, 1, 1; size(rsl(:, :, :, 1))], ...
                opts.hold, itrf * inv(trp(:, :, ic)) * opts.trfv2(:, :, 1), opts.interpr);

            % progress bar
            if ~isempty(pbar)
                pbar.Progress(pbmin + pbdif * (5 + mcef * nvol + ic) / maxn, ...
                    sprintf('%scorrected motion in vol %d/%d...', pbstr, ic, nvol));
            end
        end
    end
end

% clear progress bar
if ~isempty(pbar) && ...
    isempty(opts.pbar)
    closebar(pbar);
end

% close plot
if opts.rtplot
    delete(rtf);
end

function xyz2 = coords(xyz, p, m1, m2)
xyz2 = xyz;
m  = (inv(m2) * inv(spmtrf(p(1:3), p(4:6))) * m1);
xyz2(:, 1) = m(1) * xyz(:, 1) + m(5) * xyz(:, 2) + m(9) * xyz(:, 3) + m(13);
xyz2(:, 2) = m(2) * xyz(:, 1) + m(6) * xyz(:, 2) + m(10) * xyz(:, 3) + m(14);
xyz2(:, 3) = m(3) * xyz(:, 1) + m(7) * xyz(:, 2) + m(11) * xyz(:, 3) + m(15);
