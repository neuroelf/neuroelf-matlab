% FUNCTION ne_rm_compute: compute RFX mediation
function ne_rm_compute(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:39 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% global variable
global ne_gcfg;
cf = ne_gcfg.h.RM.RMFig;
ch = ne_gcfg.h.RM.h;
ci = ne_gcfg.c.ini.Mediation;

% get GLM
glm = cf.UserData.lastglm;
if ~isxff(glm, 'glm')
    return;
end
rtv = glm.RunTimeVars;
rtvc = rtv.CovariatesData;

% get options
mopt = struct( ...
    'bootmeth', ci.BootSE, ...
    'bootn',    ci.BootNumSmp, ...
    'bootre',   logical(ci.BootReuse), ...
    'robust',   logical(ci.Robust), ...
    'sabmeth',  ci.Strategy, ...
    'ztrans',   logical(ci.ZTrans));

% get selection
ss = ch.Subjects.Value;
xt = ch.XType.Value;
xv = ch.XList.Value;
yt = ch.YType.Value;
yv = ch.YList.Value;
mn = ch.MCons.Value;
mv = ch.MCovs.Value;
cn = ch.CCons.Value;
cv = ch.CCovs.Value;

% check mediator!
if isempty(mn) && ...
    isempty(mv)
    uiwait(warndlg('No mediator (M) selected!', 'NeuroElf - info', 'modal'));
    return;
end
if ~isempty(intersect(mn, cn)) || ...
   ~isempty(intersect(mv, cv))
    uiwait(warndlg('Mediator (M) and covariates (C) overlap!', ...
        'NeuroElf - info', 'modal'));
    return;
end

% remove Infs/NaNs from regression
if xt == 2
    ss(any(isinf(rtvc(ss, xv)) | isnan(rtvc(ss, xv)), 2)) = [];
end
if yt == 2
    ss(any(isinf(rtvc(ss, yv)) | isnan(rtvc(ss, yv)), 2)) = [];
end
if ~isempty(mv)
    ss(any(isinf(rtvc(ss, mv)) | isnan(rtvc(ss, mv)), 2)) = [];
end
if ~isempty(cv)
    ss(any(isinf(rtvc(ss, cv)) | isnan(rtvc(ss, cv)), 2)) = [];
end

% compute D.F. (for t-values)
ns = numel(ss);
df = ns - (2 + numel(mn) + numel(mv) + numel(cn) + numel(cv));
dfp = numel(mn) + numel(mv);
if df < 2
    uiwait(warndlg('Too few subjects or too many variables selected!', ...
        'NeuroElf - info', 'modal'));
    return;
end

% set pointer shape
cf.Pointer = 'watch';
drawnow;

% get data for X and Y
if xt == 1
    if any(mn == xv)
        cf.Pointer = 'arrow';
        uiwait(warndlg('Mediator (M) must not contain X.', 'NeuroElf - info', ...
            'modal'));
        return;
    end
    if any(cn == xv)
        cf.Pointer = 'arrow';
        uiwait(warndlg('Covariate (C) must not contain X.', 'NeuroElf - info', ...
            'modal'));
        return;
    end
    x = glm.RFX_conmaps(rtv.Contrasts{xv, 2}, struct('subsel', ss));
else
    if any(mv == xv)
        cf.Pointer = 'arrow';
        uiwait(warndlg('Mediator (M) must not contain X.', 'NeuroElf - info', ...
            'modal'));
        return;
    end
    if any(cv == xv)
        cf.Pointer = 'arrow';
        uiwait(warndlg('Covariate (C) must not contain X.', 'NeuroElf - info', ...
            'modal'));
        return;
    end
    x = rtvc(ss, xv);
end
if yt == 1
    if any(mn == yv)
        cf.Pointer = 'arrow';
        uiwait(warndlg('Mediator (M) must not contain Y.', 'NeuroElf - info', ...
            'modal'));
        return;
    end
    if any(cn == yv)
        cf.Pointer = 'arrow';
        uiwait(warndlg('Covariate (C) must not contain Y.', 'NeuroElf - info', ...
            'modal'));
        return;
    end
    y = glm.RFX_conmaps(rtv.Contrasts{yv, 2}, struct('subsel', ss));
else
    if any(mv == yv)
        cf.Pointer = 'arrow';
        uiwait(warndlg('Mediator (M) must not contain Y.', 'NeuroElf - info', ...
            'modal'));
        return;
    end
    if any(cv == yv)
        cf.Pointer = 'arrow';
        uiwait(warndlg('Covariate (C) must not contain Y.', 'NeuroElf - info', ...
            'modal'));
        return;
    end
    y = rtvc(ss, yv);
end

% compile mediator and covariate
m = cell(numel(mn) + numel(mv), 1);
mc = 1;
c = cell(numel(cn) + numel(cv), 1);
cc = 1;
for ic = 1:numel(mn)
    m{mc} = glm.RFX_conmaps(rtv.Contrasts{mn(ic), 2}, struct('subsel', ss));
    mc = mc + 1;
end
if ~isempty(mn)
    ms = size(m{1});
    ms(end) = [];
    rs = [ones(1, numel(ms)), ns];
    if numel(ms) < 2
        ms(2) = 1;
    end
end
for ic = 1:numel(mv)
    m{mc} = rtvc(ss, mv(ic));
    if ~isempty(mn)
        m{mc} = repmat(reshape(m{mc}, rs), ms);
    end
    mc = mc + 1;
end
for ic = 1:numel(cn)
    c{cc} = glm.RFX_conmaps(rtv.Contrasts{cn(ic), 2}, struct('subsel', ss));
    cc = cc + 1;
end
if ~isempty(cn)
    ms = size(c{1});
    ms(end) = [];
    rs = [ones(1, numel(ms)), ns];
    if numel(ms) < 2
        ms(2) = 1;
    end
end
for ic = 1:numel(cv)
    c{cc} = rtvc(ss, cv(ic));
    if ~isempty(cn)
        c{cc} = repmat(reshape(c{cc}, rs), ms);
    end
    cc = cc + 1;
end

% concatenate M and C
if ~isempty(mn)
    m = cat(ndims(m{1}) + 1, m{:});
else
    m = cat(2, m{:});
end
if ~isempty(cn)
    c = cat(ndims(c{1}) + 1, c{:});
elseif ~isempty(cv)
    c = cat(2, c{:});
else
    c = [];
end

% for each variable, check validity
msk = true;
if numel(x) ~= max(length(x))
    md = ndims(x);
    msk = msk & (sum(isinf(x) | isnan(x) | x == 0, md) <= (0.5 * ns));
end
if numel(y) ~= max(length(y))
    md = ndims(y);
    msk = msk & (sum(isinf(y) | isnan(y) | y == 0, md) <= (0.5 * ns));
end
if numel(m) ~= (ns * (numel(mn) + numel(mv)))
    if (numel(mn) + numel(mv)) > 1
        md = ndims(m) - 1;
    else
        md = ndims(m);
    end
    bv = any(isinf(m) | isnan(m), md + 1) | all(m == 0, md + 1);
    msk = msk & (sum(bv, md) <= (0.5 * ns));
end
if ~isempty(cn)
    if (numel(cn) + numel(cv)) > 1
        md = ndims(c) - 1;
        bv = any(isinf(c) | isnan(c), md + 1) | all(c == 0, md + 1);
    else
        md = ndims(c);
        bv = (isinf(c) | isnan(c) | c == 0);
    end
    msk = msk & (sum(bv, md) <= (0.5 * ns));
end

% check mask
mskn = numel(msk);
msks = sum(msk(:));
if mskn > 0 && ...
    msks ~= mskn
    usemsk = true;

    % mask data accordingly
    if numel(x) ~= max(length(x))
        x = reshape(x, mskn, ns);
        x = x(msk(:), :);
    end
    if numel(y) ~= max(length(y))
        y = reshape(y, mskn, ns);
        y = y(msk(:), :);
    end
    if numel(m) ~= (ns * (numel(mn) + numel(mv)))
        m = reshape(m, mskn, ns * (numel(mn) + numel(mv)));
        m = m(msk(:), :);
        m = reshape(m, [msks, ns, numel(mn) + numel(mv)]);
    end
    if ~isempty(cn)
        c = reshape(c, mskn, ns * (numel(cn) + numel(cv)));
        c = c(msk(:), :);
        c = reshape(c, [msks, ns, numel(cn) + numel(cv)]);
    end

% no mask needed
else
    usemsk = false;
end

% close dialog
ne_rm_closeui;
drawnow;

% get bounding box
bb = glm.BoundingBox;

% compute mediation results
[p, se, t] = mediationpset(x, m, y, c, mopt);

% mask used
if usemsk

    % get number of maps
    np = size(p, ndims(p));
    msks = size(msk);

    % and then copy to temp variable and set in mask voxels
    mv = p;
    p = zeros(mskn, np);
    p(msk(:), :) = mv;
    p = reshape(p, [msks, np]);
    mv = t;
    t = zeros(mskn, np);
    t(msk(:), :) = mv;
    t = reshape(t, [msks, np]);
end

% single value(s)
if numel(p) == max(size(p))

    % display in output
    fprintf('Path coefficients:\n');
    fprintf('      c-path:   %.7f (se: %.7f)\n', p(1), se(1));
    if numel(p) == 5
        fprintf('      a-path:   %.7f (se: %.7f)\n', p(2), se(2));
        fprintf('      b-path:   %.7f (se: %.7f)\n', p(3), se(3));
        fprintf('    a*b-path:   %.7f (se: %.7f)\n', p(4), se(4));
    else
        for mc = 2:3:numel(p)
            fprintf('   a(%d)-path:   %.7f (se: %.7f)\n', floor(mc/3), p(mc), se(mc));
            fprintf('   b(%d)-path:   %.7f (se: %.7f)\n', floor(mc/3), p(mc+1), se(mc+1));
            fprintf(' a*b(%d)-path:   %.7f (se: %.7f)\n', floor(mc/3), p(mc+2), se(mc+2));
        end
    end
    fprintf('     c''-path:   %.7f (se: %.7f)\n\n', p(end), se(end));
    tp = 2 * sdist('tcdf', -abs(t), [df + dfp; ones(numel(t) - 1, 1) .* df]);
    fprintf('t-statistics:\n');
    fprintf('      c-path:   %.7f (p = %.7f)\n', t(1), tp(1));
    if numel(p) == 5
        fprintf('      a-path:   %.7f (p = %.7f)\n', t(2), tp(2));
        fprintf('      b-path:   %.7f (p = %.7f)\n', t(3), tp(3));
        fprintf('    a*b-path:   %.7f (p = %.7f)\n', t(4), tp(4));
    else
        for mc = 2:3:numel(p)
            fprintf('   a(%d)-path:   %.7f (p = %.7f)\n', floor(mc/3), t(mc), tp(mc));
            fprintf('   b(%d)-path:   %.7f (p = %.7f)\n', floor(mc/3), t(mc+1), tp(mc+1));
            fprintf(' a*b(%d)-path:   %.7f (p = %.7f)\n', floor(mc/3), t(mc+2), tp(mc+2));
        end
    end
    fprintf('     c''-path:   %.7f (p = %.7f)\n\n', t(end), tp(end));

% or maps
else

    % get map size
    ms = size(t);
    nm = ms(end);

    % create VMP
    vmp = newnatresvmp(bb.BBox, bb.ResXYZ(1), ones(2 * nm, 1));
    vmp.RunTimeVars.TrfPlus = rtv.TrfPlus;

    % put maps into VMP
    vmp.Map(1).VMPData = single(p(:, :, :, 1));
    vmp.Map(1).Type = 15;
    vmp.Map(1).Name = 'Mediation c-path';
    vmp.Map(nm + 1).VMPData = single(t(:, :, :, 1));
    vmp.Map(nm + 1).Name = 'Mediation c-path (t-map)';
    vmp.Map(nm + 1).LowerThreshold = -sdist('tinv', 0.005, df + dfp);
    vmp.Map(nm + 1).UpperThreshold = -sdist('tinv', 0.0001, df + dfp);
    vmp.Map(nm + 1).DF1 = df + dfp;
    for mc = 2:3:(nm-1)
        vmp.Map(mc).VMPData = single(p(:, :, :, mc));
        vmp.Map(mc).Type = 15;
        vmp.Map(mc).Name = 'Mediation a-path';
        vmp.Map(nm + mc).VMPData = single(t(:, :, :, mc));
        vmp.Map(nm + mc).Name = 'Mediation a-path (t-map)';
        vmp.Map(nm + mc).LowerThreshold = -sdist('tinv', 0.005, df);
        vmp.Map(nm + mc).UpperThreshold = -sdist('tinv', 0.0001, df);
        vmp.Map(nm + mc).DF1 = df;
        vmp.Map(mc + 1).VMPData = single(p(:, :, :, mc + 1));
        vmp.Map(mc + 1).Type = 15;
        vmp.Map(mc + 1).Name = 'Mediation b-path';
        vmp.Map(nm + mc + 1).VMPData = single(t(:, :, :, mc + 1));
        vmp.Map(nm + mc + 1).Name = 'Mediation b-path (t-map)';
        vmp.Map(nm + mc + 1).LowerThreshold = -sdist('tinv', 0.005, df);
        vmp.Map(nm + mc + 1).UpperThreshold = -sdist('tinv', 0.0001, df);
        vmp.Map(nm + mc + 1).DF1 = df;
        vmp.Map(mc + 2).VMPData = single(p(:, :, :, mc + 2));
        vmp.Map(mc + 2).Type = 15;
        vmp.Map(mc + 2).Name = 'Mediation a*b-path';
        vmp.Map(nm + mc + 2).VMPData = single(t(:, :, :, mc + 2));
        vmp.Map(nm + mc + 2).Name = 'Mediation a*b-path (t-map)';
        vmp.Map(nm + mc + 2).LowerThreshold = -sdist('tinv', 0.005, df);
        vmp.Map(nm + mc + 2).UpperThreshold = -sdist('tinv', 0.0001, df);
        vmp.Map(nm + mc + 2).DF1 = df;
    end
    vmp.Map(nm).VMPData = single(p(:, :, :, nm));
    vmp.Map(nm).Type = 15;
    vmp.Map(nm).Name = 'Mediation c''-path';
    vmp.Map(2 * nm).VMPData = single(t(:, :, :, nm));
    vmp.Map(2 * nm).Name = 'Mediation c''-path (t-map)';
    vmp.Map(2 * nm).LowerThreshold = -sdist('tinv', 0.005, df);
    vmp.Map(2 * nm).UpperThreshold = -sdist('tinv', 0.0001, df);
    vmp.Map(2 * nm).DF1 = df;

    % show/update correct VMP
    ne_openfile(0, 0, vmp, true);
end
