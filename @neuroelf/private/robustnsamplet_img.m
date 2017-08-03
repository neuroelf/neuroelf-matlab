function [tv, wv, ot, b] = robustnsamplet_img(v, g, opts)
% robustnsamplet_img  - apply robust N-sample test to set of images
%
% FORMAT:       [tv, wv, ot, bv] = robustnsamplet_img(v, g [, opts])
%
% Input fields:
%
%       v           Vx1 spm_vol structure or N-D data array, XxYxZxV
%       g           Vx1 group assignment (unique, ordered)
%       opts        optional settings
%        .olsalso   flag, do OLS regression as well (default: false)
%        .olsonly   only OLS (for other functions, default: false)
%        .olstpat   OLS t-vol file pattern (def.: [pwd/G%d-G%d_ot.img])
%        .tvolpat   t-volume filename pattern (def.: [pwd/G%d-G%d_t.img])
%        .usenulls  flag, accept 0 values as real (default: false);
%        .wvols     create weight volumes
%        .wvolsinv  invert weight volumes (so that data is 1 - w)
%        .wvolspat  w-volume filename pattern (def.: [pwd/sub%03d_w.img])
%
% Output fields:
%
%       tv          either set of pairwise t-test volumes or XxYxZxP data
%       wv          either set of weighting volumes or XxYxZxV data
%       ot          either set of pairwise OLS t-test vols or data
%       bv          beta-estimate volumes from robust regression
%
% Note: t-stats will be re-computed to nominal d.f. !

% Version:  v0.9d
% Build:    14070313
% Date:     Jul-03 2014, 1:10 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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
   ((~isstruct(v) || ...
     numel(g) ~= numel(v)) && ...
    ((~isa(v, 'double') && ...
      ~isa(v, 'single')) || ...
     ndims(v) ~= 4 || ...
     numel(g) ~= size(v, 4))) || ...
    isempty(v) || ...
   ~isa(g, 'double') || ...
    numel(g) ~= length(g) || ...
    any(isinf(g) | isnan(g))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
fvols = false;
if isstruct(v)
    try
        fvols = true;
        vo = v;
        v = spm_read_vols(vo);
    catch ne_eo;
        error( ...
            'neuroelf:SPMError', ...
            'Error reading from input volumes: ''%s''.', ...
            ne_eo.message ...
        );
    end
end
if ~isa(v, 'double')
    v = double(v);
end
sv = size(v);
sv3 = sv(1:3);
ns = sv(4);
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'olsalso') || ...
   ~islogical(opts.olsalso) || ...
    numel(opts.olsalso) ~= 1
    opts.olsalso = false;
end
if ~isfield(opts, 'olsonly') || ...
   ~islogical(opts.olsonly) || ...
    numel(opts.olsonly) ~= 1
    opts.olsonly = false;
end
if opts.olsonly
    opts.wvols = false;
end
if ~isfield(opts, 'olstpat') || ...
   ~ischar(opts.oltspat) || ...
    numel(opts.olstpat) ~= length(opts.olstpat) || ...
    sum(opts.olstpat == '%') ~= 2
    opts.olstpat = [pwd '/G%d-G%d_ot.img'];
else
    opts.olstpat = opts.olstpat(:)';
end
if ~isfield(opts, 'tvolpat') || ...
   ~ischar(opts.tvolpat) || ...
    numel(opts.tvolpat) ~= length(opts.tvolpat) || ...
    sum(opts.tvolpat == '%') ~= 2
    opts.tvolpat = [pwd '/G%d-G%d_t.img'];
else
    opts.tvolpat = opts.tvolpat(:)';
end
if ~isfield(opts, 'usenulls') || ...
   ~islogical(opts.usenulls) || ...
    numel(opts.usenulls) ~= 1
    opts.usenulls = false;
end
if ~isfield(opts, 'wvols') || ...
   ~islogical(opts.wvols) || ...
    numel(opts.wvols) ~= 1
    opts.wvols = false;
end
if ~isfield(opts, 'wvolsinv') || ...
   ~islogical(opts.wvolsinv) || ...
    numel(opts.wvolsinv) ~= 1
    opts.wvolsinv = false;
end
if ~isfield(opts, 'wvolspat') || ...
   ~ischar(opts.wvolspat) || ...
    numel(opts.wvolspat) ~= length(opts.wvolspat) || ...
    sum(opts.wvolspat == '%') ~= 1
    opts.wvolspat = [pwd '/sub%03d_w.img'];
else
    opts.wvolspat = opts.wvolspat(:)';
end

% get group indices
g = g(:);
ug = unique(g);
ng = numel(ug);
df = numel(g) - ng;

% one-sample t-test
if ng < 2

    % robust
    if ~opts.olsonly
        [b, wv] = fitrobustbisquare_img(ones(ns, 1), v);
        tv = robustt(ones(ns, 1), v, b, wv, 1);
        tv(isinf(tv) | isnan(tv)) = 0;
    end

    % OLS
    if opts.olsalso || ...
        opts.olsonly
        ot = ((1 / sqrt(ns)) .* sum(v, 4)) ./ std(v, [], 4);
        ot(isinf(ot) | isnan(ot)) = 0;

        % OLS only
        if opts.olsonly
            tv = ot;
        end
    end

    % return early
    return;
end

% prepare for multi-sample t-test
gi = zeros(ng, ng);
gx = cell(1, ng);
gn = zeros(1, ng);
gc = 1;
for gc1 = 1:ng
    gx{gc1} = (g == ug(gc1));
    gn(gc1) = sum(gx{gc1});
    for gc2 = (gc1 + 1):ng
        gi(gc1, gc2) = gc;
        gc = gc + 1;
    end
end
[gi1, gi2] = ind2sub(size(gi), find(gi(:)));
np = numel(gi1);

% initialize output
tv = zeros([sv3, np]);
wv = [];
ot = [];
if opts.olsalso
    ot = tv;
end

% perform computation
nv3 = prod(sv3);
v = reshape(v, [nv3, ns]);
if ~opts.usenulls
    v(v == 0) = NaN;
end
if opts.wvols
    [b, tvl, dfvl, wv] = robustnsamplet(v, g);
    if opts.wvolsinv
        wv = 1 - wv;
    end
    wv = reshape(wv, sv);
elseif ~opts.olsonly
    [b, tvl, dfvl] = robustnsamplet(v, g);
end
if nargout > 3
    b = reshape(b, [sv3, size(b, 2)]);
end

% correct for bad values
if ~opts.olsonly
    badv = find(isinf(tvl) | isnan(tvl) | isnan(dfvl) | dfvl < 1);
    tvl(badv) = 0;
    dfvl(badv) = 1;

    % pack volumes
    for pc = 1:np
        tv(:, :, :, pc) = reshape(sdist('tinv', sdist('tcdf', ...
            tvl(:, gi1(pc), gi2(pc)), dfvl(:, gi1(pc), gi2(pc))), df), sv3);
    end
end

% perform OLS as well
if opts.olsalso || ...
    opts.olsonly
    v(isinf(v) | isnan(v)) = 0;
    sv = sum(v ~= 0, 2);
    mv = repmat(mean(v, 2) ./ (sv / size(v, 2)), [1, size(v, 2)]);
    v(v ~= 0) = v(v ~= 0) - mv(v ~= 0);
    gc = 1;
    for gc1 = 1:ng
        m1 = mean(v(:, gx{gc1}), 2);
        s1 = var(v(:, gx{gc1}), [], 2) ./ gn(gc1);
        for gc2 = (gc1 + 1):ng
            m2 = mean(v(:, gx{gc2}), 2);
            s2 = var(v(:, gx{gc2}), [], 2) ./ gn(gc2);
            t12 = (m1 - m2) ./ sqrt(s1 + s2);
            df12 = ((s1 + s2) .^ 2) ./ ...
                (s1 .* s1 ./ (gn(gc1) - 1) + s2 .* s2 ./ (gn(gc2) - 1));
            badv = find(isinf(t12) | isnan(t12) | isnan(df12) | df12 < 1);
            t12(badv) = 0;
            df12(badv) = 1;
            ot(:, :, :, gc) = reshape(-sign(t12) .* sdist('tinv', ...
                sdist('tcdf', -abs(t12), df12), df), sv3);
            gc = gc + 1;
        end
    end
    if opts.olsonly
        tv = ot;
        opts.olsalso = false;
    end
end

% create volumes and replace output arguments
if fvols
    vo = repmat(vo(1), [np, 1]);
    if opts.olsalso
        vo2 = vo;
    end
    for pc = 1:np
        vo(pc).fname = sprintf(opts.tvolpat, gi1(pc), gi2(pc));
        vo(pc) = spm_write_vol(vo(pc), tv(:, :, :, pc));
        if opts.olsalso
            vo2(pc).fname = sprintf(opts.olsvpat, gi1(pc), gi2(pc));
            vo2(pc) = spm_write_vol(vo2(pc), ot(:, :, :, pc));
        end
    end
    tv = vo;
    if opts.olsalso
        ot = vo2;
    end
    if opts.wvols
        vo = repmat(vo(1), [ns, 1]);
        for wc = 1:ns
            vo(wc).fname = sprintf(opts.wvolspat, wc);
            vo(wc) = spm_write_vol(vo(wc), wv(:, :, :, wc));
        end
        wv = vo;
    end
end
