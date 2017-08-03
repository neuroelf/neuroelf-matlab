function [b, t, df, wv, stat] = robustnsamplet(s, g)
% robustnsamplet  - performs robust fit and computes t scores
%
% FORMAT:       [b, t, df, w, stat] = robustnsamplet(s, g)
%
% Input fields:
%
%       s           N-dim double data, last, non-sigleton dim G
%       g           Gx1 group assignment (e.g. [1;1;1;1;1;2;2;2;2;2;2])
%
% Output fields:
%
%       b           (N-1)-dim double data with regression weights
%       t           (N-1)-dim double data with t-scores
%       df          (N-1)-dim double data with effective d.f.
%       w           N-dim weighting data
%       stat        (N-1)-dim stat structure with output of robustfit
%
% Note: Inf's and NaN's will be removed (correctly reflected in df!)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-08 2011, 9:16 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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
   ~isa(s, 'double') || ...
    isempty(s) || ...
   ~isnumeric(g) || ...
    isempty(g) || ...
   ~any(size(s) == numel(g)) || ...
    numel(g) ~= length(g) || ...
    any(isinf(g) | isnan(g))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% get correct dimesions
g = g(:);
ug = unique(g);
ns = numel(g);
ng = numel(ug);
sz = size(s);
dm = find(sz == numel(g));
dm = dm(end);
p1 = [setdiff(1:numel(sz), dm), dm];
dp = false;
if ~isequal(p1, 1:numel(sz))
    dp = true;
    p2 = 1:(numel(sz) - 1);
    p2 = [p2(1:(dm - 1)), numel(sz), p2(dm:end)];
    s = permute(s, p1);
    ts = sz(p1(1:end-1));
else
    ts = sz(1:end-1);
end
if numel(ts) < 2
    ts(2) = 1;
end

% build regression matrix X
X = zeros(numel(g), ng);
for uc = 1:ng
    X(g == ug(uc), uc) = 1;
end

% build group indices
gi = cell(1, ng);
gs = zeros(1, ng);
for g1c = 1:ng
    gi{g1c} = (X(:, g1c) == 1);
    gs(g1c) = sum(X(:, g1c));
end

% create output matrices
b = zeros([ts, ng]);
if nargout > 1
    t = zeros([ts, ng, ng]);
    if nargout > 2
        df = zeros([ts, ng, ng]);
        if nargout > 3
            wv = zeros(size(s));
            if nargout > 4
                [bnu, rnu, wnu, sstat] = fitrobustbisquare(X, randn(numel(g), 1));
                fns = fieldnames(sstat);
                for fnc = 1:numel(fns)
                    sstat.(fns{fnc}) = [];
                end
                stat = repmat(sstat, ts);
            end
        end
    end
end

% test xprogress
tf = prod(ts);
try
    pbar = xprogress;
    xprogress(pbar, 'setposition', [80, 200, 640, 36]);
    xprogress(pbar, 'settitle', sprintf('Performing %d-sample robust regression...', ng));
    xprogress(pbar, 0, sprintf('Running %d voxels...', tf), 'visible', 0, tf);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    pbar = [];
end

% loop
tfs = tf * ns;
tfg = tf * ng;
dgt = ns / 2;
for dc = 1:tf
    dti = dc:tf:tfs;
    dt = s(dti);
    dt = dt(:);
    dtg = (~isinf(dt) & ~isnan(dt));
    if sum(dtg) >= dgt
        try
            [bnu, r, w] = fitrobustbisquare(X(dtg, :), dt(dtg));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            continue;
        end
        b(dc:tf:end) = bnu;
        if nargout > 1
            w = sqrt(w);
            for g1c = 1:ng
                g1w = w(gi{g1c}(dtg));
                g1ws = sum(g1w .* g1w);
                g1w = g1w * (gs(g1c) / g1ws);
                g1r = r(gi{g1c}(dtg)) .* g1w;
                g1rv = var(g1r) * gs(g1c) / (g1ws * g1ws);
                for g2c = (g1c+1):ng
                    g2w = w(gi{g2c}(dtg));
                    g2ws = sum(g2w .* g2w);
                    g2w = g2w * (gs(g2c) / g2ws);
                    g2r = r(gi{g2c}(dtg)) .* g2w;
                    g2rv = var(g2r) * gs(g2c) / (g2ws * g2ws);
                    tstat = (bnu(g1c) - bnu(g2c)) / sqrt(g1rv + g2rv);
                    t(dc+(g1c-1)*tf+(g2c-1)*tfg) = tstat;
                    t(dc+(g2c-1)*tf+(g1c-1)*tfg) = tstat;
                    if nargout > 2
                        dfstat = (g1rv + g2rv) .^ 2 / ...
                            (g1rv * g1rv / (g1ws - 1) + g2rv * g2rv / (g2ws - 1));
                        df(dc+(g1c-1)*tf+(g2c-1)*tfg) = dfstat;
                        df(dc+(g2c-1)*tf+(g1c-1)*tfg) = dfstat;
                    end
                end
            end
        end
        if nargout > 3
            wv(dti(dtg)) = w .* w;
            if nargout > 4
                stat(dc) = sstat;
            end
        end
    end
    if mod(dc, 100) == 0 && ...
       ~isempty(pbar)
        xprogress(pbar, dc);
    end
end

% re-permute?
if dp
    b = squeeze(permute(b, [p2, numel(p2) + 1]));
    if nargout > 1
        t = squeeze(permute(t, [p2, (numel(p2) + 1):(numel(p2) + 2)]));
        if nargout > 2
            df = squeeze(permute(df, [p2, (numel(p2) + 1):(numel(p2) + 2)]));
            if nargout > 3
                wv = squeeze(permute(wv, [p2, numel(g)]));
                if nargout > 4
                    stat = permute(stat, p2);
                end
            end
        end
    end
else
    b = squeeze(b);
    if nargout > 1
        t = squeeze(t);
        if nargout > 2
            df = squeeze(df);
            if nargout > 3
                wv = squeeze(wv);
                if nargout > 4
                    stat = squeeze(stat);
                end
            end
        end
    end
end

% clear progress bar
if ~isempty(pbar)
    closebar(pbar);
end
