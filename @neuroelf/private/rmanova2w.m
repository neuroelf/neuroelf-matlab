function [stats, dfs, idx] = rmanova2w(y, s, f1, f2, fn)
% rm_anova2  - compute repeated-measures ANOVA with two-within factors
%
% FORMAT:       [stats, dfs, idx] = rmanova2w(y, s, f1, f2, fn)
%
% Input fields:
%
%       y           numeric input data (dependent variable, last dim)
%       s           subject-grouping variable (Lx1 vector)
%       f1, f2      within-factors grouping variable (Lx1 vectors)
%       fn          factor names (used only if y is a Lx1 vector)
%
% Output fields:
%
%       stats       if y is a Lx1 vector (single univariate case), original
%                   output is produced; otherwise a ND output is produced,
%                   where as the order of values in regression dim is
%                   - mean-sums-of-squares (MSS) within-factor 1
%                   - mean-sums-of-squares (MSS) subject-by-factor 1
%                   - mean-sums-of-squares (MSS) within-factor 2
%                   - mean-sums-of-squares (MSS) subject-by-factor 2
%                   - mean-sums-of-squares (MSS) within-factors interaction
%                   - mean-sums-of-squares (MSS) subject-by-factor intera.
%                   such that the three F-tests can be directly computed
%       dfs         degrees of freedom (in order of MSS outputs)
%       idx         F2-by-F1-by-S indices into input y (along ANOVA dim)
%
% Note: y should either be a 1-d column vector with all of numeric data, or
%       an ND-array with the last dimension being the subject-by-within
%       folded dimension (e.g. if you have a 3D brain data per subject
%       and within-factor combination with 18 subjects, and 3x2 conditions
%       the data would initially be a X-by-Y-by-Z-by-18-by-3-by-2 array,
%       which should be reshaped to a X-by-Y-by-Z-by-108 array!).
%       the grouping variables should be 1D numeric vector, each with same
%       length as y in the ANOVA dimension; each entry in each of the
%       grouping vectors indicates the subject or level number of the
%       corresponding entry in y; alternatively s, f1, and f2 can be
%       singleton values, in which case the assumption is that the order
%       of values is s-slowest and f2-fastest changing, such as in:
%       [1, 1, 1; 1, 1, 2; 1, 2, 1; 1, 2, 2; 1, 3, 1; 1, 3, 2; 2, 1, 1; ...]

% Version:  v0.9c
% Build:    13032315
% Date:     May-31 2011, 5:53 PM EST
% Author:   Aaron Schurger (2005.02.04)
%   Derived from Keppel & Wickens (2004) "Design and Analysis" ch. 18
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2004, 2005, 2011, 2013, Aaron Schurger
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
if nargin < 4 || ...
   ~isa(y, 'double') || ...
    isempty(y) || ...
   ~isa(s, 'double') || ...
   (~any(numel(s) == size(y)) && ...
    numel(s) ~= 1) || ...
    any(isinf(s(:)) | isnan(s(:)) | s(:) < 1 | s(:) ~= round(s)) || ...
   ~isa(f1, 'double') || ...
   (~any(numel(f1) == size(y)) && ...
    numel(f1) ~= 1) || ...
    any(isinf(f1(:)) | isnan(f1(:)) | f1(:) < 1 | f1(:) ~= round(f1)) || ...
   ~isa(f2, 'double') || ...
   (~any(numel(f2) == size(y)) && ...
    numel(f2) ~= 1) || ...
    any(isinf(f2(:)) | isnan(f2(:)) | f2(:) < 1 | f2(:) ~= round(f2))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
s = s(:);
if all(s < 3)
    error( ...
        'neuroelf:BadArgument', ...
        'Less than 3 subjects not supported.' ...
    );
end
f1 = f1(:);
if all(f1 < 2)
    error( ...
        'neuroelf:BadArgument', ...
        'Factor 1 must have at least 2 levels.' ...
    );
end
f2 = f2(:);
if all(f2 < 2)
    error( ...
        'neuroelf:BadArgument', ...
        'Factor 2 must have at least 2 levels.' ...
    );
end
if nargin < 5 || ...
   ~iscell(fn) || ...
    numel(fn) ~= 2
    fn = {'Factor1'; 'Factor2'};
else
    fn = fn(:);
end
if ~ischar(fn{1}) || ...
    isempty(fn{1})
    fn{1} = 'Factor1';
else
    fn{1} = fn{1}(:)';
end
if ~ischar(fn{2}) || ...
    isempty(fn{2})
    fn{2} = 'Factor2';
else
    fn{2} = fn{2}(:)';
end

% single univariate case
if size(y, 1) == numel(y)

    % store as 1xL array
    suv = true;
    y = y(:)';
    st = [1, 1];
    ny = numel(y);

% mass-univariate case
else
    suv = false;

    % get size and reshape to 2D
    sy = size(y);
    st = [sy(1, 1:end-1), 6];
    y = reshape(y, prod(st(1:end-1)), sy(end));
    ny = size(y, 2);
end
nv = size(y, 1);

% deal with singleton cases
if numel(s) == 1
    if mod(ny, s) ~= 0
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid number of maps for s * f1f2 combination.' ...
        );
    end
    s = lsqueeze(repmat(1:s, round(ny / s), 1));
end
if numel(s) ~= ny
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid number of maps or s assignments.' ...
    );
end
n = max(s);
if n ~= numel(unique(s)) || ...
    any(histcount(s, 1, n, 1) ~= round(ny / n))
    error( ...
        'neuroelf:BadArgument', ...
        'Missing s assignments or imbalanced design unsupported.' ...
    );
end

% number of within level combinations
nw = round(ny / n);
if numel(f1) == 1
    if mod(nw, f1) ~= 0
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid number of factor levels for f1.' ...
        );
    end
    f1 = repmat(lsqueeze(repmat(1:f1, round(nw / f1), 1)), n, 1);
end
if numel(f1) ~= ny
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid number of maps or f1 assignments.' ...
    );
end
a = max(f1);
if a ~= numel(unique(f1)) || ...
    any(histcount(f1, 1, a, 1) ~= round(ny / a))
    error( ...
        'neuroelf:BadArgument', ...
        'Missing f1 assignments or imbalanced design unsupported.' ...
    );
end
if mod(nw, a) ~= 0
    error( ...
        'neuroelf:BadArgument', ...
        'Factor crossing invalid.' ...
    );
end
if numel(f2) == 1
    if f2 ~= round(nw / a)
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid number of levels for f2.' ...
        );
    end
    f2 = repmat((1:f2)', n * a, 1);
    if numel(f2) ~= ny
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid number of levels for f2.' ...
        );
    end
end
b = max(f2);
if b ~= numel(unique(f2)) || ...
    any(histcount(f2, 1, b, 1) ~= round(ny / b))
    error( ...
        'neuroelf:BadArgument', ...
        'Missing f2 assignments or imbalanced design unsupported.' ...
    );
end

% test unique combinations
if size(unique([s, f1, f2], 'rows'), 1) ~= ny
    error( ...
        'neuroelf:BadArgument', ...
        's-by-f1-by-f2 combinations must be fully crossed.' ...
    );
end

% hold the index for each combination
idx = zeros([b, a, n]);
for ia = 1:a
    am = f1 == ia;
    for ib = 1:b
        bm = f2 == ib;
        for in = 1:n
            idx(ib, ia, in) = find(am & bm & s == in);
        end
    end
end

% reorder data
y = y(:, idx(:)');

% reshape accordingly
y = reshape(y, [nv, size(idx)]);

% make tables (see table 18.1, p. 402)
ab = reshape(sum(y, 4), [nv, b, a]); % across subjects
as = reshape(sum(y, 2), [nv, a, n]); % across factor 2
bs = reshape(sum(y, 3), [nv, b, n]); % across factor 1

% sum again
A = sum(ab, 2);
B = sum(ab, 3);
s = sum(as, 2);
t = sum(A, 3);

% degrees of freedom
dfa = a - 1;
dfb = b - 1;
dfab = (a - 1) * (b - 1);
dfas = (a - 1) * (n - 1);
dfbs = (b - 1) * (n - 1);
dfabs = (a - 1) * (b - 1) * (n - 1);

% bracket terms (expected value)
expa = sum(A .* A, 3) ./ (b * n);
expb = sum(B .* B, 2) ./ (a * n);
expab = sum(sum(ab .* ab, 2), 3) ./ n;
exps = sum(s .* s, 3) ./ (a * b);
expas = sum(sum(as .* as, 2), 3) ./ b;
expbs = sum(sum(bs .* bs, 2), 3) ./ a;
expy = sum(sum(sum(y .* y, 2), 3), 4);
expt = (t .* t) ./ (a * b * n);

% sums of squares
ssa = expa - expt;
ssb = expb - expt;
ssab = expab - expa - expb + expt;
ssas = expas - expa - exps + expt;
ssbs = expbs - expb - exps + expt;
ssabs = expy - expab - expas - expbs + expa + expb + exps - expt;

% mean squares
msa = ssa ./ dfa;
msb = ssb ./ dfb;
msab = ssab ./ dfab;
msas = ssas ./ dfas;
msbs = ssbs ./ dfbs;
msabs = ssabs ./ dfabs;

% return values
if suv

    % f statistic
    fa = msa ./ msas;
    fb = msb ./ msbs;
    fab = msab ./ msabs;

    % p values
    pa = sdist('fcdf', fa, dfa, dfas, true);
    pb = sdist('fcdf', fb, dfb, dfbs, true);
    pab = sdist('fcdf', fab, dfab, dfabs, true);
    dfs = [dfa, dfas, dfb, dfbs, dfab, dfabs];

    stats = {'Source','SS','df','MS','F','p';...
             fn{1}, ssa, dfa, msa, fa, pa;...
             fn{2}, ssb, dfb, msb, fb, pb;...
             [fn{1} ' x ' fn{2}], ssab, dfab, msab, fab, pab;...
             [fn{1} ' x Subj'], ssas, dfas, msas, [], [];...
             [fn{2} ' x Subj'], ssbs, dfbs, msbs, [], [];...
             [fn{1} ' x ' fn{2} ' x Subj'], ssabs, dfabs, msabs, [], []};
else
    stats = reshape([msa, msas, msb, msbs, msab, msabs], st);
end
