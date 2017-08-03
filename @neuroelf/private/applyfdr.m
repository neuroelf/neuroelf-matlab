function th = applyfdr(stats, stype, levels, df1, df2, lnmeth, fdrfac)
% applyfdr  - apply FDR thresholding to given statistic
%
% FORMAT:       th = applyfdr(stats, stype, levels, df1 [, df2 [, lnmeth]])
%
% Input fields:
%
%       stats       statistical values (p/r/t/F)
%       stype       either of {'p'}, 'r', 't', 'F'
%       levels      q(FDR) levels
%       df1         d.f. 1
%       df2         d.f. 2
%       lnmeth      method to use, default: 3
%                   0 :  c(V) = 1
%                   1 :  c(V) = ln(V) + E
%                   2 :  c(V) = 1, use last crossing
%                   3 :  c(V) = ln(V) + E, use last crossing
%
% Output fields:
%
%       th          thresholds in original statistic
%
% See also fdr_thresholds

% Version:  v0.9c
% Build:    12011212
% Date:     Jan-05 2012, 11:29 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2012, Jochen Weber
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
if nargin < 3 || ...
   (~isa(stats, 'double') && ...
    ~isa(stats, 'single')) || ...
   ~ischar(stype) || ...
    numel(stype) ~= 1 || ...
   ~any('fprt' == lower(stype)) || ...
   (~isa(levels, 'double') && ...
    ~isa(levels, 'single')) || ...
    isempty(levels) || ...
    any(isinf(levels(:)) | isnan(levels(:)) | levels(:) <= 0 | levels(:) >= 1) || ...
   (lower(stype) ~= 'p' && ...
    (nargin < 4 || ...
     (~isa(df1, 'double') && ...
      ~isa(df1, 'single')) || ...
     numel(df1) ~= 1 || ...
     isnan(df1) || ...
     df1 < 1 || ...
     (lower(stype) ~= 'f' && ...
      df1 < 2) || ...
     (lower(stype) == 'f' && ...
      (nargin < 5 || ...
       (~isa(df2, 'double') && ...
        ~isa(df2, 'single')) || ...
       numel(df2) ~= 1 || ...
       isnan(df2) || ...
       df2 < 1))))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 6 || ...
    ~isa(lnmeth, 'double') || ...
    numel(lnmeth) ~= 1 || ...
    isinf(lnmeth) || ...
    isnan(lnmeth) || ...
   ~any((0:3) == lnmeth)
    lnmeth = 3;
end
if nargin < 7
    fdrfac = [];
end
nstats = numel(stats);
levels = levels(:);
stats = stats(:);
stats(isnan(stats) | stats == 0) = [];

% for empty stats, set to Bonferroni level - eps
if isempty(stats)

    blevels = (1 / nstats) * levels - eps;

    % depends on type
    switch (lower(stype))
        case {'f'}
            th = sdist('finv', 1 - blevels, df1, df2);
        case {'p'}
            th = blevels;
        case {'r'}
            th = correlinvtstat(-sdist('tinv', 0.5 * blevels, df1), df1 + 2);
        case {'t'}
            th = -sdist('tinv', 0.5 * blevels, df1);
    end
    if lnmeth
        th = th(:, [1, 1]);
    end
    return;
end

% depends on type
switch (lower(stype))

    % F stats
    case {'f'}
        th = sdist('finv', 1 - fdr_thresholds( ...
            1 - sdist('fcdf', abs(double(stats(:))), df1, df2), ...
            levels(:), lnmeth, fdrfac), df1, df2);

    % already p values
    case {'p'}
        th = fdr_thresholds(stats(:), levels(:), lnmeth, fdrfac);

    % correlation r
    case {'r'}
        if all(sign(stats(:)) == sign(stats(1)))
            th = correlinvtstat(-sdist('tinv', fdr_thresholds( ...
                1 - sdist('tcdf', abs(correltstat(stats(:), df1 + 2)), df1), ...
                levels(:), lnmeth, fdrfac), df1), df1 + 2);
        else
            th = correlinvtstat(-sdist('tinv', fdr_thresholds( ...
                1 - sdist('tcdf', abs(correltstat(stats(:), df1 + 2)), df1), ...
                levels(:) / 2, lnmeth, fdrfac), df1), df1 + 2);
        end

    % t stats
    case {'t'}
        if all(sign(stats(:)) == sign(stats(1)))
            th = -sdist('tinv', fdr_thresholds( ...
                1 - sdist('tcdf', abs(double(stats(:))), df1), ...
                levels(:), lnmeth, fdrfac), df1);
        else
            th = -sdist('tinv', fdr_thresholds( ...
                1 - sdist('tcdf', abs(double(stats(:))), df1), ...
                levels(:) / 2, lnmeth, fdrfac), df1);
        end
end
