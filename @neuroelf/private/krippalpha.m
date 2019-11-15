function [alpha, cm] = krippalpha(d)
%KRIPPALPHA  Compute Krippendorff's alpha coefficient of agreement.
%   ALPHA = KRIPPALPHA(D) computes the alpha from data matrix D, where
%   the matrix can be given as either a Rx2 rating matrix, in which each
%   row contains a unit ID in the first column and a rating in the second
%   column or a UxR matrix, in which each row contains the ratings of R
%   raters, with NaN being missing values.
%
%   [ALPHA, CM] = KRIPPALPHA(D) also returns the coincidence matrix, CM.
%
%   Example:
%
%   % data matrix
%   D = [1, 2, 1; 1, 1, 1; NaN, 3, 3; 2, 2, 2; 2, 2, 2; ...
%        1, 1, 1; 4, 3, 4; 4, 4, NaN; 1, 1, 1; 1, 1, 1];
%   ALPHA = KRIPPALPHA(D)

% Version:  v1.1
% Build:    19111316
% Date:     Nov-13 2019, 4:04 PM EST
% Author:   Jochen Weber, Memorial Sloan Kettering Cancer Center, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2019, Jochen Weber
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

% check arguments
if nargin < 1 || ~isa(d, 'double') || ndims(d) > 2 || size(d, 2) < 2
    error('neuroelf:krippalpha:badArgument', 'Bad or missing argument.')
end

% decide between two raters and unit ID + rating
if size(d, 2) == 2 && ~any(isnan(d(:)))
    ud1 = unique(d(~isnan(d(:, 1)), 1));
    ud2 = unique(d(~isnan(d(:, 2)), 2));
    if numel(ud1) > numel(ud2) && numel(ud1) < (0.5 * size(d, 1)) && ...
        all(ud1 == round(ud1))
        mind = min(1, min(ud1));
        if mind < 1
            ud1 = ud1 - (mind - 1);
            d(:, 1) = d(:, 1) - (mind - 1);
        end
        idd = zeros(max(ud1), 1);
        idd(ud1) = 1:numel(ud1);
        d(:, 1) = idd(d(:, 1));
        h = histcount(d(:, 1), 1, max(ud1));
        mh = max(h);
        dd = NaN .* zeros(numel(ud1), mh);
        di = ones(numel(ud1), 1);
        for dc = 1:size(d, 1)
            dv = d(dc, 1);
            dix = di(dv);
            dd(dv, dix) = d(dc, 2);
            di(dv) = dix + 1;
        end
        d = dd;
    end
end

% remove invalid rows
d(sum(isnan(d), 2) >= (size(d, 2) - 1), :) = [];

% map values?
if any(~isnan(d(:)) & d(:) ~= round(d(:)))
end

% compute coincidence matrix
ud = unique(d(~isnan(d(:))));
nud = numel(ud);
cm = zeros(nud, nud);
for rc = 1:size(d, 1)
    r = d(rc, :);
    r = r(~isnan(r));
    nr = numel(r);
    for uc1 = 1:nud
        r1 = (r == ud(uc1));
        s1 = sum(r1);
        if s1 == 0
            continue;
        end
        if s1 > 1
            cm(uc1, uc1) = cm(uc1, uc1) + (s1 * (s1 - 1) / (nr - 1));
        end
        for uc2 = (uc1+1):nud
            r2 = (r == ud(uc2));
            s2 = sum(r2);
            if s2 > 0
                cm(uc1, uc2) = cm(uc1, uc2) + s1 * s2 / (nr - 1);
            end
        end
    end
end

% fill opposite triangle
for uc1 = 1:nud
    for uc2 = (uc1+1):nud
        cm(uc2, uc1) = cm(uc1, uc2);
    end
end

% compute a few intermediate variables
n = sum(cm(:));
nc = sum(cm);
ncc = sum(nc .* (nc - 1));
ncc = ncc / (n * (n - 1));

% and compute alpha
alpha = (sum(diag(cm)) / n - ncc) / (1 - ncc);
