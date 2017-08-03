function [tacc, cacc, sacc, cw] = slmdnclassify(data, conds)
%SLMDCLASSIFY  Mahalanobis-classification of searchlight data
%   TACC = SLMDCLASSIFY(DATA, CONDS) runs a classification with the matrix
%   in DATA, which must be samples-by-features, and CONDS must specify the
%   binary condition labels with values 1 and 2. The default is to run
%   a series of leave-one-out classifications (one of each condition), if
%   the number of combinations exceeds 5000, 5000 random pairs will be
%   picked. TACC then contains both the average accuracy (for the pair).
%
%   [TACC, CACC, SACC, SW] = SLMDCLASSIFY(DATA, CONDS) also returns
%   accuracy class in CACC, the accuracy per sample in SACC [0..1], and
%   the samples weight (if enough samples were given per class).

% Version:  v1.1
% Build:    16032921
% Date:     Mar-29 2016, 9:32 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2016, Jochen Weber
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

% arguments
if nargin < 2 || ~isa(data, 'double') || size(data, 1) < 4 || ...
   ~isa(conds, 'double') || size(data, 1) ~= numel(conds) || ...
    any(isinf(data(:)) | isnan(data(:))) || any(isinf(conds(:)) | isnan(conds(:)))
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
conds = conds(:);

% labels
labels = unique(conds);
if any(isinf(labels) | isnan(labels)) || numel(labels) < 2
    error('neuroelf:general:badArgument', 'Less than two or invalid labels.');
end

% some numbers we'll need, number of samples, features, and labels
nc = numel(conds);
nf = size(data, 2);
onf = ones(1, nf);
nl = numel(labels);

% storage for information, indices for each label, number of samples/label
li = cell(1, nl);
ln = zeros(1, nl);

% label-based means and total (class-weighted) mean
lmean = zeros(nl, nf);
tmean = zeros(1, nf);

% label-based and total (co)variances
lcov = cell(1, nl);
tcov = zeros(nf, nf);

% label-based weights and sum of weights
lw = cell(1, nl);
lws = zeros(1, nl);
tw = 0;
cw = zeros(1, nc);

% iterate over labels
for lc = 1:nl

    % find indices and number of samples/label
    li{lc} = lsqueeze(find(conds == labels(lc)));
    ln(lc) = numel(li{lc});

    % use robust estimate
    uri = (ln(lc) > nf);
    if uri
        try
            [lcov{lc}, lw{lc}, fmd, lmean(lc, :)] = robwcov(data(li{lc}, :));
            lws(lc) = sum(lw{lc});
        catch
            uri = false;
        end
    end
    
    % unweighted if too many features/samples/label (or error)
    if ~uri
        
        % compute mean (within class)
        lmean(lc, :) = (1 / ln(lc)) .* sum(data(li{lc}, :), 1);

        % remove mean
        databar = data(li{lc}, :) - ones(ln(lc), 1) * lmean(lc, :);

        % (co)variance
        lcov{lc} = (1 / (ln(lc) - 1)) .* databar' * databar;

        % full weights
        lw{lc} = ones(ln(lc), 1);
        lws(lc) = ln(lc);
    end

    % add to total (co)variance
    tcov = tcov + lws(lc) .* lcov{lc};
    tmean = tmean + lmean(lc, :);
    tw = tw + lws(lc);
    cw(li{lc}) = lw{lc};
end

% compute re-weighted (co)variance (of mean-removed samples)
tcov = (1 / tw) .* tcov;
tmean = (1 / nl) .* tmean;

% compute inverse (safely)
try
    uri = ((tw - nl) >= nf);
    if uri
        tcovi = inv(tcov);
    else
        tcovi = pinv(tcov);
    end
catch
    tacc = NaN;
    cacc = NaN .* ones(nl, 1);
    sacc = NaN .* ones(nc, 1);
    try
        if uri
            tcovi = pinv(tcov);
        end
    catch
        return;
    end
end

% remove mean, decorrelate
data = (data - tmean(ones(nc, 1), :)) * tcovi;
% for lc = 1:nl
%     lmean(lc, :) = (1 / lws(lc)) .* sum(lw{lc}(:, onf) .* data(li{lc}, :), 1);
% end

% compute inner product according to weighted distance, see
% http://jmlr.org/proceedings/papers/v2/ye07a/ye07a.pdf (pg. 2++)
p = permute(pairs((1:nl)'), [2, 3, 1]);
np = size(p, 1);
pd = zeros(nc, np);
for pc = 1:np
    p1 = p(pc, 1);
    p2 = p(pc, 2);
    l1 = labels(p1);
    l2 = labels(p2);
    l12 = 0.5 * (l1 + l2);
    s12 = 0.5 * (l1 - l2);
    bdiff = lmean(p(pc, 1), :) - lmean(p(pc, 2), :);
    pd(:, pc) = l12 + s12 .* sign(data * bdiff');
end

% accurate label
sacc = double(all(pd == conds(:, ones(1, np)), 2));

% total accuracy
tacc = (cw * sacc) / tw;

% class based accuracy
if nargout > 1
    cacc = zeros(1, nl);
    for lc = 1:nl
        cacc(lc) = sum(lw{lc} .* sacc(li{lc})) / lws(lc);
    end
end
