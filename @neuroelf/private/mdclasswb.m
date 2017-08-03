function [icov, lmean] = mdclasswb(data, conds)
%MDCLASSWB  Mahalanobis-classification inverse covariance and means vectors
%   [ICOV, LMEAN] = SLMDCLASSIFY(DATA, LABELS) returns the (class-mean
%   removed) inverse of the covariance matrix and the class-based means.

% Version:  v1.1
% Build:    16032921
% Date:     Mar-29 2016, 9:48 PM EST
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
if any(isinf(labels) | isnan(labels))
    error('neuroelf:general:badArgument', 'Invalid labels.');
end

% some numbers we'll need, number of samples, features, and labels
nc = numel(conds);
nf = size(data, 2);
nl = numel(labels);

% storage for information, indices for each label, number of samples/label
li = cell(1, nl);
ln = zeros(1, nl);

% label-based means and total (class-weighted) mean
lmean = zeros(nl, nf);

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
    tw = tw + lws(lc);
    cw(li{lc}) = lw{lc};
end

% compute re-weighted (co)variance (of mean-removed samples)
icov = (1 / tw) .* tcov;

% compute inverse (safely)
try
    uri = ((tw - nl) >= nf);
    if uri
        icov = inv(tcov);
    else
        icov = pinv(tcov);
    end
catch
    try
        if uri
            icov = pinv(tcov);
        end
    catch
        return;
    end
end
