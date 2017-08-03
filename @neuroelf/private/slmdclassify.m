function [tacc, cacc, sacc] = slmdclassify(data, conds, robflag)
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
% Build:    16032922
% Date:     Mar-29 2016, 10:26 PM EST
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
nc = numel(conds);
nf = size(data, 2);
nl = numel(labels);

% robust
if nargin < 3 || ~islogical(robflag) || numel(robflag) ~= 1
    robflag = false;
end

% compute full means and inverse covs
ficov = cell(1, nl);
fmean = cell(1, nl);
li = cell(1, nl);
ln = zeros(1, nl);
lo = cell(1, nl);
for lc = 1:nl
    li{lc} = find(conds == labels(lc));
    ln(lc) = numel(li{lc});
    lo{lc} = ones(ln(lc), 1);
    rankok = (ln(lc) > (1.4 * nf));
    if robflag
        [ficov{lc}, fmean{lc}] = mdclasswb(data(li{lc}, :), lo{lc});
    else
        rmdata = data(li{lc}, :);
        fmean{lc} = (1 / ln(lc)) .* sum(rmdata, 1);
        rmdata = rmdata - lo{lc} * fmean{lc};
        ficov{lc} = (1 / (ln(lc) - 1)) .* (rmdata' * rmdata);
        if rankok
            try
                ficov{lc} = inv(ficov{lc});
            catch
                try
                    ficov{lc} = pinv(ficov{lc});
                catch
                    tacc = NaN;
                    if nargout > 1
                        cacc = NaN .* ones(nl, 1);
                        if nargout > 2
                            sacc = NaN .* ones(nc, 1);
                        end
                    end
                    return;
                end
            end
        else
            try
                ficov{lc} = pinv(ficov{lc});
            catch
                tacc = NaN;
                if nargout > 1
                    cacc = NaN .* ones(nl, 1);
                    if nargout > 2
                        sacc = NaN .* ones(nc, 1);
                    end
                end
                return;
            end
        end
    end
end

% create overall accuracy output
sacc = zeros(nc, 1);
cacc = zeros(nl, 1);

% for only one pair of labels
for lc = 1:nl

    % accuracy for this class (samples)
    clacc = zeros(ln(lc), nl);
    rankok = (ln(lc) > (1.4 * nf));

    % compute test for "other" class for class 1
    for lc2 = 1:nl
        if lc == lc2
            looi = true(ln(lc), 1);
            lnm1 = ln(lc) - 1;
            loon = ones(lnm1, 1);
            lnm1 = 1 / lnm1;
            rmdata = data(li{lc}, :);
            for sc = 1:ln(lc)
                looi(sc) = false;
                rmidata = rmdata(looi, :);
                rmsdata = rmdata(sc, :);
                if robflag
                    [icovp, mp] = mdclasswb(rmidata, loon);
                else
                    mp = lnm1 .* sum(rmidata, 1);
                    rmidata = rmidata - loon * mp;
                    icovp = lnm1 .* (rmidata' * rmidata);
                    if rankok
                        try
                            icovp = inv(icovp);
                        catch
                            try
                                icovp = pinv(icovp);
                            catch
                                tacc = NaN;
                                if nargout > 1
                                    cacc = NaN .* ones(nl, 1);
                                    if nargout > 2
                                        sacc = NaN .* ones(nc, 1);
                                    end
                                end
                                return;
                            end
                        end
                    else
                        try
                            icovp = pinv(icovp);
                        catch
                            tacc = NaN;
                            if nargout > 1
                                cacc = NaN .* ones(nl, 1);
                                if nargout > 2
                                    sacc = NaN .* ones(nc, 1);
                                end
                            end
                            return;
                        end
                    end
                end
                rmsdata = rmsdata - mp;
                clacc(sc, lc) = sum((rmsdata * icovp) .* rmsdata);
                looi(sc) = true;
            end
        else
            rmdata = data(li{lc}, :) - lo{lc} * fmean{lc2};
            clacc(:, lc2) = sum((rmdata * ficov{lc2}) .* rmdata, 2);
        end
    end

    % compute trial accuracy
    ctacc = double(min(clacc, [], 2) == clacc(:, lc));
    cacc(lc) = (1 / ln(lc)) * sum(ctacc);

    % store
    sacc(li{lc}) = ctacc;
end

% total accuracy
tacc = sum(sacc) / nc;
