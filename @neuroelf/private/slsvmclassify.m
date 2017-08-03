function [tacc, cacc, sacc, asl, acv] = slsvmclassify(data, conds, opts)
%SLSVMCLASSIFY  SVM-based classification of searchlight data
%   TACC = SLSVMCLASSIFY(DATA, CONDS) runs a classification with the matrix
%   in DATA, which must be samples-by-features, and CONDS must specify the
%   binary condition labels with values 1 and 2. The default is to run
%   a series of leave-one-out classifications (one of each condition), if
%   the number of combinations exceeds 5000, 5000 random pairs will be
%   picked. TACC then contains both the average accuracy (for the pair).
%
%   [TACC, CACC, SACC] = SLCSVMLASSIFY(DATA, CONDS) also returns the mean
%   average accuracy per class and per sample.
%
%   [TACC, CACC, SACC, SI, DV] = SLCSVMLASSIFY(DATA, CONDS) also returns
%   the sample indices and (label-corrected) decision values (> 0 correct)
%   for each pairing.
%
%   TACC = SLSVMCLASSIFY(DATA, CONDS, OPTS) to set optional parameters,
%   whereas the following fields can be set in opts as a 1x1 struct:
%    .niter  - 1x1 double, number of iterations (default: 5000)
%    .pairs  - 1x1 logical, perform leave-one-pair out classification
%    .popts  - 1x1 struct, ne_svmpredict options (default: opts)
%    .topts  - 1x1 struct, ne_svmtrain options (default: opts)

% Version:  v1.1
% Build:    16033112
% Date:     Mar-31 2016, 12:24 PM EST
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
labels = unique(conds);
if any(isinf(labels) | isnan(labels)) || numel(labels) < 2
    error('neuroelf:general:badArgument', 'Less than two or invalid labels.');
end
nc = numel(conds);
nl = numel(labels);
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'niter') || ~isa(opts.niter, 'double') || numel(opts.niter) ~= 1 || ...
    isinf(opts.niter) || isnan(opts.niter) || opts.niter < 20
    niter = 5000;
else
    niter = min(1000000, round(opts.niter));
end
if ~isfield(opts, 'pairs') || ~islogical(opts.pairs) || numel(opts.pairs) ~= 1
    opts.pairs = true;
end
if ~isfield(opts, 'popts') || ~isstruct(opts.popts) || numel(opts.popts) ~= 1
    popts = opts;
else
    popts = opts.popts;
end
if ~isfield(opts, 'topts') || ~isstruct(opts.topts) || numel(opts.topts) ~= 1
    topts = opts;
else
    topts = opts.topts;
end

% pair-wise
if opts.pairs

    % find indices
    li = cell(1, nl);
    for lc = 1:nl
        li{lc} = lsqueeze(find(conds == labels(lc)));
    end
    ln = cellfun('prodofsize', li);
    if any(ln < 2)
        error('neuroelf:general:badArgument', 'Each class much at least have two samples.');
    end
    
    % maximum number of pairs
    tpairs = sum(prod(pairs(ln(:)), 3));

    % complete test
    if tpairs <= niter

        % create accuracy data
        acc = zeros(tpairs, 1);
        acl = zeros(tpairs, 2);
        acv = zeros(tpairs, 2);
        asa = zeros(tpairs, 2);
        asl = zeros(tpairs, 2);

        % select from classes
        iter = 1;
        for c1c = 1:(nl-1)
            c1i = true(1, ln(c1c));
            l1 = ones(ln(c1c) - 1, 1);
            for c1n = 1:ln(c1c)
                c1i(c1n) = false;
                for c2c = (c1c+1):nl
                    c2i = true(1, ln(c2c));
                    l2 = 2 .* ones(ln(c2c) - 1, 1);
                    l12 = [l1; l2];
                    for c2n = 1:ln(c2c)
                        c2i(c2n) = false;

                        % train model
                        model = ne_svmtrain(l12, data([li{c1c}(c1i); li{c2c}(c2i)], :), topts);

                        % predict
                        lix = [li{c1c}(c1n), li{c2c}(c2n)];
                        [mp1, mp2, mp3] = ne_svmpredict([1; 2], ...
                            data(lix, :), model, popts);

                        % store values
                        acc(iter) = mp2(1);
                        acl(iter, :) = conds(lix);
                        acv(iter, :) = [1; -1] .* mp3;
                        asa(iter, :) = double(mp1 == [1; 2]);
                        asl(iter, :) = lix';

                        % increase counter
                        iter = iter + 1;

                        % reset index
                        c2i(c2n) = true;
                    end
                end
                c1i(c1n) = true;
            end
        end

    % random sampling
    else

        % create accuracy data
        acc = zeros(niter, 1);
        acl = zeros(niter, 2);
        acv = zeros(niter, 2);
        asa = zeros(niter, 2);
        asl = zeros(niter, 2);

        % iterate
        iter = 1;
        while iter <= niter

            % pick a random pair
            rorder = randperm(nc);
            cr1 = conds(rorder(1));
            cr2 = conds(rorder(2));
            if cr1 == cr2
                continue;
            end

            % train with the rest of those two classes
            tval = true(nc, 1);
            i12 = rorder(1:2);
            tval(i12) = false;
            clab = [cr1; cr2];
            tidx = [li{cr1}(tval(li{cr1})); li{cr2}(tval(li{cr2}))];
            l12 = conds(tidx);
            model = ne_svmtrain(l12, data(tidx, :), topts);

            % predict
            [mp1, mp2, mp3] = ne_svmpredict(clab, data(i12, :), model, popts);

            % store values
            acc(iter) = mp2(1);
            acl(iter, :) = clab;
            acv(iter, :) = [1; -1] .* mp3;
            asa(iter, :) = double(mp1 == clab);
            asl(iter, :) = i12;

            % increase counter
            iter = iter + 1;

        end
    end
    
% k-folds
else
    error('neuroelf:general:notYetImplemented', 'Not yet implemented.');
end

% scale
tacc = 0.01 * mean(acc);

% additional outputs
if nargout > 1
    cacc = histcountn(acl(:), 1, nl + 1, 1, asa(:)) ./ histcountn(acl(:), 1, nl + 1, 1);
    if nargout > 2
        sacc = histcountn(asl(:), 1, nc + 1, 1, asa(:)) ./ histcountn(asl(:), 1, nc + 1, 1);
    end
end
