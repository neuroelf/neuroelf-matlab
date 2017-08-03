% mkda_plp_cluster  - attempt coordinate based clustering using K-means

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:28 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2013, 2016, Jochen Weber
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

% maximum number of clusters, number of iterations to find clusters
maxclust = 50;
stabiter = 50;
lkmrepls = 25;

% load PLP and VMP
oplp = xff('*.plp');
vmp = xff('*.vmp');

% apply mask and selection
plp = oplp.ReduceToMKDAMask(vmp);

% remove other objects
oplp.ClearObject;
vmp.ClearObject;

% get points
p = plp.Points(:, 1:3);
np = size(p, 1);

% iterate to find clustering whereas 95% of points are stable
stab = zeros(maxclust, 1);
stab(1) = 1;
stac = cell(1, maxclust);
stac{1} = ones(np, 1);

% apply clustering
for cc = 2:maxclust

    % show stability estimate for previous clustering
    dispp(sprintf(' -> stability for %d classes = %.7g', cc - 1, stab(cc - 1)));

    % stability iterations
    stc = zeros(np, stabiter);
    for sc = 1:stabiter

        % apply litekmeans clustering
        stc(:, sc) = litekmeans(p, cc, 'Replicates', lkmrepls);
    end

    % figure out stability measure
    stat = zeros(np, 1);
    statx = zeros(round(0.5 * stabiter * (stabiter - 1)), 3);
    staty = 1;
    statd = 0;
    statt = 0;
    for sc = 1:stabiter
        for sc2 = (sc+1):stabiter

            % get unique combination of classes
            ucc = unique(stc(:, [sc, sc2]), 'rows');

            % for each class
            for tcc = 1:maxclust

                % which rows in stat
                stati = find(stc(:, sc) == tcc);

                % find possible matches
                mcs = ucc(ucc(:, 1) == tcc, 2);

                % more than one match
                if numel(mcs) > 1

                    % find highest portion
                    mcn = mcs;
                    mcp = stc(stati, sc2);
                    for mcc = 1:numel(mcs)
                        mcn(mcc) = sum(mcp == mcs(mcc));
                    end

                    % stability is fractional (only if second class >= 0.5)
                    if any(mcn >= (0.5 * numel(stati)))

                        % which second class
                        mcm = find(mcn >= (0.5 * numel(stati)));
                        stata = 1 / numel(stati);
                        stati(stc(stati, sc2) ~= mcs(mcm)) = [];
                        stat(stati) = stat(stati) + stata * numel(stati);
                    end
                else

                    % stability is 1 for this class
                    stat(stati) = stat(stati) + 1;
                end
            end
            statd = statd + 1;

            % store difference
            statm = mean(stat);
            statx(staty, :) = [statm - statt, sc, sc2];
            statt = statm;
            staty = staty + 1;
        end
    end
    stat = stat ./ statd;

    % store stability
    stab(cc) = robustmean(stat);

    % try to determine classes
    if stab(cc) > 0.95

        % find pairs that were "best" (most stable between them)
        stats = sort(statx(:, 1));
        goodp = find(statx(:, 1) >= stats(1 + floor(0.99 * numel(stats))));

        % best pairing and unique index
        bestp = find(statx(:, 1) == stats(end));
        bestp = stc(:, statx(bestp(1), 2:3));
        bestc = bestp(:, 1);
        bestm = bestc + (maxclust + 1) .* bestp(:, 2);

        % absolute match
        if numel(unique(bestm)) == maxclust

            % store and continue
            stac{cc} = bestc;
            continue;
        end

        besth = histcount(bestm, 1, max(bestm));
        besti = find(besth > 0);
        for bc = 1:numel(besti)
            bestj = besti(bc);
            besth(bestj) = besth(bestj) / harmmean([ ...
                sum(bestp(:, 1) == mod(bestj, maxclust + 1)); ...
                sum(bestp(:, 2) == floor(bestj / (maxclust + 1)))]);
        end

        % all pairs
        if numel(goodp) == numel(stats)

            % it doesn't matter!
            stac{cc} = stc(:, 1);
            continue;

        % only one good pair
        elseif numel(goodp) == 1

            % choose the classing from that pair
            stc = stc(:, statx(goodp, 2:3));

        % determine from good pairs
        else
            goods = unique(lsqueeze(statx(goodp, 2:3)));

            % get classes from those
            stc = stc(:, goods);
        end

        % and unique splits
        ustc = unique(stc, 'rows');

        % for each "pairing", determine frequency (absolute and relative)
        pfreq = zeros(size(ustc, 1), 1);
        rfreq = zeros(size(ustc));
        for rc = 1:numel(pfreq)
            for ic = 1:size(stc, 2)
                rfreq(rc, ic) = sum(stc(:, ic) == ustc(rc, ic));
            end
            freqhm = harmmean(rfreq(rc, :), 2);
            pfreq(rc) = sum(all(stc == ones(np, 1) * ustc(rc, :), 2)) / freqhm;
            rfreq(rc, :) = rfreq(rc, :) ./ freqhm;
        end

        % determine classes
        stacs = zeros(np, 1);

        % for each point
        for pc = 1:np

            % match in best pair
            if besth(bestm(pc)) >= 0.9
                stacs(pc) = bestc(pc);
                continue;
            end
        end

        % all filled?
        if all(stacs > 0)
            continue;
        end

        % for each class
        for tcc = 1:maxclust

            % which rows in stat
            stati = find(stc(:, 1) == tcc);

            % find possible matches
            mcs = ucc(ucc(:, 1) == tcc, 2);

            % more than one match
            if numel(mcs) > 1

                % find highest portion
                mcn = mcs;
                mcp = stc(stati, sc2);
                for mcc = 1:numel(mcs)
                    mcn(mcc) = sum(mcp == mcs(mcc));
                end

                % stability is fractional (only if second class >= 0.5)
                if any(mcn >= (0.5 * numel(stati)))

                    % which second class
                    mcm = find(mcn >= (0.5 * numel(stati)));
                    stata = 1 / numel(stati);
                    stati(stc(stati, sc2) ~= mcs(mcm)) = [];
                    stat(stati) = stat(stati) + stata * numel(stati);
                end
            else

                % stability is 1 for this class
                stat(stati) = stat(stati) + 1;
            end
        end
    end
end
figure;
axes;
hold on;
s = cell(max(lkc), 1);
for sc = 1:numel(s)
    s{sc} = scatter3(p(lkc == sc, 1), p(lkc == sc, 2), p(lkc == sc, 3));
end
