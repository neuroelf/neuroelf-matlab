function [neis, neimtx, neidist] = srf_NeighborsNDegree(xo, nd, opts)
% SRF::NeighborsNDegree  - retrieve neighbors up to a degree
%
% FORMAT:       [nei, neimtx] = srf.NeighborsNDegree(nd [, opts])
%
% Input fields:
%
%       nd          1x1 degree to which neighbors are retrieved
%       opts        optional settings
%        .maxdist   maximum distance (in mm, default: Inf)
%        .notself   do not include vertex index itself (default: false)
%
% Output fields:
%
%       nei         cell array with list of neighbors
%       neimtx      matrix version (for faster access)

% Version:  v1.1
% Build:    16040417
% Date:     Apr-04 2016, 5:59 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2013 - 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid object for SRF::NeighborsNDegree');
end
bc = xo.C;
if nargin < 2 || ~isa(nd, 'double') || numel(nd) ~= 1 || isinf(nd) || isnan(nd) || nd < 1
    nd = 1;
else
    nd = min(round(nd), 20);
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'maxdist') || ~isa(opts.maxdist, 'double') || numel(opts.maxdist) ~= 1 || ...
    isnan(opts.maxdist) || opts.maxdist <= 0
    opts.maxdist = Inf;
end
if ~isfield(opts, 'notself') || ~islogical(opts.notself) || numel(opts.notself) ~= 1
    opts.notself = false;
end

% neighbors
crd = bc.VertexCoordinate;
oneis = bc.Neighbors;
neis = oneis;
nnei = size(neis, 1);
if nargout > 2
    neidist = neis(:, 2);
    for nc = 1:nnei
        neidist{nc} = sqrt(sum((crd(nc .* ones(neis{nc, 1}, 1), :) - crd(neis{nc, 2}, :)) .^ 2, 2))';
    end
end

% iterate (if necessary)
ond = nd;
if nd > 1
    if nargout < 3
        while nd > 1
            for nc = 1:nnei
                neis{nc, 2} = cat(2, oneis{neis{nc, 2}, 2});
            end
            nd = nd - 1;
        end
    else
        while nd > 1
            for nc = 1:nnei
                orgcrd = crd(nc, :);
                newnei = setdiff(cat(2, oneis{neis{nc, 2}, 2}), neis{nc, 2});
                newdist = zeros(1, numel(newnei));
                for nnc = 1:numel(newnei)
                    midcrd = mean(crd(intersect(oneis{newnei(nnc), 2}, neis{nc, 2}), :), 1);
                    newdist(nnc) = sqrt(sum((orgcrd - midcrd) .^ 2, 2)) + ...
                        sqrt(sum((midcrd - crd(newnei(nnc), :)) .^ 2, 2));
                end
                neis{nc, 2} = [neis{nc, 2}, newnei];
                neidist{nc} = [neidist{nc}, newdist];
            end
            nd = nd - 1;
        end
    end
elseif ~opts.notself
    if nargout < 3
        for nc = 1:nnei
            neis(nc, :) = {neis{nc, 1} + 1, [nc, neis{nc, 2}]};
        end
    else
        for nc = 1:nnei
            neis(nc, :) = {neis{nc, 1} + 1, [nc, neis{nc, 2}]};
            neidist{nc} = [0, neidist{nc}];
        end
    end
end

% unique
if ond > 1
    if opts.notself
        if nargout < 3
            for nc = 1:nnei
                neis{nc, 2} = unique(neis{nc, 2});
                neis{nc, 2}(neis{nc, 2} == nc) = [];
                neis{nc, 1} = numel(neis{nc, 2});
            end
        else
            for nc = 1:nnei
                [neis{nc, 2}, ui] = unique(neis{nc, 2});
                selfi = (neis{nc, 2} == nc);
                ui(selfi) = [];
                neis{nc, 2}(selfi) = [];
                [neidist{nc}, ui] = sort(neidist{nc}(ui));
                neis{nc, 2} = neis{nc, 2}(1, ui);
                neis{nc, 1} = numel(neis{nc, 2});
            end
        end
    else
        if nargout < 3
            for nc = 1:nnei
                neis{nc, 2} = unique([nc, neis{nc, 2}]);
                neis{nc, 1} = numel(neis{nc, 2});
            end
        else
            for nc = 1:nnei
                [neis{nc, 2}, ui] = unique(neis{nc, 2});
                [neidist{nc}, ui] = sort(neidist{nc}(ui));
                neis{nc, 2} = neis{nc, 2}(1, ui);
                neis{nc, 1} = numel(neis{nc, 2});
            end
        end
    end
end

% neighbor-to-neighbor matrix
if nargout > 1

    % total number of neighbors
    tnei = sum(cellfun('prodofsize', neis(:, 2)));
    sj = cat(2, neis{:, 2});
    si = zeros(tnei, 1);
    ss = ones(tnei, 1);
    ti = 1;
    for nc = 1:nnei
        si(ti:ti+neis{nc, 1}-1) = nc;
        ti = ti + neis{nc, 1};
    end
    neimtx = sparse(si, sj, ss, nnei, nnei, tnei);
end
