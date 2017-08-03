function [cs, cv, cl, cc] = splitclustercoords(c, v, k, d)
% splitclustercoords  - split coords of one cluster to subclusters
%
% FORMAT:       [cs, cv, cl, cc] = splitclustercoords(c, v [, k [, d]])
%
% Input fields:
%
%       c           Vx3 coordinates of values
%       v           Vx1 values (at coordinates)
%       k           size threshold for sub-clusters (default: 3)
%       d           optional 1x3 minimum distance [default: [2, 2, 2]]
%
% Output fields:
%
%       cs          list of cluster sizes
%       cv          Cx1 cell array with lists of values
%       cl          Vx4 lister of cluster voxels
%       cc          Cx1 cell array with lists of coordinates
%
% Note: the distance (d) must be given in voxel units (given that the
%       function is ignorant of the spatial resolution of the data).

% Version:  v1.1
% Build:    16030622
% Date:     Mar-06 2016, 10:59 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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
if nargin < 2 || ...
   ~isa(c, 'double') || ...
    ndims(c) ~= 2 || ...
    size(c, 2) ~= 3 || ...
    size(c, 1) < 4 || ...
    any(isinf(c(:)) | isnan(c(:)) | c(:) < 1 | c(:) > 2048 | c(:) ~= round(c(:))) || ...
    prod(max(c)) > 2e8 || ...
   ~isa(v, 'double') || ...
    numel(v) ~= size(c, 1) || ...
    numel(v) ~= max(size(v)) || ...
    any(isinf(v) | isnan(v))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
nc = size(c, 1);

% linearize given values
v = v(:);

% make sure coordinates are as low as possible (for clustering volume)
mnc = min(c) - 1;
if any(mnc > 0)
    c = c - mnc(ones(1, nc), :);
end
mxc = max(2, max(c));

% get minimum value and patch values to be >= 1
mnv = min(v);
mv = (1 - mnv) + v;

% create volume filled with a value that is definitely lower than anything
clvn = -(abs(mnv) * 2 + nc);
clv = zeros(mxc) + clvn;

% fill with patched values
clv(sub2ind(mxc, c(:, 1), c(:, 2), c(:, 3))) = mv;

% check there is only one cluster
numclusters = clustercoordsc(clv > 0, 3, 1);
if numel(numclusters) > 1
    error( ...
        'neuroelf:InternalError', ...
        'Subclusters only valid for single cluster.' ...
    );
end

% what threshold for subclusters
if nargin < 3 || ...
   ~isa(k, 'double') || ...
    numel(k) ~= 1 || ...
    isinf(k) || ...
    isnan(k) || ...
    k < 1
    k = 3;
elseif k >= (nc ^ (2 / 3))
    k = round(nc ^ (2 / 3));
else
    k = round(k);
end

% what distance %% TO BE IMPLEMENTED %%
if nargin < 4 || ...
   ~isa(d, 'double') || ...
    numel(d) ~= 3 || ...
    any(isinf(d(:)') | isnan(d(:)') | d(:)' < 1.5 | d(:)' > (mxc / 2))
    d = [2, 2, 2];
end
dmid = 1 + ceil(d);
dsiz = dmid + ceil(d);
dmsk = false(dsiz);
for xc = 1:dsiz(1)
    for yc = 1:dsiz(2)
        for zc =1:dsiz(3)
            dd = ([xc, yc, zc] - dmid) ./ d;
            if sum(dd .* dd) < 1
                dmsk(xc, yc, zc) = true;
            end
        end
    end
end
[dmsi{1:3}] = ind2sub(dsiz, find(dmsk(:)));
dmsi{1} = dmsi{1} - dmid(1);
dmsi{2} = dmsi{2} - dmid(2);
dmsi{3} = dmsi{3} - dmid(3);

% prepare other arrays
ci = zeros(1, nc);

% sort values and coordinates
[mv, mi] = sort(mv, 'descend');
v = v(mi);
c = c(mi, :);

% set max peak + distance to -1
ci(1) = 1;
cs = zeros(1, nc);
tcrd = validcoords( ...
    [dmsi{1} + c(1, 1), dmsi{2} + c(1, 2), dmsi{3} + c(1, 3)], mxc);
clv(tcrd) = -1;

% create peaks matrix (assuming there cannot be more than ceil(nc/2))
peaks = zeros(ceil(0.5 * nc), 3);
peaks(1, :) = c(1, :);

% iterate over remaining values and
numclus = 1;
for rc = 2:nc

    % get coordinate
    dc = c(rc, :);

    % check for other clusters
    df = max(dc - 1, 1);
    dt = min(dc + 1, mxc);
    cltouch = lsqueeze(clv(df(1):dt(1), df(2):dt(2), df(3):dt(3)));
    cltouch = cltouch(cltouch > clvn & cltouch < 0);

    % touches cluster -> mark to one with closest peak
    if ~isempty(cltouch)

        % peaks of touching clusters
        cltouch = unique(-cltouch);
        if numel(cltouch) > 1
            tpeaks = peaks(cltouch, :);
            clid = cltouch(minpos(sum((tpeaks - dc(ones(1, numel(cltouch)), :)) .^ 2, 2)));
        else
            clid = cltouch;
        end
        ci(rc) = clid;
        clv(dc(1), dc(2), dc(3)) = -clid;
        cs(clid) = cs(clid) + 1;

    % does not touch
    else

        % increase number of identified sub-clusters
        numclus = numclus + 1;
        cs(numclus) = 1;
        tcrd = validcoords( ...
            [dmsi{1} + dc(1), dmsi{2} + dc(2), dmsi{3} + dc(3)], mxc);
        clv(tcrd(clv(tcrd) >= 0)) = -numclus;
        ci(rc) = numclus;

        % record peak
        peaks(numclus, :) = dc;
    end
end
cs = cs(1:numclus);

% apply sub-clustering k-threshold
for rc = numclus:-1:1

    % smaller ?
    if cs(rc) < k

        % for each voxel find adjacent cluster with highest number
        cci = find(ci == rc);
        ccoords = c(cci, :);
        for scc = cs(rc):-1:1
            % get coordinate
            dc = ccoords(scc, :);

            % find other clusters
            df = max(dc - 1, 1);
            dt = min(dc + 1, mxc);
            cltouch = lsqueeze(clv(df(1):dt(1), df(2):dt(2), df(3):dt(3)));
            cltouch(cltouch == -rc | cltouch < -nc | cltouch >= 0) = [];
            if isempty(cltouch)
                clid = 0;
            else
                clid = -max(cltouch(cltouch < 0));
                if clid <= numel(cs)
                    cs(clid) = cs(clid) + 1;
                else
                    clid = 0;
                end
            end

            % reasign
            clv(dc(1), dc(2), dc(3)) = -clid;
            ci(cci(scc)) = clid;
        end

        % then remove cluster
        cs(rc) = [];
        numclus = numclus - 1;
    end
end

% renumber clusters
cn = find(ci == 0);
c(cn, :) = [];
ci(cn) = [];
v(cn) = [];
mi(cn) = [];
clnum = unique(ci);
clnid = zeros(1, max(clnum));
clnid(clnum) = 1:numel(clnum);
ci = clnid(ci);

% add coordinate
if any(mnc > 0)
    c = c + mnc(ones(1, size(c, 1)), :);
end

% create output
cv = cell(numclus, 1);
cl = [c, ci(:)];
cl(:, end) = ci(:);
cc = cell(numclus, 1);
for rc = 1:numclus
    cv{rc} = lsqueeze(v(ci == rc));
    cc{rc} = c(ci == rc, :);
end

% re-sort cluster coordinate to original order
[mv, mi] = sort(mi);
cl = cl(mi, :);


%%% sub functions

% validcoords
function c = validcoords(c, m)
c(any(c < 1, 2) | any(c > m(ones(1, size(c, 1)), :), 2), :) = [];
c = sub2ind(m, c(:, 1), c(:, 2), c(:, 3));
