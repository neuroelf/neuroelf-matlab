function tsm = srf_SphereToSphereMapTri(xo, xo2, opts)
% SRF::SphereToSphereMapTri  - triangle-based Sphere-To-Sphere mapping
%
% FORMAT:       tsm = srf.SphereToSphereMapTri(srf2 [, opts])
%
% Input fields:
%
%       srf         target sphere
%       srf2        source sphere
%       opts        optional settings
%        .comp      try compiled version
%        .fwithin   force within triangle
%
% Output fields:
%
%       tsm         TSM file with mapping target <- sphere
%
% Using: findfirst, mesh_trianglestoneighbors, mesh_trimapmesh, normvecs,
%        spherecoords.

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:46 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% neuroelf library
global ne_methods;

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf') || ...
    numel(xo2) ~= 1 || ~xffisobject(xo2, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% options
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'comp') || ~islogical(opts.comp) || numel(opts.comp) ~= 1
    opts.comp = true;
end
if ~isfield(opts, 'fwithin') || ~islogical(opts.fwithin) || numel(opts.fwithin) ~= 1
    opts.fwithin = false;
end
fwi = opts.fwithin;

% get contents
bc1 = xo.C;
bc2 = xo2.C;
c1 = bc1.VertexCoordinate;
c2 = bc2.VertexCoordinate;
nc1 = size(c1, 1);
nc2 = size(c2, 1);
t2 = bc2.TriangleVertex;
[n2, bn2, triref2] = ne_methods.mesh_trianglestoneighbors(nc2, t2);

% compiled version
if opts.comp
    [t, vl] = ne_methods.mesh_trimapmesh( ...
        ne_methods.normvecs(c1 - repmat(bc1.MeshCenter, [nc1, 1])), ...
        ne_methods.normvecs(c2 - repmat(bc2.MeshCenter, [nc2, 1])), t2, triref2);
    ntf = find(t == 0);
    if isempty(ntf)
        tsm = xff('new:tsm');
        tsc = tsm.C;
        tsc.NrOfTargetVertices = nc1;
        tsc.NrOfSourceVertices = nc2;
        tsc.NrOfSourceTriangles = size(t2, 1);
        tsc.SourceTriangleOfTarget = t;
        tsc.TriangleEdgeLengths = vl;
        tsm.C = tsc;
        return;
    end
end

% compute a few things first
c1(ntf, 1) = c1(ntf, 1) - (min(c1(:, 1)) + max(c1(:, 1))) / 2;
c1(ntf, 2) = c1(ntf, 2) - (min(c1(:, 2)) + max(c1(:, 2))) / 2;
c1(ntf, 3) = c1(ntf, 3) - (min(c1(:, 3)) + max(c1(:, 3))) / 2;
c1 = c1(ntf, :);
c2(:, 1) = c2(:, 1) - (min(c2(:, 1)) + max(c2(:, 1))) / 2;
c2(:, 2) = c2(:, 2) - (min(c2(:, 2)) + max(c2(:, 2))) / 2;
c2(:, 3) = c2(:, 3) - (min(c2(:, 3)) + max(c2(:, 3))) / 2;
r1 = sqrt(sum(c1 .* c1, 2));
r2 = sqrt(sum(c2 .* c2, 2));
sp1 = ne_methods.spherecoords(c1);
sp2 = ne_methods.spherecoords(c2);

% put the data into quad areas, first estimate a good area size
qs = [180, 120, 90, 60, 48, 36, 30, 24, 20, 16, 12, 10, 8, 6, 4, 2];
qs = qs(ne_methods.findfirst(floor(nc2 ./ (qs .* qs)) >= 16));

% if even 2 is too much (only 4!)
if isempty(qs)
    error('neuroelf:xff:badArgument', 'Too few vertices for matching algorithm.');
end

% get some interim numbers
qsh = qs / 2;
qss = qs * qs;
qsm = qs - 1;
qssm = qss - qs;

% and now put the data in there
q1 = 1 + min(qssm, qs * floor(qsh + qsh * cos(sp1(:, 2)))) + ...
    min(qsm, qsh + floor(qsh * sp1(:, 3) / pi));
q2 = 1 + min(qssm, qs * floor(qsh + qsh * cos(sp2(:, 2)))) + ...
    min(qsm, qsh + floor(qsh * sp2(:, 3) / pi));
[q2s, q2i] = sort(q2);
q2d = [find([1, diff(q2s(:)')]), numel(q2) + 1];
if numel(q2d) ~= (qss + 1)
    if qs > 2
        warning('neuroelf:xff:internal', ...
            'Vertices are highly unevenly spread, matching slowed down.');
        qs = qs / 2;
        qsh = qs / 2;
        qss = qs * qs;
        qsm = qs - 1;
        qssm = qss - qs;

        % and now put the data in there
        q1 = 1 + min(qssm, qs * floor(qsh + qsh * cos(sp1(:, 2)))) + ...
            min(qsm, qsh + floor(qsh * sp1(:, 3) / pi));
        q2 = 1 + min(qssm, qs * floor(qsh + qsh * cos(sp2(:, 2)))) + ...
            min(qsm, qsh + floor(qsh * sp2(:, 3) / pi));
        [q2s, q2i] = sort(q2);
        q2d = [find([1, diff(q2s(:)')]), numel(q2) + 1];
    end
    if numel(q2d) ~= (qss + 1)
        error('neuroelf:xff:badArgument', ...
            'Emtpy patch(es) in spherical surface, distribution too uneven.');
    end
end
q2c = cell(1, qss);
for qc = 1:qss
    q2c{qc} = q2i(q2d(qc):q2d(qc+1)-1);
end

% check spherical shape
if (mean(r1) / std(r1)) < 60 || (mean(r2) / std(r2)) < 60
    error('neuroelf:xff:badArgument', 'Shapes not spherical enough for matching.');
end

% make spherical with r = 1
c1 = c1 ./ r1(:, [1, 1, 1]);
c2 = c2 ./ r2(:, [1, 1, 1]);

% get neighbors list
tf = nc1;
n2 = n2(:, 2);

% init found array
mf = zeros(tf, 1);
nf = zeros(tf, 1);
nl = zeros(tf, 2);

% iterate until all vertices mapped
for m1 = 1:tf

    % get coordinate of target vertex
    cm1 = c1(m1, :);

    % get the vertices in the same quad area
    m2x = q2c{q1(m1)};

    % compute distance
    dst = sqrt(sum((cm1(ones(1, numel(m2x)), :) - c2(m2x, :)) .^ 2, 2));

    % use the one with the shortest distance
    [mval, mpos] = min(dst);
    m2 = m2x(mpos);

    % iterate until match for next vertex found
    while true

        % compute distance between m1 and m2 incl. neighbors
        m2x = cat(2, n2{n2{m2}});
        dst = sqrt(sum((cm1(ones(1, numel(m2x)), :) - c2(m2x, :)) .^ 2, 2));

        % find shortest distance
        [mval, mpos] = min(dst);

        % break if at candidate
        if m2x(mpos) == m2
            break;
        end

        % otherwise set to shortest element
        m2 = m2x(mpos);
    end

    % get directional vectors
    cm2 = c2(m2, :);
    tbn = n2{m2};
    tbs = numel(tbn);
    cvec = c2([tbn, tbn(1)], :) - cm2(ones(1, tbs + 1), :);
    cvcs = sum(cvec .* cvec, 2);
    cvcp = sum(cvec(1:end-1, :) .* cvec(2:end, :), 2);
    tvec = cm1 - cm2;

    % try to find two vectors best fitting target
    tuv = zeros(2, tbs);
    for nc = 1:tbs
        tuv(:, nc) = (1 / (cvcs(nc) * cvcs(nc+1) - cvcp(nc) * cvcp(nc))) * ...
            [cvcs(nc+1), -cvcp(nc); -cvcp(nc), cvcs(nc)] * cvec(nc:nc+1, :) * tvec';
    end

    % find entry/ies with all >= 0
    tnc = find(all(tuv >= 0));

    % if only one match and sum of that < 1
    if numel(tnc) == 1 && sum(tuv(:, tnc)) < 1

        % store and continue
        mf(m1) = m2;
        nf(m1) = tnc;
        nl(m1, :) = tuv(:, tnc)';
        continue;
    end

    % quickly test whether there is an entry with a value just below zero
    tnc = find(all(tuv >= -0.001));
    if numel(tnc) == 1 && sum(tuv(:, tnc)) < 1

        % store and continue
        mf(m1) = m2;
        nf(m1) = tnc;
        nl(m1, :) = tuv(:, tnc)';
        continue;
    end

    % check the wider neighborhood for a better match
    q2xy = q2(m2) - 1;
    qsx = mod(q2xy, qs);
    qsy = (q2xy - qsx) / qs;
    qsx = mod(qs - 1 + [qsx - 1, qsx, qsx + 1], qs);
    qsy = mod(qs - 1 + [qsy - 1, qsy, qsy + 1]', qs);
    qsxy = 1 + qsx([1, 1, 1], :) + qs * qsy(:, [1, 1, 1]);
    m2x = cat(1, q2c{qsxy(:)});

    % compute distances
    dst = sqrt(sum((cm1(ones(1, numel(m2x)), :) - c2(m2x, :)) .^ 2, 2));

    % remove all where distance is too great
    m2x = m2x(dst <= (3 * mval));
    m2xn = numel(m2x);
    if m2xn > 0
        ntun = cell(1, m2xn);
        ntuv = cell(1, m2xn);
        for m2c = 1:m2xn
            stbn = [n2{m2x(m2c)} n2{m2x(m2c)}(1)];
            ntun{m2c} = stbn;
            stbs = numel(stbn);
            cvec = c2(stbn, :) - c2(m2x(m2c * ones(1, stbs)), :);
            tvec = c1(m1, :) - c2(m2x(m2c), :);
            stbs = stbs - 1;
            ntuv{m2c} = zeros(2, stbs);
            for nc = 1:stbs
                cvm = cvec(nc:nc+1, :);
                cvo = cvm * cvm';
                mpi = (1 / (cvo(1) * cvo(4) - cvo(2) * cvo(3))) * [cvo(4), -cvo(2); -cvo(3), cvo(1)];
                ntuv{m2c}(:, nc) = mpi * cvm * tvec';
            end
        end
        for m2c = m2xn:-1:1
            stnc = find(all(ntuv{m2c} >= 0, 1) & sum(ntuv{m2c}, 1) <= (2 / 3));
            if numel(stnc) ~= 1
                ntuv(m2c) = [];
                ntun(m2c) = [];
                m2x(m2c) = [];
            else
                ntuv{m2c} = ntuv{m2c}(:, stnc);
                ntun{m2c} = ntun{m2c}(stnc);
            end
        end
        if ~isempty(ntuv)
            mntusum = sum(ntuv{1});
            mntupos = 1;
            for m2c = 2:numel(m2x)
                if sum(ntuv{m2c}) < mntusum
                    mntusum = sum(ntuv{m2c});
                    mntupos = m2c;
                end
            end
            mf(m1) = m2x(mntupos);
            nf(m1) = find(n2{m2x(mntupos)} == ntun{mntupos});
            nl(m1, :) = ntuv{mntupos}';
            continue;
        end
    end

    tnc = find(all(tuv >= -0.001, 1) & sum(abs(tuv), 1) < 3);
    if isempty(tnc)
        tuv = [0; 0];
        tnc = 1;
    end
    [tuvm1, tuvm2] = min(sum(tuv(:, tnc)));
    nc = tnc(tuvm2);

    % force into triangle ?
    if fwi
        tuv(:, nc) = tuv(:, nc) ./ sum(abs(tuv(:, nc)));
    end

    % put into arrays
    mf(m1) = m2;
    nf(m1) = nc;
    nl(m1, :) = tuv(:, nc)';
end

% create TSM object
tsm = xff('new:tsm');
tsc = tsm.C;
tsc.NrOfTargetVertices = tf;
tsc.NrOfSourceVertices = size(c2, 1);
tsc.SourceOfTarget = mf;
tsc.NeighborVectors = nf;
tsc.NeighborLengths = nl;
tsm.C = tsc;
