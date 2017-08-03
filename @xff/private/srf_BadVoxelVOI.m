function voi = srf_BadVoxelVOI(xo, opts)
% SRF::BadVoxelVOI  - find bad segmented voxels
%
% FORMAT:       voi = srf.BadVoxelVOI([opts]);
%
% Input fields:
%
%       opts        optional settings
%        .bridges   try to find bridges
%
% Output fields:
%
%       voi         VOI object
%
% Using: mesh_morph, mesh_trianglestoneighbors.

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:01 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
mesh_morph = ne_methods.mesh_morph;

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'bridges') || numel(opts.bridges) ~= 1 || ...
   (~islogical(opts.bridges) && ~isa(opts.bridges, 'double'))
    opts.bridges = false;
else
    if ~islogical(opts.bridges)
        opts.bridges = (opts.bridges ~= 0);
    end
end

% get coordinates
bc = xo.C;
c = bc.VertexCoordinate;
n = bc.Neighbors;
t = bc.TriangleVertex;
nc = size(c, 1);
[an, bn] = ne_methods.mesh_trianglestoneighbors(nc,t);
n2 = an(:, 2);

% build neighbor list
nl1 = zeros(12 * nc, 1);
nl2 = zeros(12 * nc, 1);
nli = 0;
for nlc = 1:nc
    tn = n2{nlc};
    tnn = numel(tn);
    nl1(nli + 1:nli + tnn) = nlc;
    nl2(nli + 1:nli + tnn) = tn;
    nli = nli + tnn;
end
nl1 = nl1(1:nli);
nl2 = nl2(1:nli);

% create voi object
voi = xff('new:voi');
voic = voi.C;
voic.VOI(:) = [];
vc = 1;

% check bad neighbors list first
if ~isempty(bn)

    % loop over bad neighbors
    for bnc = 1:numel(bn)

        % get bad coordinates
        bnl = unique([bn(bnc), n2{bn(bnc)}, n{bn(bnc), 2}]);
        badc = 128 - bc.VertexCoordinate(bnl, :);
        badx = badc(:, 3);
        bady = badc(:, 1);
        badz = badc(:, 2);

        % put into voi
        badv = unique(floor([ ...
            badx    , bady    , badz    ; ...
            badx + 1, bady    , badz    ; ...
            badx    , bady + 1, badz    ; ...
            badx + 1, bady + 1, badz    ; ...
            badx    , bady    , badz + 1; ...
            badx + 1, bady    , badz + 1; ...
            badx    , bady + 1, badz + 1; ...
            badx + 1, bady + 1, badz + 1; ...
            ]), 'rows');

        % fill VOI
        voic.VOI(vc).NrOfVoxels = size(badv, 1);
        voic.VOI(vc).Voxels = badv;
        voic.VOI(vc).Name = sprintf('neighbors_%d', bn(bnc));
        voic.VOI(vc).Color = [255, 0, 0];
        vc = vc + 1;
    end

    % put back
    voic.NrOfVOIs = vc - 1;
    voi.C = voic;
    return;
end

% perform smoothing
c = mesh_morph(c, n, t, struct('force', 0.04, 'niter', 250));

% compute average and minimum vertex distances
md = sqrt(sum((c(nl1, :) - c(nl2, :)) .^ 2, 2));

% find entries that are too close
badi = [];
badf = find(md < 0.01);

% bad coordinates after smoothing
if ~isempty(badf)

    % get bad coordinates
    badi = unique(nl1(badf));
    badc = 128 - bc.VertexCoordinate(badi, :);
    badx = badc(:, 3);
    bady = badc(:, 1);
    badz = badc(:, 2);

    % put into voi
    badv = unique(floor([ ...
        badx    , bady    , badz    ; ...
        badx + 1, bady    , badz    ; ...
        badx    , bady + 1, badz    ; ...
        badx + 1, bady + 1, badz    ; ...
        badx    , bady    , badz + 1; ...
        badx + 1, bady    , badz + 1; ...
        badx    , bady + 1, badz + 1; ...
        badx + 1, bady + 1, badz + 1; ...
        ]), 'rows');
    voic.VOI(1).NrOfVoxels = size(badv, 1);
    voic.VOI.Voxels = badv;

    % put back into VOI and return!
    voic.VOI.Name = 'smoothing';
    voic.VOI.Color = [192, 0, 0];
    voi.C = voic;
    return;
end

% no bad voxels after smoothing, apply some inflation
ns = 0;
niter = 25;
for ic = 1:12
    try
        tc = mesh_morph(c, n, t, struct('force', 0.5, 'niter', niter, 'areac', 1, 'distc', 2.5));
        if any(isnan(tc(:)))
            error('NaNERROR');
        end
        c = tc;
        ns = ns + niter;
    catch xfferror
        if niter == 1
            warning('neuroelf:xff:MEXError', ...
                'Error inflating mesh (after %d pre-steps): %s.', ns, xfferror.message);
            delete(voi);
            return;
        end
        niter = round(niter / 2);
        continue;
    end
    md = sqrt(sum((c(nl1, :) - c(nl2, :)) .^ 2, 2));
    if any(md < 0.001)
        badi = unique(nl1(md < 0.01));
        break;
    end
end

% any bad vertices (left)
while ~isempty(badi)

    % get first bad vertex information
    badv = badi(1);
    chkv = badv;

    % check neighbors
    chkvo = chkv;
    while ~isempty(chkv)
        chkvi = chkv;
        chkv = [];
        for cc = 1:numel(chkvi)
            cn = intersect(badi, n2{chkvi(cc)});
            chkv = [chkv, setdiff(cn, chkvo)];
        end
        chkv = unique(chkv);
        chkvo = union(chkvo, chkv);
    end

    % remove from index list
    badi = setdiff(badi, chkvo);

    % get bad coordinates
    badc = 128 - bc.VertexCoordinate(chkvo, :);
    badx = badc(:, 3);
    bady = badc(:, 1);
    badz = badc(:, 2);

    % put into voi
    badv = unique(floor([ ...
        badx    , bady    , badz    ; ...
        badx + 1, bady    , badz    ; ...
        badx    , bady + 1, badz    ; ...
        badx + 1, bady + 1, badz    ; ...
        badx    , bady    , badz + 1; ...
        badx + 1, bady    , badz + 1; ...
        badx    , bady + 1, badz + 1; ...
        badx + 1, bady + 1, badz + 1; ...
        ]), 'rows');
    voic.VOI(vc).Name = sprintf('Vertex_%d', chkvo(1));
    voic.VOI(vc).NrOfVoxels = size(badv, 1);
    voic.VOI(vc).Voxels = badv;
    voic.VOI(vc).Color = [255, 0, 64];
    vc = vc + 1;
end

% check for bridges
if vc < 2 && opts.bridges
    tc = mesh_morph(c, n, t, struct('force', 0.8, 'niter', 2000, 'areac', 1, 'distc', 2.5));
end

% put back into VOI and return!
voic.NrOfVOIs = vc - 1;
voi.C = voic;
