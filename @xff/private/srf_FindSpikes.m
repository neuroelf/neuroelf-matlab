function [spikes, npmin, vdmin, vdmean] = srf_FindSpikes(xo, cutoff)
% SRF::FindSpikes  - find spiked vertices (1-based indices)
%
% FORMAT:       [spikes, npmin, vdmean] = srf.FindSpikes([cutoff]);
%
% Input fields:
%
%       cutoff      cutoff value for critical sqrt(abs(n*n')),
%                   default: 1/3
%
% Output fields:
%
%       spikes      Vx1 vertex list of found spikes
%       npmin       Vx1 minimal product list per vertex
%       vdmin       Vx1 min vertex-to-neighbors distance
%       vdmean      Vx1 mean vertex-to-neighbors distance
%
% Note: the srf object remains unchanged, often it is more useful
%       to apply e.g. coloring to an unsmoothed mesh

% Version:  v1.1
% Build:    16021119
% Date:     Feb-11 2016, 7:49 PM EST
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

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isa(cutoff, 'double') || numel(cutoff) ~= 1 || ...
    isinf(cutoff) || isnan(cutoff) || cutoff <= 0 || cutoff >= 1
    cutoff = 1/3;
else
    cutoff = real(cutoff);
end

% get vertices, normals, and neighbors
vx = bc.VertexCoordinate;
no = bc.VertexNormal;
nb = bc.Neighbors;
nv = size(no, 1);

% initialize array
npmin = ones(nv, 1);
vdmin = zeros(nv, 1);
vdmean = zeros(nv, 1);

% iterate over vertex normals
for nc = 1:nv

    % get minimal product of this vertex' neighbors with normal
    tnb = nb{nc, 2};
    npmin(nc) = min(sqrt(abs(no(nc, :) * no(tnb, :)')));

    % calc mean vertex-neighbor distance
    vd = sqrt(sum((repmat(vx(nc, :), [numel(tnb), 1]) - vx(tnb, :)) .^ 2, 2));
    vdmin(nc) = min(vd);
    vdmean(nc) = mean(vd);
end

% find spikes
spikes = find((npmin .* (vdmean / mean(vdmean)) .* (vdmin / mean(vdmin))) <= cutoff);
