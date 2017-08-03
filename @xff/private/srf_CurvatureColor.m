function xo = srf_CurvatureColor(xo, csmooth)
% SRF::CurvatureColor  - color according to curvature
%
% FORMAT:       srf.CurvatuteColor([csmooth])
%
% Input fields:
%
%       csmooth     amount of smoothing (nr of neighbors to mean)
%
% No output fields.

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
if nargin < 2 || ~isnumeric(csmooth) || isempty(csmooth) || isinf(csmooth(1)) || csmooth(1) > 5
    csmooth = 0;
end

% build curvature
p = bc.VertexCoordinate;
n = bc.VertexNormal;
nei = bc.Neighbors;
numv = bc.NrOfVertices;
curve = zeros(numv, 1);

% iterate over vertices
for vc = 1:numv

    % set initial curvature estimate to 0
    cp = 0;

    % get reference point
    rp = p(vc, :);

    % get neighbors
    pcnx = nei{vc, 2};

    % iterate over neighbors
    for pc = pcnx(:)'

        % get vector and normal
        tv = p(pc, :) - rp;
        nv = n(pc, :);

        % add to curvature estimate
        cp = cp + acos((tv * nv') / sqrt((tv * tv') * (nv * nv')));
    end

    % set into array
    curve(vc) = cp;
end

% apply some smoothing
for cc = 1:csmooth

    % reduce by mean
    rcurve = curve - mean(curve);

    % fill array anew
    curve = zeros(numv, 1);
    for vc = 1:numv

        % get curvature estimate
        cp = rcurve(vc);

        % get neighbors
        pcnx = nei{vc, 2};

        % add curvature of neighbors
        for pc1 = pcnx(:)'
            cp = cp + rcurve(vc);
        end

        % and store back
        curve(vc) = cp;
    end

    % limit data
    cstd  = std(curve);
    cmean = mean(curve);
    curve(curve < (cmean - cstd)) = cmean - cstd;
    curve(curve > (cmean + cstd)) = cmean + cstd;
    curve = fix(0.5 + (curve - min(curve)) / (max(curve) - min(curve) + eps));
end

% split into 0 and 1
mc = mean(curve);
curve = double(curve(:) >= mc);

% create color indexing
curve(:, 2:4) = 0;

% set into file
bc.VertexColor = curve;
xo.C = bc;
