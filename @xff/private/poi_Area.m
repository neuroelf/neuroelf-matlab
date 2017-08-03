function [pas, pa] = poi_Area(xo, srf, poin)
% POI::Area  - calcuate area of POI
%
% FORMAT:       [pas, pa] = poi.Area(srf [, poin])
%
% Input fields:
%
%       srf         surface file object (for coordinates)
%       poin        number of POI (if not given set to 1)
%
% Output fields:
%
%       pas         sum of POI area
%       pa          area of within-POI triangles

% Version:  v1.1
% Build:    16020314
% Date:     Feb-03 2016, 2:36 PM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'poi') || ...
    numel(srf) ~= 1 || ~xffisobject(srf, true, 'srf')
    error('neuroelf:xff:badArguments', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if numel(bc.POI) < 1
    error('neuroelf:xff:invalidObject', ...
        'POI object must contain at least one POI definition.');
end
if nargin < 3 || ~isa(poin, 'double') || numel(poin) ~= 1 || ...
    isnan(poin) || poin < 1 || poin ~= fix(poin) || poin > numel(bc.POI)
    poin = 1;
end
srfc = srf.C;

% get POI vertices
poiv = bc.POI(poin).Vertices(:);

% find triangles that match
tv = srfc.TriangleVertex;
for tvc = 1:numel(poiv)
    tv(tv == poiv(tvc)) = Inf;
end
ti = find(sum(isinf(tv),2) == 3);
if numel(ti) == 0
    error('neuroelf:xff:internalError', 'No triangles in POI.');
end

% get vertices of triangles and replace with numbers from 1 to numel(poiv)
tv = srfc.TriangleVertex(ti, :);
for tvc = 1:numel(poiv)
    tv(tv == poiv(tvc)) = tvc;
end

% get coordinates
tc = srfc.VertexCoordinate(poiv, :);

% calculate sides for Heron / Archimedes formula
% (see http://www.mste.uiuc.edu/dildine/heron/triarea.html for details)
ts = [ ...
    sqrt((tc(tv(:,1),1) - tc(tv(:,2),1)) .^ 2 + ...
         (tc(tv(:,1),2) - tc(tv(:,2),2)) .^ 2 + ...
         (tc(tv(:,1),3) - tc(tv(:,2),3)) .^ 2), ...
    sqrt((tc(tv(:,2),1) - tc(tv(:,3),1)) .^ 2 + ...
         (tc(tv(:,2),2) - tc(tv(:,3),2)) .^ 2 + ...
         (tc(tv(:,2),3) - tc(tv(:,3),3)) .^ 2), ...
    sqrt((tc(tv(:,3),1) - tc(tv(:,1),1)) .^ 2 + ...
         (tc(tv(:,3),2) - tc(tv(:,1),2)) .^ 2 + ...
         (tc(tv(:,3),3) - tc(tv(:,1),3)) .^ 2)];

% use Heron / Archimedes formula
tcs = sum(ts, 2) / 2;
pa = sqrt(tcs .* (tcs-ts(:,1)) .* (tcs-ts(:,2)) .* (tcs-ts(:,3)));

% global area
pas = sum(pa);
