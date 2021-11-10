function [m, cut, im] = tom_MarkSpot(xo, crd, csz, msz)
% TOM::MarkSpot  - mark and extract a spot around a coordinate
%
% FORMAT:       [mspec, cut, im] = tom.MarkSpot(coord [, cutsize]);
%
% Input fields:
%
%       coord       3-element coordinate around which to cut (x, y, z)
%       cutsize     optional cut size (default: 256 pixels)
%       marksize    optional mark size (default: 0 = no marking)
%
% Output fields:
%
%       mspec       texture map specification (tnum, cx, cy)
%       cut         cut-out piece
%       im          full texture image

% Version:  v1.1
% Build:    21102012
% Date:     Oct-20 2021, 12:50 PM EST
% Author:   Jochen Weber, NeuroElf.net, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2021, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'tom')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if ~isa(crd, 'double') || numel(crd) ~= 3 || any(isinf(crd) | isnan(crd))
    error('neuroelf:xff:badArgument', 'Bad or missing COORD argument.');
end
if nargin < 3 || ~isa(csz, 'double') || numel(csz) ~= 1 || isinf(csz) || isnan(csz)
    csz = 256;
else
    csz = max(40, abs(round(csz)));
end
if nargin < 4 || ~isa(msz, 'double') || numel(msz) ~= 1 || isinf(msz) || isnan(msz)
    msz = 0;
else
    msz = min(round(0.8 * csz), abs(round(msz)));
end

% extract from object
t = xo.C;
v = t.VertexCoordinate;
tr = t.TriangleVertex;
tm = t.CornerTexVtxMap;

% find closest vertex (dist in mm)
d = sqrt(sum((v - ones(size(v, 1), 1) * crd(:)') .^ 2, 2));
[mv, mpos] = min(d);
if mv > 7.5
    error('neuroelf:xff:noMatchFound', 'No such coordinate in dataset.')
end

% find relevant triangles, precompute vectors
vpos = v(mpos, :);
vcrd = crd - vpos;
tc = find(any(tr == mpos, 2));
nt = numel(tc);
tcs = tr(tc, :);
tc2 = (tcs(:, 2) == mpos);
tc3 = (tcs(:, 3) == mpos);
tcs(tc2, :) = tcs(tc2, [2, 3, 1]);
tcs(tc3, :) = tcs(tc3, [3, 1, 2]);

% which triangle does the point fall into
if mv > 0
    v2 = v(tcs(:, 2), :) - v(tcs(:, 1), :);
    l2 = v2 ./ (sqrt(sum(v2 .* v2, 2)) * [1, 1, 1]);
    v3 = v(tcs(:, 3), :) - v(tcs(:, 1), :);
    l3 = v3 ./ (sqrt(sum(v3 .* v3, 2)) * [1, 1, 1]);
    vn = cross(l2, l3, 2);
    vd = sqrt(sum(vn .* vn, 2));
    vn = vn ./ (vd * [1, 1, 1]);
    vcd = sum(vn .* vcrd(ones(nt, 1), :), 2);
    pcrd = crd(ones(nt, 1), :) - (vcd * [1, 1, 1]) .* vn;
    pdst = pcrd - vpos(ones(nt, 1), :);
    insidx = -ones(nt, 1);
    ins1 = zeros(nt, 1);
    ins2 = zeros(nt, 1);
    for tcc = 1:nt
        tv = [v2(tcc, :)', v3(tcc, :)'] \ pdst(tcc, :)';
        if all(tv >= -1e-5)
            insidx(tcc) = min(tv);
            ins1(tcc) = tv(1);
            ins2(tcc) = tv(2);
        end
    end
    [insmax, insidx] = max(insidx);
    if insmax == -1
        error('neuroelf:xff:triangleMismatch', 'Could not locate point in triangles.');
    end
    ins1 = ins1(insidx);
    ins2 = ins2(insidx);
else
    insidx = 1;
    ins1 = 0;
    ins2 = 0;
end
tr = tc(insidx);
ttr = t.CornerTexVtxMap(tm(:, 1) == tr, [2, 3]);
[~, tti] = sort(ttr(:, 1));
ttr = ttr(tti, :);
tvc = double(t.TexVertACoord(ttr(:, 2), :));
tvm = t.TexVertAMap(ttr(:, 2));
msel = tvm(1);
if ~all(tvm == msel)
    error('neuroelf:xff:triangleMismatch', 'Texture map mismatch between vertices.');
end
if tc2(insidx)
    tvc = tvc([2, 3, 1], :);
elseif tc3(insidx)
    tvc = tvc([3, 1, 2], :);
end
tv1 = tvc(2, :) - tvc(1, :);
tv2 = tvc(3, :) - tvc(1, :);
mcrd = tvc(1, :) + ins1 * tv1 + ins2 * tv2;

% combine
m = [msel, mcrd];

% pass on to secondary function
if nargout > 1
    [cut, im] = tom_ExtractSpot(xo, m, csz, msz);
end
