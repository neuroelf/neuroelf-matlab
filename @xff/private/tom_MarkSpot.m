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
%
% Using: catstruct.

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

% neuroelf library
global ne_methods;

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
v = double(t.Vertices);
tr = t.Triangles;
tm = t.CornerTexVtxMap;

% find closest vertex (dist in mm)
d = sqrt(sum((v - ones(size(v, 1), 1) * crd(:)') .^ 2, 2));
[mv, mpos] = min(d);
if mv > 7.5
    error('neuroelf:xff:noMatchFound', 'No such coordinate in dataset.')
end

% find relevant triangles
tc = find(any(tr == mpos, 2));

% and seek corresponding texture vertices
ti = false(numel(t.TexVertAMap), 1);
for c = 1:numel(tc)
    trow = (tm(:, 1) == tc(c));
    ti(unique(tm(trow, 3))) = true;
end

% extract corresponding maps and coordinates
maps = t.TexVertAMap(ti);
mcrd = double(t.TexVertACoord(ti, :));

% remove irrelevant maps
if any(maps ~= maps(1))
    msel = mode(maps);
    mcrd = mcrd(maps == msel, :);
else
    msel = maps(1);
end

% find "central" coordinate
mcrd = mean(mcrd);
mcrd(2) = 1 - mcrd(2);
m = [msel, mcrd];

% pass on to secondary function
if nargout > 1
    [cut, im] = tom_ExtractSpot(xo, m, csz, msz);
end
