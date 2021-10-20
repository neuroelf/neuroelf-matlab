function [m, cut, im] = tom_MarkSpot(xo, crd, csz)
% TOM::MarkSpot  - mark and extract a spot around a coordinate
%
% FORMAT:       [mspec, cut, im] = tom.MarkSpot(coord [, cutsize]);
%
% Input fields:
%
%       coord       3-element coordinate around which to cut
%       cutsize     optional cut size (default: 256 pixels)
%
% Output fields:
%
%       mspec       texture map specification (index, x, y)
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
if nargin < 2 || ~isa(crd, 'double') || numel(crd) ~= 3
    error('neuroelf:xff:badArgument', 'Bad or missing COORD argument.');
end
if nargin < 3 || ~isa(csz, 'double') || numel(csz) ~= 1 || isinf(csz) || isnan(csz)
    csz = 256;
else
    csz = max(40, abs(round(csz)));
end

% extract from object
t = xo.C;
v = double(t.Vertices);
tr = t.Triangles;
tm = t.CornerTexVtxMap;

% find closest vertex (dist in mm)
d = sqrt(sum((v - ones(size(v, 1), 1) * crd(:)') .^ 2, 2));
[mv, mpos] = min(d);
if mv > 5
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
m = [msel, mcrd];

% read corresponding image
f = t.Field;
fn = {f.ContentName};
txs = find(strcmpi(fn, 'txtrjpg_') | strcmpi(fn, 'txtrjpga'));
txc = t.Field(txs(msel));
tname = tempname;
fid = fopen(tname, 'wb');
if fid < 1
    error('Cannot write temp file.');
end
fwrite(fid, txc.Content, 'uint8');
fclose(fid);
im = imread(tname);
delete(tname);

% mark and cut spot
isz = size(im);
row = 1 + floor(isz(1) * (1 - mcrd(2)));
col = 1 + floor(isz(2) * mcrd(1));
chalf = ceil(0.5 * csz);
if mcrd(1) < 0.5
    fcol = max(1, col - chalf);
    tcol = fcol + csz - 1;
else
    tcol = min(isz(2), col + chalf);
    fcol = tcol + 1 - csz;
end
if mcrd(2) < 0.5
    trow = min(isz(1), row + chalf);
    frow = trow + 1 - csz;
else
    frow = max(1, row - chalf);
    trow = frow + csz - 1;
end
cut = im(frow:trow, fcol:tcol, :);
im([frow:frow+3, trow-3:trow], fcol:tcol, :) = 255;
im(frow:trow, [fcol:fcol+3, tcol-3:tcol], :) = 255;
