function xo = srf_AddBitmap(xo, bmpfile, cc, xdir, ydir)
% SRF::AddBitmap  - add a 2d bitmap to a SRF object
%
% FORMAT:       [srf] = srf.AddBitmap(bmpfile, cc, xdir, ydir)
%
% Input fields:
%
%       bmpfile     either a XxYx3 image matrix or filename for imread
%       cc          center coordinate (BV internal coordinate)
%       xdir        direction alogn columns
%       ydir        direction along rows
%
% Output fields:
%
%       srf         altered SRF object

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:33 PM EST
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

% check arguments
if nargin < 5 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf') || ...
   ((~ischar(bmpfile) || isempty(bmpfile)) && (~isnumeric(bmpfile) || size(bmpfile, 3) ~= 3)) || ...
   ~isa(cc, 'double') || numel(cc) ~= 3 || any(isinf(cc) | isnan(cc) | cc < -256 | cc > 256) || ...
   ~isa(xdir, 'double') || numel(xdir) ~= 3 || any(isinf(xdir) | isnan(xdir) | xdir < -32 | xdir > 32) || ...
    all(xdir == 0) || ~isa(ydir, 'double') || numel(ydir) ~= 3 || ...
    any(isinf(ydir) | isnan(ydir) | ydir < -32 | ydir > 32) || ...
    all(ydir(:) == 0 | ydir(:) == xdir(:))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
cc = cc(:)';
xdir = xdir(:)';
ydir = ydir(:)';
xnrm = xdir ./ sqrt(sum(xdir .* xdir));
ynrm = ydir ./ sqrt(sum(ydir .* ydir));
znrm = cross(xnrm, ynrm);

% check bmpfile
if ischar(bmpfile)
    try
        bmpfile = imread(bmpfile(:)');
        if size(bmpfile, 3) ~= 3
            error('BAD_IMAGE');
        end
    catch xfferror
        neuroelf_lasterr(xfferror);
        error('neuroelf:xff:badArgument', 'Invalid image filename.');
    end
end

% current SRF settings
nv = size(bc.VertexCoordinate, 1);

% generate coordinates, etc.
xs = size(bmpfile, 2);
ys = size(bmpfile, 1);
fc = cc + ((1 - xs) / 2) * xdir + ((1 - ys) / 2) * ydir;
[xg, yg] = meshgrid(0:(xs - 1), 0:(ys - 1));
ic = [fc(1) + xdir(1) * xg(:) + ydir(1) * yg(:), ...
    fc(2) + xdir(2) * xg(:) + ydir(2) * yg(:), ...
    fc(3) + xdir(3) * xg(:) + ydir(3) * yg(:)];
nrm = repmat(znrm, [size(ic, 1), 1]);
ccodes = [nan(size(ic, 1), 1), double(reshape(bmpfile, [size(ic, 1), 3]))];
nei = cell(nv, 2);
ym = ys - 1;
sp = [-1, -ys, -ym, 1, ys, ym] + nv;
for xc = 1:xs
    for yc = 1:ys
        ni = yc + (xc - 1) * ys;
        if xc == 1
            if yc == 1
                nei(ni, :) = {2, [ni + 1, ni + ys] + nv};
            elseif yc == ys
                nei(ni, :) = {3, [ni + ys, ni + ym, ni - 1] + nv};
            else
                nei(ni, :) = {4, [ni + 1, ni + ys, ni + ym, ni - 1] + nv};
            end
        elseif xc == xs
            if yc == 1
                nei(ni, :) = {3, [ni + 1, ni - ys, ni - ym] + nv};
            elseif yc == ys
                nei(ni, :) = {2, [ni - 1, ni - ys] + nv};
            else
                nei(ni, :) = {4, [ni - 1, ni - ys, ni - ym, ni + 1] + nv};
            end
        else
            if yc == 1
                nei(ni, :) = {4, [ni - ys, ni - ym, ni + 1, ni + ys] + nv};
            elseif yc == ys
                nei(ni, :) = {4, [ni + ys, ni + ym, ni - 1, ni - ys] + nv};
            else
                nei(ni, :) = {6, sp + ni};
            end
        end
    end
end
ctri = [reshape([1:(ys - 1);(1 + ys):(2 * ys - 1)], [2 * (ys - 1), 1]), ...
    reshape(floor(2:0.5:(ys + 0.5)), [2 * (ys - 1), 1]), ...
    reshape(floor((ys + 1.5):0.5:(2 * ys)), [2 * (ys - 1), 1])];
tri = zeros(2 * (xs - 1) * (ys - 1), 3);
colc = nv;
for k = 1:(2 * (ys - 1)):size(tri, 1)
    tri(k:(k - 1 + size(ctri, 1)), :) = ctri + colc;
    colc = colc + ys;
end

% add BMP
bc.VertexCoordinate = [bc.VertexCoordinate; ic];
bc.VertexNormal = [bc.VertexNormal; nrm];
bc.VertexColor = [bc.VertexColor; ccodes];
bc.Neighbors = [bc.Neighbors; nei];
bc.TriangleVertex = [bc.TriangleVertex; tri];
bc.NrOfVertices = size(bc.VertexCoordinate, 1);
bc.NrOfTriangles = size(bc.TriangleVertex, 1);
xo.C = bc;
