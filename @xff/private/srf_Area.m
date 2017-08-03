function [a, ta, ca] = srf_Area(xo)
% SRF::Area  - compute surface area
%
% FORMAT:       [a, ta, ca] = srf.Area;
%
% No input fields.
%
% Output fields:
%
%       a           total area
%       ta          triangle area
%       ca          area estimate per vertex
%
% Note: this method requires the MEX file meshmorph to be compiled!

% Version:  v1.1
% Build:    16031616
% Date:     Mar-16 2016, 4:36 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'fsbf', 'srf'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end

% get content
bc = xo.C;

% get coordinates and triangle vertices
crd = bc.VertexCoordinate;
tri = bc.TriangleVertex;

% get triangle points
c1 = crd(tri(:,1), :);
c2 = crd(tri(:,2), :);
c3 = crd(tri(:,3), :);

% compute sides
c12 = c1 - c2;
c23 = c2 - c3;
c31 = c3 - c1;

% compute lengths of sides
c1 = sqrt(sum(c12 .* c12, 2));
c2 = sqrt(sum(c23 .* c23, 2));
c3 = sqrt(sum(c31 .* c31, 2));

% Heron formula
s = 0.5 .* (c1 + c2 + c3);
ta = sqrt(s .* (s - c1) .* (s - c2) .* (s - c3));

% total area is the sum
a = sum(ta);

% area estimate per vertex
if nargout > 2

    % distribute
    ca = zeros(size(crd, 1), 1);
    for c = 1:size(tri, 1)
        t = tri(c, :);
        ca(t) = ca(t) + ta(c);
    end

    % weight
    ca = (1 / 3) .* ca;
end
