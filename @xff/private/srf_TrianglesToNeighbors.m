function [n, t] = srf_TrianglesToNeighbors(xo, triref)
% SRF::TrianglesToNeighbors  - get full 1. order neighbors lists
%
% FORMAT:       [nei, tri] = srf.TrianglesToNeighbors([triref]);
%
% Input fields:
%
%       triref      optional boolean flag, create triangle ref list as well
%
% Output fields:
%
%       nei         Vx2 cell array, neighbors lists (for Neighbors property)
%       tri         Vx1 cell array, triangle back-reference lists
%
% Using: mesh_trianglestoneighbors.

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:39 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get contents
bc = xo.C;
nc = size(bc.VertexCoordinate, 1);

% which call
if nargin > 1 && islogical(triref) && numel(triref) == 1 && triref
    [n, bd, t] = ne_methods.mesh_trianglestoneighbors(nc, bc.TriangleVertex);
else
    n = ne_methods.mesh_trianglestoneighbors(nc, bc.TriangleVertex);
end
