function xo = srf_UpdateNeighbors(xo)
% SRF::UpdateNeighbors  - update Neighbors from triangles and colors
%
% FORMAT:       [srf = ] srf.UpdateNeighbors
%
% No input/output fields.
%
% Using: mesh_trianglestoneighbors.

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:36 PM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end

% get content
bc = xo.C;

% get deletion status of vertices
del = (bc.VertexColor(:, 1) > 2^31);

% create new triangles list
try
    n = ne_methods.mesh_trianglestoneighbors(size(bc.VertexCoordinate, 1), bc.TriangleVertex);
catch xfferror
    rethrow(xfferror);
end

% any deleted
if any(del)

    % remove those from all neighbor lists
    for vc = 1:size(n, 1)
        if any(del(n{vc, 2}))
            n{vc, 2}(del(n{vc, 2})) = [];
            n{vc, 1} = numel(n{vc, 2});
            if n{vc, 1} == 0
                n{vc, 1} = 1;
                n{vc, 2} = vc;
            end
        end
    end
end

% set back into object
bc.Neighbors = n;
xo.C = bc;
