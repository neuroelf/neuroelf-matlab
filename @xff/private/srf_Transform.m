function xo = srf_Transform(xo, mat44)
% SRF::Transform  - apply 4x4 quaternion to coordinates/normals
%
% FORMAT:       [srf] = srf.Transform(mat44)
%
% Input fields:
%
%       mat44       4x4 affine transformation matrix
%
% Output fields:
%
%       srf         altered object

% Version:  v1.1
% Build:    16021121
% Date:     Feb-11 2016, 9:25 PM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf') || ...
   ~isa(mat44, 'double') || numel(size(mat44)) ~= 2 || any(size(mat44) ~= 4) || ...
    any(isinf(mat44(:)) | isnan(mat44(:)))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;

% apply transformations
v = bc.VertexCoordinate;
v = v - repmat(bc.MeshCenter, [size(v, 1), 1]);
v(:, 4) = 1;
n = bc.VertexNormal;
v = (mat44 * v')';
n = (mat44(1:3, 1:3) * n')';
n = n ./ repmat(sqrt(sum(n .* n, 2)), [1, 3]);
n(any(isinf(n), 2) | any(isnan(n), 2), :) = sqrt(3);

% put back into object
bc.VertexCoordinate = v(:, 1:3) + repmat(bc.MeshCenter, [size(v, 1), 1]);
bc.VertexNormal = n;
xo.C = bc;
