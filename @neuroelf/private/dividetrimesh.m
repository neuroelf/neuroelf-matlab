function [t, v] = dividetrimesh(t, v)
% dividetrimesh -  divide triangles of a mesh (sub-tesselation)
%
% FORMAT:       [t, v] = dividetrimesh(t, v)
%
% Input/output fields:
%
%       t           triangles
%       v           vertices

% Version:  v0.9d
% Build:    14062814
% Date:     Jun-28 2014, 2:10 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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
if nargin < 1 || ...
   ~isa(t, 'double') || ...
    ndims(t) > 2 || ...
    size(t, 2) ~= 3 || ...
    any(isinf(t(:)) | isnan(t(:)) | t(:) < 1 | t(:) ~= fix(t(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
nv = max(t(:));
rc = true;
if nargin < 2
    v = zeros(nv, 3);
    rc = false;
elseif ~isa(v, 'double') || ...
   ~isequal(size(v), [nv, 3]) || ...
    any(isinf(v(:)) | isnan(v(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad v argument.' ...
    );
end

% copy triangles
ot = t;
nt = size(t, 1);

% prepare larger array
t(4 * nt, 3) = 0;

% create lookup arrays -> for vertex recomputation
if rc
    vl = zeros(4 * nt, 2);
end

% -> and for splits (up to 24 direct neighbors supported)
sp = zeros(nv, 24, 2);
spi = ones(nv, 1);

% initialize counter
nnv = nv + 1;

% for each triangle
for c = 1:nt

    % get the three vertices
    tri = ot(c, :);

    % first pair (P1, P2) not already split?
    if ~any(sp(tri(1), :, 1) == tri(2))

        % add split coordinate
        sp1 = nnv;
        nnv = nnv + 1;

        % add to splits
        sp(tri(1), spi(tri(1)), :) = [tri(2), sp1];
        spi(tri(1)) = spi(tri(1)) + 1;
        sp(tri(2), spi(tri(2)), :) = [tri(1), sp1];
        spi(tri(2)) = spi(tri(2)) + 1;

    % or already split
    else
        sp1 = sp(tri(1), sp(tri(1), :, 1) == tri(2), 2);
    end

    % repeat for second pair
    if ~any(sp(tri(2), :, 1) == tri(3))
        sp2 = nnv;
        nnv = nnv + 1;
        sp(tri(2), spi(tri(2)), :) = [tri(3), sp2];
        spi(tri(2)) = spi(tri(2)) + 1;
        sp(tri(3), spi(tri(3)), :) = [tri(2), sp2];
        spi(tri(3)) = spi(tri(3)) + 1;
    else
        sp2 = sp(tri(2), sp(tri(2), :, 1) == tri(3), 2);
    end

    % and for third pair also
    if ~any(sp(tri(3), :, 1) == tri(1))
        sp3 = nnv;
        nnv = nnv + 1;
        sp(tri(3), spi(tri(3)), :) = [tri(1), sp3];
        spi(tri(3)) = spi(tri(3)) + 1;
        sp(tri(1), spi(tri(1)), :) = [tri(3), sp3];
        spi(tri(1)) = spi(tri(1)) + 1;
    else
        sp3 = sp(tri(3), sp(tri(3), :, 1) == tri(1), 2);
    end

    % replace original triangle with new one
    t(c, :) = [tri(1), sp1, sp3];

    % and add three new triangles
    t(nt+3*c-2:nt+3*c, :) = [tri(2), sp2, sp1; tri(3), sp3, sp2; sp1, sp2, sp3];

    % add recomputation vectors
    if rc
        vl([sp1, sp2, sp3], :) = [tri(1:2); tri(2:3); tri([3, 1])];
    end
end

% recompute vertices
if rc
    nnv = nnv - 1;
    v(nv+1:nnv, :) = 0.5 * (v(vl(nv+1:nnv, 1), :) + v(vl(nv+1:nnv, 2), :));
end
