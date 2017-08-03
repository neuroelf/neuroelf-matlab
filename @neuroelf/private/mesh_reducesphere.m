function [t, v] = mesh_reducesphere(t, v)
%MESH_REDUCESPHERE  reduce triangle faces of a regular icosahedron surface
%   T = MESH_REDUCESPHERE(T, V) reduces the triangles from an icosahedron
%   surface such that one regular triangle division is removed again.
%
%   See also: MESH_TRIDIVIDE.

% Version:  v1.1
% Build:    16031812
% Date:     Mar-18 2016, 12:13 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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
if nargin ~= 2 || ~isa(t, 'double') || ~isa(v, 'double') || ndims(t) > 2 || ndims(v) > 2 || ...
    isempty(t) || isempty(v) || size(t, 2) ~= 3 || size(v, 2) ~= 3 || ...
   ~any(size(t, 1) == (20 .* (4 .^ (1:9)))) || ~any(size(v, 1) ~= (2 + 10 .* (4 .^ (1:9)))) || ...
    round(size(t, 1) / size(v, 1)) ~= 2 || any(isinf(v(:)) | isnan(v(:))) || ...
    any(isinf(t(:)) | isnan(t(:)) | t(:) < 1 | t(:) > size(v, 1) | t(:) ~= round(t(:))) || ...
    any(histcount(t(:), 1, size(v, 1), 1) > 6)
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end

% generate vertex to triangle mapping
nv = size(v, 1);
vtm = sparse(repmat((1:size(t, 1))', 3, 1), t(:), ones(numel(t), 1), size(t, 1), nv, numel(t));

% start with the first to-keep vertex
nextv = find(histcount(t(:), 1, nv, 1) == 5);
nextr = [];

% repeat until done
while ~isempty(nextv)
    
    % iterate over vertices
    for vc = 1:numel(nextv)

        % find all triangles with this vertex
        v1 = nextv(vc);
        vtri = find(vtm(:, v1));

        % iterate over those
        for rc = 1:numel(vtri)

            % still needs collapsing
            if t(vtri(rc), 1) > 0

                % collapse: get triangles for other two vertices
                tri = vtri(rc);
                tvs = t(tri, :);
                if tvs(1) == v1
                    ps = 1;
                    v2 = tvs(2);
                    v3 = tvs(3);
                elseif tvs(2) == v1
                    ps = 2;
                    v2 = tvs(3);
                    v3 = tvs(1);
                else
                    ps = 3;
                    v2 = tvs(1);
                    v3 = tvs(2);
                end
                v2tri = find(vtm(:, v2));
                v3tri = find(vtm(:, v3));

                % intersection (without tri) is the triangle opposite
                otri = intersect(v2tri, v3tri);
                otri(otri == tri) = [];

                % get opposite vertex (to original vertex)
                v4 = t(otri, :);
                v4(v4 == v2 | v4 == v3) = [];

                % repeat with otri
                v4tri = find(vtm(:, v4));
                v5tri = intersect(v2tri, v4tri);
                v5tri(v5tri == otri) = [];
                v6tri = intersect(v3tri, v4tri);
                v6tri(v6tri == otri) = [];

                % get correct vertices
                cv2 = t(v5tri, :);
                cv2(cv2 == v2 | cv2 == v4) = [];
                cv3 = t(v6tri, :);
                cv3(cv3 == v3 | cv3 == v4) = [];

                % replace
                if ps == 1
                    t(tri, :) = [ -v1, -cv2, -cv3];
                elseif ps == 2
                    t(tri, :) = [-cv3,  -v1, -cv2];
                else
                    t(tri, :) = [-cv2, -cv3,  -v1];
                end

                % invalidate remaining triangles
                t([otri, v5tri, v6tri], :) = 0;

                % what next
                n = [cv2, cv3];
                nextr(end+1:end+numel(n)) = n;
            end
        end
    end
    nextv = nextr;
    nextr = [];
end

% remove empty triangles from list
t(all(t == 0, 2), :) = [];

% find all remaining vertices
rv = -t(t < 0);
mrv = max(rv);
urv = unique(rv);

% requires compressing
if numel(rv) ~= numel(t) || numel(urv) ~= mrv
    rv = abs(t(t ~= 0));
    mrv = max(rv);
    urv = unique(rv);
    mmat = zeros(1, mrv);
    mmat(urv) = 1:numel(urv);
    v = v(urv, :);
    t = reshape(t(mmat(abs(t))), size(t));
else
    t = -t;
    v = v(1:mrv, :);
end



% sub-function (for recursive calling
function [t, n] = collapse(t, v1, tri, vtm)

