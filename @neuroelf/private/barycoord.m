function [u, v] = barycoord(t1, t2, t3, v)
% barycoord  - returns barycentric coordinates of point vs. triangle
%
% FORMAT:       [u, v] = barycoord(t1, t2, t3, v)
%
% Input fields:
%
%       t1, t2, t3  coordinates of triangle points
%       v           point to find in triangle with (t2-t1) and (t3-t1)
%
% Output fields:
%
%       u, v        units to multiply vectors with to get from t1 to v

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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
if nargin < 4 || ...
   ~isa(t1, 'double') || ...
    isempty(t1) || ...
   ~isa(t2, 'double') || ...
    isempty(t2) || ...
   ~isa(t3, 'double') || ...
    isempty(t3) || ...
   ~isa(v, 'double') || ...
    isempty(v)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument' ...
    );
end

% adapt t1 and v if needed
if size(t1, 1) == 1 && ...
    size(t2, 1) > 1
    t1 = t1(ones(1, size(t2, 1)), :);
end
if size(v, 1) == 1 && ...
    size(t1, 1) > 1
    v = v(ones(1, size(t1, 1)), :);
end

% compute helping vectors
va = t2(:, 1) - t1(:, 1);
vb = t3(:, 1) - t1(:, 1);
vc = t1(:, 1) -  v(:, 1);
vd = t2(:, 2) - t1(:, 2);
ve = t3(:, 2) - t1(:, 2);
vf = t1(:, 2) -  v(:, 2);
vg = t2(:, 3) - t1(:, 3);
vh = t3(:, 3) - t1(:, 3);
vi = t1(:, 3) -  v(:, 3);

% compute u and v
vdg = vd + vg;
veh = ve + vh;
vfi = vf + vi;
u = (vb .* (vfi) - vc .* veh) ./ (va .* veh - vb .* vdg);
v = (va .* (vfi) - vc .* vdg) ./ (vb .* vdg - va .* veh);

% for invalid entries
r = isinf(u) | isnan(u) | isinf(v) | isnan(v);
if any(r)
    r = find(r);
    vag = vd(r) + vg(r);
    vbh = ve(r) + vh(r);
    vci = vf(r) + vi(r);
    u(r) = (ve(r) .* vci - vf(r) .* vbh) ./ (vd(r) .* vbh - ve(r) .* vag);
    v(r) = (vd(r) .* vci - vf(r) .* vag) ./ (ve(r) .* vag - vd(r) .* vbh);
end
