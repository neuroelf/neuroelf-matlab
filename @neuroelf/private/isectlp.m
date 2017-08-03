function [t, u, v] = isectlp(pl, vl, pp, v1p, v2p)
% isectlp  - intersection of line and plain
%
% FORMAT:       [t, u, v] = isectlp(pl, vl, pp, np [, v2p])
%
% Input fields:
%
%       pl          point on line (p0 of line)
%       vl          line vector (l = pl + t * vl)
%       pp          point on plain (p0 of plain)
%       v1p, v2p    in-plain vectors of plain (E = pp + u * v1p + v * v2p)
%
% Output fields:
%
%       t, u, v     distances from p0 of line and plain in units of vectors

% Version:  v0.9d
% Build:    14072317
% Date:     Jul-23 2014, 5:43 PM EST
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

% check arguments
if nargin ~= 5 || ...
   ~isa(pl, 'double') || ...
   ~isa(vl, 'double') || ...
   ~isa(pp, 'double') || ...
   ~isa(v1p, 'double') || ...
   ~isa(v2p, 'double') || ...
   ~any(size(pl) == 3) || ...
   ~isequal(size(pl), size(vl)) || ...
   ~isequal(size(pl), size(pp)) || ...
   ~isequal(size(pl), size(v1p)) || ...
   ~isequal(size(pl), size(v2p))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% vector dim
vs = size(pl);
vd = findfirst(vs == 3);
nv = numel(pl) / 3;
rsa = [3, 1, nv];
vs(vd) = [];
if numel(vs) < 2
    vs(2) = 1;
end

% reshape and permute if necessary
if vd == 1 || ...
    prod(vs(1:vd-1)) == 1
    nop = true;
    pl = reshape(pl, rsa);
    vl = reshape(vl, rsa);
    pp = reshape(pp, rsa);
    v1p = reshape(v1p, rsa);
    v2p = reshape(v2p, rsa);
else
    nop = false;
    rsb = [prod(vs(1:vd-1)), 3, prod(vs(vd:end))];
    pl = reshape(permute(reshape(pl, rsb), [2, 1, 3]), rsa);
    vl = reshape(permute(reshape(vl, rsb), [2, 1, 3]), rsa);
    pp = reshape(permute(reshape(pp, rsb), [2, 1, 3]), rsa);
    v1p = reshape(permute(reshape(v1p, rsb), [2, 1, 3]), rsa);
    v2p = reshape(permute(reshape(v2p, rsb), [2, 1, 3]), rsa);
end

% create matrix and vectors for required solution
t = cat(2, -vl, v1p, v2p);
v = pl - pp;

% find t, u, v
v = transmul(invnd(t), v);

% separate into different arrays
t = v(1, 1, :);
u = v(2, 1, :);
v = v(3, 1, :);

% reshape
if nop
    t = reshape(t, vs);
    u = reshape(u, vs);
    v = reshape(v, vs);
else
    rsb = [prod(vs(vd:end)), prod(vs(1:vd-1))];
    t = reshape(t, rsb);
    u = reshape(u, rsb);
    v = reshape(v, rsb);
end
