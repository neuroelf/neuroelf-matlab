function [v, f] = mesh3d(meshsize, opts, varargin)
% mesh3d  - generate a flat mesh grid for patch in 3D coords
%
% FORMAT:       [v, f] = mesh3d(meshsize, opts)
%
% Input fields:
%
%       meshsize    1x2 mesh size (e.g. [80, 64])
%       opts        optional settings
%        .orient    either of 1, 2, or {3}, determining the flat axis
%        .trf       4x4 transformation matrix (will be applied at the end)
%        .xline     1x2 vector with [mx, bx] for x = mx * x0 + bx
%        .yline     1x2 vector with [my, by] for y = mx * y0 + by
%        .zvalue    fixed values of 3rd axis (default: 0)
%
% Output fields:
%
%       v           Vx3 vertices
%       f           Fx3 face vertex indices (1-based, triangles)

% Version:  v0.9d
% Build:    14071115
% Date:     Jul-11 2014, 3:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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
if nargin > 1 && ...
    isa(opts, 'double') && ...
    numel(opts) == 1 && ...
    isa(meshsize, 'double') && ...
    numel(meshsize) == 1
    meshsize = [meshsize, opts];
    if nargin > 2
        opts = varargin{1};
    end
end
if nargin < 1 || ...
   ~isa(meshsize, 'double') || ...
    numel(meshsize) ~= 2 || ...
    any(isinf(meshsize) | isnan(meshsize) | meshsize < 0 | meshsize ~= round(meshsize))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'orient') || ...
   ~isa(opts.orient, 'double') || ...
    numel(opts.orient) ~= 1 || ...
    isinf(opts.orient) || ...
    isnan(opts.orient) || ...
   ~any(opts.orient == (1:6))
    opts.orient = 3;
end
if ~isfield(opts, 'trf') || ...
   ~isa(opts.trf, 'double') || ...
   ~isequal(size(opts.trf), [4, 4]) || ...
    any(isinf(opts.trf(:)) | isnan(opts.trf(:))) || ...
    any(opts.trf(4, :) ~= [0, 0, 0, 1])
    opts.trf = [];
end
if ~isfield(opts, 'xline') || ...
   ~isa(opts.xline, 'double') || ...
    numel(opts.xline) ~= 2 || ...
    any(isinf(opts.xline) | isnan(opts.xline)) || ...
    opts.xline(1) == 0
    opts.xline = [1, 0];
end
if ~isfield(opts, 'yline') || ...
   ~isa(opts.yline, 'double') || ...
    numel(opts.yline) ~= 2 || ...
    any(isinf(opts.yline) | isnan(opts.yline)) || ...
    opts.yline(1) == 0
    opts.yline = [1, 0];
end
if ~isfield(opts, 'zvalue') || ...
   ~isa(opts.zvalue, 'double') || ...
    numel(opts.zvalue) ~= 1 || ...
    isinf(opts.zvalue) || ...
    isnan(opts.zvalue)
    opts.zvalue = 0;
end

% generate grid
[vx, vy] = ndgrid(0:meshsize(1), 0:meshsize(2));
if opts.zvalue == 0
    vz = zeros(numel(vx), 1);
else
    vz = opts.zvalue .* ones(numel(vx), 1);
end

% linearize for computations
vx = vx(:);
vy = vy(:);

% apply lines for x and y direction
if opts.xline(1) ~= 1
    vx = opts.xline(1) .* vx + opts.xline(2);
elseif opts.xline(2) ~= 0
    vx = vx + opts.xline(2);
end
if opts.yline(1) ~= 1
    vy = opts.yline(1) .* vy + opts.yline(2);
elseif opts.yline(2) ~= 0
    vy = vy + opts.yline(2);
end

% put together to 3D coordinates according to orientation
switch (opts.orient)
    case {1}
        v = [vz, vx, vy];
    case {2}
        v = [vx, vz, vy];
    case {3}
        v = [vx, vy, vz];
    case {4}
        v = [vz, vy, vx];
    case {5}
        v = [vy, vz, vx];
    case {6}
        v = [vy, vx, vz];
end

% generate first pair of faces (as column vector)
f = [1; 2; 3 + meshsize(1); 3 + meshsize(1); 2 + meshsize(1); 1];

% generate the necessary faces for first column of data
f = f * ones(1, meshsize(1)) + ones(6, 1) * (0:(meshsize(1)-1));

% linearize for next dimension
f = f(:);

% then replicate and add per-column offsets
f = f * ones(1, meshsize(2)) + ...
    ones(numel(f), 1) * (0:(meshsize(1)+1):((meshsize(1)+1)*(meshsize(2)-1)));

% reshape
f = reshape(f, 3, 2 * prod(meshsize));

% apply transformation
if ~isempty(opts.trf)

    % add necessary ones
    f(4, :) = 1;

    % multiply
    f = opts.trf * f;

    % remove ones
    f(4, :) = [];
end

% transpose
f = f';
