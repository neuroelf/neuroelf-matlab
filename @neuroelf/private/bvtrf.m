function m44 = bvtrf(t, r, s, sh, o)
% bvtrf  - compile a transformation matrix as in the 3D volume tools
%
% FORMAT:       trf = bvtrf(trans, rot [, scale [,shear [, orig]]])
%
% Input fields:
%
%       trans       1x3 translation values (X,Y,Z)
%       rot         1x3 rotation values (around X,Y,Z, in GRAD)
%       scale       1x3 scaling values (256: no scaling, in X,Y,Z dir)
%       shear       1x3 shearing vectors (XY, XZ, YZ, not in BV)
%       orig        1x3 origin vector (default: 0, 0, 0)
%
% Output fields:
%
%       trf         4x4 transformation matrix to apply on system coords
%
% Note: the TFMatrix store in TRF files is the *inverse* of the result

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
if nargin < 1 || ...
   ~isa(t, 'double') || ...
    numel(t) ~= 3 || ...
    any(isinf(t) | isnan(t) | t < -512 | t > 512)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid translation values.' ...
    );
end
t = t(:)';
if nargin < 2 || ...
   ~isa(r, 'double') || ...
    numel(r) ~= 3 || ...
    any(isinf(r) | isnan(r) | r < -360 | r > 360)
    warning( ...
        'neuroelf:BadArgument', ...
        'Invalid rotation values; assuming no rotation.' ...
    );
    r = [0, 0, 0];
end
r = r(:)';
if nargin < 3 || ...
    isempty(s)
    s = [256, 256, 256];
elseif ~isa(s, 'double') || ...
    numel(s) ~= 3 || ...
    any(isinf(s) | isnan(s) | abs(s) < 64 | abs(s) > 1024)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid scaling argument.' ...
    );
end
s = s(:)';
if nargin < 4 || ...
    isempty(sh)
    sh = [0, 0, 0];
elseif ~isa(sh, 'double') || ...
    numel(sh) ~= 3 || ...
    any(isinf(sh) | isnan(sh) | sh < -256 | sh > 256)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid shearing argument.' ...
    );
end
sh = sh(:)';
if nargin < 5 || ...
    isempty(o)
    o = [0, 0, 0];
elseif ~isa(o, 'double') || ...
    numel(o) ~= 3 || ...
    any(isinf(o) | isnan(o) | o < -256 | o > 256)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid shearing argument.' ...
    );
end
o = o(:);

% origin appliance
o44 = eye(4);
b44 = eye(4);
if any(o ~= 0)
    o44(1:3, 4) = -o;
    b44(1:3, 4) = o;
end

% translation
t44 = eye(4);
t44(1:3, 4) = -t([2, 3, 1]);

% rotation
r44 = tfmatrix(struct( ...
    'type', {'r', 'r', 'r'}, ...
    'xyz' , {[-r(2) * pi / 180, 0, 0], ...
             [0,  r(3) * pi / 180, 0], ...
             [0, 0, -r(1) * pi / 180]}))';

% scaling
s44 = diag([s([2, 3, 1]) ./ 256, 1]);

% shearing
h44 = eye(4);
h44(1, 2) = sh(2);
h44(1, 3) = sh(3);
h44(2, 3) = sh(1);

% complete matrix
m44 = t44 * o44 * r44 * s44 * b44 * h44;
