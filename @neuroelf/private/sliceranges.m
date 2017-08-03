function xyzfromto = sliceranges(v, nelem)
% sliceranges  - compute x/y/z from and to ranges for slicing directions
%
% FORMAT:       xyzfromto = sliceranges(v [, nelem])
%
% Input fields:
%
%       v           3D numeric volume
%       nelem       optional neutral element (default: 0)
%
% Output fields:
%
%       xyzfromto   struct with fields
%        .xfromy    X-by-1 list with Y-indices in X-slice direction from
%                   data is different from nelem
%        .xtoy      X-by-1 list with Y-indices in X-slice direction to ...
%        .xfromz
%        .xtoz      same (and also for other directions)

% Version:  v0.9c
% Build:    14022820
% Date:     Feb-21 2014, 5:41 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://www.mathworks.com/matlabcentral/fileexchange/21993-viewer3d

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
if nargin < 1 || ...
   ~isnumeric(v) || ...
    isempty(v) || ...
    ndims(v) < 3
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
v = v(:, :, :, 1);
if nargin < 2 || ...
   ~isnumeric(nelem) || ...
    numel(nelem) ~= 1
    nelem = 0;
end
if istransio(nelem)
    nelem = nelem(:);
end
if isa(nelem, 'double') || ...
    isa(nelem, 'single')
    if isinf(nelem)
        sn = sign(nelem);
        nelem = randn(1, 1);
        if isa(v, 'double') || ...
            isa(v, 'single')
            ii = isinf(v);
            if any(ii(:))
                v(ii) = (sn * nelem) .* sign(v(ii));
            end
        end
    elseif isnan(nelem)
        nelem = randn(1, 1);
        if isa(v, 'double') || ...
            isa(v, 'single')
            v(isnan(v)) = nelem;
        end
    end
end

% output
s = size(v);
if numel(s) < 3
    s(3) = 1;
end
x = NaN .* ones(s(1), 1);
y = NaN .* ones(s(2), 1);
z = NaN .* ones(s(3), 1);
xyzfromto = struct( ...
    'xfromy', x, ...
    'xtoy',   x, ...
    'xfromz', x, ...
    'xtoz',   x, ...
    'yfromx', y, ...
    'ytox',   y, ...
    'yfromz', y, ...
    'ytoz',   y, ...
    'zfromx', z, ...
    'ztox',   z, ...
    'zfromy', z, ...
    'ztoy',   z);

% iterate over slices
for c = 1:s(1)
    sl = reshape(v(c, :, :), s(2), s(3));
    d1 = any(sl ~= nelem, 2);
    d2 = any(sl ~= nelem, 1);
    if any(d1)
        xyzfromto.xfromy(c) = findfirst(d1);
        xyzfromto.xtoy(c) = findfirst(d1, -1);
        xyzfromto.xfromz(c) = findfirst(d2);
        xyzfromto.xtoz(c) = findfirst(d2, -1);
    end
end
for c = 1:s(2)
    sl = reshape(v(:, c, :), s(1), s(3));
    d1 = any(sl ~= nelem, 2);
    d2 = any(sl ~= nelem, 1);
    if any(d1)
        xyzfromto.yfromx(c) = findfirst(d1);
        xyzfromto.ytox(c) = findfirst(d1, -1);
        xyzfromto.yfromz(c) = findfirst(d2);
        xyzfromto.ytoz(c) = findfirst(d2, -1);
    end
end
for c = 1:s(3)
    sl = v(:, :, c);
    d1 = any(sl ~= nelem, 2);
    d2 = any(sl ~= nelem, 1);
    if any(d1)
        xyzfromto.zfromx(c) = findfirst(d1);
        xyzfromto.ztox(c) = findfirst(d1, -1);
        xyzfromto.zfromy(c) = findfirst(d2);
        xyzfromto.ztoy(c) = findfirst(d2, -1);
    end
end
