function vol = dilate3d(vol, dop, r, res)
% dilate3d  - perform 3D dilation with operator
%
% FORMAT:       vol = dilate3d(vol [, dop, r, res])
%
% Input fields:
%
%       vol         3-D boolean volume
%       dop         dilation operator with odd size ! (default 3x3x3 cube)
%       r           1x1 double radius, weighted dilation
%       res         1x3 radius resolution (default distance 1)
%
% Output fields:
%
%       vol         3-D boolean volume of same size as input

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
    ndims(vol) ~= 3 || ...
   ~any(strcmpi(class(vol), {'logical', 'uint8', 'uint16'}))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if isempty(vol)
    vol = false(size(vol));
    return;
end
if nargin < 2 || ...
    ndims(dop) > 3 || ...
    any(size(dop) ~= (1 + 2 .* floor(size(dop) ./ 2))) || ...
   ~any(strcmpi(class(dop), {'logical', 'uint8', 'uint16'}))
    dop = uint16(true(3, 3, 3));
end

% convert to uint16
vol = uint16(vol);
dop = uint16(dop);

% get sizes
vsz = size(vol);
dsz = size(dop);
asz = vsz + dsz - 1;
vsx = vsz(1) - 1;
vsy = vsz(2) - 1;
vsz = vsz(3) - 1;

% special case
if nargin > 2 && ...
    isa(r, 'double') && ...
    numel(r) == 1 && ...
   ~isnan(r) && ...
    r > 0 && ...
    r <= (1 + sqrt(sum(dsz .* dsz)))

    % resolution argument
    if nargin < 4 || ...
       ~isa(res, 'double') || ...
        numel(res) ~= 3 || ...
        any(isinf(res) | isnan(res) | res < 0)
        res = ones(1, 3);
    else
        res = res(:)';
    end

    dsr = ceil(dsz / 2);
    for xc = 1:dsz(1)
        for yc = 1:dsz(2)
            for zc = 1:dsz(3)
                dop(xc, yc, zc) = ...
                    round(3 * max(0, r - sqrt(sum( ...
                    (res .* (dsr - [xc, yc, zc])) .^ 2))));
            end
        end
    end
end

% create erosion array
dil = uint16([]);
dil(1:asz(1), 1:asz(2), 1:asz(3)) = uint16(0);

% perform calculus
for x = 1:dsz(1)
    for y = 1:dsz(2)
        for z = 1:dsz(3)
            if dop(x,y,z) > 0
                dil(x:(x + vsx), y:(y + vsy), z:(z + vsz)) = ...
                    dil(x:(x + vsx), y:(y + vsy), z:(z + vsz)) + ...
                    dop(x, y, z) * vol;
            end
        end
    end
end

% erosion
esx = ceil(dsz(1) / 2);
esy = ceil(dsz(2) / 2);
esz = ceil(dsz(3) / 2);
vol(vol == 0) = 1;
vol = (dil(esx:(esx + vsx), esy:(esy + vsy), esz:(esz + vsz)) >= ...
    ceil(single(vol) .* single(mean(dop(:)))));
