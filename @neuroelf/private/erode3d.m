function vol = erode3d(vol, eop, r, res)
% erode  - perform 3D erosion with operator
%
% FORMAT:       vol = erode3d(vol [, eop, r, res])
%
% Input fields:
%
%       vol         3-D boolean volume
%       eop         erosion operator with odd size ! (default 3x3x3 cube)
%       r           1x1 double radius, if given, weighted erosion
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
    ndims(eop) > 3 || ...
    any(size(eop) ~= (1 + 2 .* floor(size(eop) ./ 2))) || ...
   ~any(strcmpi(class(eop), {'logical', 'uint8', 'uint16'}))
    eop = uint16(true(3, 3, 3));
end

% convert to uint16
vol = uint16(vol);
eop = uint16(eop);

% get sizes
vsz = size(vol);
esz = size(eop);
asz = vsz + esz - 1;
vsx = vsz(1) - 1;
vsy = vsz(2) - 1;
vsz = vsz(3) - 1;

% special case
if nargin > 2 && ...
    isa(r, 'double') && ...
    numel(r) == 1 && ...
   ~isnan(r) && ...
    r > 0 && ...
    r <= (1 + sqrt(sum(esz .* esz)))

    % resolution argument
    if nargin < 4 || ...
       ~isa(res, 'double') || ...
        numel(res) ~= 3 || ...
        any(isinf(res) | isnan(res) | res < 0)
        res = ones(1, 3);
    else
        res = res(:)';
    end

    % calculus of operator
    esr = ceil(esz / 2);
    for xc = 1:esz(1)
        for yc = 1:esz(2)
            for zc = 1:esz(3)
                eop(xc, yc, zc) = ...
                    round(3 * max(0, r - sqrt(sum( ...
                    (res .* (esr - [xc, yc, zc])) .^ 2))));
            end
        end
    end
end

% create erosion array
ero = uint16([]);
ero(1:asz(1), 1:asz(2), 1:asz(3)) = uint16(0);

% perform calculus
for x = 1:esz(1)
    for y = 1:esz(2)
        for z = 1:esz(3)
            if eop(x,y,z) > 0
                ero(x:(x + vsx), y:(y + vsy), z:(z + vsz)) = ...
                    ero(x:(x + vsx), y:(y + vsy), z:(z + vsz)) + ...
                    eop(x, y, z) * vol;
            end
        end
    end
end

% erosion
esx = ceil(esz(1) / 2);
esy = ceil(esz(2) / 2);
esz = ceil(esz(3) / 2);
vol(vol == 0) = 1;
vol = (ero(esx:(esx + vsx), esy:(esy + vsy), esz:(esz + vsz)) >= ...
    (vol .* sum(eop(:))));
