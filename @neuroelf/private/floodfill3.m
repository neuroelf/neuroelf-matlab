function [vol, marked] = floodfill3(vol, sx, sy, sz, meth)
% floodfill3  - keep adjacent values from seed point
%
% FORMAT:       [vol, marked] = floodfill3(vol, sx, sy, sz [, meth])
%
% Input fields:
%
%       vol         XxYxZ logical/uint8 array, true/1 voxels are linked
%       sx, sy, sz  start index (1-based)
%       meth        either of {'face'}, 'edge', 'vertex', 'xyface'
%
% Output fields:
%
%       vol         same as input array with only one cluster
%       marked      number of marked elements
%
% Note: if available uses the MEX function floodfill3c.

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 2:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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
   (~islogical(vol) && ...
    ~strcmpi(class(vol), 'uint8')) || ...
    isempty(vol) || ...
   ~isa(sx, 'double') || ...
    numel(sx) ~= 1 || ...
    isnan(sx) || ...
   ~isa(sy, 'double') || ...
    numel(sy) ~= 1 || ...
    isnan(sy) || ...
   ~isa(sz, 'double') || ...
    numel(sz) ~= 1 || ...
    isnan(sz) || ...
    sx < 1 || ...
    sx > size(vol, 1) || ...
    sy < 1 || ...
    sy > size(vol, 2) || ...
    sz < 1 || ...
    sz > size(vol, 3) || ...
    any([sx, sy, sz] ~= fix([sx, sy, sz]))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 5 || ...
   ~ischar(meth) || ...
    isempty(meth)
    meth = 'f';
else
    meth = lower(meth(1));
end
if ~any('efvx4567' == meth)
    meth = 'f';
end

switch (meth)
    case {'e'}
        md = 9;
        meth = 2;
    case {'f'}
        md = 3;
        meth = 1;
    case {'v'}
        md = 13;
        meth = 3;
    case {'x', '4'}
        md = 2;
        meth = 4;
    case {'5'}
        md = 2;
        meth = 5;
    case {'6'}
        md = 4;
        meth = 6;
    case {'7'}
        md = 4;
        meth = 7;
end

% C-function available
if exist('floodfill3c', 'file') == 3

    % call floodfill3c and return
    try
        [vol, marked] = floodfill3c(vol, [sx, sy, sz], meth);
    catch ne_eo;
        rethrow(ne_eo);
    end
    return;
end

% otherwise procede ...
st = size(vol);
vx = st(1);
vy = st(2);
if numel(st) > 2
    vz = st(3);
else
    vz = 1;
    st(3) = vz;
end

% check whether sx, sy, sz are in volume
if ~vol(sx, sy, sz) && ...
   ~any([5, 7] == meth)
    vol(:) = false;
    marked = 0;
    return;
end

% convert vol to uint8
vol = uint8(vol);
if ~any([5, 7] == meth)
    vol(sx, sy, sz) = 2;
    marked = 1;
else
    marked = 0;
    for dc = 1:size(vol, 3)
        if vol(sx, sy, dc) == 1
            vol(sx, sy, dc) = 2;
            marked = marked + 1;
        end
    end
    sz = find(squeeze(vol(sx, sy, :)) > 1);
    sx = sx * osz(sy);
    sy = sy * osz(sz);
end

% and list of directions (1-3 face, 1-9 edge, 1-13 vertices)
%  1: (  +x  ,  -x  )  2: (  +y  ,  -y  )  3: (  +z  ,  -z  )
%  4: ( +x+y , -x-y )  5: ( +x-y , -x+y )
%  6: ( +x+z , -x-z )  7: ( +x-z , -x+z )
%  8: ( +y+z , -y-z )  9: ( +y-z , -y+z )
% 10: (+x+y+z,-x-y-z) 11: (+x+y-z,-x-y+z)
% 12: (+x-y+z,-x+y-z) 13: (+x-y-z,-x+y+z)
dl = [...
     1,  0,  0;  0,  1,  0;  0,  0,  1; ...
     1,  1,  0;  1, -1,  0; ...
     1,  0,  1;  1,  0, -1; ...
     0,  1,  1;  0,  1, -1; ...
     1,  1,  1;  1,  1, -1; ...
     1, -1,  1;  1, -1, -1  ...
    ];
if (md == 4)
    dl(3, :) = [];
end

% setup hook point arrays
hx = zeros(1, md * ceil(numel(vol) ^ 0.66));
hy = zeros(1, md * ceil(numel(vol) ^ 0.66));
hz = zeros(1, md * ceil(numel(vol) ^ 0.66));

% iterate while new points found
while ~isempty(sx)

    % get indices of to checking voxels
    hc = 1;
    for dc = 1:md
        dr = dl(dc, :);
        dx = dr(1);
        dy = dr(2);
        dz = dr(3);
        if dx ~= 0
            tx = [sx - dx, sx + dx];
        else
            tx = [sx, sx];
        end
        if dy ~= 0
            ty = [sy - dy, sy + dy];
        else
            ty = [sy, sy];
        end
        if dz ~= 0
            tz = [sz - dz, sz + dz];
        else
            tz = [sz, sz];
        end
        gt = find(tx > 0 & tx <= vx & ty > 0 & ty <= vy & tz > 0 & tz <= vz);
        gs = numel(gt);
        hx(hc:(hc + gs - 1)) = tx(gt);
        hy(hc:(hc + gs - 1)) = ty(gt);
        hz(hc:(hc + gs - 1)) = tz(gt);
        hc = hc + gs;
    end

    % only accept where further coordinates
    hc = hc - 1;
    tc = sub2ind(st, hx(1:hc), hy(1:hc), hz(1:hc));
    tc = unique(tc(vol(tc) == 1));

    % set these to 2 as well
    vol(tc) = 2;
    marked = marked + numel(tc);

    % then continue with these
    [sx, sy, sz] = ind2sub(st, tc);
end

% convert back
vol = (vol == 2);
