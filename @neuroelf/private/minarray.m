function [d, off, sz] = minarray(d, minthresh, maxthresh, fringe)
% minarray  - minimize array size
%
% FORMAT:       [d, off, sz] = minarray(d, minthresh, maxthresh, fringe)
%
% Input fields:
%
%       d           2d or 3d numeric data array
%       minthresh   minimum value (default: 0, > comparison)
%       maxthresh   maximum value (default: Inf, <= comparison)
%       fringe      added after border find (default: 0)
%
% Output fields:
%
%       d           minimized array
%       off         offset (1-based)
%       sz          new array size

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 11:08 PM EST
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
if nargin < 1 || ...
   (~isnumeric(d) && ...
    ~islogical(d))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
sz = size(d);
if isempty(d)
    off = [1, 1, 1];
    return;
end
if nargin < 2 || ...
   ~isnumeric(minthresh) || ...
    isempty(minthresh)
    minthresh = 0;
else
    minthresh = minthresh(1);
end
try
    minthresh = eval([class(d) '(minthresh)']);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end
if nargin < 3 || ...
   ~isnumeric(maxthresh) || ...
    isempty(maxthresh)
    maxthresh = Inf;
else
    maxthresh = maxthresh(1);
end
try
    maxthresh = eval([class(d) '(maxthresh)']);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end
if nargin < 4 || ...
   ~isnumeric(fringe) || ...
    numel(fringe) ~= 1 || ...
    isnan(fringe) || ...
    fringe < 0 || ...
    fringe > 512 || ...
    fringe ~= fix(fringe)
    fringe = 0;
end

% first find X border
nxc = [sz(2) * sz(3), 1];
for x1 = 1:sz(1)
    if any(reshape(d(x1, :, :), nxc) > minthresh) && ...
        any(reshape(d(x1, :, :), nxc) <= maxthresh)
        break;
    end
end
for xe = sz(1):-1:x1
    if any(reshape(d(xe, :, :), nxc) > minthresh) && ...
        any(reshape(d(xe, :, :), nxc) <= maxthresh)
        break;
    end
end
sxc = x1:xe;

% find Y border
nyc = [numel(sxc) * sz(3), 1];
for y1 = 1:sz(2)
    if any(reshape(d(sxc, y1, :), nyc) > minthresh) && ...
        any(reshape(d(sxc, y1, :), nyc) <= maxthresh)
        break;
    end
end
for ye = sz(2):-1:y1
    if any(reshape(d(sxc, ye, :), nyc) > minthresh) && ...
        any(reshape(d(sxc, ye, :), nyc) <= maxthresh)
        break;
    end
end
syc = y1:ye;

% find Z border
nzc = [numel(sxc) * numel(syc), 1];
for z1 = 1:sz(3)
    if any(reshape(d(sxc, syc, z1), nzc) > minthresh) && ...
        any(reshape(d(sxc, syc, z1), nzc) <= maxthresh)
        break;
    end
end
for ze = sz(3):-1:z1
    if any(reshape(d(sxc, syc, ze), nzc) > minthresh) && ...
        any(reshape(d(sxc, syc, ze), nzc) <= maxthresh)
        break;
    end
end

% fringe
x1 = max(1, x1 - fringe);
xe = min(sz(1), xe + fringe);
y1 = max(1, y1 - fringe);
ye = min(sz(2), ye + fringe);
z1 = max(1, z1 - fringe);
ze = min(sz(3), ze + fringe);

% reframe data
d = d(x1:xe, y1:ye, z1:ze);
off = [x1, y1, z1];
sz = size(d);
