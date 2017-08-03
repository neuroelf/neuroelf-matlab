function mview = renderv3dvm(r, s, t)
% renderv3dvm  - construct a 4x4 transformation matrix, R * S * T
%
% Format:       mview = renderv3dvm(r, s, t)
%
% Input fields:
%
%       r           rotation vector [rx, ry, rz]
%       s           (re-) size vector [sx, sy, sz]
%       t           translation vector [tx, ty, tz]
%
% Output fields:
%
%       mview       4x4 transformation matrix
%
% Function is written by D.Kroon University of Twente (October 2008)

% Version:  v0.9c
% Build:    14022011
% Date:     Feb-20 2014, 11:54 AM EST
% Author:   Dirk-Jan Kroon, University of Twente, Enschede, NL
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://www.mathworks.com/matlabcentral/fileexchange/21993-viewer3d

% Copyright (c) 2008, Dirk-Jan Kroon
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
   ~isa(r, 'double') || ...
    numel(r) ~= 3 || ...
    any(isinf(r) | isnan(r))
    r = [0, 0, 0];
end
R = RotationMatrix(r);
if nargin < 2 || ...
   ~isa(s, 'double') || ...
    numel(s) ~= 3 || ...
    any(isinf(s) | isnan(s) | s == 0)
    s = [1, 1, 1];
end
S = ResizeMatrix(s);
if nargin < 3 || ...
   ~isa(t, 'double') || ...
    numel(t) ~= 3 || ...
    any(isinf(t) | isnan(t))
    t = [0, 0, 0];
end
T = TranslateMatrix(t);

% compute output
mview = R * S * T;

% sub-functions
function R = RotationMatrix(r)
Rx =[1 0 0 0; 0 cosd(r(1)) -sind(r(1)) 0; 0 sind(r(1)) cosd(r(1)) 0; 0 0 0 1];
Ry =[cosd(r(2)) 0 sind(r(2)) 0; 0 1 0 0; -sind(r(2)) 0 cosd(r(2)) 0; 0 0 0 1];
Rz =[cosd(r(3)) -sind(r(3)) 0 0; sind(r(3)) cosd(r(3)) 0 0; 0 0 1 0; 0 0 0 1];
R = Rx * Ry * Rz;

function S = ResizeMatrix(s)
S = [1/s(1) 0 0 0;
     0 1/s(2) 0 0;
     0 0 1/s(3) 0;
     0 0 0 1];

function T = TranslateMatrix(t)
T = [1 0 0 -t(1);
     0 1 0 -t(2);
     0 0 1 -t(3);
     0 0 0 1];
