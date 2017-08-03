function xo = vtc_MinFrame(xo)
% VTC::MinFrame  - cut away all-0 data in a VTC
%
% FORMAT:       [vtc = ] vtc.MinFrame;
%
% No input fields.
%
% Output fields:
%
%       vtc         minimized VTC
%
% Using: findfirst.

% Version:  v1.1
% Build:    16021320
% Date:     Feb-13 2016, 8:46 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2011, 2014, 2016, Jochen Weber
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

% neuroelf library
global ne_methods;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;

% this only works if VTCData is not transio
if istransio(bc.VTCData)
    bc.VTCData = resolve(bc.VTCData);
end

% for each spatial dimension, create a logical array
s1 = true(1, size(bc.VTCData, 2));
s2 = true(1, size(bc.VTCData, 3));
s3 = true(1, size(bc.VTCData, 4));

% then look at data
for sc = 1:numel(s1)
    if all(lsqz(bc.VTCData(:, sc, :, :)) == 0)
        s1(sc) = false;
    else
        break;
    end
end
for sc = numel(s1):-1:1
    if s1(sc) && all(lsqz(bc.VTCData(:, sc, :, :)) == 0)
        s1(sc) = false;
    else
        break;
    end
end
for sc = 1:numel(s2)
    if all(lsqz(bc.VTCData(:, :, sc, :)) == 0)
        s2(sc) = false;
    else
        break;
    end
end
for sc = numel(s2):-1:1
    if s2(sc) && all(lsqz(bc.VTCData(:, :, sc, :)) == 0)
        s2(sc) = false;
    else
        break;
    end
end
for sc = 1:numel(s3)
    if all(lsqz(bc.VTCData(:, :, :, sc)) == 0)
        s3(sc) = false;
    else
        break;
    end
end
for sc = numel(s3):-1:1
    if s3(sc) && all(lsqz(bc.VTCData(:, :, :, sc)) == 0)
        s3(sc) = false;
    else
        break;
    end
end

% no changes
if all(s1) && all(s2) && all(s3)
    return;
end

% empty dataset?
if ~any(s1)
    bc.VTCData = reshape(bc.VTCData([]), [size(bc.VTCData, 1), 0, 0, 0]);
    bc.XEnd = bc.XStart;
    bc.YEnd = bc.YStart;
    bc.ZEnd = bc.ZStart;
    xo.C = bc;
    return;
end

% patch data
bc.VTCData = bc.VTCData(:, s1, s2, s3);

% find first indices (shift)
f1 = ne_methods.findfirst(s1) - 1;
f2 = ne_methods.findfirst(s2) - 1;
f3 = ne_methods.findfirst(s3) - 1;

% patch start/end
bc.XStart = bc.XStart + bc.Resolution * f1;
bc.YStart = bc.YStart + bc.Resolution * f2;
bc.ZStart = bc.ZStart + bc.Resolution * f3;
bc.XEnd = bc.XStart + bc.Resolution * sum(s1);
bc.YEnd = bc.YStart + bc.Resolution * sum(s2);
bc.ZEnd = bc.ZStart + bc.Resolution * sum(s3);

% set back
xo.C = bc;


% local implementation of ne_methods.lsqueeze
function v = lsqz(v)
v = v(:);
