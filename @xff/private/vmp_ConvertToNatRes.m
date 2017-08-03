function xo = vmp_ConvertToNatRes(xo, res)
% VMP::ConvertToNatRes  - convert to native resolution map
%
% FORMAT:       [vmp] = vmp.ConvertToNatRes([res]);
%
% Input fields:
%
%       res         resolution, if not given, try first 3 then 2
%
% Note: the changes are in-place, so, no new object is being created

% Version:  v1.1
% Build:    16021320
% Date:     Feb-13 2016, 8:24 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
if bc.NativeResolutionFile
    return;
end

% guess resolution
bbox = aft_BoundingBox(xo);
osz = diff(bbox.BBox) - 1;
if nargin < 2 || numel(res) ~= 1 || ~isnumeric(res) || ~any((2:12) == res)

    % test 3
    if all(mod(osz, 3) == 0)
        res = 3;
    else
        res = 2;
    end
end
if any(mod(osz, res) > 0)
    error('neuroelf:xff:badArgument', 'VMP cannot be resampled, irregular grid size.');
end

% make changes
ofv = bc.FileVersion;
bc.Resolution = res;
bc.NativeResolutionFile = 1;
bc.FileVersion = 5;
xo.C = bc;
vmp_Update(xo, 'FileVersion', struct('type', '.', 'subs', 'FileVersion'), ofv);
