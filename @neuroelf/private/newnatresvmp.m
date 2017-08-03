function vmp = newnatresvmp(bbox, res, mtype)
% newnatresvmp  - create natural resolution VMP object
%
% FORMAT:       vmp = newnatresvmp([bbox, res [, mtype]])
%
% Input fields:
%
%       bbox        bounding box (from BoundingBox)
%       res         resolution (1x1 double, from BoundingBox)
%       mtype       if given 1xN numeric map type (default [1])
%
% Output fields
%
%       vmp         NR-VMP object with requested properties

% Version:  v1.1
% Build:    16061423
% Date:     Jun-14 2016, 11:54 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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
if nargin < 1 || ...
   ~isa(bbox, 'double') || ...
    ndims(bbox) > 2 || ...
    any(size(bbox) ~= [2, 3]) || ...
    any(bbox(:) < 0 | bbox(:) > 255 | isnan(bbox(:))) || ...
    any(bbox(1, :) >= bbox(2, :))
    bbox = [44, 38, 44; 242, 194, 212];
end
if nargin < 2 || ...
   ~isa(res, 'double') || ...
    numel(res) ~= 1 || ...
   ~any((1:12) == res)
    res = 3;
end
if nargin < 3 || ...
   ~isa(mtype, 'double') || ...
    isempty(mtype) || ...
    any(isinf(mtype(:)) | isnan(mtype(:)) | mtype(:) < 0 | mtype(:) ~= fix(mtype(:)))
    mtype = 1;
else
    mtype = mtype(:);
end
if res > 1
    bbox(2, :) = bbox(2, :) + 0.5;
end
vsz = round(diff(bbox) ./ res);
if res > 1
    bbox(2, :) = bbox(1, :) + res .* vsz;
end
if any(vsz == 0)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% create object
xffroot = xff();
vmp = xff('new:vmp');
vmpu = xffroot.UpdateState('vmp', false);
vmp.XStart = bbox(1, 1);
vmp.XEnd = bbox(2, 1);
vmp.YStart = bbox(1, 2);
vmp.YEnd = bbox(2, 2);
vmp.ZStart = bbox(1, 3);
vmp.ZEnd = bbox(2, 3);
vmp.Resolution = res;
vmp.Map(1).VMPData = single(zeros(diff(bbox) ./ res));
ofv = vmp.FileVersion;
if res > 1
    vmp.NativeResolutionFile = 1;
    vmp.FileVersion = 6;
    vmp.Update('FileVersion', struct('type', '.', 'subs', 'FileVersion'), ofv);
else
    vmp.NativeResolutionFile = 0;
    vmp.FileVersion = 4;
    vsz = vsz + 1;
end
vmp.Map.VMPData = single(zeros(vsz));
vmp.Map.Type = 1;
vmp.Map.DF2 = 0;
vmp.Map.BonferroniValue = prod(vsz);
vmp.Map.Name = '';
vmp.Map(2:numel(mtype)) = vmp.Map(1);
vmp.NrOfMaps = numel(vmp.Map);

% set type(s)
for mc = 1:numel(mtype)
    vmp.Map(mc).Type = mtype(mc);
end

% make setting as before
xffroot.UpdateState('vmp', vmpu);
