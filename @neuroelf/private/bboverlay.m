function [o1, o2, sz] = bboverlay(v1, v2)
% bboverlay  - bounding box overlay for same-resolution objects
%
% FORMAT:       [o1, o2, sz] = bboverlay(v1, v2)
%
% Input fields:
%
%       v1, v2      xff voxel-based objects of same type and resolution
%
% Output fields:
%
%       o1, o2      offset additions into arrays of v1 and v2
%       sz          size of array that overlays

% Version:  v0.9b
% Build:    10080400
% Date:     Aug-03 2010, 11:43 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% check arguments
if nargin < 2 || ...
    numel(v1) ~= 1 || ...
   ~isxff(v1, {'ava', 'cmp', 'glm', 'hdr', 'msk', 'vdw', 'vmp', 'vmr', 'vtc'}) || ...
    numel(v2) ~= 1 || ...
   ~isxff(v2, {'ava', 'cmp', 'glm', 'hdr', 'msk', 'vdw', 'vmp', 'vmr', 'vtc'})
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% depending on filetype
switch (lower(v1.Filetype))

    % for Analyze/NIftI files
    case {'hdr'}
        if ~strcmpi(v2.Filetype, 'hdr') || ...
            any(any(v1.CoordinateFrame.Trf(1:3, 1:3) ~= v2.CoordinateFrame.Trf(1:3, 1:3)))
            error( ...
                'neuroelf:BadArgument', ...
                'Objects do not match in type, resolution and/or orientation.' ...
            );
        end

    % for BVQX data files
    case {'ava', 'cmp', 'glm', 'msk', 'vdw', 'vmp', 'vtc'}
        if any(strcmpi(v2.Filetype, {'hdr', 'vmr'})) || ...
            v1.Resolution ~= v2.Resolution
            error( ...
                'neuroelf:BadArgument', ...
                'Objects do not match in resolution.' ...
            );
        end

    % for BVQX anatomical files
    case {'vmr'}
        if ~strcmpi(v2.Filetype, 'vmr') || ...
            any([v1.FramingCube, v1.VoxResX, v1.VoxResY, v1.VoxResZ] ~= ...
                [v2.FramingCube, v2.VoxResX, v2.VoxResY, v2.VoxResZ])
            error( ...
                'neuroelf:BadArgument', ...
                'Objects do not match in type and/or resolution.' ...
            );
        end

        % compute maximum offset and minimum offset
        on1 = [v1.OffsetX, v1.OffsetY, v1.OffsetZ];
        on2 = [v2.OffsetX, v2.OffsetY, v2.OffsetZ];
        sz1 = size(v1.VMRData);
        sz2 = size(v2.VMRData);
        mon = max(on1, on2);
        moff = min(on1 + sz1, on2 + sz2);
        o1 = (mon - on1);
        o2 = (mon - on2);
        sz = moff - mon;
end

% control for empty arrays
if any(sz < 1)
    sz = [0, 0, 0];
end
