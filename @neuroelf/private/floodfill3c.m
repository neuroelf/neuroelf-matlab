function [v, n] = floodfill3c(v, s, m)
% floodfill3c  - MEX implementation of flood-filling algorithm
%
% FORMAT:       [v, n] = floodfill3c(v, s, m);
%
% Input fields:
%
%       v           input/output volume
%       s           start index (1-based)
%       m           method:
%                   1 face-connectivity (3D)
%                   2 edge-connectivity (3D)
%                   3 vertex-connectivity (3D)
%                   4 face-connectivity (XY-slice)
%                   5 face-connectivity (XY-slice, zstart = :)
%                   6 edge-connectivity (XY-slice)
%                   7 edge-connectivity (XY-slice, zstart = :)
%
% Output fields:
%
%       v           flooded data (only voxels connected to start)
%       n           number of voxels

% Version:  v0.9a
% Build:    10073109
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
if nargin ~= 3 || ...
   (~islogical(v) && ...
    ~strcmpi(class(v), 'uint8')) || ...
    isempty(v) || ...
    ndims(v) ~= 3 || ...
   ~isa(s, 'double') || ...
    numel(s) ~= 3 || ...
    any(isinf(s) | isnan(s) | s < 1 | s > size(v) | s ~= fix(s)) || ...
   ~isa(m, 'double') || ...
    numel(m) ~= 1 || ...
    isinf(m) || ...
    isnan(m) || ...
   ~any(m == (1:7))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% bail out for MEX
error( ...
    'neuroelf:MEXMissing', ...
    'This is a compiled function, but the MEX file is missing.' ...
);
