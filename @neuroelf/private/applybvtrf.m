function tc = applybvtrf(c, t, it, bb)
% applybvtrf  - apply BV transformation matrix to system coords
%
% FORMAT:       tc = applybvtrf(c, t [, it [, bb]])
%
% Input fields:
%
%       c           Cx3 coordinate list
%       t           4x4 transformation matrix (from TRF file, type 2)
%       it          inverse transformation flag, default: true (!)
%       bb          optional .BoundingBox (default: VMR/VTC space 1mm)
%
% Output fields:
%
%       tc          transformed coordinates
%
% Note: since BV stores the transformation as to sample the intensities
%       the matrix is actually already an inverse of the configured
%       matrix. To keep this logic, the it parameter must be set to false

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
if nargin < 2 || ...
   ~isa(c, 'double') || ...
   ~isa(t, 'double') || ...
    size(c, 2) ~= 3 || ...
    numel(t) ~= 16 || ...
    any(size(t) ~= 4) || ...
    any(isinf(t(:)) | isnan(t(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing input argument.' ...
    );
end
if nargin < 3 || ...
   ~islogical(it) || ...
    isempty(it)
    it = true;
end
if nargin < 4
    bb = struct( ...
        'BBox',   [0, 0, 0; 255, 255, 255], ...
        'DimXYZ', [256, 256, 256], ...
        'FCube',  [256, 256, 256], ...
        'ResXYZ', [1, 1, 1]);
elseif ~isstruct(bb) || ...
    numel(bb) ~= 1 || ...
   ~isfield(bb, 'BBox') || ...
   ~isfield(bb, 'DimXYZ') || ...
   ~isfield(bb, 'FCube') || ...
   ~isfield(bb, 'ResXYZ')
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid .BoundingBox in bb argument.' ...
    );
end

% swap and append columns, subtract center coordinate
if any(bb.ResXYZ < 1)
    cc = 0.5 .* (bb.ResXYZ .* (bb.FCube - 1));
else
    cc = 127.5;
end
if all(cc == cc(1))
    tc = [c(:, [2, 3, 1]) - cc(1), ones(size(c, 1), 1)];
else
    tc = [c(:, [2, 3, 1]) - cc(ones(size(c, 1), 1), :), ones(size(c, 1), 1)];
end

% apply transformation
if it(1)
    tc = tc * inv(t)';
else
    tc = tc * t';
end

% swap and remove columns again and add center again
if all(cc == cc(1))
    tc = cc(1) + tc(:, [3, 1, 2]);
else
    tc = cc(ones(size(c, 1), 1), :) + tc(:, [3, 1, 2]);
end
