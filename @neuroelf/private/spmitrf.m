function p = spmitrf(m)
% spmitrf  - uncompile an affine transformation into paramters
%
% FORMAT:       p = spmitrf(m)
%
% Input fields:
%
%       m           4x4 affine transformation matrix
%
% Output fields:
%
%       p           1x4 cell array for use with spmtrf

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
if nargin ~= 1 || ...
   ~isa(m, 'double') || ...
   ~isequal(size(m), [4, 4]) || ...
    any(isinf(m(:)) | isnan(m(:))) || ...
    any(m(4, :) ~= [0, 0, 0, 1])
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing transformation matrix supplied.' ...
    );
end
p = cell(1, 4);

% translation
p{1} = m(1:3, 4)';

% decompose rotation matrix into zooms
p{2} = [0, 0, 0];
r = m(1:3, 1:3);
rc = chol(r' * r);
p{3} = diag(rc)';
if det(r) < 0
    p{3}(1) = -p{3}(1);
end

% shear
rc = diag(diag(rc)) \ rc;
p{4} = rc([4, 7, 8]);

% create untranslated for rotation
r0 = spmtrf([0, 0, 0], p{2:4});
r0 = r0(1:3, 1:3);
r0 = r / r0;

% compute rotation parameters
p{2}(2)  = asin(rang(r0(1, 3)));
if (abs(p{2}(2)) - pi / 2) ^ 2 < 1e-9
	p{2}(1) = 0;
	p{2}(3) = atan2(-rang(r0(2, 1)), rang(-r0(3, 1) / r0(1, 3)));
else
	c = cos(p{2}(2));
	p{2}(1) = atan2(rang(r0(2, 3) / c), rang(r0(3, 3) / c));
	p{2}(3) = atan2(rang(r0(1, 2) / c), rang(r0(1, 1) / c));
end

% There may be slight rounding errors making b>1 or b<-1.
function a = rang(b)
a = min(max(b, -1), 1);
