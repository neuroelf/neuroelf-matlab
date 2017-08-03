function g = maxgrad(d)
% maxgrad  - return maximum direct neighbor gradient (diff) per voxel
%
% FORMAT:       mg = maxgrad(d)
%
% Input fields:
%
%       d           3D data
%
% Output fields:
%
%       mg          max. gradient

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
    ndims(d) ~= 3 || ...
   ~isnumeric(d)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument in call.' ...
    );
end

% make double
if ~isa(d, 'double')
    d = double(d);
end

% create three convolution operators
cop = zeros(3, 3, 3);
cop(2, 2, 2) =  1;

% compute first direct gradient
cop(1, 2, 2) = -1;
g = abs(conv3d(d, cop, 2));
cop(1, 2, 2) = 0;

% then max with others
cop(3, 2, 2) = -1;
g = max(g, abs(conv3d(d, cop, 2)));
cop(3, 2, 2) = 0;
cop(2, 1, 2) = -1;
g = max(g, abs(conv3d(d, cop, 2)));
cop(2, 1, 2) = 0;
cop(2, 3, 2) = -1;
g = max(g, abs(conv3d(d, cop, 2)));
cop(2, 3, 2) = 0;
cop(2, 2, 1) = -1;
g = max(g, abs(conv3d(d, cop, 2)));
cop(2, 2, 1) = 0;
cop(2, 2, 3) = -1;
g = max(g, abs(conv3d(d, cop, 2)));
