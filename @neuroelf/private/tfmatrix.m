function tfm = tfmatrix(tf)
% tfmatrix  - build 4x4 transformation matrix for different steps
%
% FORMAT:       tfm = tfmatrix(steps)
%
% Input fields:
%
%       steps       1xN struct array with steps (in order of performance)
%        .type      one of
%                   - 't' / 'translate'  - translation
%                   - 'r' / 'rotate'     - rotation
%                   - 's' / 'scale'      - scaling
%                   - 'h' / 'shear'      - shearing
%        .xyz       1x3 double for command
%
% Note: as with SPM, if a multi-rotation is given (rotation around
%       more than one axis), rotation is performed in reverse order,
%       first around Z axis, then around Y axis, and finally around X.
%       angles must be specified in radiens.

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
if nargin < 1 || ...
   ~isstruct(tf) || ...
   ~isfield(tf, 'type') || ...
   ~isfield(tf, 'xyz')
    if nargin > 0 && ...
        isa(tf, 'double')
        tf = tf(:)';
        if numel(tf) < 7
            tf(7:9) = 1;
        end
        tfi = [tf zeros(1, 3)];
        tf    = struct('type', 'h', 'xyz', tfi(10:12));
        tf(2) = struct('type', 's', 'xyz', tfi( 7: 9));
        tf(3) = struct('type', 'r', 'xyz', tfi( 4: 6));
        tf(4) = struct('type', 't', 'xyz', tfi( 1: 3));
    else
        error( ...
            'neuroelf:BadArgument', ...
            'Bad input argument supplied.' ...
        );
    end
end
tf  = tf(:);
tfm = eye(4);
for c = 1:length(tf)
    if ~ischar(tf(c).type) || ...
       ~isa(tf(c).xyz, 'double') || ...
        numel(tf(c).xyz) ~= 3 || ...
        any(isnan(tf(c).xyz) | isinf(tf(c).xyz))
        error( ...
            'neuroelf:BadArgument', ...
            'Bad input argument supplied.' ...
        );
    end
    switch lower(tf(c).type(:)'), case {'t', 'translate'}
        tfm = [[eye(3);0 0 0], [tf(c).xyz(:);1]] * tfm;
    case {'r', 'rotate'}
        r = tf(c).xyz(:);
        if r(3) ~= 0
            tfm = [cos(r(3)) sin(r(3)) 0 0; -sin(r(3)) cos(r(3)) 0 0; 0 0 1 0; 0 0 0 1] * tfm;
        end
        if r(2) ~= 0
            tfm = [cos(r(2)) 0 sin(r(2)) 0; 0 1 0 0; -sin(r(2)) 0 cos(r(2)) 0; 0 0 0 1] * tfm;
        end
        if r(1) ~= 0
            tfm = [1 0 0 0; 0 cos(r(1)) sin(r(1)) 0; 0 -sin(r(1)) cos(r(1)) 0; 0 0 0 1] * tfm;
        end
    case {'s', 'scale'}
        tfm = diag([tf(c).xyz(:)', 1]) * tfm;
    case {'h', 'shear'}
        h = tf(c).xyz(:)';
        tfm = [1 h(1:2) 0; 0 1 h(3) 0; 0 0 1 0; 0 0 0 1] * tfm;
    otherwise
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid transformation type: %s', ...
            tf(c).type(:)' ...
        );
    end
end
