function c = rangegrid(r)
% rangegrid  - convert a flex method range into the corresponding grid
%
% FORMAT:       coords = rangegrid(rangespec)
%
% Input fields:
%
%       rangespec   4xA specification with all(isinf(1, :)) = true
%                   and range defined as values in rows 2:3:4
%
% Output fields:
%
%       coords      coordinates returned by ndgrid (double)

% Version:  v0.9b
% Build:    10081312
% Date:     Aug-13 2010, 12:27 PM EST
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

% argument check
if nargin ~= 1 || ...
   ~isa(r, 'double') || ...
    isempty(r) || ...
    ndims(r) ~= 2 || ...
    size(r, 1) ~= 4 || ...
    size(r, 2) > 4 || ...
   ~all(isinf(r(1, :))) || ...
    any(any(isinf(r(2:4, :))) | any(isnan(r(2:4, :))))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing rangespec argument.' ...
    );
end

% how many dimensiosn
switch (size(r, 2))

    % 1D
    case {1}

        % make column vector
        c = (r(2):r(3):r(4))';

        % nothing else needed
        return;

    % 2D
    case {2}

        % ndgrid with two arguments
        [c{1}, c{2}] = ndgrid( ...
            r(2):r(3):r(4), r(2, 2):r(3, 2):r(4, 2));

    % 3D
    case {3}

        % ndgrid with three arguments
        [c{1}, c{2}, c{3}] = ndgrid( ...
            r(2):r(3):r(4), r(2, 2):r(3, 2):r(4, 2), r(2, 3):r(3, 3):r(4, 3));

    % 4D
    case {4}

        % ndgrid with four arguments
        [c{1}, c{2}, c{3}, c{4}] = ndgrid( ...
            r(2):r(3):r(4), r(2, 2):r(3, 2):r(4, 2), r(2, 3):r(3, 3):r(4, 3), ...
            r(2, 4):r(3, 4):r(4, 4));
end

% concatenate
for cc = 1:numel(c)
    c{cc} = c{cc}(:);
end
c = cat(2, c{:});
