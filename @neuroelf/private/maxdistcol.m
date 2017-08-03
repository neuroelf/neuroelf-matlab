function mdc = maxdistcol(c1, c2)
% maxdistcol  - compute color with maximal color distance to both (in HSV)
%
% FORMAT:       mdc = maxdistcol(c1, c2)
%
% Input fields:
%
%       c1, c2      [0 .. 255] RGB-based or [0 .. 1] HSV-based colors
%
% Output fields:
%
%       mdc         maximally distant color

% Version:  v0.9c
% Build:    14012213
% Date:     Jan-22 2014, 1:06 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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
if nargin ~= 2 || ...
   ~strcmp(class(c1), class(c2)) || ...
   ~isnumeric(c1) || ...
   ~isequal(size(c1), size(c2)) || ...
    size(c1, ndims(c1)) ~= 3 || ...
    any(isinf(c1(:)) | isnan(c1(:)) | c1(:) < 0 | c1(:) >= 256 | isinf(c2(:)) | isnan(c2(:)) | c2(:) < 0 | c2(:) >= 256)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% size
csize = size(c1);
cdim = numel(csize);
c1 = reshape(c1, prod(csize(1:cdim-1)), csize(cdim));
c2 = reshape(c2, prod(csize(1:cdim-1)), csize(cdim));

% convert
c1rgb = false;
if any(c1(:) > 1)
    c1rgb = true;
    if ~isa(c1, 'uint8')
        c1 = hsvconv(uint8(round(c1)), 2);
    else
        c1 = hsvconv(c1, 2);
    end
end
if any(c2(:) > 1)
    if ~isa(c2, 'uint8')
        c2 = hsvconv(uint8(round(c2)), 2);
    else
        c2 = hsvconv(c2, 2);
    end
end

% compute average
mdc = 0.5 .* (c1 + c2);

% compute max. distance to component
maxd = max(min(abs(c1(:, 1) - mdc(:, 1)), abs((c1(:, 1) + 1) - mdc(:, 1))), ...
    min(abs(c2(:, 1) - mdc(:, 1)), abs((c2(:, 1) + 1) - mdc(:, 1))));

% for smaller < 0.25
rmd = (maxd < 0.25);
mdc(rmd, 1) = mdc(rmd, 1) + 0.5;
mdc(mdc(:, 1) > 1, 1) = mdc(mdc(:, 1) > 1, 1) - 1;

% convert back
if c1rgb
    mdc = hsvconv(mdc, 1);
end

% reshape
mdc = reshape(mdc, csize);
