function m = montagemix(m1, m2, mixval)
% montagemix  - mix two RGB images depending on pixel values
%
% FORMAT:       mixed = montagemix(c1, c2, mixval)
%
% Input fields:
%
%       c1, c2      XxYx3 RGB images (uint8 or single)
%       mixval      value between 0 and 6 where
%                   - 0 any pixel in m1 > 0 wins
%                   - 1 pixel are weighted with max(m1)
%                   - 2 pixel are weighted max(m1) / (max(m1) + max(m2))
%                   - 3 50% mixing
%                   - 4 pixel are weighted max(m2) / (max(m1) + max(m2))
%                   - 5 pixel are weighted with max(m1)
%                   - 6 any pixel in m2 > 0 wins
%
% Output fields:
%
%       mixed       mixed RGB image

% Version:  v0.9c
% Build:    11042819
% Date:     Apr-28 2011, 5:48 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
if nargin < 3 || ...
   (~isa(m1, 'uint8') && ...
    ~isa(m1, 'single')) || ...
   (~isa(m2, 'uint8') && ...
    ~isa(m2, 'single')) || ...
    ndims(m1) > 3 || ...
    ndims(m2) > 3 || ...
    size(m1, 1) ~= size(m2, 1) || ...
    size(m1, 2) ~= size(m2, 2) || ...
   ~any(size(m1, 3) == [1, 3]) || ...
   ~any(size(m2, 3) == [1, 3]) || ...
    numel(mixval) ~= 1 || ...
   ~isa(mixval, 'double') || ...
    isinf(mixval) || ...
    isnan(mixval)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% ensure that all sizes are RGB compliant
if size(m1, 3) < 3
    m1 = m1(:, :, [1, 1, 1]);
end

% handle empty cases
if isempty(m1)
    m = uint8(m1);
    return;
end

% ensure size of m2
if size(m2, 3) < 3
    m2 = m2(:, :, [1, 1, 1]);
end

% for maximal visibility of overlay
if mixval >= 6

    % mixing is a on/off function depending on m2
    mix = single(~any(m2 > 0, 3));

% for next tier of values ([5 ... 6[)
elseif mixval >= 5

    % start with the lower boundary value
    mix = single(1 / 255) .* single(255 - max(m2, [], 3));

    % if between 5 and 6
    if mixval > 5

        % submix with higher value
        mixh = single(~any(m2 > 0, 3));
        mix = (mixval - 5) .* mixh + (6 - mixval) .* mix;
    end

% for next tier of values ([3 ... 5[)
elseif mixval > 3

    % compute required terms
    mix1 = single(1 / 255) .* single(max(m1, [], 3));
    mix2 = single(1 / 255) .* single(max(m2, [], 3));
    mix = limitrangec(1 - mix2 ./ (mix1 + mix2), 0, 1, 0);

    % sub-tiers
    if mixval > 4
        mix = (mixval - 4) .* (1 - mix2) + (5 - mixval) .* mix;
    elseif mixval < 4
        mix = (mixval - 3) .* mix + (4 - mixval) .* single(0.5);
    end

% for mixval == 3
elseif mixval == 3

    % 50/50 mix
    mix = single(0.5);

% for next tier of values (]1 ... 3[)
elseif mixval > 1

    % compute required terms
    mix1 = single(1 / 255) .* single(max(m1, [], 3));
    mix2 = single(1 / 255) .* single(max(m2, [], 3));
    mix = limitrangec(mix1 ./ (mix1 + mix2), 0, 1, 0);

    % sub-tiers
    if mixval > 2
        mix = (mixval - 2) .* single(0.5) + (3 - mixval) .* mix;
    elseif mixval < 2
        mix = (mixval - 1) .* mix + (2 - mixval) .* mix1;
    end

% for mixval == 1
elseif mixval > 0

    % start with the higher boundary value
    mix = single(1 / 255) .* single(max(m1, [], 3));

    % if between 5 and 6
    if mixval < 1

        % submix with higher value
        mix = mixval .* mix + (1 - mixval) .* single(any(m1 > 0, 3));
    end

% for maximal visibility of m1
else
    mix = single(any(m1 > 0, 3));
end

% perform mixing
if numel(mix) > 1
    mix = mix(:, :, [1, 1, 1]);
end
m = uint8(mix .* single(m1) + single(1 - mix) .* single(m2));
