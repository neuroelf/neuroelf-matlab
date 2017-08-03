%SPRINTFBX  Convert single or double to char.
%   HXSTRING = SPRINTFBX(VALUE) converts the double value (or single value
%   or array) into a string using either SPRINTF('%bx', VALUE) or a more
%   manual conversion.
%
%   See also SPRINTF.

% Version:  v1.1
% Build:    16031110
% Date:     Mar-11 2016, 10:04 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 -2016, Jochen Weber
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

% function definition
function a = sprintfbx(d)

% not valid if not either double or single
if nargin < 1 || (~isa(d, 'double') && ~isa(d, 'single'))
    error('neuroelf:sprintfbx:invalidInputType', 'Invalid input type: %s.', class(d));
end

% empty, create empty string
if isempty(d)
    a = char(zeros(1, 0));
    return;
end

% linearize
d = d(:);

% for doubles
if isa(d, 'double')

    % pass on to sprintf (fast!)
    a = sprintf('%bx', real(d));

% for singles
else

    % number of values
    nd = numel(d);

    % create output (as char array)
    a = repmat(' ', nd, 8);

    % keep track of Inf and NaN values
    di = isinf(d);
    dn = isnan(d);

    % record sign as well
    s = sign(d);
    n = (s < 0);
    if any(di)
        p = (s > 0);
    end

    % take absolute value
    d = abs(double(d));

    % compute exponent
    x = floor(log2(d));
    e = 127 + x;

    % set to 0 where necessary
    e(di | dn | s == 0) = 0;

    % correct negatives
    e = e + 256 .* n;

    % determines first two hex places
    a(:, 1:2) = reshape(sprintf('%02x', floor(0.5 * e)), 2, nd)';

    % correct numbers
    d = d .* (2 .^ (-x)) - 1;
    d(s == 0) = 0;
    d = (2 ^ 23) .* d;
    e = 8388608 .* mod(e, 2);

    % store rest
    a(:, 3:8) = reshape(sprintf('%06x', e + floor(d)), 6, nd)';

    % infinities
    if any(di)
        a(di & n, :) = repmat('ff800000', sum(di & n), 1);
        a(di & p, :) = repmat('7f800000', sum(di & p), 1);
    end

    % NaNs
    if any(dn)
        a(dn, :) = repmat('ffc00000', sum(dn), 1);
    end

    % reshape as single string
    a = reshape(a', 1, 8 * nd);
end
