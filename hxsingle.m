%HXSINGLE  Transform a single into a 1x8 (or longer) ASCII string.
%   HXSTRING = HXSINGLE(VALUE) converts the number(s) in VALUE to their
%   hexadecimal representation, such that they can be stored and later
%   converted back to numbers without loss of precision.
%
%   VALUE = HXSINGLE(HXSTRING) converts the ASCII representation in
%   HXSTRING back to its single representation.
%
%   Importantly, for more than one value (vector, matrix, or array),
%   VALUE will be linearized with (:) first; hence for proper conversion
%   to the original VALUE, it is necessary to separately store the SIZE.
%
%   See also ANY2ASCII, SIZE, EVAL, HXDOUBLE.

% Version:  v1.1
% Build:    16031110
% Date:     Mar-11 2016, 10:06 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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
function hxs = hxsingle(inp)

% persistent lookup values
persistent hxs_hxr;
if isempty(hxs_hxr)
    hxs_hxr = [ ...
        nan * ones(1, 48) ...
          0   1   2   3   4   5   6   7   8   9 nan nan nan nan nan nan ...
        nan  10  11  12  13  14  15 nan nan nan nan nan nan nan nan nan ...
        nan * ones(1, 16) ...
        nan  10  11  12  13  14  15 nan nan nan nan nan nan nan nan nan ...
        nan * ones(1,144)];
end

% number of arguments
if nargin ~= 1
    error('neuroelf:general:badNumberOfInputs', 'Invalid number of inputs.');
end

% for a double value/array
if isa(inp, 'single')

    % pass on to sprintfbx
    hxs = sprintfbx(inp);

% for a char
elseif ischar(inp)

    % linearize
    inp = inp(:)';

    % number of double values
    ns  = floor(numel(inp) / 8);

    % prepare output
    hxs = nan .* zeros(1, ns);

    % iterate
    for hxc = 1:ns

        % get (lower-char) string
        ipp = lower(inp((hxc-1)*8+1:hxc*8));

        % NaN?
        if strcmp(ipp, 'ffc00000')
            continue;
        end

        % only allow ASCII, and add one for indexing into lookup
        ipi = hxs_hxr(min(double(ipp), 255) + 1);
        if any(isnan(ipi))
            continue;

        % all 0
        elseif all(ipi == 0)
            hxs(hxc) = 0;
            continue;
        end

        % compute exponent
        iex = 32 * ipi(1) + 2*ipi(2) + fix(ipi(3) / 8);

        % get sign correct
        if iex > 255
            isg = -1;
            iex = iex - 256;
        else
            isg = 1;
        end

        % compute mantissa
        imt = 1 + bitand(7, ipi(3)) / 8 + hex2dec(ipp(4:8)) / (2 .^ 23);
        hxs(hxc) = isg * imt * 2 .^ (iex - 127);
    end

    % convert at the end
    hxs = single(hxs);

% unsupported class
else
    error('neuroelf:hxsingle:invalidInputType', 'Invalid input type: %s.', class(inp));
end
