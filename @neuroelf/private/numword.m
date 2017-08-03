function w = numword(n)
% numword  - convert a number into English words
%
% FORMAT:       words = numword(number)
%
% Input fields:
%
%       number      1x1 double number
%
% Output fields:
%
%       words       1xC string representing the number

% Version:  v0.9c
% Build:    13020216
% Date:     Feb-02 2013, 4:15 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, 2013, Jochen Weber
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

% persistent variable
persistent nww;
if isempty(nww)
    nww = struct( ...
        'tw', {{'one', 'two', 'three', 'four', 'five', 'six', 'seven', ...
            'eight', 'nine', 'ten', 'eleven', 'twelve', 'thirteen', ...
            'fourteen', 'fifteen', 'sixteen', 'seventeen', 'eighteen', ...
            'nineteen'}}, ...
        'ten', {{'twenty', 'thirty', 'forty', 'fifty', 'sixty', 'seventy', ...
            'eighty', 'ninety'}});
end

% argument check
if nargin ~= 1 || ...
   ~isa(n, 'double') || ...
    numel(n) ~= 1 || ...
   (~isinf(n) && ...
    abs(n) >= 1e9)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% special cases
if isinf(n)
    if n < 0
        w = 'minus infinity';
    else
        w = 'infinity';
    end
    return;
elseif isnan(n)
    w = 'not-a-number';
    return;

% negative numbers
elseif n < 0
    w = ['minus ' numword(-n)];
    return;
end

% round if discrepancy is small
if abs(n - fix(n)) < sqrt(eps)
    n = fix(n);
end

% fractional values (alone)
if abs(n) < 1
    w = 'point';
    fdigs = 8;
    while abs(n) > sqrt(eps) && ...
        fdigs > 0
        n = 10 * n;
        nf = floor(n);
        n = n - nf;
        if nf > 0
            w = [w ' ' nww.tw{nf}];
        else
            w = [w ' zero'];
        end
        fdigs = fdigs - 1;
    end
    return;

% fractional value (with full value)
elseif n ~= fix(n)
    fn = floor(n);
    w = [numword(fn) ' ' numword(n - fn)];
    return;
end

% small numbers
if n == 0
    w = 'zero';
elseif n < 20
    w = nww.tw{n};
elseif n < 100
    if mod(n, 10) == 0
        w = nww.ten{0.1 * n - 1};
    else
        w = sprintf('%s-%s', nww.ten{floor(0.1 * n) - 1}, nww.tw{mod(n, 10)});
    end
else
    w = '';
    if n >= 1000000
        w = [w numword(floor(n / 1000000)) ' million '];
        n = mod(n, 1000000);
    end
    if n >= 1000
        w = [w numword(floor(n / 1000)) ' thousand '];
        n = mod(n, 1000);
    end
    if n >= 100
        w = [w numword(floor(n / 100)) ' hundred '];
        n = mod(n, 100);
        if n > 0
            w = [w 'and '];
        end
    end
    if n > 0
        w = [w numword(n)];
    end
    if w(end) == ' '
        w(end) = [];
    end
end
