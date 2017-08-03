function hd = hexdump(av)
% hexdump  - dump hexadecimal values into string
%
% FORMAT:       hd = hexdump(asciivals)
%
% Input fields:
%
%       asciivals   1xN double or char array
%
% Output fields:
%
%       hd          hexdumped string
%
% See also bitdump

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
   (~ischar(av) && ...
    ~isa(av, 'double'))
    error( ...
        'neuroelf:BadArgument', ...
        'Missing or wrong argument.' ...
    );
end

% get lenth and size
l = length(av);
p = 1;

% prepare output
hd = '';

% iterate until end
while (p < l)

    % get stretch
    pa = av(p:min(p+15, l));

    % print into string
    hd = [hd sprintf('%08X:  %-50s %-16s\n', p - 1, ...
        sprintf('%02X ', double(pa)), ...
        saveascii(pa))];

    % go on
    p = p + 16;
end


% internal function to remove ascii characters
function sa = saveascii(av)
    av(av < 32 | av > 127) = 46;
    sa = sprintf('%c', av);
% end of function sa = saveascii(av)
