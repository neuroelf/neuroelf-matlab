function bitdump(buffer)
% bitdump  - dump bits in a uint buffer
%
% FORMAT:       bitdump(buffer)
%
% Input fields:
%
%       buffer      integer array buffer (only uint8/16/32 !)
%
% No output fields. (will display in console window)

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
   ~isnumeric(buffer) || ...
   ~any(strcmpi(class(buffer), {'uint16', 'uint32', 'uint8'}))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% check class and size
cl = lower(class(buffer));
ps = numel(buffer);
m1 = uint32(255 * 256 * 256);
m2 = uint32(255 * 256);
m3 = uint32(255);
m4 = uint16(255);
switch (cl)
    case {'uint8'}
        newbuffer = buffer(:); % already good
    case {'uint16'}
        % convert to uint8
        newbuffer = uint8(0);
        newbuffer(2 * ps) = 0;
        for bc = 1:ps
            newbuffer(2 * bc - 1) = uint8(bitshift(buffer(bc), -8));
            newbuffer(2 * bc) = uint8(bitand(buffer(bc), m4));
        end
    case {'uint32'}
        % convert to uint8
        newbuffer = uint8(0);
        newbuffer(4 * ps) = 0;
        for bc = 1:ps
            newbuffer(4 * bc - 3) = uint8(bitshift(buffer(bc), -24));
            newbuffer(4 * bc - 2) = uint8(bitshift(bitand(buffer(bc), m1), -16));
            newbuffer(4 * bc - 1) = uint8(bitshift(bitand(buffer(bc), m2), -8));
            newbuffer(4 * bc) = uint8(bitand(buffer(bc), m3));
        end
end
ps = numel(newbuffer);

b(1) = uint8(128);
b(2) = uint8(64);
b(3) = uint8(32);
b(4) = uint8(16);
b(5) = uint8(8);
b(6) = uint8(4);
b(7) = uint8(2);
b(8) = uint8(1);

for bc = 1:ps
    disp(sprintf('%08X: %-8s', bc - 1, ...
        sprintf('%c', 46 - 4 * (double(bitand(b, newbuffer(bc))) ~= 0))));
end
