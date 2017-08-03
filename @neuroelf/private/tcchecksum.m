function tcc = tcchecksum(src)
% tcchecksum  - checksum using textcrypt (not reliable for binary data)
%
% FORMAT:       tcc = tcchecksum(src)
%
% Input fields:
%
%       src         text-based input
%
% Output fields:
%
%       tcc         checksum string (32-digit hex)

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
   ~ischar(src)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% create long-enough size argument to attach to src to get unique hash
srcs = size(src);
srcs(end+1:64) = 1;
srcs = sprintf('%d-', srcs);
srcs(end) = [];

% run through textcrypt with password 'tcchecksum'
tcc = textcrypt([src(:)', srcs], 'tcchecksum');

% make sure it's a multiple of 16 characters
if mod(numel(tcc), 16) ~= 0
   tcc(end+1:16*ceil(0.0625 * numel(tcc))) = ' ';
end

% reshape to a 16-by-... array, build sum and keep modulus over 256
tcc = mod(sum(reshape(tcc, 16, 0.0625 * numel(tcc)), 2), 256);

% then create hash string
tcc = sprintf('%02x', tcc');
