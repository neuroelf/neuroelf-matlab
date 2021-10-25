function decoded = t2k_DecodeByteStream(xo)
% T2K::DecodeByteStream  - decode the byte stream in a T2K object
%
% FORMAT:       decoded = t2k.DecodeByteStream();
%
% No onput fields.
%
% Output fields:
%
%       decoded     decoded byte stream

% Version:  v1.1
% Build:    21102517
% Date:     Oct-25 2021, 5:29 PM EST
% Author:   Jochen Weber, NeuroElf.net, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2021, Jochen Weber
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
if nargin < 1 || numel(xo) ~= 1 || ...
   ~xffisobject(xo, true, 't2k')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
encoded = bc.ByteStream;
le = numel(encoded);
if mod(le, 4) ~= 0
    encoded = [encoded, uint8(zeros(1, 4 - mod(le, 4)))];
end
encoded = typecast(encoded, 'uint32');
decoded = encoded;
nv = numel(encoded);

% initialize key
key = uint32(2779);
ba = uint32(31);
for vc = 1:nv
    dv = bitxor(encoded(vc), key);
    ks = double(bitand(dv, ba));
    if ks ~= 0
        key = bitor(bitshift(key, ks), bitshift(key, ks-32));
    end
    decoded(vc) = dv;
end
decoded = typecast(decoded, 'uint8');
if numel(decoded) > le
    decoded = decoded(1:le);
end
