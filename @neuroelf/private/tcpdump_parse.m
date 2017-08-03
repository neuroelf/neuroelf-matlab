function packet = tcpdump_parse(filename, out)
% tcpdump_parse  - parse tcpdump outfile (-w) file
%
% FORMAT:       packet = tcpdump_parse(filename [, out])
%
% Input fields:
%
%       filename    name of dumped file
%       out         if given, dump output
%
% Output fields:
%
%       packet      1xN struct with packet fields
%
% Note: only a few protocols and types are currently supported

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-08 2011, 9:16 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% try reading file
try
    filec = double(binread(filename));
    if length(filec) < 64
        error('Too few bytes in file.');
    end
catch ne_eo;
    error( ...
        'neuroelf:ErrorReadingFile', ...
        'Error reading file (%s).', ...
        ne_eo.message ...
    );
end

% some constant
kd = datenum('01-Jan-1970');

% filecheck
magic = sprintf('%8X', dword(filec(1:4)));
if ~strcmp(magic, 'A1B2C3D4')
    error( ...
        'neuroelf:BadFileContent', ...
        'No/bad tcpdump file specified.' ...
    );
end
snaplen = dword(filec(17:20));

% get real content
mpack = struct;
mpack.Type           = 'MASTER';
mpack.Magic          = magic;
mpack.MagicMajorVer  = word(filec(5:6));
mpack.MagicMinorVer  = word(filec(7:8));
mpack.TimeZoneOffset = dword(filec(9:12));
mpack.TimeStampAccur = dword(filec(13:16));
mpack.SnapshotLength = snaplen;
mpack.LinkLayerType  = dword(filec(21:24));
packet = struct( ...
    'Date',     '', ...
    'FrameLen', 0, ...
    'SavedLen', 0, ...
    'Version',  [], ...
    'IHL',      [], ...
    'TOS',      [], ...
    'Length',   [], ...
    'Ident',    [], ...
    'FlagRSV',  [], ...
    'FlagDFR',  [], ...
    'FlagMFR',  [], ...
    'FragOff',  [], ...
    'TTL',      [], ...
    'Protocol', [], ...
    'CheckSum', 0, ...
    'Source',   '', ...
    'Target',   '', ...
    'PROTO',    mpack, ...
    'PDATA',    uint8([]) ...
);
l = length(filec);
c = 2;
p = 25;
while ((p + 30) < l)
    pack = struct;
    pack.Date = datestr(kd + dword(filec(p:p + 3)) / 86400);
    pack.FrameLen = dword(filec(p + 8:p + 11));
    pack.SavedLen = dword(filec(p + 12:p + 15));
    p = p + 30;
    pack.Version  = floor(filec(p) ./ 16);
    pack.IHL      = mod(filec(p), 16);
    if pack.IHL < 5
        error( ...
            'neuroelf:BadFileContent', ...
            'Invalid IP header length in packet.' ...
        );
    end
    pack.TOS      = filec(p + 1);
    pack.Length   = nword(filec(p + 2:p + 3));
    pack.Ident    = nword(filec(p + 4:p + 5));
    pack.FlagRSV  = filec(p + 6) > 127;
    pack.FlagDFR  = mod(filec(p + 6), 128) > 63;
    pack.FlagMFR  = mod(filec(p + 6), 64) > 31;
    pack.FragOff  = mod(nword(filec(p + 6:p + 7)), 32 * 256);
    pack.TTL      = filec(p + 8);
    pack.Protocol = filec(p + 9);
    pack.CheckSum = sprintf('%4X', nword(filec(p + 10:p + 11)));
    pack.Source   = sprintf('%d.%d.%d.%d', filec(p + 12:p + 15));
    pack.Target   = sprintf('%d.%d.%d.%d', filec(p + 16:p + 19));

    sp = p + 20;
    % skip remaining IP header
    if pack.IHL > 5
        sp = sp + 4 * (pack.IHL - 5);
    end

    % what protocol
    phead = struct;
    pdata = [];
    switch (pack.Protocol)
        case {  6} % TCP
            phead.Type     = 'TCP';
            phead.SrcPort  = nword(filec(sp:sp + 1));
            phead.DestPort = nword(filec(sp + 2:sp + 3));
        case { 17} % UDP
            phead.Type = 'UDP';
            phead.SrcPort  = nword(filec(sp:sp + 1));
            phead.DestPort = nword(filec(sp + 2:sp + 3));
            phead.Length   = nword(filec(sp + 4:sp + 5));
            phead.CheckSum = nword(filec(sp + 6:sp + 7));
            pdata = uint8(filec(sp + 8:sp + phead.Length - 1));
        otherwise
            phead.Type     = '???';
            phead.SrcPort  = 0;
            phead.DestPort = 0;
    end

    % store data
    pack.PROTO = phead;
    pack.PDATA = pdata;

    packet(c) = pack;
    c = c + 1;
    p = p + pack.Length;
end

if nargin < 2 || ...
    isempty(out) || ...
    ~out(1)
    return;
end

for c = 2:length(packet)
    disp(sprintf('%s  %s.%d -> %s.%d: %s, length %d', ...
        packet(c).Date, ...
        packet(c).Source, ...
        packet(c).PROTO.SrcPort, ...
        packet(c).Target, ...
        packet(c).PROTO.DestPort, ...
        packet(c).PROTO.Type, ...
        length(packet(c).PDATA)));
    disp(hexdump(packet(c).PDATA));
end



function hd = hexdump(av)
    l = length(av);
    p = 1;
    hd = '';
    while (p < l)
        pa = av(p:min(p+15, l));
        hd = [hd sprintf('    %08x:  %-50s %-16s\n', p - 1, ...
            sprintf('%02x ', pa), ...
            saveascii(pa))];
        p = p + 16;
    end
% end of function hd = hexdump(av)

function sa = saveascii(av)
    av(av < 32 | av > 127) = 46;
    sa = sprintf('%c', av);
% end of function sa = saveascii(av)

%function dwval = ndword(av)
%    dwval = av(4) + av(3) * 256 + av(2) * 65536 + av(1) * 16777216;
% end of function dwval = ndword(av)

function dwval = nword(av)
    dwval = av(2) + av(1) * 256;
% end of function wval = nword(av)

function dwval = dword(av)
    dwval = av(1) + av(2) * 256 + av(3) * 65536 + av(4) * 16777216;
% end of function dwval = dword(av)

function wval = word(av)
    wval = av(1) + av(2) * 256;
% end of function wval = word(av)
