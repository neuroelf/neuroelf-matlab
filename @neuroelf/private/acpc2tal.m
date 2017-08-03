function acpcc = acpc2tal(acpcc, tal, ital, hires)
% acpc2tal  - convert AC-PC coordinates into TAL coordinates
%
% FORMAT:       talc = acpc2tal(acpcc, tal [, invtal [, hires]])
%
% Input fields:
%
%       acpcc       Cx3 AC-PC coordinates
%       tal         either 8x3 TAL coordinates or TAL xff object
%       invtal      boolean flag, perform inverse operation (default: false)
%       hires       compute transformation for 0.5mm (default: false)
%
% Output fields:
%
%       talc        Talairach coordinates (or AC-PC for inverse)
%
% Note: this function works with BV system (TAL axes order but reverse
%       orientation, so that 128 - BVSys = TAL) coordinates. If any of
%       the coordinates is below zero, the function assumes TAL instead
%       and equally returns TAL coordinates

% Version:  v1.1
% Build:    16051016
% Date:     May-10 2016, 4:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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

% persistent storage
persistent acpc2tal_istdtal;
if isempty(acpc2tal_istdtal)
    talobj = xff('new:tal');
    acpc2tal_istdtal = getcont(talobj);
    delete(talobj);
end

% argument check
if nargin < 2 || ~isa(acpcc, 'double') || ndims(acpcc) ~= 2 || size(acpcc, 2) ~= 3 || ...
    any(isinf(acpcc(:)) | isnan(acpcc(:)) | acpcc(:) < -128 | acpcc(:) > 768) || ...
   (~isxff(tal, 'tal') && (~isa(tal, 'double') || ndims(tal) ~= 2 || size(tal, 2) ~= 3 || ...
    any(isinf(tal(:)) | isnan(tal(:)) | tal(:) < 0 | tal(:) > 256)))
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
if nargin < 3 || ~islogical(ital) || numel(ital) ~= 1
    ital = false;
end
if nargin < 4 || ~islogical(hires) || numel(hires) ~= 1
    hires = false;
end

% if xff TAL object then extract points
if isxff(tal)
    tal = [tal.AC; tal.PC; tal.AP; tal.PP; tal.SP; tal.IP; tal.RP; tal.LP];
end

% standard TAL coordinates
stt = acpc2tal_istdtal;
stt = [stt.AC; stt.PC; stt.AP; stt.PP; stt.SP; stt.IP; stt.RP; stt.LP];

% already standard ?
if all(tal(:) == stt(:))
    return;
end

% hires
if hires
    if ital
        tal = 0.5 .* tal;
    else
        stt = 2 .* stt;
    end
end

% for backwards transform, reverse tal system and std system
if ital
    [tal, stt] = deal(stt, tal);
end

% get 12 subvolume coords
txc = [tal(7, 3), tal(1, 3), tal(8, 3)];
tyc = [tal(3, 1), tal(1, 1), tal(2, 1), tal(4, 1)];
tzc = [tal(5, 2), tal(1, 2), tal(6, 2)];
txd = diff(txc);
tyd = diff(tyc);
tzd = diff(tzc);
sxc = [stt(7, 3), stt(1, 3), stt(8, 3)];
syc = [stt(3, 1), stt(1, 1), stt(2, 1), stt(4, 1)];
szc = [stt(5, 2), stt(1, 2), stt(6, 2)];
sxd = diff(sxc);
syd = diff(syc);
szd = diff(szc);
xd = sxc - txc;
yd = syc - tyc;
zd = szc - tzc;

% transform X
if xd(1) ~= 0
    b1 = (acpcc(:, 1) <  txc(1));
end
if any(xd(1:2) ~= 0)
    b2 = (acpcc(:, 1) >= txc(1) & acpcc(:, 1) <= txc(2));
end
if any(xd(2:3) ~= 0)
    b3 = (acpcc(:, 1) >  txc(2) & acpcc(:, 1) <= txc(3));
end
if xd(3) ~= 0
    b4 = (acpcc(:, 1) >  txc(3));
end
if xd(1) ~= 0
    acpcc(b1, 1) = xd(1) + acpcc(b1, 1);
end
if any(xd(1:2) ~= 0)
    acpcc(b2, 1) = sxc(1) + (sxd(1) / txd(1)) * (acpcc(b2, 1) - txc(1));
end
if any(xd(2:3) ~= 0)
    acpcc(b3, 1) = sxc(2) + (sxd(2) / txd(2)) * (acpcc(b3, 1) - txc(2));
end
if xd(3) ~= 0
    acpcc(b4, 1) = xd(3) + acpcc(b4, 1);
end

% transform Y
if yd(1) ~= 0
    b1 = (acpcc(:, 2) <  tyc(1));
end
if any(yd(1:2) ~= 0)
    b2 = (acpcc(:, 2) >= tyc(1) & acpcc(:, 2) <  tyc(2));
end
if any(yd(2:3) ~= 0)
    b3 = (acpcc(:, 2) >=  tyc(2) & acpcc(:, 2) <= tyc(3));
end
if any(yd(3:4) ~= 0)
    b4 = (acpcc(:, 2) >  tyc(3) & acpcc(:, 2) <= tyc(4));
end
if yd(4) ~= 0
    b5 = (acpcc(:, 2) >  tyc(4));
end
if yd(1) ~= 0
    acpcc(b1, 2) = yd(1) + acpcc(b1, 2);
end
if any(yd(1:2) ~= 0)
    acpcc(b2, 2) = syc(1) + (syd(1) / tyd(1)) * (acpcc(b2, 2) - tyc(1));
end
if any(yd(2:3) ~= 0)
    acpcc(b3, 2) = syc(2) + (syd(2) / tyd(2)) * (acpcc(b3, 2) - tyc(2));
end
if any(yd(3:4) ~= 0)
    acpcc(b4, 2) = syc(3) + (syd(3) / tyd(3)) * (acpcc(b4, 2) - tyc(3));
end
if yd(4) ~= 0
    acpcc(b5, 2) = yd(4) + acpcc(b5, 2);
end

% transform Z
if zd(1) ~= 0
    b1 = (acpcc(:, 3) <  tzc(1));
end
if any(zd(1:2) ~= 0)
    b2 = (acpcc(:, 3) >= tzc(1) & acpcc(:, 3) <= tzc(2));
end
if any(zd(2:3) ~= 0)
    b3 = (acpcc(:, 3) >  tzc(2) & acpcc(:, 3) <= tzc(3));
end
if zd(3) ~= 0
    b4 = (acpcc(:, 3) >  tzc(3));
end
if zd(1) ~= 0
    acpcc(b1, 3) = zd(1) + acpcc(b1, 3);
end
if any(zd(1:2) ~= 0)
    acpcc(b2, 3) = szc(1) + (szd(1) / tzd(1)) * (acpcc(b2, 3) - tzc(1));
end
if any(zd(2:3) ~= 0)
    acpcc(b3, 3) = szc(2) + (szd(2) / tzd(2)) * (acpcc(b3, 3) - tzc(2));
end
if zd(3) ~= 0
    acpcc(b4, 3) = zd(3) + acpcc(b4, 3);
end
