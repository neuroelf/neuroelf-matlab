function xo = srf_SpreadColor(xo, v)
% SRF::SpreadColor  - spread colored rings from one vertex
%
% FORMAT:       [srf] = srf.SpreadColor([vertex])
%
% Input fields:
%
%       vertex      optional vertex number (by default 1)
%
% Output fields:
%
%       srf         colored surface

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:42 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
cc = bc.NrOfVertices;
if nargin < 2 || ~isa(v, 'double') || numel(v) ~= 1 || isnan(v) || v < 1 || v > cc || v ~= fix(v)
    v = 1;
end

% create color
col = repmat([nan, 255, 0, 0], [cc, 1]);

% get neighbors
nei = bc.Neighbors;
neir = nei(:, 2);

% prepare lists of done and to-do vertices
neidone = uint16(zeros(1, cc));
neicons = false(1, cc);
ncolorb = {uint32([1, zeros(1, 12)])};
ncolorb(2:cc) = ncolorb;

% mark first vertex
neicons(v) = true;

% step counters
st = 0;

% continue while still vertices found
while any(neicons)

    % mark colors
	col(neicons, 3) = mod(10 * st, 256);
	col(neicons, 4) = mod(10 * floor(st / 10), 256);
	st = st + 1;

    % set in done
    neici = find(neicons);
    neidone(neici) = st;
    consnei = neir(neici);
    neicons(:) = false;
    for c = 1:numel(neici)
        tnei = consnei{c};
        fnei = tnei(neidone(tnei) == 0);
        neicons(tnei) = true;
        for tc = 1:numel(fnei)
            tcc = fnei(tc);
            ncolorb{tcc}(1) = ncolorb{tcc}(1) + 1;
            ncolorb{tcc}(ncolorb{tcc}(1)) = neici(c);
        end
    end
    neicons(neidone > 0) = false;
end

% check by whom vertices were colored
for c = 1:cc
    if ncolorb{c}(1) > 5
        neicons(c) = true;
    end
end
col(neicons, 2) = 64;

% back into storage
bc.VertexColor = col;
xo.C = bc;
