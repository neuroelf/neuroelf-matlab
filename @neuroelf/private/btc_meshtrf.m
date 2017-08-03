% FUNCTION tm = btc_meshtrf: compute hgtransform/surface transform matrix
function tm = btc_meshtrf(cc)

% Version:  v1.1
% Build:    16031900
% Date:     Mar-19 2016, 12:07 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% for hgtransform, expect struct
if isstruct(cc)
    ang = (pi / 180) .* [cc.anglex, cc.angley];
    cang = cos(ang);
    sang = sin(ang);
    zoom = cc.zoom;

    % compute global transformation matrices
    atrf = [cang(2), 0, sang(2), 0; 0, 1, 0, 0; -sang(2), 0, cang(2), 0; 0, 0, 0, 1] * ...
           [cang(1), sang(1), 0, 0; -sang(1), cang(1), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
    tm = [1, 0, 0, cc.trans(1); 0, 1, 0, cc.trans(2); 0, 0, 1, cc.trans(3); 0, 0, 0, 1] * ...
          atrf * [zoom, 0, 0, 0; 0, zoom, 0, 0; 0, 0, zoom, 0; 0, 0, 0, 1];

% for surface, expect cell
else
    iptrf = cc{1};
    cscc = cos(cc{2});
    sscc = sin(cc{2});
    ipzoom = cc{3};
    if numel(ipzoom) == 1
        ipzoom = ipzoom([1, 1, 1]);
    end
    iatrf = [1, 0, 0, 0; 0, cscc(1), sscc(1), 0; 0, -sscc(1), cscc(1), 0; 0, 0, 0, 1] * ...
            [cscc(2), 0, sscc(2), 0; 0, 1, 0, 0; -sscc(2), 0, cscc(2), 0; 0, 0, 0, 1] * ...
            [cscc(3), sscc(3), 0, 0; -sscc(3), cscc(3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
    tm = [1, 0, 0, iptrf(1); 0, 1, 0, iptrf(2); 0, 0, 1, iptrf(3); 0, 0, 0, 1] * ...
          iatrf * [ipzoom(1), 0, 0, 0; 0, ipzoom(2), 0, 0; 0, 0, ipzoom(3), 0; 0, 0, 0, 1];
end
