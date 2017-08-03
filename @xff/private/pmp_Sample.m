function pmpv = pmp_Sample(xo, coords, mnum)
% PMP::Sample  - sample PMP at given coordinates
%
% FORMAT:       pmpv = pmp.Sample(coords [, mnum]);
%
% Input fields:
%
%       coords      Cx3 cartesian coordinates to sample PMP at (centered!)
%       mnum        optional map number (default: 1)
%
% Output fields:
%
%       pmpv        Cx1 PMP values at coords
%
% Using: flexinterpn_methods, spherecoords.

% Version:  v1.1
% Build:    16020314
% Date:     Feb-03 2016, 2:35 PM EST
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

% neuroelf library
global ne_methods;

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'pmp') || ...
   ~isa(coords, 'double') || ndims(coords) > 2 || size(coords, 2) ~= 3
    error('neuroelf:xff:badArguments', 'Invalid call to %s.', mfilename);
end

% get file contents
bc = xo.C;
r1 = bc.ThetaResolution;
r2 = bc.PhiResolution;
if nargin < 3 || ~isa(mnum, 'double') || numel(mnum) ~= 1 || ...
    isinf(mnum) || isnan(mnum) || mnum < 1 || mnum > numel(bc.Map)
    mnum = 1;
else
    mnum = round(mnum);
end

% get spherical coordinates
scoords = ne_methods.spherecoords(coords);

% compute indices
crd1 = 1 + mod(r1 + (r1 / (2 * pi)) * scoords(:, 3), r1);
crd2 = 1 + mod(r2 + (r2 / pi) * scoords(:, 2), r2);

% build map for interpolation
imap = bc.Map(mnum).PMPData([1:r1, 1], [1:r2, r2]);

% interpolate
pmpv = ne_methods.flexinterpn_method(imap, [crd1, crd2], 0, 'linear');
