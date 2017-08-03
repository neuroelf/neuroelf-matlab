function [varargout] = fsbf_Combine(xo, xo2, cbopt)
% FSBF::Combine  - combines two FSBFs into one
%
% FORMAT:       [cfsbf [, ...]] = fsbf1.Combine(fsbf2, cbopt);
%
% Input fields:
%
%       fsbf2       surface to use to combine given FSBF with
%       cbopt       struct with optional settings
%       .type       one of
%                   'backtoback' - rotate one mesh 180 in XY plane
%                   'custom'     - build custom scene, see below
%                   'gapped'     - join contents with a gap
%                   'wholebrain' - simply join contents (default)
%       .color1     1x4 double for .Color1 field in SRF options
%       .color2     1x4 double for .Color2 field in SRF options
%       .filename   store combined under new filename
%       .gap        1x1 double, mm to insert between two meshes
%                   (applied for all types accordingly, defaults:
%                   backtoback: 25, gapped: 100, outandin: 20,
%                   outintb: 25, patched: 25, spm2: 25, wholebrain: 0)
%       .linkedsrf  filename of linked SRF, set to empty if not given
%       .mtc1       xff MTC object for the first FSBF
%       .mtc2       xff MTC object for the second FSBF
%       .smp1       xff SMP object for the first FSBF
%       .smp2       xff SMP object for the second FSBF
%       .ssm1       xff SSM object for the first FSBF
%       .ssm2       xff SSM object for the second FSBF
%       .transform  1x2 cell array with 4x4 double transformation matrices,
%                   needed for custom scenaries
%
% Output fields:
%
%       cfsbf       combined FSBF
%       ...         combined MTC/SMP/SSM (if given, in that order!)
%
% Using: tfmatrix.

% Version:  v1.1
% Build:    16031616
% Date:     Mar-16 2016, 4:52 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || numel(xo2) ~= 1 || ...
   ~xffisobject(xo, true, 'fsbf') || ~xffisobject(xo2, true, 'fsbf')
    error('neuroelf:xff:badArguments', 'Invalid call to %s.', mfilename);
end
if nargin < 3 || ~isstruct(cbopt) || numel(cbopt) ~= 1
    cbopt = struct;
end

% pass on
try
    [varargout{1:max(1, nargout)}] = srf_Combine(xo, xo2, cbopt);
catch xfferror
    rethrow(xfferror);
end
