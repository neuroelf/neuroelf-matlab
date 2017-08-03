function nsph = fsbf_ApplySSM(xo, ssm, sphere, smiter)
% FSBF::ApplySSM  - apply SSM file on vertex coordinates
%
% FORMAT:       sph = fsbf.ApplySSM(ssm, sphere [, smiter])
%
% Input fields:
%
%       fsbf        folded mesh (e.g. RECOSM mesh)
%       ssm         SSM mapping to sphere mesh
%       sphere      sphere mesh SRF used for mapping
%       smiter      number of smoothing steps (default = 20, force 0.03)
%
% Output fields:
%
%       sph         SPH mesh with NrOfVertices of sphere, folded as fsbf
%
% Using: mesh_morph.

% Version:  v1.1
% Build:    16031813
% Date:     Mar-18 2016, 1:11 PM EST
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
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'fsbf') || ...
    numel(ssm) ~= 1 || ~xffisobject(ssm, true, 'ssm') || ...
    numel(sphere) ~= 1 || ~xffisobject(sphere, true, {'fsbf', 'srf'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 4 || numel(smiter) ~= 1 || ~isa(smiter, 'double') || ...
    isinf(smiter) || isnan(smiter) || smiter < 0 || smiter > 3000
    smiter = 20;
else
    smiter = round(smiter);
end

% pass on
try
    nsph = srf_ApplySSM(xo, tsm, sphere, smiter);
catch xfferror
    rethrow(xfferror);
end
