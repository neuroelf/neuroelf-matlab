function nxo = mtc_ApplySSM(xo, ssm)
% MTC::ApplySSM  - resample MTC with SSM mapping
%
% FORMAT:       newmtc = mtc.ApplySSM(ssm)
%
% Input fields:
%
%       ssm         Sphere-to-Sphere-Mapping (SSM) object
%
% Output fields:
%
%       newmtc      MTC with resampled timecourses

% Version:  v1.1
% Build:    16020917
% Date:     Feb-09 2016, 5:06 PM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'mtc') || ...
    numel(ssm) ~= 1 || ~xffisobject(ssm, true, 'ssm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
ssmc = ssm.C;
if size(bc.MTCData, 2) ~= ssmc.NrOfSourceVertices
    error('neuroelf:xff:badArgument', ...
        'NrOfVertices must match with NrOfSourceVertices in SSM object.');
end

% copy object
nxo = aft_CopyObject(xo);
nbc = nxo.C;

% resample data and update fields
nbc.MTCData = bc.MTCData(:, ssmc.SourceOfTarget);
nbc.NrOfVertices = size(nbc.MTCData, 2);
nxo.C = nbc;
