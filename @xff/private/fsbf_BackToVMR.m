function xo = fsbf_BackToVMR(xo, opts)
% FSBF::BackToVMR  - back-project vertices to VMR space
%
% FORMAT:       [fsbf = ] fsbf.BackToVMR([opts])
%
% Input fields:
%
%       opts        optional settings
%        .fillmode  either of {'nearest'} or 'linear'
%        .fillvmax  maximal filling value (combination, default: 200)
%        .nfrom     Vertex + from * Normal (default: -0.25)
%        .nstep     step along normal(default: 0.5)
%        .nto       Vertex + to * Normal (default: nfrom + nstep)
%        .res       VMR resolution, either of {1} or 0.5
%        .smp       SMP map object (if different marking is requested)
%        .smpmap    sub-map of SMP (default: 1)
%        .tcode     target color code (default: 235)
%        .triovsmp  triangular oversampling factor (default: 2)
%        .vmr       insert data into VMR (after sampling)
%
% Output fields:
%
%       fsbf        FSBF object with VertexVMRData set accordingly

% Version:  v1.1
% Build:    16031614
% Date:     Mar-16 2016, 2:05 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'fsbf')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end

% pass on
try
    srf_BackToVMR(xo, opts);
catch xfferror
    rethrow(xfferror);
end
