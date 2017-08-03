function [varargout] = fsbf_Morph(xo, varargin)
% FSBF::Morph  - apply vertex morphing
%
% FORMAT:       [fsbf, densm] = fsbf.Morph([niter, force [, type [, opts]])
%
% Input fields:
%
%       niter       number of iterations (default 1)
%       force       smoothing force applied to morphing (default: 0.07)
%       type        either of 'even', {'smooth'}
%       opts        optional settings
%        .areac     if given and true, keep area constant
%        .areaw     area-weighted smoothing (higher precedence, false)
%        .distc     additional distortion correction force (>0 ... 3)
%        .distw     distance-weighted smoothing (lower precedence, false)
%        .distwsq   square-of-distance weighting (also sets distw, false)
%        .norm      force along normal vector (default: 0)
%        .normramp  ramp-up normal force from 0 to .norm (default: false)
%        .pbar      optional progress bar object (xfigure/xprogress, [])
%        .show      show during morphing (default: true)
%        .sphere    additionally applied to-sphere force
%        .title     display title for progress bar
%
% Output fields:
%
%       fsbf        altered object
%       densm       density map SMP
%
% Using: mesh_morph.

% Version:  v1.1
% Build:    16031616
% Date:     Mar-16 2016, 4:43 PM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'fsbf')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end

% pass on
try
    [varargout{1:max(1, nargout)}] = srf_Morph(xo, varargin{:});
catch xfferror
    rethrow(xfferror);
end
