function [varargout] = fsbf_Smooth(xo, varargin)
% FSBF::Smooth  - apply smoothing morph
%
% FORMAT:       [fsbf, densm] = fsbf.Smooth([niter [, force, opts]])
%
% Input fields:
%
%       niter       number of iterations, default 150
%       force       smoothing force, default 0.07
%       opts        optional settings
%        .pbar      progress bar flag (default: from config file)
%        .show      show during smoothing (default: true)
%
% Output fields:
%
%       fsbf        smoothed surface
%       densm       density SMP
%
% Note: this method simply passes to FSBF::Morph.

% Version:  v1.1
% Build:    16031616
% Date:     Mar-16 2016, 4:40 PM EST
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
    [varargout{1:max(1, nargout)}] = srf_Smooth(xo, varargin{:});
catch xfferror
    rethrow(xfferror);
end
