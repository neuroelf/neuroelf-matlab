function [xo, densm] = srf_Inflate(xo, niter, force, opts)
% SRF::Inflate  - apply inflation morph
%
% FORMAT:       [srf, densm] = srf.Inflate([niter [, force]])
%
% Input fields:
%
%       niter       number of iterations, default 3000
%       force       smoothing force, default 0.7
%
% Output fields:
%
%       srf         inflated surface
%       densm       density SMP
%
% Note: this method simply passes to SRF::Morph. It automatically set
%       area and distortion correction (to 6 * smoothing)

% Version:  v1.1
% Build:    16031617
% Date:     Mar-16 2016, 5:12 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'fsbf', 'srf'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isa(niter, 'double') || numel(niter) ~= 1 || ...
    isnan(niter) || niter < 1
    niter = 3000;
else
    niter = min(25000, floor(niter));
end
if nargin < 3 || ~isa(force, 'double') || numel(force) ~= 1 || ...
    isnan(force) || force <= 0 || force >= 1
    force = 0.7;
end
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
opts.areac = 1;
opts.distwsq = true;
opts.distc = max(1.5, min(3 - sqrt(eps), 6 * force));
opts.title = 'Inflating';

% depending on nargout
if nargout > 1
    [xo, densm] = srf_Morph(xo, niter, force, 'smooth', opts);
else
    srf_Morph(xo, niter, force, 'smooth', opts);
end
