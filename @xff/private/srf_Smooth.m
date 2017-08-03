function [xo, densm] = srf_Smooth(xo, niter, force, opts)
% SRF::Smooth  - apply smoothing morph
%
% FORMAT:       [srf, densm] = srf.Smooth([niter [, force, opts]])
%
% Input fields:
%
%       niter       number of iterations, default 150
%       force       smoothing force, default 0.07
%       opts        optional settings
%        .pbar      progress bar flag (default: from config file)
%        .show      show during smoothing (default: true, SRF only)
%
% Output fields:
%
%       srf         smoothed surface
%       densm       density SMP
%
% Note: this method simply passes to SRF::Morph.

% Version:  v1.1
% Build:    16031616
% Date:     Mar-16 2016, 4:45 PM EST
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

% global variable for config
global xffsngl;
stc = xffsngl.CONF.settings;
if ~isfield(stc, 'Morphing')
    stc = struct('Morphing', struct);
end
stc = stc.Morphing;
if isfield(stc, 'Smooth')
    stc = stc.Smooth;
else
    stc = struct('Force', 0.07, 'NrOfIterations', 150);
end

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'fsbf', 'srf'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isa(niter, 'double') || numel(niter) ~= 1 || ...
    isnan(niter) || niter < 0 || niter > 1e5
    niter = stc.NrOfIterations;
else
    niter = floor(niter);
end
if nargin < 3 || ~isa(force, 'double') || numel(force) ~= 1 || ...
    isnan(force) || force < -16 || force > 16 || force == 0
    force = stc.Force;
end
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'pbar')
    opts.pbar = [];
end
if ~isfield(opts, 'distw') || numel(opts.distw) ~= 1 || ...
   (~isa(opts.distw, 'double') && ~islogical(opts.distw))
    opts.distw = 0;
else
    opts.distw = logical(double(opts.distw) ~= 0);
end
if ~isfield(opts, 'distwl') || numel(opts.distwl) ~= 1 || ...
   (~isa(opts.distwl, 'double') && ~islogical(opts.distwl))
    opts.distwl = 0;
else
    opts.distwl = logical(double(opts.distwl) ~= 0);
end
if ~isfield(opts, 'distwsq') || numel(opts.distwsq) ~= 1 || ...
   (~isa(opts.distwsq, 'double') && ~islogical(opts.distwsq))
    opts.distwsq = 0;
else
    opts.distwsq = logical(double(opts.distwsq) ~= 0);
end
if ~isfield(opts, 'show') || ~islogical(opts.show) || numel(opts.show) ~= 1
    opts.show = true;
end
if ~isfield(opts, 'stepsize') || ~isa(opts.stepsize, 'double') || numel(opts.stepsize) ~= 1 || ...
    isinf(opts.stepsize) || isnan(opts.stepsize) || opts.stepsize < 1
    opts.stepsize = min(250, max(50, ceil(niter ^ (2 / 3))));
else
    opts.stepsize = min(niter, floor(opts.stepsize));
end

% depending on nargout
if nargout > 1
    [xo, densm] = srf_Morph(xo, niter, force, 'smooth', struct('pbar', opts.pbar, ...
        'distw', opts.distw, 'distwl', opts.distwl, 'distwsq', opts.distwsq, ...
        'show', opts.show, 'stepsize', opts.stepsize));
else
    srf_Morph(xo, niter, force, 'smooth', struct('pbar', opts.pbar, ...
        'distw', opts.distw, 'distwl', opts.distwl, 'distwsq', opts.distwsq, ...
        'show', opts.show, 'stepsize', opts.stepsize));
end
