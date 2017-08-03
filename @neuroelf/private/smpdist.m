function s = smpdist(d, sz, opts)
% smpdist  - sample values from a distribution
%
% FORMAT:       s = smpdist(d, sz [, opts])
%
% Input fields:
%
%       d           1D vector defining the distribution (e.g. a residual)
%       sz          required size of the sampled values
%       opts        optional settings
%        .smooth    smoothing of the distribution, one of
%                   {'none'}  - no smoothing (simply draw with replacement)
%                   'norm'    - add a normally distributed error
%                   'uniform' - add a uniform distributed error
%        .smoothk   smoothing kernel, default: sqrt(std(d) / numel(d))
%
% Output fields:
%
%       s           sampled values

% Version:  v0.9b
% Build:    10082916
% Date:     Aug-24 2010, 12:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the
%       distribution.
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
if nargin < 2 || ...
   ~isa(d, 'double') || ...
    numel(d) ~= max(size(d)) || ...
    any(isinf(d) | isnan(d)) || ...
   ~isnumeric(sz) || ...
    numel(sz) ~= size(sz, 2) || ...
    any(isinf(sz) | isnan(sz) | sz < 0 | sz ~= fix(sz)) || ...
    prod(sz) > 1e8
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'smooth') || ...
   ~ischar(opts.smooth) || ...
   ~any(strcmpi(opts.smooth(:)', {'none', 'norm', 'uniform'}))
    opts.smooth = 'none';
else
    opts.smooth = lower(opts.smooth(:)');
end

% no smoothing
if strcmp(opts.smooth, 'none')

    % create sampling
    s = d(ceil(numel(d) .* rand(sz)));

    % return
    return;
end

% smoothing kernel
if ~isfield(opts, 'smoothk') || ...
   ~isa(opts.smoothk, 'double') || ...
    numel(opts.smoothk) ~= 1 || ...
    isinf(opts.smoothk) || ...
    isnan(opts.smoothk) || ...
    opts.smoothk < 0
    opts.smoothk = sqrt(std(d) / numel(d));
end

% adding noise from the normal distribution
if strcmp(opts.smooth, 'norm')

    % sample
    s = d(ceil(numel(d) .* rand(sz))) + opts.smoothk .* randn(sz);

% adding noise from the uniform distribution
else

    % sample
    s = d(ceil(numel(d) .* rand(sz))) + ...
        (sqrt(12) * opts.smoothk) .* (-0.5 + rand(sz));
end
