function xo = mat_Resample(xo, sfreq, nfreq, opts)
% MAT::Resample  - resample channel data to new frequency
%
% FORMAT:       [mat = ] mat.Resample(sfreq, freq [, opts])
%
% Input fields:
%
%       sfreq       original frequency
%       nfreq       new frequency
%       opts        optional settings
%        .cubic     either 1x1 or 1xC boolean flag, use cubic/gaussian
%                   interpolation (default, 1x1 true; for stim channels
%                   you would want to turn this off)
%
%
% Output fields:
%
%       mat         altered object
%
% Using: resampleaa.

% Version:  v1.1
% Build:    16020515
% Date:     Feb-05 2016, 3:55 PM EST
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

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'mat') || ...
   ~isa(sfreq, 'double') || numel(sfreq) ~= 1 || isinf(sfreq) || isnan(sfreq) || ...
    sfreq <= 0 || sfreq > 10000 || ...
   ~isa(nfreq, 'double') || numel(nfreq) ~= 1 || isinf(nfreq) || isnan(nfreq) || ...
    nfreq <= 0 || nfreq > 10000
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'cubic') || ~islogical(opts.cubic) || ...
   ~any([1, size(bc.Data, 2)] == numel(opts.cubic))
    opts.cubic = true(1, size(bc.Data, 2));
elseif numel(opts.cubic) == 1
    opts.cubic = opts.cubic(1, ones(1, size(bc.Data, 2)));
end

% set kernel to []
k = [];

% compute factor
fac = sfreq / nfreq;

% do nothing?
if fac == 1
    return;
end
if fac < 1
    ub = size(bc.Data, 1) + 1 - fac;
elseif fac > 1
    ub = size(bc.Data, 1) + (fac - 1) / fac;
end
newdata = zeros(floor((ub - 1) / fac) + 1, size(bc.Data, 2));

% iterate over channels
for cc = 1:size(newdata, 2)

    % with cubic/gaussian interpolation
    if opts.cubic(cc)
        [newdata(:, cc), k] = ...
            ne_methods.resampleaa(mat_ChannelData(xo, cc), fac, 1, 0, k);

    % without
    else
        olddata = mat_ChannelData(xo, cc);
        newdata(:, cc) = olddata(round(1:fac:ub));
    end
end

% set back in memory
bc.Data = newdata;
xo.C = bc;
