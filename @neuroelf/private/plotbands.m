function h = plotbands(y, opts)
%PLOTBANDS  Plots the data in separate bands.
%   PLOTBANDS(Y) plots the data in the columns in Y into separate bands.
%
%   PLOTBANDS(Y, OPTS) allows to set options in a 1x1 struct:
%
%     .spread   distance between bands (default: 2 times the median STD)
%     .trans    transform the column of Y, one of {'none'}, 'psc', 'z'
%     .x        give the X ordinate values for Y (must be size(Y, 1)-by-1)

% Version:  v1.1
% Build:    16060616
% Date:     Jun-06 2016, 4:00 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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
if nargin < 1 || ~isnumeric(y) || size(y, 1) < 2 || ndims(y) > 2
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
y = double(y);
szy = size(y);
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'spread') || ~isa(opts.spread, 'double') || numel(opts.spread) ~= 1 || ...
    isinf(opts.spread) || isnan(opts.spread) || opts.spread <= 0
    opts.spread = [];
end
if ~isfield(opts, 'trans') || ~ischar(opts.trans) || isempty(opts.trans) || ~any(lower(opts.trans(1)) == 'npz')
    opts.trans = 'n';
end
if ~isfield(opts, 'x') || ~isa(opts.x, 'double') || numel(opts.x) ~= szy(1)
    opts.x = (0:szy(1)-1)';
end

% transform
if opts.trans == 'p'
    y = psctrans(y) - 100;
elseif opts.trans == 'z'
    y = ztrans(y);
else
    y = y - repmat(mean(y), szy(1), 1);
end

% spread needed
if isempty(opts.spread)
    opts.spread = 2 * median(std(y));
end

% plot
h = plot(opts.x, y + ones(szy(1), 1) * (0:opts.spread:(szy(2)-0.9)*opts.spread));
