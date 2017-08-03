function ro = ranktrans(rd, dim, opts)
% ranktrans  - return rank transform of data
%
% FORMAT:       ranks = ranktrans(data [, dim [, opts]])
%
% Input fields:
%
%       data        numeric data
%       dim         dimension along to transform (default: last)
%       ppts        optional settings
%       .meancenter flag, center around mean
%       .nozero     flag, do not "rank" zero samples
%
% Output fields:
%
%       ranks       1:size(data, dim) rank-transformed data
%
% Note: ties will get the mean rank.

% Version:  v0.9c
% Build:    12121415
% Date:     Dec-14 2012, 3:14 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2012, Jochen Weber
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
if nargin < 1 || ...
   ~isnumeric(rd) || ...
    numel(rd) < 1
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing first argument.' ...
    );
end
if ~isa(rd, 'double')
    rd = double(rd);
end
if nargin < 2 || ...
   ~isa(dim, 'double') || ...
    numel(dim) ~= 1 || ...
    isinf(dim) || ...
    isnan(dim) || ...
    dim < 1 || ...
    dim > ndims(rd)
    dim = ndims(rd);
    if size(rd, dim) == 1
        dim = dim - 1;
    end
end
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'meancenter') || ...
   ~islogical(opts.meancenter) || ...
    numel(opts.meancenter) ~= 1
    opts.meancenter = false;
end
if ~isfield(opts, 'nozero') || ...
   ~islogical(opts.nozero) || ...
    numel(opts.nozero) ~= 1
    opts.nozero = false;
end

% get elements that are Inf, NaN
bv = isinf(rd) | isnan(rd);
if opts.nozero
    bv = bv | (rd == 0);
end
mmm = minmaxmean(rd, 4);
rd(bv) = mmm(2) + 1;

% first pass, sort forwards
[rs, ri] = sort(rd, dim);
[ri, ro] = sort(ri, dim);

% to resolve ties (to mean rank) sort reversed array
sr = {':'};
sr = sr(ones(1, ndims(rd)));
sr = struct('type', '()', 'subs', {sr});
sr.subs{dim} = size(rd, dim):-1:1;
[rs, ri] = sort(subsref(rd, sr), dim);
[ri, rs] = sort(ri, dim);

% then build average
ro = 0.5 .* ro + 0.5 .* subsref(rs, sr);

% mean center?
if opts.meancenter
    rma = ones(1, ndims(ro));
    rma(dim) = size(ro, dim);
    ro = ro - repmat(sum(ro .* (~bv), dim) ./ sum(~bv, dim, 'double'), rma);
end

% remove bad values
if opts.nozero
    ro(ro == 0) = eps;
end
bv = bv | isinf(ro) | isnan(ro);
ro(bv) = 0;
