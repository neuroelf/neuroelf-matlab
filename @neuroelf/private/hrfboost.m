function [hb, hw] = hrfboost(hd, opts)
% hrfboost  - assume that best fit is the (root of) sum of squares
%
% FORMAT:       [hb, hw] = hrfboost(hd [, opts])
%
% Input fields:
%
%       hd          HRF with derivates (as last non-singleton dimension)
%       opts        optional settings
%        .bf        TxB basis functions (used to construct summed shape)
%        .comp      computation, one of {'auc'}, 'boost', 'max', 'posauc'
%        .wcutoff   cut-off for output to be set to NaN / 0 (default: 0.5)
%        .wtype     weight function type, either of 'R2' or {'varfract'}
%
% Output fields:
%
%       hb          boosted HRF estimate
%       hw          weight for confidence of HRF estimate
%
% Note: if no basis function set is given, the computation is forced to
%       boost (sign of first beta * sqrt or sum of squares)

% Version:  v0.9c
% Build:    12042615
% Date:     Apr-09 2012, 12:04 PM EST
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
   ~isnumeric(hd)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'bf') || ...
   ~isa(opts.bf, 'double') || ...
    any(isinf(opts.bf(:)) | isnan(opts.bf(:))) || ...
    size(opts.bf, 1) < 2 || ...
   ~any([numel(hd), size(hd, ndims(hd))] == size(opts.bf, 2))
    opts.bf = [];
    opts.comp = 'b';
end
if ~isfield(opts, 'comp') || ...
   ~ischar(opts.comp) || ...
    isempty(opts.comp) || ...
   ~any('abchmp' == lower(opts.comp(1)))
    opts.comp = 'a';
else
    opts.comp = lower(opts.comp(1));
    if any('ch' == opts.comp)
        opts.comp = 'b';
    end
end
if ~isfield(opts, 'wcutoff') || ...
   ~isa(opts.wcutoff, 'double') || ...
    numel(opts.wcutoff) ~= 1 || ...
    isinf(opts.wcutoff) || ...
    isnan(opts.wcutoff) || ...
    opts.wcutoff < 0
    opts.wcutoff = 0.5;
else
    opts.wcutoff = min(1, opts.wcutoff);
end
if ~isfield(opts, 'wtype') || ...
   ~ischar(opts.wtype) || ...
    isempty(opts.wtype) || ...
   ~any('rv' == lower(opts.wtype(1)))
    opts.wtype = 'v';
else
    opts.wtype = lower(opts.wtype(1));
end

% get size and dimensions
sz = size(hd);
nd = ndims(hd);

% special case for empty input
if isempty(hd)
    if nd == 2 && ...
       (sz(2) == 1 || ...
        sz(2) == 0)
        sz = [0, 0];
    else
        sz(end) = 1;
    end
    hb = reshape(hd, sz);
    if nargout > 1
        hw = hb;
    end
    return;
end

% for 2-dim data check if sz(2) == 1
if nd == 2 && ...
    sz(2) == 1

    % switch dims
    hd = hd(:)';
    sz = size(hd);
end

% reshape input
ts = [sz(1:nd-1), 1];
hd = reshape(hd, [prod(ts), size(hd, nd)]);

% we need weights
if nargout > 1 || ...
    opts.comp ~= 'b'

    % produce product (sum across basis functions)
    if isempty(opts.bf)
        opts.bf = hrf('twogamma', 0.1, 5, 15, 6, 0, 1, 1, [], size(hd, 2) - 1);
    end

    % compute sum of scaled basis functions
    hb = hd * opts.bf';

    % normalization factor, depending on computation
    switch (opts.comp)

        % AUC (signed)
        case {'a'}
            nfac = 1 / sum(opts.bf(:, 1));

        % MAX
        case {'m'}
            nfac = 1 / max(max(opts.bf(:, 1)), -min(opts.bf(:, 1)));

        % positive-only AUC (unsigned positive)
        case {'p'}
            nfac = 1 / sum(opts.bf(opts.bf(:, 1) > 0, 1));
    end

    % compute weights
    if nargout > 1

        % R2 (weight is the correlation ^ 2 between sum(sBF) ./. main BF)
        if opts.wtype == 'r'

            % compute w
            [null, hw] = cov_nd(hb, hd(:, 1) * opts.bf(:, 1)');

            % square (and limit to [0, 1] range)
            hw = limitrangec(hw .* hw, 0, 1, 0);

        % variance portion
        else

            % compute variance portion (and limit to [0, 1] interval)
            hw = limitrangec(var(hd(:, 1) * opts.bf(:, 1)', [], 2) ./ var(hb, [], 2), ...
                0, 1, 0);
        end

        % reshape
        hw = reshape(hw, ts);
    end
end

% computation type
switch (opts.comp)

    % AUC (normalized, signed)
    case {'a'}
        hb = nfac .* reshape(sum(hb, 2), ts);

    % classic (Calhoun) HRF boost
    case {'b'}
        hb = reshape(sign(hd(:, 1)) .* sqrt(sum(hd .* hd, 2)), ts);

    % MAX
    case {'m'}

        % get max/min
        maxb = max(hb, [], 2);
        minb = min(hb, [], 2);

        % mask
        mskb = (abs(minb) > abs(maxb)) & (minb < 0);

        % replace values
        maxb(mskb) = minb(mskb);

        % computation
        hb = nfac .* reshape(maxb, ts);

    % positive-only AUC (unsigned)
    case {'p'}

        % replace <0 values
        hb(hb < 0) = 0;

        % computation
        hb = nfac .* reshape(sum(hb, 2), ts);
end
