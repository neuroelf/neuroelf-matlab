function [p, h] = prctilend(l, v, opts)
% prctilend  - N-dim percentile lookup
%
% FORMAT:       [p, h] = prctilend(l, v [, opts])
%
% Input fields:
%
%       l           lookup field (2d histogram or raw data)
%       v           variate (N-d, used to lookup percentile)
%       opts        optional settings
%        .cumsum    boolean flag, compute cumsum over lookup field (true)
%        .dim       dim along to lookup (default: first that matches)
%        .interp    interpolation, either of 'linear', {'nearest'}
%        .ltype     either of 'hist' or {'raw'}
%        .mask      mask to apply to v before lookup
%        .normalize normalization value, 1x2 double ([0, 1])
%        .rmax      range maximum (maximally reachable/reached value, 0)
%        .rmin      range maximum (maximally reachable/reached value, 1)
%        .rstep     range step size (default: 0.001)
%        .tails     tails of test, one of {'both'}, 'negative', 'positive'
%
% Output fields:
%
%       p           percentile of each variate value
%       h           histogram of l between rmin:rstep:rmax (only for data)
%
% Note: either of the dimensions of l must be equal to numel(v)
%
%       the .cumsum, .interp, and .rmax/.rmin flags are only used
%       if .ltype is 'hist'
%
%       the output will be scaled between .normalize(1) and .normalize(2)

% Version:  v0.9c
% Build:    11052414
% Date:     May-19 2011, 12:20 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
if nargin < 2 || ...
   ~isnumeric(l) || ...
    ndims(l) > 2 || ...
    numel(l) < 2 || ...
   ~isnumeric(v)
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
if ~isfield(opts, 'cumsum') || ...
   ~islogical(opts.cumsum) || ...
    numel(opts.cumsum) ~= 1
    opts.cumsum = true;
end
if ~isfield(opts, 'dim') || ...
   ~isa(opts.dim, 'double') || ...
    numel(opts.dim) ~= 1 || ...
   (~isequal(opts.dim, 1) && ...
    ~isequal(opts.dim, 2))
    opts.dim = [];
end
if ~isfield(opts, 'interp') || ...
   ~ischar(opts.interp) || ...
    isempty(opts.interp) || ...
   ~any(lower(opts.interp(1)) == 'eln')
    opts.interp = 'n';
else
    opts.interp = lower(opts.interp(1));
end
if ~isfield(opts, 'ltype') || ...
   ~ischar(opts.ltype) || ...
    isempty(opts.ltype) || ...
   ~any(lower(opts.ltype(1)) == 'dhnr')
    opts.ltype = 'r';
else
    opts.ltype = lower(opts.ltype(1));
end
if ~isfield(opts, 'mask') || ...
   ~islogical(opts.mask) || ...
    numel(opts.mask) ~= numel(v)
    mask = [];
else
    mask = opts.mask(:);
    sv = size(v);
end
if ~isfield(opts, 'normalize') || ...
   ~isa(opts.normalize, 'double') || ...
    numel(opts.normalize) ~= 2 || ...
    any(isinf(opts.normalize) | isnan(opts.normalize)) || ...
    opts.normalize(1) == opts.normalize(2)
    opts.normalize = [0, 1];
end
if ~isfield(opts, 'rmax') || ...
   ~isa(opts.rmax, 'double') || ...
    numel(opts.rmax) ~= 1 || ...
    isinf(opts.rmax) || ...
    isnan(opts.rmax)
    opts.rmax = 1;
end
if ~isfield(opts, 'rmin') || ...
   ~isa(opts.rmin, 'double') || ...
    numel(opts.rmin) ~= 1 || ...
    isinf(opts.rmin) || ...
    isnan(opts.rmin)
    opts.rmin = 0;
end
if opts.rmax <= opts.rmin
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid combination of .rmax/.rmin options.' ...
    );
end
fr = opts.rmax - opts.rmin;
if ~isfield(opts, 'rstep') || ...
   ~isa(opts.rstep, 'double') || ...
    numel(opts.rstep) ~= 1 || ...
    isinf(opts.rstep) || ...
    isnan(opts.rstep) || ...
    opts.rstep <= 0
    opts.rstep = 0.001 * fr;
end
if ~isfield(opts, 'tails') || ...
   ~ischar(opts.tails) || ...
    isempty(opts.tails) || ...
   ~any(lower(opts.tails(1)) == '2bnp')
    opts.tails = 'b';
else
    opts.tails = lower(opts.tails(1));
end

% get/check dim
sdim = false;
if isempty(opts.dim)
    if isempty(mask)
        fd = 3 - findfirst(size(l) == numel(v));
    else
        fd = 3 - findfirst(size(l) == sum(mask));
    end
    if isempty(fd)
        sdim = true;
        fd = findfirst(size(l) > 1);
    end
else
    fd = opts.dim;
    if (isempty(mask) && ...
        size(l, 3 - fd) ~= numel(v)) || ...
       (~isempty(mask) && ...
        size(l, 3 - fd) ~= sum(mask))
        error( ...
            'neuroelf:BadArgument', ...
            'Dimensions mismatch.' ...
        );
    end
end

% for histogram data
if opts.ltype == 'h'

    % convert to double
    l = double(l);

    % get size (-1!)
    if ~isempty(fd)
        ld = size(l, fd) - 1;
    end

    % cumsum?
    if opts.cumsum
        l = cumsum(l, fd);
    end

    % compute lookup value
    if isempty(mask)
        v = 1 + (ld / fr) .* (v(:) - opts.rmin);
    else
        v = 1 + (ld / fr) .* (lsqueeze(v(mask)) - opts.rmin);
    end
    vi = (1:numel(v))';

    % all from a single distribution
    if sdim

        % simply look up
        if opts.interp ~= 'l'
            p = reshape(l(round(limitrangec(v(:), 1, ld + 1))), size(v));

        % or interpolate
        else
            p = reshape(flexinterpn(l, (limitrangec(v(:), 1, ld + 1))), size(v));
        end

    % separate distributions
    else

        % along first dim
        if fd == 1

            % nearest
            if opts.interp ~= 'l'

                % use indexarray
                p = reshape(indexarray(l, ...
                    [limitrangec(v, 1, ld + 1), vi]), size(v));

            % linear
            else

                % use flexinterpn (with default kernel)
                p = reshape(flexinterpn(l, ...
                    [limitrangec(v, 1, ld + 1), vi]), size(v));
            end

            % set values before first edge to 0
            p(v < 1) = 0;

            % apply normalization
            p = p ./ reshape(l(end, :), size(v));

        % along second dim
        else

            % nearest
            if opts.interp ~= 'l'

                % use indexarray
                p = reshape(indexarray(l, ...
                    [vi, limitrangec(v, 1, ld + 1)]), size(v));

            % linear
            else

                % use flexinterpn (with default kernel)
                p = reshape(flexinterpn(l, ...
                    [vi, limitrangec(v, 1, ld + 1)]), size(v));
            end

            % set values before first edge to 0
            p(v < 1) = 0;

            % apply normalization
            p = p ./ reshape(l(:, end), size(v));
        end
    end

% data
else

    % mask data
    if ~isempty(mask)
        v = lsqueeze(v(mask));
    end

    % start with 0 assignment
    p = zeros(size(v));
    ld = size(l, fd);

    % also create histogram
    if nargout > 1
        h = zeros(1, numel(opts.rmin:opts.rstep:opts.rmax));
        harg = {opts.rmin, opts.rmax, opts.rstep};
        hdo = true;
    else
        hdo = false;
    end

    % along first dim
    if fd == 1

        % transpose v
        v = reshape(v, 1, numel(v));

        % iterate along diml
        for lc = 1:ld

            % get lookup part
            lp = l(lc, :);

            % compare value to data
            pa = (v > lp);

            % then add
            p(pa) = p(pa) + 1;

            % and also add .5 for equal values
            pa = (v == lp);
            p(pa) = p(pa) + 0.5;

            % also add to histogram
            if hdo
                h = h + histcount(lp, harg{:});
            end
        end

    % along second dim
    else

        % iterate along dim2
        for lc = 1:ld

            % get lookup part
            lp = l(:, lc);

            % compare value to data
            pa = (v > lp);

            % then add
            p(pa) = p(pa) + 1;

            % and also add .5 for equal values
            pa = (v == lp);
            p(pa) = p(pa) + 0.5;

            % also add to histogram
            if hdo
                h = h + histcount(lp, harg{:});
            end
        end
    end
end

% compute percentile
p = p ./ ld;

% handle tails
if opts.tails ~= 'n'
    p(p > 0.5) = p(p > 0.5) - (1 / ld);
end
if opts.tails ~= 'p'
    p(p < 0.5) = p(p < 0.5) + (1 / ld);
end

% apply normalization
nm = opts.normalize(1);
nr = opts.normalize(2) - nm;
if nr ~= 1
    p = nr .* p;
end
if nm ~= 0
    p = p + nm;
end

% mask
if ~isempty(mask)
    po = p;
    p = zeros(sv);
    p(mask) = po;
    mask = ~mask;

    % correct "tailing"
    if opts.tails == 'n'
        p(mask) = opts.normalize(2);
    elseif opts.tails ~= 'p'
        p(mask) = nm + 0.5 * nr;
    end
end
