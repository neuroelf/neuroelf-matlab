function hs = histscale(h, sod, opts)
% histscale  - scale existing histogram to new data size
%
% FORMAT:       hc = histscale(h, sod [, opts]);
%
% Input fields:
%
%       h           existing histogram
%       sod         new size (1x1 double) or data to match histogram with
%       opts        optional settings
%        .bins      bin defition (default: 0:numel(h)), given as
%                   1x2 double: from:to
%                   1x3 double: from:step:to
%                   1xH double: actual bins (H = numel(h) + 1)
%        .dist      generate distribution (default: false)
%        .distiv    distribution with inter-value ranges (default: false)
%        .smooth    kernel to smooth histogram (default: 0)
%
% Output fields:
%
%       hs          output histogram (or data)

% Version:  v0.9c
% Build:    12112316
% Date:     Nov-23 2012, 10:41 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2012, Jochen Weber
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
   (~isa(h, 'double') && ...
    ~isa(h, 'uint32')) || ...
    isempty(h) || ...
    numel(h) ~= max(size(h)) || ...
    any(isinf(h) | isnan(h) | h < 0) || ...
   ~isnumeric(sod) || ...
    isempty(sod) || ...
   (numel(sod) == 1 && ...
    ~isa(sod, 'double'))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
nh = numel(h);
h = double(h(:));
h = h ./ sum(h);
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'bins') || ...
   ~isa(opts.bins, 'double') || ...
   ~any(numel(opts.bins) == [2, 3, (numel(h) + 1)]) || ...
    any(isinf(opts.bins(:)) | isnan(opts.bins(:)));
    opts.bins = 0:numel(h);
elseif numel(opts.bins) ~= (numel(h) + 1)
    if numel(opts.bins) == 2
        opts.bins = opts.bins(1):opts.bins(2);
    else
        opts.bins = opts.bins(1):opts.bins(2):opts.bins(3);
    end
else
    opts.bins = sort(opts.bins(:))';
end
if ~isfield(opts, 'dist') || ...
   ~islogical(opts.dist) || ...
    numel(opts.dist) ~= 1
    opts.dist = false;
end
if ~isfield(opts, 'distiv') || ...
   ~islogical(opts.distiv) || ...
    numel(opts.distiv) ~= 1
    opts.distiv = false;
end
if ~isfield(opts, 'smooth') || ...
   ~isa(opts.smooth, 'double') || ...
    numel(opts.smooth) ~= 1 || ...
    isinf(opts.smooth) || ...
    isnan(opts.smooth) || ...
    opts.smooth < 0
    opts.smooth = 0;
else
    opts.smooth = min(max(16, sqrt(numel(h))), opts.smooth);
end

% smooth histogram
if opts.smooth > 0
    h = flexinterpn(h, [Inf; 1; 1; numel(h)], smoothkern(opts.smooth, 0), 1);
    h = h ./ sum(h);
end

% single input size
if numel(sod) == 1
    nsod = floor(sod);
else
    nsod = numel(sod);
    opts.dist = true;
end

% scale histogram
hs = round(nsod .* h);

% make sure numbers match
if sum(hs) ~= nsod

    % compute difference
    d = nsod - sum(hs);

    % more than size of h
    if abs(d) > (0.5 * nh)
        hs = hs + round(d / nh);
    end

    % compute rounding errors
    r = (nsod .* h) - hs;

    % sort by error
    [sr, si] = sort(r);

    % missing elements
    if d > 0

        % add to elements with largest positive error
        hs(si(end+1-d:end)) = hs(si(end+1-d:end)) + 1;

    % too many elements
    else

        % subtract from elements with largest negative error
        d = abs(d);
        hs(si(1:d)) = hs(si(1:d)) - 1;
    end
end

% generate distribution
if opts.dist

    % input
    hi = hs;

    % output
    if numel(sod) == 1
        hs = zeros(sum(hi), 1);
        ti = [];
    else
        [hs, ti] = sort(sod(:));
        hs = sod;
    end

    % iterate
    tix = 1;
    for sc = 1:nh

        % elements
        ne = hi(sc);

        % get binrange
        bfrom = opts.bins(sc);
        bto = opts.bins(sc + 1);
        bstep = (bto - bfrom) / ne;

        % intervalues
        if opts.distiv

            % target indices
            if ~isempty(ti)

                % generate and target-assign
                hs(ti(tix:tix+ne-1)) = bfrom:bstep:(bto - 0.5 * bstep);

            % no target indices
            else

                % generate and assign
                hs(tix:tix+ne-1) = bfrom:bstep:(bto - 0.5 * bstep);
            end

        % bin centers (close to!)
        else

            % target indices
            if ~isempty(ti)

                % target-assign value
                hs(ti(tix:tix+ne-1)) = 0.5 * (bfrom + bto - bstep);

            % no target indices
            else

                % assign value
                hs(tix:tix+ne-1) = 0.5 * (bfrom + bto - bstep);
            end
        end

        % increase counter
        tix = tix + ne;
    end
end
