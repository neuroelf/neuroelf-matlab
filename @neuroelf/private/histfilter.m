function fd = histfilter(d, r, opts)
% histfilter  - histogram based filter
%
% FORMAT:       d = histfilter(d, r [, opts])
%
% Input fields:
%
%       d           2D/3D data
%       r           radius (shape, see opts)
%       opts        optional settings
%        .histbins  number of histogram bins, default: 64
%        .histlims  1x2 limits for histogram (default: auto-detect)
%        .histsm    histogram smoothing, default: sqrt(histbins)
%        .mask      mask volume (only operate on those voxels/pixels)
%        .method    either of 'median', {'mode'}, 'rmean'
%        .mrange    1x2 range used for rmean method, default [0.25, 0.75]
%        .pbar      progress bar object
%        .prange    progress range
%        .shape     either of 'box' or {'sphere'} (equals square vs. circle)
%        .smooth    smoothing factor (default: 1, must be 0 < .smooth <= 1)
%
% Output fields:
%
%       d           filtered data

% Version:  v0.9d
% Build:    14052519
% Date:     May-25 2014, 7:33 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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

% parse options
if nargin < 2 || ...
   ~isnumeric(d) || ...
    isempty(d) || ...
    ndims(d) > 3 || ...
   ~isa(r, 'double') || ...
    numel(r) ~= 1 || ...
    isinf(r) || ...
    isnan(r) || ...
    r < 1
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if istransio(d)
    d = resolve(d);
end
if isinteger(d)
    d = single(d);
end
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'histbins') || ...
   ~isa(opts.histbins, 'double') || ...
    numel(opts.histbins) ~= 1 || ...
    isinf(opts.histbins) || ...
    isnan(opts.histbins) || ...
    opts.histbins < 4
    opts.histbins = 64;
else
    opts.histbins = min(1024, ceil(opts.histbins));
end
if ~isfield(opts, 'histlims') || ...
   ~isa(opts.histlims, 'double') || ...
    numel(opts.histlims) ~= 1 || ...
    any(isinf(opts.histlims) | isnan(opts.histlims)) || ...
    opts.histlims(1) <= opts.histlims(2)
    if isa(d, 'double') || ...
        isa(d, 'single')
        opts.histlims = minmaxmean(d, 4);
    else
        opts.histlims = minmaxmean(d);
    end
    opts.histlims(3:end) = [];
end
histfrom = opts.histlims(1);
histto = opts.histlims(2);
histstep = (histto - histfrom) / opts.histbins;
histto = histto + 1.5 * histstep;
if ~isfield(opts, 'histsm') || ...
   ~isa(opts.histsm, 'double') || ...
    numel(opts.histsm) ~= 1 || ...
    isinf(opts.histsm) || ...
    isnan(opts.histsm) || ...
    opts.histsm < 0
    opts.histsm = sqrt(opts.histbins);
else
    opts.histsm = min(2 * sqrt(opts.histbins), opts.histsm);
end
histsmk = {[0; 1; 0], smoothkern(opts.histsm)};
histsmk{2}(histsmk{2} < (0.001 * max(histsmk{2}))) = [];
histsmf = {1, 1};
if ~isfield(opts, 'mask') || ...
   ~islogical(opts.mask) || ...
   ~isequal(size(opts.mask), size(d))
    opts.mask = false(0, 0);
end
if ~isfield(opts, 'method') || ...
   ~ischar(opts.method) || ...
    isempty(opts.method) || ...
   ~any(strcmpi(opts.method(:)', {'median', 'mode', 'rmean'}))
    opts.method = 'mode';
else
    opts.method = lower(opts.method(:)');
end
switch (opts.method)
    case {'median'}
        mmeth = 0;
    case {'mode'}
        mmeth = 1;
    case {'rmean'}
        mmeth = 2;
end
if ~isfield(opts, 'mrange') || ...
   ~isa(opts.mrange, 'double') || ...
    numel(opts.mrange) ~= 2 || ...
    any(isinf(opts.mrange) | isnan(opts.mrange) | opts.mrange < 0 | opts.mrange > 1) || ...
    opts.mrange(1) >= opts.mrange(2)
    opts.mrange = [0.25, 0.75];
end
mrfrom = opts.mrange(1);
mrto = opts.mrange(2);
if ~isfield(opts, 'pbar') || ...
    numel(opts.pbar) ~= 1 || ...
   (~isa(opts.pbar, 'xfigure') && ...
    ~isa(opts.pbar, 'xprogress'))
    opts.pbar = [];
else
    pbarvis = opts.pbar.Visible;
end
pbar = opts.pbar;
if ~isfield(opts, 'prange') || ...
   ~isa(opts.prange, 'double') || ...
    numel(opts.prange) ~= 2 || ...
    any(isinf(opts.prange) | isnan(opts.prange) | opts.prange < 0 | opts.prange > 1) || ...
    opts.prange(1) <= opts.prange(2)
    opts.prange = [0, 1];
end
pmin = opts.prange(1);
pdiff = max(0.001, opts.prange(2) - pmin);
if ~isfield(opts, 'shape') || ...
   ~ischar(opts.shape) || ...
    isempty(opts.shape) || ...
   ~any(lower(opts.shape(1)) == 'bs')
    opts.shape = 's';
else
    opts.shape = lower(opts.shape(1));
end
if opts.shape == 's'
    cr = ceil(r);
else
    cr = floor(r);
end
if ~isfield(opts, 'smooth') || ...
   ~isa(opts.smooth, 'double') || ...
    numel(opts.smooth) ~= 1 || ...
    isinf(opts.smooth) || ...
    isnan(opts.smooth) || ...
    opts.smooth <= 0 || ...
    opts.smooth > 1
    opts.smooth = 1;
end
smf = opts.smooth;
smfm = 1 - smf;

% data size and elements
sd = size(d);
dd = numel(sd);

% empty mask
if isempty(opts.mask)
    mi = [];
    nd = numel(d);
else
    mi = find(opts.mask(:));
    nd = numel(mi);
    d(~opts.mask) = NaN;
end

% array to hold coordinates of shape
s = cell(1, dd);

% span box (up to required size)
b = repmat({-cr:cr}, 1, 3);
[s{:}] = ndgrid(b{1:dd});

% linearize coordinates (to be stored in numeric array s)
for dc = 1:dd
    s{dc} = s{dc}(:);
end

% concatenate
s = cat(2, s{:});

% compute length
sl = sqrt(sum(s .* s, 2));

% sort
[sl, sli] = sort(sl);
s = s(sli, :);

% for spherical (circular) shape
if opts.shape == 's'

    % remove coordinates > radius from 0
    s(sl > r, :) = [];
end

% get number of elements
ns = size(s, 1);

% process no more than 10M elements at once
if mmeth ~= 1
    stepsize = ceil(1e7 / max(ns, opts.histbins));
else
    stepsize = ceil(1e7 / max(ns, 4 * opts.histbins));
end

% generate data array holding the sampled data
did = d(1);
did(1:stepsize, 1:ns) = 0;

% histogram sampling
if mmeth ~= 1
    hmsi = [Inf, Inf; 1, 1; 1, 1; stepsize, opts.histbins + 1];
    osh = ones(1, ns);
else
    hmsi = [Inf, Inf; 1, 1; 1, 0.25; stepsize, opts.histbins + 1];
end

% process
fd = d;
idx = cell(1, dd);
os = ones(stepsize, 1);
oh = ones(1, opts.histbins + 2);
hinv = histto + 3 * histstep;
for v = 1:stepsize:nd

    % progress
    if ~isempty(pbar)
        pval = (v - 1) / nd;
        pbar.Progress(pmin + pdiff * pval, ...
            sprintf('Histogram filtering (%.1f%% done)', 100 * pval));
        if v == 1
            pbar.Visible = 'on';
        end
    end

    % generate indices
    vto = v + stepsize - 1;
    if vto > nd
        remi = vto - nd;
        vto = nd;

        % compute new stepsize
        stepsize = stepsize - remi;
        os = ones(stepsize, 1);
        hmsi(4, 1) = stepsize;

        % and shrink required array
        did(stepsize+1:end, :) = [];
    end
    if isempty(mi)
        [idx{:}] = ind2sub(sd, lsqueeze(v:vto));
    else
        [idx{:}] = ind2sub(sd, lsqueeze(mi(v:vto)));
    end
    for dc = 1:dd
        idx{dc} = idx{dc}(:);
    end
    idxa = cat(2, idx{:});

    % sample
    did(:, 1) = lsqueeze(indexarray(d, idxa));
    for sc = 2:ns
        did(:, sc) = lsqueeze(indexarraynb(d, idxa + os * s(sc, :)));
    end

    % remove bad entries
    did(did == 0 | isinf(did) | isnan(did)) = hinv;

    % compute histogram
    hm = histcount(did, histfrom, histto, histstep, 2);

    % remove last entry
    hm(:, end) = 0;

    % for all but mode
    if mmeth ~= 1

        % scale
        hm = hm ./ (sum(hm, 2) * oh);

        % cumulative sum
        chm = cumsum(hm, 2);

    % for mode
    else

        % compute smooth histogram
        hm = flexinterpn(hm, hmsi, histsmk, histsmf);
    end

    % depending on method -> median
    if mmeth == 0

        % greater than 0.5
        chm = histfrom + histstep .* (max(1, sum(chm <= 0.5, 2)) - 1);

    % -> mode
    elseif mmeth == 1

        % find peak
        chm = histfrom + (0.25 * histstep) .* (maxpos(hm, [], 2) - 1);

    % -> range
    else

        % determine ranges from and to
        cfrom = histfrom + histstep .* sum(chm <= mrfrom, 2);
        cto = histfrom + histstep .* (max(1, sum(chm <= mrto, 2)) - 1);

        % set to NaN in data to mean
        did(did < (cfrom(:) * osh) | did > (cto(:) * osh)) = NaN;

        % final data
        chm = meannoinfnan(did, 2);
    end

    % set initial value to 0
    chm(chm == histfrom) = 0;

    % store
    if isempty(mi)
        if smf == 1
            fd(v:vto) = chm;
        else
            fd(v:vto) = smf .* chm + smfm .* lsqueeze(d(v:vto));
        end
    else
        if smf == 1
            fd(mi(v:vto)) = chm;
        else
            fd(mi(v:vto)) = smf .* chm + smfm .* lsqueeze(d(mi(v:vto)));
        end
    end
end

% adapt overall contrast
if isempty(mi)
    fd = fd .* (meannoinfnan(d(:), 1, true) / meannoinfnan(fd(:), 1, true));
else
    fd(opts.mask) = fd(opts.mask) .* ...
        (meannoinfnan(d(opts.mask), 1, true) / meannoinfnan(fd(opts.mask), 1, true));
end

% visible flag of progress
if ~isempty(pbar)
    pbar.Visible = pbarvis;
end
