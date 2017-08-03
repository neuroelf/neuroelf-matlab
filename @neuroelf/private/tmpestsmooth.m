function [fwhm, fi, cr] = tmpestsmooth(t, tr, opts)
% resestsmooth  - estimate temporal smoothness in signal
%
% FORMAT:       [fwhm, fi] = tmpestsmooth(t, tr [, opts])
%
% Input fields:
%
%       t           temporal data
%       tr          TR (sec/sample)
%       opts        optional settings struct
%        .detrend   detrend data first (default: false)
%        .filter    filter data first (default: true)
%        .ftcutoff  filter cutoff (in seconds, default: 128)
%        .ftgsig    filter global signal (default: false)
%        .ftype     filter type: {'dct'}, 'fourier', 'poly'
%        .maxlag    maximum lag to consider (default: 1)
%        .tdim      temporal dimension (default: last)
%
% Output fields:
%
%       fwhm        [X, Y, Z] overall FWHM estimate
%       fi          smoothness image
%       cr          auto-correlation images (bias corrected)

% Version:  v0.9c
% Build:    12121717
% Date:     Dec-12 2012, 2:28 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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
if nargin < 1 || ...
   ~isnumeric(t) || ...
    numel(t) < 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isa(tr, 'double') || ...
    numel(tr) ~= 1 || ...
    isinf(tr) || ...
    isnan(tr) || ...
    tr <= 0
    tr = 1;
end
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'tdim') || ...
   ~isa(opts.tdim, 'double') || ...
    numel(opts.tdim) ~= 1 || ...
    isinf(opts.tdim) || ...
    isnan(opts.tdim) || ...
   ~any(opts.tdim == 1:ndims(t))
    opts.tdim = ndims(t);
    if size(t, opts.tdim) == 1
        opts.tdim = 1;
    end
end
nt = size(t, opts.tdim);
if ~isfield(opts, 'maxlag') || ...
   ~isa(opts.maxlag, 'double') || ...
    numel(opts.maxlag) ~= 1 || ...
    isinf(opts.maxlag) || ...
    isnan(opts.maxlag)
    opts.maxlag = 1;
else
    opts.maxlag = round(min(floor(0.5 * nt), max(1, opts.maxlag)));
end
if ~isfield(opts, 'detrend') || ...
   ~islogical(opts.detrend) || ...
    numel(opts.detrend) ~= 1
    opts.detrend = false;
end
if ~isfield(opts, 'filter') || ...
   ~islogical(opts.filter) || ...
    numel(opts.filter) ~= 1
    opts.filter = true;
end
if ~isfield(opts, 'ftcutoff') || ...
   ~isa(opts.ftcutoff, 'double') || ...
    numel(opts.ftcutoff) ~= 1 || ...
    isnan(opts.ftcutoff) || ...
    opts.ftcutoff <= (2 * tr)
    opts.ftcutoff = 128;
end
if ~isfield(opts, 'ftgsig') || ...
   ~islogical(opts.ftgsig) || ...
    numel(opts.ftgsig) ~= 1
    opts.ftgsig = false;
end
if ~isfield(opts, 'ftype') || ...
   ~ischar(opts.ftype) || ...
    isempty(opts.ftype) || ...
   ~any(strcmpi(opts.ftype(:)', {'d', 'dct', 'f', 'fourier', 'p', 'poly'}))
    opts.ftype = 'd';
else
    opts.ftype = lower(opts.ftype(1));
end

% resolve data
if ~isa(t, 'double')
    t = double(t);
end

% correct dim for cov_nd
if opts.tdim ~= ndims(t)
    t = permute(t, [setdiff(1:ndims(t), opts.tdim), opts.tdim]);
end
tdim = ndims(t);
nl = opts.maxlag;

% get size
st = size(t);
st(tdim) = [];
if numel(st) == 1
    st(2) = 1;
end
sts = [prod(st), nt];

% filter
if opts.detrend
    if opts.filter
        X = ztrans((0:(1/(nt-1)):1)');
    else
        t = tempfilter(t, struct('tdim', tdim, 'tempdt', true));
    end
else
    X = zeros(nt, 0);
end
if opts.filter
    if opts.ftgsig
        X = [X, ztrans(meannoinfnan(reshape(t, sts), 1)')];
    end
    tfopts = struct('tdim', tdim, 'nuisreg', X);
    if opts.ftype == 'd' || ...
        opts.ftype == 'p'
        tfopts.tempdct = opts.ftcutoff / tr;
    else
        tfopts.tempsc = max(1, floor((nt * tr) / opts.ftcutoff));
    end
    t = tempfilter(t, tfopts);
end

% compute correlation maps
cr = zeros([prod(st), nl]);

% perform computation
for lc = 1:nl

    % auto-correlation
    [cnv, cnc] = cov_nd(t, t, lc);

    % bias-correction
    cr(:, lc) = min(1, cnc(:) + 0.6 / (nt - 1));
end

% compute FWHM estimate
ffac = (nt * tr * sqrt(8 * log(2)) / (nt - 1));
fwhm = ffac .* sqrt(-1 / (4 * log(max(0, ...
    fisherr2z(meannoinfnan(fisherr2z(cr(:, 1)), 1, true), true)))));
cr(cr == 1) = NaN;
if nargin > 1
    fi = ffac .* sqrt(-1 ./ (4 .* log(max(0, cr(:, 1)))));
    fi = reshape(fi, st);
    if nargout > 2
        if st(end) == 1 && ...
            nl > 1
            st(end) = [];
        end
        cr = reshape(cr, [st, nl]);
    end
end
