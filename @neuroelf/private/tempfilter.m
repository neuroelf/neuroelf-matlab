function [fd, X, w] = tempfilter(d, opts)
% tempfilter  - filter temporal data
%
% FORMAT:       [fd, X, w] = tempfilter(d, opts)
%
% Input fields:
%
%       d           data to filter
%       opts        mandatory struct but with optional fields
%        .nuisreg   nuisance regressors, also added to filter matrix
%        .orthpoly  orthogonalize polynomials (faster computation later)
%        .pbar      handle to xfigure::Progress or xprogress object
%        .prange    range to display progress over (default: [0, 1])
%        .spat      enable spatial filtering (default: false)
%        .spkern    smoothing kernel in sampling units (default: [2, 2, 2])
%        .tdim      temporal filter dimension (default: 1)
%        .temp      enable temporal filtering (default: true)
%        .tempdct   DCT-based filtering (min. wavelength, default: Inf)
%        .tempdt    detrend (default: true, is overriden by dct/sc)
%        .temphp    temporal highpass (inv. smoothing) in units (def: 0)
%        .temphpr   regress highpass instead of subtract (default: false)
%        .templp    temporal lowpass (smoothing) kernel in units (def: 0)
%        .temppoly  set of orthogonal polynomials (default: 0)
%        .tempsc    sin/cos set of frequencies (number of pairs, def: 0)
%        .trobust   perform temporal filtering robustly (default: false)
%
% Output fields:
%
%       fd          filtered data (in input datatype, scaled if needed)
%       X           filtering matrix
%       w           weights of robust regression
%
% Note: low-pass filtering now uses flexinterpn (speed gain) and the
%       smoothing kernel is determined using a R=8 sinc window interpolator

% Version:  v0.9d
% Build:    14072317
% Date:     Jul-23 2014, 5:44 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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
   (~isnumeric(d) && ...
    ~istransio(d)) || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument in call.' ...
    );
end
if ~isfield(opts, 'nuisreg') || ...
   ~isa(opts.nuisreg, 'double') || ...
    ndims(opts.nuisreg) ~= 2 || ...
    any(isinf(opts.nuisreg(:)) | isnan(opts.nuisreg(:)))
    opts.nuisreg = [];
end
if ~isfield(opts, 'orthpoly') || ...
   (~islogical(opts.orthpoly) && ...
    ~isnumeric(opts.orthpoly)) || ...
    numel(opts.orthpoly) ~= 1
    opts.orthpoly = false;
else
    opts.orthpoly = (double(opts.orthpoly) ~= 0);
end
if ~isfield(opts, 'pbar') || ...
    numel(opts.pbar) ~= 1 || ...
   ~any(strcmpi(class(opts.pbar), {'xfigure', 'xprogress'}))
    opts.pbar = [];
end
if ~isfield(opts, 'prange') || ...
   ~isa(opts.prange, 'double') || ...
    numel(opts.prange) ~= 2 || ...
    any(isinf(opts.prange) | isnan(opts.prange) | opts.prange < 0 | opts.prange > 1)
    opts.prange = [0, 1];
else
    opts.prange = sort(opts.prange(:))';
    if opts.prange(2) == opts.prange(1)
        opts.prange(2) = opts.prange(1) + sqrt(eps);
    end
end
if ~isfield(opts, 'tdim') || ...
   ~isa(opts.tdim, 'double') || ...
    numel(opts.tdim) ~= 1 || ...
   ~any(opts.tdim == (1:ndims(d)))
    opts.tdim = 1;
end
if numel(d) == size(d, opts.tdim)
    opts.spat = false;
end
if ~isfield(opts, 'spat') || ...
   (~islogical(opts.spat) && ...
    ~isnumeric(opts.spat)) || ...
    numel(opts.spat) ~= 1
    opts.spat = false;
else
    opts.spat = (opts.spat ~= 0);
end
if ~isfield(opts, 'spkern') || ...
   ~isa(opts.spkern, 'double') || ...
    numel(opts.spkern) ~= 3 || ...
    any(isinf(opts.spkern) | isnan(opts.spkern) | opts.spkern < 0 | opts.spkern > 10)
    opts.spkern = [2, 2, 2];
else
    opts.spkern = opts.spkern(:)';
end
if all(opts.spkern) < 0.5
    opts.spat = false;
end
if ~isfield(opts, 'temp') || ...
   (~islogical(opts.temp) && ...
    ~isnumeric(opts.temp)) || ...
    numel(opts.temp) ~= 1
    opts.temp = true;
else
    opts.temp = (opts.temp ~= 0);
end
if ~isfield(opts, 'tempdct') || ...
   ~isa(opts.tempdct, 'double') || ...
    numel(opts.tempdct) ~= 1 || ...
    isinf(opts.tempdct) || ...
    isnan(opts.tempdct) || ...
    opts.tempdct < 4
    opts.tempdct = Inf;
end
if ~isfield(opts, 'tempdt') || ...
   ~islogical(opts.tempdt) || ...
    numel(opts.tempdt) ~= 1
    opts.tempdt = true;
end
if ~isfield(opts, 'temphp') || ...
   ~isa(opts.temphp, 'double') || ...
    numel(opts.temphp) ~= 1 || ...
    isinf(opts.temphp) || ...
    isnan(opts.temphp) || ...
    opts.temphp < 0
    opts.temphp = 0;
elseif opts.temphp > 128
    warning(...
        'neuroelf:LongComputation', ...
        'Overlarge high-pass filter. Long computation, try using different method.' ...
    );
end
if ~isfield(opts, 'temphpr') || ...
   ~islogical(opts.temphpr) || ...
    numel(opts.temphpr) ~= 1
    opts.temphpr = false;
end
if ~isfield(opts, 'templp') || ...
   ~isa(opts.templp, 'double') || ...
    numel(opts.templp) ~= 1 || ...
    isinf(opts.templp) || ...
    isnan(opts.templp) || ...
    opts.templp <= 0
    opts.templp = 0;
elseif opts.templp > 128
    warning(...
        'neuroelf:LongComputation', ...
        'Overlarge low-pass filter. Long computation, try using different method.' ...
    );
end
if ~isfield(opts, 'temppoly') || ...
   ~isa(opts.temppoly, 'double') || ...
    numel(opts.temppoly) ~= 1 || ...
    isinf(opts.temppoly) || ...
    isnan(opts.temppoly) || ...
    opts.temppoly < 0 || ...
    opts.temppoly >= (size(d, opts.tdim) - 1)
    opts.temppoly = 0;
else
    opts.temppoly = floor(opts.temppoly);
end
if ~isfield(opts, 'tempsc') || ...
   ~isa(opts.tempsc, 'double') || ...
    numel(opts.tempsc) ~= 1 || ...
   ~any(opts.tempsc == (1:floor(0.5 * (size(d, opts.tdim)-1))))
    opts.tempsc = 0;
end
if isinf(opts.tempdct) && ...
   ~opts.tempdt && ...
    opts.temphp == 0 && ...
    opts.templp == 0 && ...
    opts.temppoly == 0 && ...
    opts.tempsc == 0
    opts.temp = false;
elseif opts.temppoly  > 0
    opts.tempdct = Inf;
    opts.tempdt = false;
    opts.tempsc = 0;
elseif opts.tempsc > 0
    opts.tempdct = Inf;
    opts.tempdt = false;
end
if ~isfield(opts, 'trobust') || ...
   ~islogical(opts.trobust) || ...
    numel(opts.trobust) ~= 1
    opts.trobust = false;
end

% return early if neither spatial/temporal filtering and no nuisreg
if ~opts.spat && ...
   ~opts.temp && ...
    isempty(opts.nuisreg)
    fd = d;
    return;
end

% initialize filtering matrix
X = zeros(size(d, opts.tdim), 0);

% get data
dt = class(d);
if istransio(d)
    if ~strcmp(dt, 'double')
        fd = single(resolve(d));
    else
        fd = resolve(d);
    end
elseif ~strcmp(dt, 'double')
    fd = single(d);
else
    fd = d;
end

% permute if necessary
di = opts.tdim;
tperm = 1:ndims(fd);
if di > 1
    tperm = [di, tperm(tperm ~= di)];
    [sperm, fperm] = sort(tperm(:));
    fd = permute(fd, tperm);
end

% get original sizing and number of time points
vs = size(fd);
if numel(vs) < numel(tperm)
    vs(end+1:numel(tperm)) = 1;
end
nv = vs(1);

% use nuisance regressors ?
if size(opts.nuisreg, 1) ~= nv
    opts.nuisreg = zeros(nv, 0);
elseif ~isempty(opts.nuisreg)

    % remove invalid entries
    opts.nuisreg(:, any(isinf(opts.nuisreg) | isnan(opts.nuisreg)) | ...
        sum(abs(diff(opts.nuisreg))) < sqrt(eps)) = [];

    % z-transform (for matrix)
    opts.nuisreg = ztrans(opts.nuisreg);

    % and force temporal filtering
    opts.temp = true;
end

% DCT overrides detrending
if opts.temp && ...
   ~isinf(opts.tempdct)
    opts.tempdct = floor(2 * nv / opts.tempdct);
    opts.tempdt = false;
end

% spatial filter first if needed
if opts.spat

    % if too few dims, make it so
    fds = size(fd);
    if numel(fds) < 4
        fd = reshape(fd, [fds(1), ones(1, 4 - numel(fds)), fds(2:end)]);
        opts.spkern(1:(4 - numel(fds))) = 0;
    end
    ssr = size(fd);
    ssf = ssr(2:end);
    ssr(1) = 1;

    % use smoothdata3 to iterate over volumes
    for vc = 1:nv
        fd(vc, :, :, :, :) = reshape( ...
            smoothdata3(reshape(fd(vc, :, :, :, :), ssf), opts.spkern), ssr);
    end

    % reshaping again?
    if numel(fds) < 4
        fd = reshape(fd, fds);
    end
end

% temporal filter
if opts.temp

    % reshape if necessary
    if numel(vs) > 2
        fd = reshape(fd, nv, prod(vs(2:end)));
    end

    % only detrend
    if opts.tempdt

        % build X
        X = [ones(nv, 1), (-1:(2/(nv-1)):1)', opts.nuisreg];

    % DCT-based filtering
    elseif ~isinf(opts.tempdct)

        % build X
        X = [ones(nv, 1), zeros(nv, opts.tempdct), opts.nuisreg];
        n = 0:(nv - 1);
        for dc = 1:opts.tempdct
            X(:, dc + 1) = cos(pi * (2 * n + 1) * dc / (2 * nv));
        end

    % (legendre) polynomial based filtering
    elseif opts.temppoly > 0

        % build X
        X = [ones(nv, 1), zeros(nv, opts.temppoly), opts.nuisreg];
        n = (1 / nv) .* ((-nv + 1):2:(nv - 1))';
        n = 0.5 .* (abs(n) + abs(n(end:-1:1)));
        n(1:(floor(0.5*nv))) = -n(1:(floor(0.5*nv)));
        if mod(nv, 2) ~= 0
            n(ceil(0.5*nv)) = 0;
        end
        X(:, 2) = n;
        for dc = 2:opts.temppoly
            X(:, dc + 1) = (1 / (dc + 1)) .* ...
                ((2 * dc + 1) .* n .* X(:, dc) - dc .* X(:, dc - 1));
        end

        % enforce orthogonality
        if opts.orthpoly
            for dc = opts.temppoly:-1:2
                X(:, dc + 1) = X(:, dc + 1) - X(:, 1:dc) * ...
                    ((X(:, 1:dc)' * X(:, 1:dc)) \ (X(:, 1:dc)' * X(:, dc + 1)));
            end
        end

        % re-scale
        X(:, 2:end) = X(:, 2:end) ./ (ones(nv, 1) * max(abs(X(:, 2:end)), [], 1));

    % sin/cos set filtering
    elseif opts.tempsc > 0

        % build X
        X = [ones(nv, 1), zeros(nv, 2 * opts.tempsc), opts.nuisreg];
        n = 0:(nv - 1);
        for pc = 1:opts.tempsc
            X(:, (pc * 2)) = sin(pi * (2 * n + 1) * pc / nv);
            X(:, (pc * 2 + 1)) = cos(pi * (2 * n + 1) * pc / nv);
        end

    % only nuisance regressors
    else

        % simply use those
        X = [ones(nv, 1), opts.nuisreg];
    end

    % perform high-pass filter first
    if opts.temphp > 0

        % for regression
        if opts.temphpr

            % perform lowpass filter
            fdl = ztrans(tempfilter(fd, struct('tempdt', false, 'templp', opts.temphp)));

            % then regress out
            b = sum((1 / nv) .* (fdl .* ztrans(fd)));
            fd = fd - repmat(b, nv, 1) .* fdl;
        else

            % subtraction (simply inverse of low-pass filter)
            fd = ones(size(fd, 1), 1) * mean(fd, 1) + ...
                (fd - tempfilter(fd, struct('tempdt', false, 'templp', opts.temphp)));
        end
    end

    % perform regression
    if size(X, 2) > 1
        if opts.trobust
            if nargout < 3
                b = fitrobustbisquare_img(X, fd', [], [], opts)';
            else
                [b, w] = fitrobustbisquare_img(X, fd', [], [], opts);
                b = b';
                w = reshape(w', vs);
                if di > 1
                    w = permute(w, fperm');
                end

            end
        else
            b = inv(X' * X);
            b(abs(b) < sqrt(eps)) = 0;
            if sum(b(:) ~= 0) <= (0.25 * numel(b))
                b = sparse(b);
            end
            b = b * (X' * double(fd));
            if nargout > 2
                if di > 1
                    w = ones(vs(fperm));
                else
                    w = ones(vs);
                end
            end
        end
    else
        b = 0;
        if nargout > 2
            if di > 1
                w = ones(vs(fperm));
            else
                w = ones(vs);
            end
        end
    end

    % keep mean
    b(1, :) = 0;

    % subtract from data
    if size(b, 1) > 1
        fd = fd - X * b;
    end

    % low-pass filter?
    if opts.templp > 0

        % get kernel
        if opts.templp > 3
            k = smoothkern(opts.templp, 1.42e-5 / sqrt(opts.templp), ...
                false, 'lanczos8');
        else
            k = smoothkern(opts.templp, 1e-6, false, 'lanczos8');
        end

        % use flexinterpn to do the smoothing
        if size(fd, 2) > 1
            ficrd = [Inf, Inf; 1, 1; 1, 1; size(fd)];
            fikern = {k, [0;1;0]};
            fiksiz = {1, 1};
        else
            ficrd = [Inf; 1; 1; numel(fd)];
            fikern = k;
            fiksiz = 1;
        end

        % use flexinterpn
         fd = flexinterpn(fd, ficrd, fikern, fiksiz);
    end

    % reshape again
    fd = reshape(fd, vs);

    % remove mean from filter matrix
    if nargout > 1
        X(:, 1) = [];
    end
end

% datatype
if ~strcmp(dt, class(fd))

    % if not single
    if ~strcmpi(dt, 'single')

        % get minmaxmean
        mmm = minmaxmean(fd);
        rng = eps + (mmm(2) - mmm(1));

        % then recompute
        switch(dt)
            case {'int16'}
                if mmm(1) < -32768 || ...
                    mmm(2) >= 32767.5
                    fd = -32768 + (65535 / rng) * (fd - mmm(1));
                end
            case {'int32'}
                if mmm(1) < -2147483648 || ...
                    mmm(2) >= 2147483647.5
                    fd = -2147483648 + (4294967295 / rng) * (fd - mmm(1));
                end
            case {'int8'}
                if mmm(1) < -128 || ...
                    mmm(2) >= 127.5
                    fd = -128 + (255 / rng) * (fd - mmm(1));
                end
            case {'uint16'}
                if mmm(1) < 0 || ...
                    mmm(2) >= 65535.5
                    fd = (65535 / rng) * (fd - mmm(1));
                end
            case {'uint32'}
                if mmm(1) < 0 || ...
                    mmm(2) >= 4294967295.5
                    fd = (4294967295 / rng) * (fd - mmm(1));
                end
            case {'uint8'}
                if mmm(1) < 0 || ...
                    mmm(2) >= 255.5
                    fd = (255 / rng) * (fd - mmm(1));
                end
        end

        % round and cast
        eval(['fd=' dt '(round(fd));']);

    % for single
    else

        % just convert
        fd = single(fd);
    end
end

% re-permute if necessary
if di > 1
    fd = permute(fd, fperm');
end
