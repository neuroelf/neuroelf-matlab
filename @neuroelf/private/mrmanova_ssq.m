function [ssq, df1, df2, msq, yw] = mrmanova_ssq(y, ri, fi, rd, fd, yw, c, opts)
% mrmanova_ssq - compute SSQ for a given set of indices (levels)
%
% FORMAT:       [ssq, df1, df2, msq, yw] = mrmanova_ssq(y, ri, fi, rd, fd [, yw [, c [, opts]]])
%
% Input fields:
%
%       y           N-D double data
%       ri          Rx1 cell array with random-indices to group
%       fi          Fx1 cell array with fixed-indices to group
%       rd          random factor dimension (must be <= ndims(y))
%       fd          fixed factor dimension (must be <= ndims(y) & ~= rd)
%       yw          optional N-D double weights (must match in size with y)
%       c           optional SxC double covariates, with S == size(y, rd)
%       opts        1x1 struct with optional settings
%        .harmmean  compute harmonic mean instead of plain mean (false)
%        .robust    force robust estimates
%
% Output fields:
%
%       ssq         sum-of-squares
%       df1         d.f. for ssq (1x1 or same size as ssq, if yw is given)
%       df2         error d.f. for ssq
%       msq         mean sum-of-squares (optional)
%       yw          weights output (useful if robust estimation is used)
%
% Note: if yw is set to NaN, robust regression will be used to determine
%       yw from y according to the selection of indices

% Version:  v0.9c
% Build:    11060214
% Date:     Jun-01 2011, 11:03 PM EST
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
if nargin < 5 || ...
   ~isa(y, 'double') || ...
    isempty(y) || ...
   ~isa(ri, 'cell') || ...
    numel(ri) ~= max(size(ri)) || ...
   ~isa(fi, 'cell') || ...
    numel(fi) ~= max(size(fi)) || ...
   ~isa(rd, 'double') || ...
    numel(rd) ~= 1 || ...
    isinf(rd) || ...
    isnan(rd) || ...
   ~any(rd == 1:ndims(y)) || ...
   ~isa(fd, 'double') || ...
    numel(fd) ~= 1 || ...
    isinf(fd) || ...
    isnan(fd) || ...
   ~any(fd == 1:ndims(y)) || ...
    fd == rd || ...
    numel(ri) > size(y, rd)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
ysz = size(y);
nd = numel(ysz);
nri = ysz(rd);
nfi = ysz(fd);
for rc = 1:numel(ri)
    if ~isa(ri{rc}, 'double') || ...
        isempty(ri{rc}) || ...
        numel(ri{rc}) ~= max(size(ri{rc})) || ...
        any(isinf(ri{rc}) | isnan(ri{rc}) | ri{rc} < 1 | ri{rc} > nri | ri{rc} ~= fix(ri{rc})) || ...
        numel(ri{rc}) ~= numel(unique(ri{rc}))
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid indexing for random factor level %d.', ...
            rc ...
        );
    end
    ri{rc} = ri{rc}(:);
end
if numel(cat(1, ri{:})) ~= nri
    error( ...
        'neuroelf:BadArgument', ...
        'Incomplete sampling of random index dim.' ...
    );
end
if numel(ri) ~= nri
    for rc = 1:numel(ri)
        if numel(ri{rc}) < 2
            error( ...
                'neuroelf:BadArgument', ...
                'Random factor groupings require at least 2 units per level.' ...
            );
        end
    end
else
    c = [];
end
for fc = 1:numel(fi)
    if ~isa(fi{fc}, 'double') || ...
        isempty(fi{fc}) || ...
        numel(fi{fc}) ~= max(size(fi{fc})) || ...
        any(isinf(fi{fc}) | isnan(fi{fc}) | fi{fc} < 1 | fi{fc} > nfi | fi{fc} ~= fix(fi{fc})) || ...
        numel(fi{fc}) ~= numel(unique(fi{fc}))
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid indexing for random factor level %d.', ...
            fc ...
        );
    end
    fi{fc} = fi{fc}(:);
end
if numel(cat(1, fi{:})) ~= nfi
    error( ...
        'neuroelf:BadArgument', ...
        'Incomplete sampling of fixed index dim.' ...
    );
end
fw = false;
if nargin < 6 || ...
    isempty(yw)
    yw = [];
elseif ~isa(yw, 'double') || ...
   (numel(yw) ~= 1 && ...
    ~isequal(size(yw), size(y))) || ...
   (numel(yw) == 1 && ...
    ~isnan(yw))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad weights argument.' ...
    );
elseif numel(yw) > 1
    yw = limitrangec(yw, 0, 1, 0);
    y(yw == 0) = 0;
elseif numel(yw) == 1
    yw = ones(size(y));
    fw = true;
end
if nargin < 7 || ...
    isempty(c)
    c = [];
elseif ~isa(c, 'double') || ...
    ndims(c) > 2 || ...
    size(c, 1) ~= size(y, rd) || ...
    size(c, 2) >= (size(y, rd) - numel(ri) * numel(fi)) || ...
    any(isinf(c(:)) | isnan(c(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad covariates argument.' ...
    );
end
if nargin < 8 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'harmmean') || ...
   ~islogical(opts.harmmean) || ...
    numel(opts.harmmean) ~= 1
    hm = false;
else
    hm = opts.harmmean;
end
if ~isfield(opts, 'robust') || ...
   ~islogical(opts.robust) || ...
    numel(opts.robust) ~= 1
    fw = false;
else
    fw = fw || opts.robust;
end

% get sizes correct
ysr = repmat({':'}, nd);
str = ysr;
ysz([rd, fd]) = 1;
td = nd + 1;
ysz(td) = numel(ri) * numel(fi);

% create output
ssq = zeros(ysz);
if numel(yw) > 1
    df1 = zeros(ysz);
    df2 = zeros(ysz);
else
    df1 = 0;
    df2 = 0;
end

% iterate over random levels
tc = 1;
for rc = 1:numel(ri)

    % set random-level source indices
    ysr{rd} = ri{rc};

    % covariates for this part
    if ~isempty(c)
        cri = [ztrans(c(ri{rc}, :)), ones(numel(ri{rc}), 1)];
    else
        cri = [];
    end

    % iterate over fixed levels
    for fc = 1:numel(fi)

        % set fixed-level source indices
        ysr{fd} = fi{fc};

        % set target index
        str{td} = tc;
        tc = tc + 1;

        % weights need to be estimated
        if fw
        end

        % deal with covariate
        if ~isempty(cri)

            % robustly?
            if fw
            else
            end
        end

        % weights are available
        if numel(yw) == numel(y)

            % compute SSQ with weights
            ssq(str{:}) = (1 / (sum(sum(yw(ysr{:}), fd), rd) .^ 2)) .* ...
                (sum(sum(yw(ysr{:}) .* y(ysr{:}), fd), rd) .^ 2);

            % compute df
            df2f = sum(sum(yw(ysr{:}), fd), rd);
            df1(str{:}) = (max(1, size(cri, 2)) / (numel(ri{rc}) * numel(fi{fc}))) .* df2f;
            df2(str{:}) = df2f - df1(str{:});

        % no weights
        else

            % compute SSQ
            ssq(str{:}) = (1 / (numel(ri{rc}) * numel(fi{fc}))) .* ...
                (sum(sum(y(ysr{:}), fd), rd) .^ 2);

            % degrees of freedom
            df1 = df1 + max(1, size(cri, 2));
            df2 = df2 + numel(ri{rc}) * numel(fi{fc}) - max(1, size(cri, 2));
        end
    end

    % harmonic means?
    if hm && ...
        numel(ri{rc}) > 1
    end
end

% compute sum
ssq = sum(ssq, td);
if numel(df1) > 1
    df1 = sum(df1, td);
    df2 = sum(df2, td);
end

% compute msq?
if nargout > 3
    msq = (1 ./ df1) .* ssq;
end
