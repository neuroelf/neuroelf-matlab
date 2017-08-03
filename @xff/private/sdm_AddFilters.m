function xo = sdm_AddFilters(xo, opts)
% SDM::AddFilters  - add filtering regressors to SDM
%
% FORMAT:       [sdm =] sdm.AddFilters(opts)
%
% Input fields:
%
%       opts        settings
%        .constant  also add constant regressor (default: false)
%        .ftype     either of 'dct', {'fourier'}, 'linear', 'none', 'poly'
%        .nuisreg   nuisance regressors (VxR matrix, default: [])
%        .number    number of filters/frequencies (cycles, 3)
%        .orthpoly  orthogonalize polynomial filters (potential speedup)
%        .range     1x2 (first and last time point, default: [1, T])
%        .sess      session number (for names, empty if not given)
%        .timepts   number of time points (valid only for empty SDM)
%
% Output fields:
%
%       sdm         altered SDM with added regressors
%
% Using: ztrans.

% Version:  v1.1
% Build:    16013123
% Date:     Jan-31 2016, 11:07 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'sdm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get design matrix right
bc = xo.C;
sdm = bc.SDMMatrix;

% options check
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'constant') || ~islogical(opts.constant) || numel(opts.constant) ~= 1
    opts.constant = false;
end
if ~isfield(opts, 'ftype') || ~ischar(opts.ftype) || ...
   ~any(strcmpi(opts.ftype(:)', {'dct', 'fourier', 'linear', 'none', 'poly'}))
    opts.ftype = 'fourier';
else
    opts.ftype = lower(opts.ftype(:)');
end
if ~isfield(opts, 'nuisreg') || ~isa(opts.nuisreg, 'double') || ndims(opts.nuisreg) ~= 2 || ...
    size(opts.nuisreg, 1) ~= size(sdm, 1) || any(isinf(opts.nuisreg(:)) | isnan(opts.nuisreg(:)))
    opts.nuisreg = [];
end
if ~isfield(opts, 'orthpoly') || ~islogical(opts.orthpoly) || numel(opts.orthpoly) ~= 1
    opts.orthpoly = false;
end
if ~isfield(opts, 'range') || ~isa(opts.range, 'double') || numel(opts.range) ~= 2 || ...
    any(isinf(opts.range) | isnan(opts.range) | opts.range < 1) || opts.range(2) < opts.range(1)
    opts.range = [];
else
    opts.range = fix(opts.range(:)');
end
if ~isfield(opts, 'sess') || ~isa(opts.sess, 'double') || numel(opts.sess) ~= 1 || ...
    isinf(opts.sess) || isnan(opts.sess) || opts.sess < 1 || opts.sess ~= fix(opts.sess)
    opts.sess = '';
else
    opts.sess = sprintf('S%d-', opts.sess);
end
if ~isempty(sdm) || ~isfield(opts, 'timepts') || ~isa(opts.timepts, 'double') || ...
    numel(opts.timepts) ~= 1 || isinf(opts.timepts) || isnan(opts.timepts) || opts.timepts < 4
    if isempty(sdm)
        error('neuroelf:xff:badArgument', 'Empty SDM: Number of timepoints required.');
    end
    opts.timepts = size(sdm, 1);
else
    opts.timepts = fix(opts.timepts);
    sdm = zeros(opts.timepts, 0);
    bc.PredictorNames(:) = [];
    bc.PredictorColors = zeros(0, 3);
end
bc.NrOfDataPoints = opts.timepts;
if isempty(opts.nuisreg)
    opts.nuisreg = zeros(opts.timepts, 0);
end
if isempty(opts.range)
    rg = [1, bc.NrOfDataPoints];
else
    rg = [opts.range(1), min(opts.range(2), bc.NrOfDataPoints)];
end
if ~isfield(opts, 'number') || ~isa(opts.number, 'double') || numel(opts.number) ~= length(opts.number) || ...
    any(isinf(opts.number) | isnan(opts.number) | opts.number < 0 | ...
    opts.number > floor(0.5 * (opts.timepts - size(sdm, 2))))
    opts.number = 3;
end
if numel(opts.number) == 1
    opts.number = 1:opts.number;
else
    opts.number = unique(opts.number);
end
opts.number = opts.number(:)';
nv = 1 + rg(2) - rg(1);

% temporal filter type
switch (opts.ftype)

    % DCT-based filtering
    case 'dct'

        % build X
        X = zeros(nv, numel(opts.number));
        Xn = cell(1, numel(opts.number));
        Xc = zeros(numel(opts.number), 3);
        n = 0:(nv - 1);
        for dc = 1:numel(opts.number)
            X(:, dc) = cos(pi * (2 * n + 1) * opts.number(dc) / (2 * nv));
            Xn{dc} = sprintf('%sDCF-f%d', opts.sess, opts.number(dc));
            Xc(dc, :) = max(0, [240, 224, 224] - 3 * opts.number(dc));
        end

    % sin/cos set filtering
    case 'fourier'

        % build X
        X = zeros(nv, 2 * numel(opts.number));
        Xn = cell(1, 2 * numel(opts.number));
        Xc = zeros(2 * numel(opts.number), 3);
        n = 0:(nv - 1);
        for pc = 1:numel(opts.number)
            X(:, 2 * pc - 1) = sin(pi * (2 * n + 1) * opts.number(pc) / nv);
            X(:, 2 * pc) = cos(pi * (2 * n + 1) * opts.number(pc) / nv);
            Xn{2 * pc - 1} = sprintf('%ssin-f%d', opts.sess, opts.number(pc));
            Xn{2 * pc} = sprintf('%scos-f%d', opts.sess, opts.number(pc));
            Xc(2 * pc - 1, :) = max(0, [224, 240, 224] - 2 * opts.number(pc));
            Xc(2 * pc, :) = max(0, [224, 224, 240] - 2 * opts.number(pc));
        end

    % only detrend
    case 'linear'

        % build X
        X = ne_methods.ztrans(1:nv);
        Xn = {[opts.sess 'linear_trend']};
        Xc = [224, 224, 224];

    % polynomial filtering
    case 'poly'

        % build X
        X = [ones(nv, 1), zeros(nv, max(opts.number))];
        Xn = cell(1, nv);
        if nv > 1
            Xc = round(repmat((224:(-192/(nv-1)):32)', 1, 3));
        else
            Xc = [128, 128, 128];
        end
        n = (1 / (nv - 1)) .* ((-nv + 1):2:(nv - 1))';
        n = 0.5 .* (abs(n) + abs(n(end:-1:1)));
        n(1:(floor(0.5*nv))) = -n(1:(floor(0.5*nv)));
        if mod(nv, 2) ~= 0
            n(ceil(0.5*nv)) = 0;
        end
        X(:, 2) = n;
        Xn{1} = 'polydeg001';
        for dc = 2:max(opts.number)
            Xn{dc} = sprintf('polydeg%03d', dc);
            X(:, dc + 1) = (1 / (dc + 1)) .* ...
                ((2 * dc + 1) .* n .* X(:, dc) - dc .* X(:, dc - 1));
        end

        % enforce orthogonality
        if opts.orthpoly
            for dc = max(opts.number):-1:2
                X(:, dc + 1) = X(:, dc + 1) - X(:, 1:dc) * ...
                    ((X(:, 1:dc)' * X(:, 1:dc)) \ (X(:, 1:dc)' * X(:, dc + 1)));
            end
        end

        % re-scale and cut away constant
        X = X(:, 2:end) ./ (ones(nv, 1) * max(abs(X(:, 2:end)), [], 1));

        % sub-select
        X = X(:, opts.number);
        Xn = Xn(opts.number);
        Xc = Xc(opts.number, :);

    % no filters
    otherwise
        X = zeros(nv, 0);
        Xn = {};
        Xc = zeros(0, 3);
end

% add nuisance regressors
if ~isempty(opts.nuisreg)
    nns = size(opts.nuisreg, 2);
    Xns = numel(Xn);
    X(:, end+1:end+nns) = opts.nuisreg;
    Xc = [Xc; floor(255.999 .* rand(nns, 3))];
    Xn = [Xn, cell(1, nns)];
    for dc = 1:nns
        Xn{Xns+dc} = sprintf('nuisreg%03d', dc);
    end
end

% add constant too?
if opts.constant && (isempty(bc.SDMMatrix) || ~any(all(bc.SDMMatrix == 1)))
    X(:, end + 1) = 1;
    Xn{end + 1} = 'Constant';
    Xc(end + 1, :) = [255, 255, 255];
end

% add to content
bc.FirstConfoundPredictor = min(bc.FirstConfoundPredictor, size(sdm, 2) + 1);
if size(X, 1) == size(sdm, 1)
    bc.SDMMatrix = [sdm, X];
else
    bc.SDMMatrix = [sdm, [zeros(rg(1)-1, size(X, 2)); X; zeros(size(X, 1)-rg(2), size(X, 2))]];
end
bc.NrOfPredictors = size(bc.SDMMatrix, 2);
bc.PredictorNames = [bc.PredictorNames(:)', Xn(:)'];
bc.PredictorColors = [bc.PredictorColors; Xc];
bc.RTCMatrix = bc.SDMMatrix;
cp = find(any(bc.SDMMatrix ~= 0, 1) & all(diff(bc.SDMMatrix) == 0, 1));
if ~isempty(cp)
    bc.RTCMatrix(:, cp) = [];
    bc.IncludesConstant = 1;
    if numel(cp) > 1
        cp(end) = [];
        bc.SDMMatrix(:, cp) = [];
        bc.PredictorNames(cp) = [];
        bc.PredictorColors(cp, :) = [];
        cp = find(any(bc.SDMMatrix ~= 0, 1) & all(diff(bc.SDMMatrix) == 0, 1));
    end
    if cp ~= size(bc.SDMMatrix, 2)
        bc.SDMMatrix(:, cp) = [];
        bc.SDMMatrix(:, end+1) = 1;
        bc.PredictorNames(cp) = [];
        bc.PredictorNames(end+1) = {'Constant'};
        bc.PredictorColors(cp, :) = [];
        bc.PredictorColors(end+1, :) = 255;
    end
else
    bc.IncludesConstant = 0;
end

% set back
xo.C = bc;
