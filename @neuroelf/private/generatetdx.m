function [X, Xh, Xu, Xs] = generatetdx(th, td, lupcols, lupcodes, cols, opts)
% generatetdx  - generate a design matrix X from tabular data with lookups
%
% FORMAT:       [X, Xh, Xu, Xs] = generatetdx(th, td, lupcols, lupcodes, cols [, opts])
%
% Input fields:
%
%       th          tabular data header (list of column names, cell array)
%       td          tabular data
%       lupcols     lookup column(s), for multiple use cell array
%       lupcodes    lookup codes (determines the size of the design matrix)
%       cols        columns to extract from the table
%       opts        optional settings
%        .splitcols columns which are used to generate a random effect (SU)
%        .tosplit   which columns to split (non-matching -> intercept)
%
% Output fields:
%
%       X           design matrix (intercept first column/s)
%       Xh          design matrix column names (headers)
%       Xu          which rows are usable (non NaN)
%       Xs          cell array with logical indices of split rows
%
% Note: This function can be used in tandem with modelcomp to determine
%       which (additional) regressors explain a significant amount of
%       variance in data.

% Version:  v0.9c
% Build:    13030716
% Date:     Feb-02 2013, 1:08 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2013, Jochen Weber
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
   ~iscell(th) || ...
    isempty(th) || ...
   ~all(cellfun(@ischar, th(:))) || ...
   ~isa(td, 'double') || ...
    isempty(td) || ...
    size(td, 2) ~= numel(th) || ...
    isempty(lupcols) || ...
   (~ischar(lupcols) && ...
    ~iscell(lupcols)) || ...
   ~isa(lupcodes, 'double') || ...
   (iscell(lupcols) && ...
    size(lupcodes, 2) ~= numel(lupcols)) || ...
   (ischar(lupcols) && ...
    size(lupcodes, 2) ~= 1) || ...
    isempty(cols) || ...
   (~ischar(cols) && ...
    ~iscell(cols))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
th = th(:);
if ~iscell(lupcols)
    lupcols = {lupcols(:)'};
else
    lupcols = lupcols(:);
    for cc = 1:numel(lupcols)
        lupcols{cc} = lupcols{cc}(:)';
    end
end
if numel(lupcols) > numel(th) || ...
   ~all(cellfun(@ischar, lupcols))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument LUPCOLS.' ...
    );
end
lupcoli = multimatch(lupcols, th);
if any(lupcoli == 0) || ...
    numel(unique(lupcoli)) ~= numel(lupcols)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid lookup column name (not found or duplicate).' ...
    );
end
if ~iscell(cols)
    cols = {cols(:)'};
else
    cols = cols(:);
    for cc = 1:numel(cols)
        cols{cc} = cols{cc}(:)';
    end
end
if numel(cols) > numel(th) || ...
   ~all(cellfun(@ischar, cols))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument COLS.' ...
    );
end
coli = multimatch(cols, th);
if any(coli == 0) || ...
    numel(unique(coli)) ~= numel(cols)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid column name (not found or duplicate).' ...
    );
end
if nargin < 6 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'splitcols') || ...
    isempty(opts.splitcols) || ...
   (~ischar(opts.splitcols) && ...
    ~iscell(opts.splitcols))
    opts.splitcols = {};
elseif ischar(opts.splitcols)
    opts.splitcols = {opts.splitcols(:)'};
else
    opts.splitcols = opts.splitcols(:);
    if ~all(cellfun(@ischar, opts.splitcols))
        error( ...
            'neuroelf:BadArgument', ...
            'Bad .splitcols field content.' ...
        );
    end
end
if ~isempty(opts.splitcols)
    splitcoli = multimatch(opts.splitcols, th);
    if any(splitcoli == 0) || ...
        numel(unique(splitcoli)) ~= numel(opts.splitcols)
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid split column name (not found or duplicate).' ...
        );
    end
else
    splitcoli = [];
end
if ~isfield(opts, 'tosplit') || ...
    isempty(opts.tosplit) || ...
   (~ischar(opts.tosplit) && ...
    ~iscell(opts.tosplit))
    opts.tosplit = {};
elseif ischar(opts.tosplit)
    opts.tosplit = {opts.tosplit(:)'};
else
    opts.tosplit = opts.tosplit(:);
    if ~all(cellfun(@ischar, opts.tosplit))
        error( ...
            'neuroelf:BadArgument', ...
            'Bad .tosplit field content.' ...
        );
    end
end
if ~isempty(opts.tosplit)
    tospliti = multimatch(opts.tosplit, cols);
    if numel(unique(tospliti)) ~= numel(opts.tosplit)
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid to-split column name (not found or duplicate).' ...
        );
    end
else
    tospliti = [];
end

% sizes
n = size(lupcodes, 1);
nc = numel(cols);
tnc = nc + 1;
on = ones(n, 1);
tn = size(td, 1);
nsp = numel(splitcoli);
nts = numel(tospliti);

% initialize to-use vector
Xu = true(n, 1);

% lookup complete
lup = zeros(n, 1);
for rc = 1:n
    lupi = findfirst(all(td(:, lupcoli) == (ones(tn, 1) * lupcodes(rc, :)), 2));
    if isempty(lupi)
        Xu(rc) = false;
    else
        lup(rc) = lupi;
    end
end

% no split columns
if nsp == 0
    splitdata = zeros(n, 1);
    usplitdata = 0;

% get data from split columns
else

    % get data to determine split
    splitdata = zeros(n, nsp);
    splitdata(Xu, :) = td(lup(Xu), splitcoli);
    usplitdata = unique(splitdata(Xu, :), 'rows');
end
nsplit = size(usplitdata, 1);
uspliti = zeros(n, 1);
for rc = 1:nsplit
    uspliti(all(splitdata == (on * usplitdata(rc, :)), 2)) = rc;
end

% generate design matrix (FFX)
X = NaN .* ones(n, (tnc - nts) + nsplit * nts);
sX = size(X, 2);
Xh = cell(sX, 1);

% FFX constant (intercept)
if nsplit == 1 || ...
   ~any(tospliti == 0)
    X(:, 1) = 1;
    Xh{1} = 'constant';
    Xi = 2;

% RFX (split) constant (intercept)
else

    % set intercepts parts to 0
    X(:, 1:nsplit) = 0;
    for rc = 1:nsplit
        X(uspliti == rc, rc) = 1;
        Xh{rc} = sprintf('constant_%03d', rc);
    end
    Xi = nsplit + 1;
end

% fill in FFX data
ffxcol = coli(setdiff(1:numel(coli), tospliti));
for cc = 1:numel(ffxcol)

    % header name
    Xh{Xi} = th{ffxcol(cc)};

    % data
    X(Xu, Xi) = td(lup(Xu), ffxcol(cc));

    % increase column counter
    Xi = Xi + 1;
end

% fill in RFX data
Xs = cell(1, nsplit);
rfxcol = coli(tospliti(tospliti > 0));
for cc = 1:numel(rfxcol)

    % iterate over data
    for rc = 1:nsplit
        Xh{Xi} = sprintf('%s_%03d', th{rfxcol(cc)}, rc);
        Xs{rc} = (Xu & uspliti == rc);
        X(:, Xi) = 0;
        X(Xs{rc}, Xi) = td(lup(Xs{rc}), rfxcol(cc));
        Xi = Xi + 1;
    end
end

% which rows to use
Xu = ~any(isnan(X), 2);
for rc = 1:nsplit
    Xs{rc} = (Xu & uspliti == rc);
end
