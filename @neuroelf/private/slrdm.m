function [varargout] = slrdm(data, cond, run, opts)
%SLRDM  Compute 1-correlation dissimilarity matrix lower-triangle elements.
%   RDM = SLRDM(DATA) returns the lower triangle elements of the
%   dissimilarity matrix between the rows (trials) in DATA (assuming all
%   come from the same condition). RDM has a size of 1-by-T*(T-1)/2.
%
%   [RDM1, RDM2, ...] = SLRDM(DATA, COND) returns the RDMs per condition.
%   the sizes of the RDMs can be different (for unbalanced designs, etc.)
%
%   [RDM1, RDM2, ...] = SLRDM(DATA, COND, RUN) sets elements to NaN which
%   come from the same run.
%
%   [RDM1, RDM2, ...] = SLRDM(DATA, COND, RUN, OPTS) allows to set options
%   in the fields of a 1x1 struct OPTS:
%
%        .ctype     computation type, {'corr'}, 'rankcorr'

% todo: add 'kendtau' to ctype (code!)

% argument check
if nargin < 1 || (~isa(data, 'double') && ~isa(data, 'single')) || isempty(data)
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
data = double(data);
nt = size(data, 1);

if nargin < 2 || isempty(cond)
    cond = 1;
elseif ~isa(cond, 'double') || numel(cond) ~= nt || any(isinf(cond(:)) | isnan(cond(:)))
    error('neuroelf:general:badArgument', 'Bad COND argument.');
else
    cond = cond(:);
end

% no run information given
if nargin < 3 || ~isa(run, 'double')
    run = [];
elseif numel(run) ~= nt || any(isinf(run(:)) | isnan(run(:)) | run(:) < 1)
    error('neuroelf:general:badArgument', 'Invalid run information.');
else
    run = run(:);
end

% get indices to return
cc1 = (1:nt)' * ones(1, nt);
cc2 = ones(nt, 1) * (1:nt);
ccm = (cc1 > cc2);

% which measure
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'comp') || ~ischar(opts.comp) || isempty(opts.comp) || ...
   ~any('kr' == lower(opts.comp(1)))
    opts.comp = 'c';
else
    opts.comp = lower(opts.comp(1));
end

% compute full correlation matrix
if opts.comp == 'c'
    data = corrcoef(data');
elseif opts.comp == 'r'
    data = corrcoef(ranktrans(data, 2)');
else
    data = kendtaur(data');
end
cm = 1 - corrcoef(data');

% run-comparison matrix (which trials to NaN out)
if ~isempty(run)
    cm((run * ones(1, nt)) == (ones(nt, 1) * run')) = NaN;
end

% generate outputs
uc = unique(cond);
nuc = numel(uc);
varargout = cell(1, nuc);

% single condition
if nuc == 1
    varargout{1} = cm(ccm);
    return;
end

% iterate over conditions
for cc = 1:nuc
    varargout{cc} = cm(ccm & cond(cc1) == uc(cc) & cond(cc2) == uc(cc));
end
