function [d, s] = slrsa(data, cond, run, opts)
%SLRSA  Perform search-light based RSA
%   [D, S] = SLRSA(SLDATA, COND) computes the (relative) greater similarity
%   within condition over cross-condition in output D(1) and the similarity
%   matrix S (CxC). D(2) is the average within condition similarity and
%   D(3) is the average across-condition similarity (and D(1) is D(2)-D(3)).
%
%   [D, S] = SLRSA(SLDATA, COND, RUN) only uses elements that do not
%   occur in the same run to lower the effect of global signal on trial
%   estimates.

% argument check
if nargin < 2 || (~isa(data, 'double') && ~isa(data, 'single')) || ~isa(cond, 'double') || ...
    numel(cond) ~= size(data, 1) || any(isinf(cond(:)) | isnan(cond(:)))
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
data = double(data);
cond = cond(:);
nt = numel(cond);

% no run information given
if nargin < 3 || ~isa(run, 'double')
    run = (1:nt)';
elseif numel(run) ~= nt || any(isinf(run(:)) | isnan(run(:)) | run(:) < 1)
    error('neuroelf:general:badArgument', 'Invalid run information.');
else
    run = run(:);
end

% get indices to compute difference from
cc1 = cond * ones(1, nt);
cc2 = ones(nt, 1) * cond';
ccm = (cc1 == cc2);

% run-comparison matrix (which trials to count)
rcm = ((run * ones(1, nt)) ~= (ones(nt, 1) * run'));

% which measure
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'comp') || ~ischar(opts.comp) || ~strcmpi(opts.comp(:)', 'cov')
    opts.comp = 'r';
else
    opts.comp = 'v';
end

% compute full correlation matrix
if opts.comp == 'r'
    cm = corrcoef(data');
else
    cm = cov(data');
end

% set invalid trials
cm(~rcm) = NaN;

% generate outputs
d = zeros(3, 1);
if opts.comp == 'r'
    sm = cm(ccm);
    d(2) = mean(fisherr2z(sm(~isnan(sm(:)))));
    sm = cm(~ccm);
    d(3) = mean(fisherr2z(sm(~isnan(sm(:)))));
    d(1) = d(2) - d(3);
else
    sm = cm(ccm);
    d(2) = mean(sm(~isnan(sm(:))));
    sm = cm(~ccm);
    d(3) = mean(sm(~isnan(sm(:))));
    d(1) = d(2) - d(3);
end

% no need to compute full condition-by-condition matrix?
if nargout < 2
    return;
end

% continue
uc = unique(cond);
nuc = numel(uc);
s = zeros(nuc, nuc);
for c1 = 1:nuc
    for c2 = c1:nuc
        sm = cm(cc1 == c1 & cc2 == c2);
        if opts.comp == 'r'
            s(c1, c2) = mean(fisherr2z(sm(~isnan(sm(:)))));
        else
            s(c1, c2) = mean(sm(~isnan(sm(:))));
        end
        if c1 ~= c2
            s(c2, c1) = s(c1, c2);
        end
    end
end
