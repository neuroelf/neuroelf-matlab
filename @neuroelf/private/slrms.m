function [trms, crms, drms] = slrms(data, cond, run, opts)
%SLRMS  Compute RMS differences between trials and condition means.
%   TRMS = SLRMS(DATA, COND, RUN) computes the RMS (root-mean-squared)
%   differences (averaged across features/voxels) for each trial and the
%   condition means of "out-of-runs" trials (for multiple runs) or all
%   remaining trials (for a single run) as a Trial-by-Condition matrix.
%
%   [TRMS, CRMS] = SLRMS(DATA, COND, RUN) also returns the average value
%   for each trial-to-condition-mean pairing in a CxC matrix.
%
%   [TRMS, CRMS, DRMS] = SLRMS(DATA, COND, RUN) finally also returns a
%   single measure of the average across-over-within, within-only, and
%   across-only condition RMS measure (3x1 array).
%
%   TRMS = SLRMS(DATA, COND, RUN, OPTS) allows to set options as fields in
%   a 1x1 OPTS struct:
%
%        .remfmean  remove mean across features (global signal, def: true)
%        .scbytrms  scale RMS by estimate of total RMS (default: true)

% argument check
if nargin < 2 || (~isa(data, 'double') && ~isa(data, 'single')) || isempty(data) || ...
   ~isa(cond, 'double') || numel(cond) ~= size(data, 1) || any(isinf(cond(:)) | isnan(cond(:)))
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end

% make sure data is double
data = double(data);

% number of features and conditions
nf = size(data, 2);
cond = cond(:);
nt = numel(cond);

% run information
if nargin < 3 || isempty(run)
    run = 1;
elseif ~isa(run, 'double') || numel(run) ~= size(data, 1) || any(isinf(run(:)) | isnan(run(:)) | run(:) < 1)
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
else
    run = floor(run(:));
end

% recode to unique conditions (in same order)
uc = unique(cond);
nuc = numel(uc);
uci = zeros(1, max(uc));
uci(uc) = 1:nuc;
cond = uci(cond);
cond = cond(:);

% options
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'remfmean') || ~islogical(opts.remfmean) || numel(opts.remfmean) ~= 1
    opts.remfmean = true;
end
if ~isfield(opts, 'scbytrms') || ~islogical(opts.scbytrms) || numel(opts.scbytrms) ~= 1
    opts.scbytrms = true;
end

% apply mean-removal (across features)?
if opts.remfmean
    data = data - ((1 / nf) .* sum(data, 2)) * ones(1, nf);
end

% compute and scale by total RMS?
if opts.scbytrms
    trms = sqrt(sum((1 / (nt * nf)) .* sum(data .* data, 2)));
    data = (1 / trms) .* data;
end

% arrays for condition indices and number of condition indices per cond.
ci = zeros(nt, nuc);
nci = zeros(nuc, 1);

% for single run
if all(run == run(1))
    
    % repmat full condition means (mdata): trials-by-features-by-condition
    for cc = 1:nuc
        ci(:, cc) = double(cond == cc);
        nci(cc) = sum(ci(:, cc));
    end
    mdata = repmat(permute(diag(1 ./ nci) * (ci' * data), [3, 2, 1]), nt, 1);
    
    % for each trial/condition
    ndci = nci ./ (nci - 1);
    for cc = 1:nuc

        % replace condition mean with unbiased version (taking trial out)
        cci = (ci(:, cc) > 0);
        mdata(cci, :, cc) = ndci(cc) .* (mdata(cci, :, cc) - (1 / nci(cc)) .* data(cci, :));
    end

% multiple runs
else

    % recode to 1:unique-runs
    ur = unique(run);
    nur = numel(ur);
    uri = zeros(1, max(ur));
    uri(ur) = 1:nur;
    run = uri(run);
    run = run(:);

    % iterate over runs and conditions
    ri = zeros(nt, nuc, nur);
    nri = zeros(nuc, nur);
    for cc = 1:nuc
        ci(:, cc) = double(cond == cc);
        nci(cc) = sum(ci(:, cc));
        for rc = 1:nur
            ri(:, cc, rc) = double(cond == cc & run ~= rc);
            nri(cc, rc) = sum(ri(:, cc, rc));
        end
    end

    % compute means that are needed (conditions-by-runs)
    mdata = diag(1 ./ nri(:)) * reshape(ri, nt, nuc * nur)' * data;

    % for each trial, pick the correct means
    mdata = permute(reshape(mdata(ones(nt, 1) * (1:nuc) + (nuc .* ((run - 1) * ones(1, nuc))), :), nt, nuc, nf), [1, 3, 2]);
end

% compute full difference set
mdata = repmat(data, [1, 1, nuc]) - mdata;

% compute RMS per trial
trms = squeeze(sqrt((1 / nf) .* sum(mdata .* mdata, 2)));

% compute means
if nargout > 1

    % compute
    crms = ((1 ./ nci) * ones(1, nuc)) .* (ci' * trms);

    % and additional simple measure
    if nargout > 2
        erms = eye(nuc) > 0;
        nrms = nci * nci';
        arms = sum(sum((1 - erms) .* nrms .* crms)) ./ sum(sum((1 - erms) .* nrms));
        wrms = sum(sum(erms .* nrms .* crms)) ./ sum(sum(erms .* nrms));
        drms = [arms - wrms; wrms; arms];
    end
end

