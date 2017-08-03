function d = sldegdet(data, cond, opts)
%SLDEGDET  Compute degree of determination of Searchlight by central voxel

if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'slremfmean') || ~islogical(opts.slremfmean) || numel(opts.slremfmean) ~= 1
    slremfmean = true;
else
    slremfmean = opts.slremfmean;
end

% remove mean
sd = size(data);
if slremfmean
    data = data - mean(data, 2) * ones(1, sd(2));
end

% regress out condition means
uc = unique(cond);
nuc = numel(uc);
cm = zeros(sd(1), nuc);
for cc = 1:nuc
    cm(cond == uc(cc), cc) = 1;
end
data = data - cm * (((cm' * cm) \ cm') * data);

% compute variances
vd = var(data);

% regress out effect of central voxel and interactions
cm = data(:, ones(1, nuc)) .* cm;
data = data - cm * (((cm' * cm) \ cm') * data);

% new variances
nvd = var(data);

% fraction of variance
d = 1 - (sum(nvd(2:end)) / sum(vd(2:end)));
