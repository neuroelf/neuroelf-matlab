function data = repbadvs(data, cutoff, idx, ovs)
%REPBADVS  Replace bad values in data along last dimension.
%   DATA = REPBADVS(DATA, CUTOFF) replaced the bad values (NaN by default)
%   in DATA with estimates from a FFT-like spectrum analysis, using
%   CUTOFF as the cut-off frequency, given in elements of data-resolution,
%   such that for a 1Hz signal a cutoff of .2Hz would be given as 5.
%
%   DATA = REPBADVS(DATA, CUTOFF, IDX) replaces data at points specified in
%   IDX, where MAX(IDX) must be smaller than SIZE(DATA, NDIMS(DATA)).
%

% data?
if nargin < 2 || (~isa(data, 'double') && ~isa(data, 'single')) || numel(cutoff) ~= 1
    error('MATLAB:general:badArgument', 'Bad or missing argument.');
end

% get data size
ds = size(data);
nv = prod(ds(1:end-1));
nt = ds(end);

% compute pairs of frequencies
np = round(nt / cutoff);
if (2 * np) > (nt - 3)
    error('MATLAB:repbadvs:invalidValue', 'Cut-off value too low.');
end

% reshape
data = reshape(data, nv, nt);

% oversampling
if nargin == 3 && numel(idx) == 1 && ~isinf(idx) && ~isnan(idx) && any(idx == (2:12))
    ovs = idx;
    idx =[];
end
if exist('ovs', 'var') == 0 || ~isa(ovs, 'double') || numel(ovs) ~= 1 || isinf(ovs) || isnan(ovs) || ovs < 2
    ovs = 1;
else
    ovs = min(12, round(ovs));
end

% bad data points
if nargin < 3 || ~isa(idx, 'double') || isempty(idx)
    dk = ceil(0.5 * ds(1:end-1));
    dk = 1 + sum([1, cumprod(ds(1:end-2))] .* (dk - 1));
    idx = find(isnan(data(dk, :)));
end

% nothing to do?
if isempty(idx)
    return;
end

% too much data?
ndata = numel(data);
ncrit = 8 * (2 ^ 30) / (8 * 3 * (ovs + ovs * (ovs > 1)));
if ndata > ncrit

    % loop
    dpc = ceil(ndata / ceil(nt * ndata / ncrit));
    for cdata = 1:dpc:nv

        % replace part
        tdata = min(nv, cdata + dpc - 1);
        data(cdata:tdata, :) = repbadvs(data(cdata:tdata, :), cutoff, idx, ovs);
    end

    % reshape and return
    data = reshape(data, ds);
    return;
end

% over sampling
if ovs > 1
    odata = flexinterpn_method(odata);
end

% good indices
gidx = setdiff((1:nt)', idx(:));

% create matrix required for filtering
w = zeros(1 + 2 * np, nt);
w(1, :) = 1;
for pc = 1:np
    w(2*pc, :) = cos((2 * pi * pc / nt) .* (0:(nt-1)));
    w(2*pc+1, :) = sin((2 * pi * pc / nt) .* (0:(nt-1)));
end
wres = w(:, gidx);
if (size(wres, 2) - 8) < size(wres, 1)
    error('MATLAB:repbadvs:tooFewDataPoints', 'Not enough good data points.');
end
wrep = w(:, idx);

% compute estimates for good data
wb = (data(:, gidx) * wres') / (wres * wres');

% replace in data
data(:, idx) = wb * wrep;

% reshape back
data = reshape(data, ds);
