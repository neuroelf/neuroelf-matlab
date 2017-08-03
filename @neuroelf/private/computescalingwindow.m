function [win, mmh] = computescalingwindow(d, win, sci, scs)

% argument check
if nargin < 1 || ...
   ~isnumeric(d) || ...
    isempty(d)
    win = [0, 1];
    return;
end
if nargin < 2 || ...
   ~isa(win, 'double') || ...
    numel(win) ~= 2 || ...
    any(isnan(win)) || ...
    win(2) < win(1)
    win = [-Inf, Inf];
else
    win = win(:)';
end
if nargin < 3 || ...
   ~isa(sci, 'double') || ...
    numel(sci) ~= 1 || ...
    isinf(sci) || ...
    isnan(sci)
    sci = 0;
end
if nargin < 4 || ...
   ~isa(scs, 'double') || ...
    numel(scs) ~= 1 || ...
    isinf(scs) || ...
    isnan(scs) || ...
    scs == 0
    scs = 1;
end

% auto window computation
if any(isinf(win))

    % compute min/max over data
    mmm = minmaxmean(d);

% otherwise use window
else
    mmm = minmaxmean(win);
end

% compute difference
mmd = mmm(2) - mmm(1);

% if difference allows a histogram to be made
if mmd > (256 * eps)

    % use histcount
    mmh = histcount(d, mmm(1), mmm(2) - mmd / 512, mmd / 256);
    mmh = mmh(:);

% no histogram available, still variability
elseif mmd > 0

    % set to flat number vector
    mmh = round(numel(d) / 256) .* ones(256, 1);

% no histogram available AT ALL
else

    % set to number of elements and a vector of zeros
    mmh = [numel(d); zeros(255, 1)];
end

% in case of a datatype with too few elements
if isinteger(d) && ...
    mmd < 256

    % smooth zero values of histogram
    hz = 1 + find(mmh(2:end-1) == 0);
    mmh(hz) = 0.5 .* (mmh(hz - 1) + mmh(hz + 1));
end

% and make sure the first and last entries are no larger than the rest
mmhm = max(mmh(2:end-1));
mmh(1) = min(mmh(1), mmhm);
mmh(end) = min(mmh(end), mmhm);

% scaling for hdr
if sci ~= 0 || ...
    scs ~= 1
    if scs ~= 0
        mmm = sci + scs .* mmm(1:2);
    else
        mmm = sci + mmm(1:2);
    end
end

% replace Inf members
win(isinf(win) & win < 0) = mmm(1);
win(isinf(win)) = mmm(2);
if win(1) == win(2)
    win(2) = win(1) + sqrt(eps);
end
