function c = slloco(d, m, w)

% input check
if nargin < 1 || ~isa(d, 'double') || isempty(d) || ndims(d) > 2
    error('neuroelf:general:badArgument', 'Bad or missing data argument.');
end
if nargin < 2 || ~ischar(m) || isempty(m) || ~any(strcmpi(m(:)', {'mean', 'median'}))
    m = 'a';
else
    m = lower(m(3));
end
mm = size(d, 1);
if nargin > 2 && (~isa(w, 'double') || size(w, 1) ~= mm || ndims(w) > 2)
    error('neuroelf:general:badArgument', 'Invalid W argument.');
elseif nargin < 3 || isempty(w)
    uw = false;
    w = 1;
else
    uw = true;
    sw = sqrt((sum(w) - 1) ./ (mm - 1));
    mn = 1:mm;
end

% compute correlations
c = zeros(1, size(w, 2));
for cc = 1:numel(c)
    
    % compute
    if uw
        wd = sparse(mn, mn, w(:, cc), mm, mm, mm);
        cm = (1 / sw(cc)) .* (d' * wd * d);
    else
        cm = d' * d;
    end

    % ensure good form
    cm(1:(nn+1):end) = NaN;
    cm(isinf(cm) | cm == 0) = NaN;

    % measure
    if m == 'a'
        cx = mean(cm(~isnan(cm)));
    else
        cx = median(cm(~isnan(cm)));
    end
    if ~isempty(cx) && ~isnan(cx)
        c(cc) = cx;
    end
end
