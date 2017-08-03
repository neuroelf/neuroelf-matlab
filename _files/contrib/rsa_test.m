% clear everything
clear all;

% ask for filename
[filename, filepath] = uigetfile({'*.mat', 'Beta-correlations MAT-file (*.mat)'}, ...
    'Please select a MAT file containing Beta-correlations data...');
if isequal(filename, 0)
    error('No file selected.');
end
if isempty(filepath)
    filepath = pwd;
end
filename = [filepath filesep filename];

% load data file
disp(filename);
v = load(filename);
if ~isfield(v, 'bcorr') || ...
   ~isfield(v, 'betas') || ...
   ~isfield(v, 'conds') || ...
   ~isfield(v, 'study') || ...
   ~isfield(v, 'tvi') || ...
   ~isfield(v, 'voic')
    error('Not a valid Beta-correlations file.');
end

% reassign in current workspace
bcorr = v.bcorr;
betas = v.betas;
conds = v.conds;
study = v.study;
tvi   = v.tvi;
voic  = v.voic;

% remove v from workspace
clear v;

% determine which run each of the conditions is in
nconds = numel(conds);
condrun = zeros(nconds, 1);
for c1 = 1:nconds
    for c2 = 1:7
        if any(strcmp(study(c2).RunTimeVars.Predictors, conds{c1}))
            condrun(c1) = c2;
            break;
        end
    end
end

% get the different conditions
uconds = strrep(conds(~cellfun('isempty', regexp(conds, '_T001'))), '_T001', '');
nuconds = numel(uconds);
ucondi = cell(nuconds, 1);
for c1 = 1:nuconds
    ucondi{c1} = find(~cellfun('isempty', regexp(conds, [uconds{c1} '_T\d+'])));
end

% for each combination, find out what the correlation matrix is
ucondp = cell(nuconds, nuconds);
for c1 = 1:nuconds
    for c2 = c1:nuconds

        % start with the full grid
        [from, to] = ndgrid(ucondi{c1}(:)', ucondi{c2}(:)');
        from = from(:);
        to = to(:);

        % determine which of them to X out (all within runs)
        cfrom = condrun(from);
        samerun = (cfrom == condrun(to));
        if any(cfrom(samerun) ~= 7)
            from(samerun) = [];
            to(samerun) = [];
        end
        ucondp{c1, c2} = [from, to];
        if c1 ~= c2
            ucondp{c2, c1} = ucondp{c1, c2};
        end
    end
end

% generate average matrices
mcorr = zeros(nuconds, nuconds, size(bcorr, 3));
mdcorr = mcorr;
mcorrt = mcorr;
for sc = 1:size(bcorr, 3)

    % access source data
    for c1 = 1:nuconds
        for c2 = c1:nuconds
            ci = ucondp{c1, c2};
            ci(:, 3) = sc;
            bcorrval = indexarray(bcorr, ci);

            % remove bad values
            bcorrval(isnan(bcorrval) | isinf(bcorrval) | bcorrval == 0 | bcorrval == 1) = [];

            % average fisherized value
            bcorrval = fisherr2z(bcorrval);
            mcorr(c1, c2, sc) = mean(bcorrval);
            mdcorr(c1, c2, sc) = median(bcorrval);

            % N for t-test
            tn = harmmean(numel(unique(ci(:, 1))), numel(unique(ci(:, 2))));
            mcorrt(c1, c2, sc) = ...
                sqrt(tn) * mcorr(c1, c2, sc) / std(bcorrval(:));
            if c1 ~= c2
                mcorr(c2, c1, sc) = mcorr(c1, c2, sc);
                mdcorr(c2, c1, sc) = mdcorr(c1, c2, sc);
                mcorrt(c2, c1, sc) = mcorrt(c1, c2, sc);
            end
        end
    end
end

% create main difference t-tests
tcorr = zeros(nuconds, nuconds);
tcorrd = zeros(nuconds, nuconds);
tcorrt = zeros(nuconds, nuconds);
for c1 = 1:nuconds
    for c2 = (c1 + 1):nuconds
        tcorr(c1, c2) = sqrt(size(mcorr, 3)) * ...
            mean(squeeze(mcorr(c1, c1, :) - mcorr(c2, c2, :))) / ...
            std(squeeze(mcorr(c1, c1, :) - mcorr(c2, c2, :)));
        tcorrd(c1, c2) = sqrt(size(mdcorr, 3)) * ...
            mean(squeeze(mdcorr(c1, c1, :) - mdcorr(c2, c2, :))) / ...
            std(squeeze(mdcorr(c1, c1, :) - mdcorr(c2, c2, :)));
        tcorrt(c1, c2) = sqrt(size(mcorrt, 3)) * ...
            mean(squeeze(mcorrt(c1, c1, :) - mcorrt(c2, c2, :))) / ...
            std(squeeze(mcorrt(c1, c1, :) - mcorrt(c2, c2, :)));
        tcorr(c2, c1) = -tcorr(c1, c2);
        tcorrd(c2, c1) = -tcorrd(c1, c2);
        tcorrt(c2, c1) = -tcorrt(c1, c2);
    end
end

% save again
clear c1 c2 cfrom ci from samerun sc to;
save(filename);
