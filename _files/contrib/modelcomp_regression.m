% script to transcode the 3 Excel (CSV) sheets of
% - onsets
% - subject-dyad/pair variables
% - single-subject variables
% into three (combined) variables (tabular datasets, td), as well as
% protocol files (for regular regression to create dyad-maps)

% switch to folder where data lives
cd /Volumes/sna/


% read onset data
t = acsvread('NetworkFaces16subjects.csv', ',', struct('asmatrix', true, 'convert', true, 'headline', ''));
nt1 = numel(t);

% get fieldnames
th1 = fieldnames(t);
pidcol = findfirst(strcmpi(th1, 'perceiverid'));
tidcol = findfirst(strcmpi(th1, 'targetid'));
pidtid = findfirst(strcmpi(th1, 'perceiveridtargetid'));
onscol = findfirst(strcmpi(th1, 'onsettime'));
if isempty(pidcol) || ...
    isempty(tidcol) || ...
    isempty(pidtid) || ...
    isempty(onscol)
    error( ...
        'SNA:BadFileContent', ...
        'Onset CSV file doesn''t contain required column.' ...
    );
end
nh1 = numel(th1);

% extend for run number information
runcol = nh1 + 1;
th1{runcol} = 'RunNumber';

% generate tabular data (with 0s)
td1 = zeros(nt1, nh1);

% for each header field
for c = 1:nh1

    % get data for that "column"
    tv = {t.(th1{c})};

    % and store in tabular data
    td1(:, c) = cat(1, tv{:});
end

% compute difference in onsets (to detect runs)
d = [1e6; diff(td1(:, onscol))];

% generate column that codes for runs
dr = zeros(size(td1, 1), 1);

% fill column at run switch with 1s
dr(d>20000|d<0) = 1;

% and generate a consecutive order list
dr = cumsum(dr);

% then store back into tabular data
td1(:, runcol) = dr;


% read dyad-rating file
t = acsvread('PerceiverTargetRelational16subjects.csv', ',', struct('asmatrix', true, 'convert', true, 'headline', ''));
nt2 = numel(t);

% get header names
th2 = fieldnames(t);
nh2 = numel(th2);

% generate tabular data (with 0s; for second table)
td2 = zeros(nt2, nh2);

% fill in with values
for c = 1:nh2
    tv = {t.(th2{c})};
    td2(:, c) = cat(1, tv{:});
end

% extend first table with header names (first three columns match!!) + data
exfrom = size(td1, 2) + 1;
exto = size(td1, 2) + nh2 - 3;
th1(exfrom:exto) = th2(4:end);
td1(:, exfrom:exto) = NaN;

% get unique pairings, but remove "ghost" faces (target ID ends in 00)
tvu = unique(td1(:, pidtid));
tvu(mod(tvu, 100) == 0) = [];

% for each unique dyad
for c=1:numel(tvu)

    % store into first table (additional columns) the data from table 2
    td1(td1(:, pidtid) == tvu(c), exfrom:exto) = ...
        repmat(td2(td2(:, pidtid) == tvu(c), 4:end), ...
        sum(td1(:, pidtid) == tvu(c)), 1);
end


% read single-subject data
t = acsvread('individualdifferences16subjects012213.csv', ',', struct('asmatrix', true, 'convert', true, 'headline', ''));
nt3 = numel(t);

% get header and data
th3 = fieldnames(t);
sidcol = findfirst(strcmpi(th3, 'subjectid'));
nh3 = numel(th3);
nnh3 = nh3 - 1;

% generate tabular data (with 0s; for second table)
td3 = zeros(nt3, nh3);

% fill in with values
for c = 1:nh3
    tv = {t.(th3{c})};
    td3(:, c) = cat(1, tv{:});
end

% extend first table with header names (first column is subject) + data
exfrom = size(td1, 2) + 1;
exto = size(td1, 2) + nnh3;
exfrom2 = size(td1, 2) + nnh3 + 1;
exto2 = size(td1, 2) + 2 * nnh3;
th1(exfrom:exto) = th3(2:end);
th1(exfrom2:exto2) = th3(2:end);
for c = exfrom:exto
    th1{c} = ['Perceiver' th1{c}];
end
for c = exfrom2:exto2
    th1{c} = ['Target' th1{c}];
end
td1(:, exfrom:exto2) = NaN;

% unique perceivers and targets
psub = unique(td1(:, pidcol));
tsub = unique(td1(:, tidcol));

% fill in data
for c = 1:numel(psub)
    if any(td3(:, sidcol) == psub(c))
        td1(td1(:, pidcol) == psub(c), exfrom:exto) = ...
            repmat(td3(td3(:, sidcol) == psub(c), 2:end), sum(td1(:, pidcol) == psub(c)), 1);
    end
end
for c = 1:numel(tsub)
    if any(td3(:, sidcol) == tsub(c))
        td1(td1(:, tidcol) == tsub(c), exfrom2:exto2) = ...
            repmat(td3(td3(:, sidcol) == tsub(c), 2:end), sum(td1(:, tidcol) == tsub(c)), 1);
    end
end

% extend second table as well
exfrom = size(td2, 2) + 1;
exto = size(td2, 2) + nnh3;
exfrom2 = size(td2, 2) + nnh3 + 1;
exto2 = size(td2, 2) + 2 * nnh3;
th2(exfrom:exto) = th3(2:end);
th2(exfrom2:exto2) = th3(2:end);
for c = exfrom:exto
    th2{c} = ['Perceiver' th2{c}];
end
for c = exfrom2:exto2
    th2{c} = ['Target' th2{c}];
end
td2(:, exfrom:exto2) = NaN;

% unique perceivers and targets
psub = unique(td2(:, pidcol));
tsub = unique(td2(:, tidcol));

% fill in data
for c = 1:numel(psub)
    if any(td3(:, sidcol) == psub(c))
        td2(td2(:, pidcol) == psub(c), exfrom:exto) = ...
            repmat(td3(td3(:, sidcol) == psub(c), 2:end), sum(td2(:, pidcol) == psub(c)), 1);
    end
end
for c = 1:numel(tsub)
    if any(td3(:, sidcol) == tsub(c))
        td2(td2(:, tidcol) == tsub(c), exfrom2:exto2) = ...
            repmat(td3(td3(:, sidcol) == tsub(c), 2:end), sum(td2(:, tidcol) == tsub(c)), 1);
    end
end


% save data
save SNA_NetworkFaces_info td1 td2 td3 th1 th2 th3


% generate protocols (for each run)
for c = 1:td1(end, runcol)

    % find rows pertaining to that run
    rr = find(td1(:, runcol) == c);

    % get data for that run
    tdr = td1(rr, :);

    % generate new PRT
    prt = xff('new:prt');

    % get target subject IDs
    tsub = unique(tdr(:, tidcol));

    % each of them is a condition!
    for cc = 1:numel(tsub)

        % get target-condition data
        tdrt = tdr(tdr(:, tidcol) == tsub(cc), :);

        % and compute onsets (relative to first onset in run!)
        tdro = tdrt(:, onscol) - tdr(1, onscol);

        % then add to the protocol
        prt.AddCond(sprintf('T%d', tsub(cc)), [tdro, tdro + 1000]);
    end

    % and save the PRT (with the correct naming)
    prt.SaveAs(sprintf('protocols/SNA%d_run%d.prt', tdr(1), 1 + mod(c+1, 2)));
    prt.ClearObject;
end

% generate MDM
mdm = xff('new:mdm');

% find VTC files
vtcs = findfiles([pwd '/subjects/*'], '*_MNI.vtc', 'depth=1');

% remove all but first two runs
vtcs = vtcs(~cellfun('isempty', regexp(vtcs, '_RUN0[12]')));

% get unique subject IDs (from VTCs!)
subids = vtcs;
for c = 1:numel(subids)
    [trash, subids{c}] = fileparts(subids{c});
end
subids = unique(regexprep(subids, '_.*$', ''));

% get number of subjects
numsubs = numel(subids);

% find protocol files
prts = findfiles([pwd '/protocols'], 'SNA*.prt');

% remove non-matching PRTs
prts(multimatch(prts, subids, true) == 0) = [];

% store in MDM
mdm.XTC_RTC = [vtcs, prts];

% find realignment parameter files
rps = findfiles([pwd '/subjects/*/func*/*networkfa*'], 'rp*.txt', 'depth=1');

% store in MDM RunTimeVars
mdm.RunTimeVars.MotionParameters = rps;
mdm.RunTimeVars.AutoSave = true;

% save MDM
mdm.SaveAs(sprintf('SNA_NetworkFaces_%dsubjects.mdm', numsubs));

% compute GLM (for those type of maps)
% - using motion parameters
% - setting the "Rest" condition (BV-lingo!) to '' (no rest defined in PRT)
% - temporal filtering with 80s HPF cut-off, using DCT filter regressors
glm = mdm.ComputeGLM(struct( ...
    'motpars',   true, ...
    'restcond',  '', ...
    'tfilter',   80, ...
    'tfilttype', 'dct'));

% save GLM
glm.SaveAs(sprintf('SNA_NetworkFaces_%dsubjects_OLS.glm', numsubs));

% browse
glm.Browse;

% get the list of subjects (Perceivers) and predictors (Targets)
s = glm.Subjects;
ns = numel(s);
sp = glm.SubjectPredictors;

% remove constant predictor
sp(end) = [];
nsp = numel(sp);

% create a cell array to contain these
mn = cell(ns * nsp,1);

% create list of map names that is "PerceiverIDTargetID" compatible
for sc = 1:ns
    for pc = 1:nsp
        mn{(sc - 1) * 27 + pc} = [s{sc}(4:6) sp{pc}(2:4)];
    end
end

% and a numeric version of that
mnn = str2double(mn);

% as well as an extract of maps
maps = single(zeros(66,52,56));
maps(1, 1, 1, ns * nsp) = 0;
for sc = 1:ns
    for pc = 1:27
        maps(:, :, : ,(sc - 1) * 27 + pc) = ...
            glm.GLMData.Subject(sc).BetaMaps(:, :, :, pc);
    end
end

% remove maps that do not match in group for Perceiver/Target dyads
kmaps = true(numel(mn), 1);
kmaps(floor(mnn / 100000) ~= floor(mod(mnn, 1000) / 100)) = false;
mn = mn(kmaps);
mnn = mnn(kmaps);
maps = maps(:, :, :, kmaps);

% save data for later use (regressions, etc.)
save SNA_NetworkFaces_DyadTrialData maps mn mnn

% find column numbers for variables based on names
pef = findfirst(strcmpi(th2, 'Perceiverineigenvectorfriendship'));
tef = findfirst(strcmpi(th2, 'Targetineigenvectorfriendship'));

% get list of "perceivers" and "targets" from available map names
perceiver = floor(mnn/1000);
uperceiver = unique(perceiver);
target = mod(mnn, 1000);
utarget = unique(target);

% create design matrix (FFX-intercept, one column)
X = ones(numel(mnn), 1);

% set second column to NaNs (in case we don't have that value for a dyad)
X(:, 2) = nan;

% for each dyad (numeric map name)
for c = 1:numel(mnn)

    % if that value is present in the second data table (td2)
    if any(td2(:, pidtid) == mnn(c))

        % set this into the second column of the "design matrix"
        X(c,2) = td2(td2(:, pidtid) == mnn(c), tef);
    end
end

% figure out which maps to use
Xn = ~any(isnan(X), 2);

% and possibly remove "self" faces
Xn = Xn & (mod(mnn, 1001) ~= 0);

% compare the two models (FFX intercept only with FFX intercept + regr.)
[f, df1, df2] = modelcomp(X(Xn, :) , X(Xn, 1), maps(:, :, :, Xn), 4);

% remove "missing values"
f(isnan(f)) = 0;

% estimate smoothness (via residual)
[b, irtc, ptc] = calcbetas(X(Xn, :), maps(:, :, :, Xn), 4);
res = maps(:, :, :, Xn) - ptc;
fwhm = resestsmooth(res, [3, 3, 3]);

% generate VMP and store results
vmp = glm.RFX_tMap([1, zeros(1, nsp)]);
vmp.Map.Type = 4; % F-test
vmp.Map.DF1 = df1;
vmp.Map.DF2 = df2;
vmp.Map.VMPData = single(f);
vmp.Map.Name = 'Target in eigenvector friendship (FFX-correlation, F-test)';
vmp.Map.RunTimeVars.FWHMResEst = fwhm;
vmp.Browse;

% generete a design matrix with RFX intercept for perceivers
X = zeros(numel(mnn), numel(uperceiver));
X(:, end+1) = NaN;
for c = 1:numel(uperceiver)
    X(perceiver == uperceiver(c), c) = 1;
end

% then add, as last regressor, FFX regressor (tef column again)
for c = 1:numel(mnn)
    if any(td2(:, pidtid) == mnn(c))
        X(c,end) = td2(td2(:, pidtid) == mnn(c), tef);
    end
end
Xn = ~any(isnan(X), 2);
Xn = Xn & (mod(mnn, 1001) ~= 0);

% compare RFX intercept only with FFX regressor added
[f, df1, df2] = modelcomp(X(Xn, :) , X(Xn, 1:end-1), maps(:, :, :, Xn), 4);
f(isnan(f)) = 0;
[b, irtc, ptc] = calcbetas(X(Xn, :), maps(:, :, :, Xn), 4);
res = maps(:, :, :, Xn) - ptc;
fwhm = resestsmooth(res, [3, 3, 3]);

% add to VMP
vmp.Map(end+1) = vmp.Map(end);
vmp.Map(end).VMPData = single(f);
vmp.Map(end).VMPDataCT = [];
vmp.Map(end).DF1 = df1;
vmp.Map(end).DF2 = df2;
vmp.Map(end).Name = 'Target in eigenvector friendship (FFX-corr, RFX-intercept, F-test)';
vmp.Map(end).RunTimeVars.FWHMResEst = fwhm;
vmp.Browse;

% create model with several regressors
selected = find(multimatch(th2, {'HowCloseAreYouTo'; 'Perceiverineigenvectorfriendship'; 'Targetineigenvectorfriendship'}) > 0);
X = ones(numel(mnn), 1);
X(:,2:1+numel(selected)) = NaN;
for c=1:numel(mnn)
    if any(td2(:, pidtid) == mnn(c))
        X(c,2:end) = td2(td2(:, pidtid) == mnn(c), selected);
    end
end
Xn = ~any(isnan(X), 2);
Xn = Xn & (mod(mnn, 1001) ~= 0);

% test 2nd variable (on top of the rest; full model vs. reduced model)
[f, df1, df2] = modelcomp(X(Xn, :) , X(Xn, [1, 3, 4]), maps(:, :, :, Xn), 4);
f(isnan(f)) = 0;
[b, irtc, ptc] = calcbetas(X(Xn, :), maps(:, :, :, Xn), 4);
res = maps(:, :, :, Xn) - ptc;
fwhm = resestsmooth(res, [3,3,3]);

% store in VMP
vmp.Map(end+1) = vmp.Map(end);
vmp.Map(end).VMPData =single(f);
vmp.Map(end).DF1 = df1;
vmp.Map(end).DF2 = df2;
vmp.Map(end).Name = 'How close are you to (above and beyond P/T in eigenvector friendship) (FFX-corr, F-test)';
vmp.Browse;
