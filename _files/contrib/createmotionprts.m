function t = createmotionprts(file, file2)
% createmotionprts  - create PRTs from EPrime output fields for MOTION

% conditions and colors
conds = { ...
    'InstLook', ...
    'InstNeg', ...
    'FoodLook', ...
    'FoodNeg', ...
    'AlcLook', ...
    'AlcNeg', ...
    'FoodRating', ...
    'AlcRating'};
condcols = [ ...
     96,  96, 224; ...
     32,  32, 160; ...
     64, 255,  64; ...
      0, 128,   0; ...
    255,  64,  64; ...
    128,   0,   0; ...
     96, 224,  96; ...
    224,  96,  96];

% read the first file
t = acsvread(file, char(9), struct('convert', 1, 'readstart', 'BehavioralData', 'headline', ''));
if nargin > 1
    t = t(1:(20 * floor(numel(t) / 20)));
    t2 = acsvread(file2, char(9), struct('convert', 1, 'readstart', 'BehavioralData', 'headline', ''));
    t = joinstructs(t, t2);
end

% replace ExperimentName with Version field
x = str2double(strrep(strrep({t.ExperimentName}, 'ROC_Version', ''), '_scanner', ''));
t(1).Version = 1;
for c = 1:numel(x)
    t(c).Version = x(c);
end
t = rmfield(t, 'ExperimentName');

% ensure the same subject number is used throughout
x = repmat({t(1).Subject}, 1, numel(t));
[t.Subject] = deal(x{:});

% and only keep the letter
x = strrep({t.Procedure_Block_}, 'RunProc', '');
[t.Procedure_Block_] = deal(x{:});

% remove un-used fields
t = rmfield(t, 'Clock_Information');
t = rmfield(t, 'Display_RefreshRate');
t = rmfield(t, 'Group');
t = rmfield(t, 'Image_DurationError');
t = rmfield(t, 'Image_OffsetDelay');
t = rmfield(t, 'Instruct_DurationError');
t = rmfield(t, 'Instruct_OffsetDelay');
t = rmfield(t, 'ISI_DurationError');
t = rmfield(t, 'ISI_OffsetDelay');
t = rmfield(t, 'ITI_DurationError');
t = rmfield(t, 'ITI_OffsetDelay');
t = rmfield(t, 'Procedure_Trial_');
t = rmfield(t, 'RandomSeed');
t = rmfield(t, 'Rating_DurationError');
t = rmfield(t, 'Rating_OffsetDelay');
if isfield(t, 'RunA')
    t = rmfield(t, 'RunA');
    t = rmfield(t, 'RunA_Cycle');
    t = rmfield(t, 'RunA_Sample');
end
if isfield(t, 'RunB')
    t = rmfield(t, 'RunB');
    t = rmfield(t, 'RunB_Cycle');
    t = rmfield(t, 'RunB_Sample');
end
if isfield(t, 'RunC')
    t = rmfield(t, 'RunC');
    t = rmfield(t, 'RunC_Cycle');
    t = rmfield(t, 'RunC_Sample');
end
if isfield(t, 'RunD')
    t = rmfield(t, 'RunD');
    t = rmfield(t, 'RunD_Cycle');
    t = rmfield(t, 'RunD_Sample');
end
t = rmfield(t, 'RunList_Cycle');
t = rmfield(t, 'RunList_Sample');
t = rmfield(t, 'Running_Block_');
t = rmfield(t, 'Running_Trial_');
t = rmfield(t, 'SessionDate');
t = rmfield(t, 'SessionTime');
t = rmfield(t, 'SessionTimeUtc');
t = rmfield(t, 'SynchWithScanner_OnsetTime');
if isfield(t, 'TrialType')
    t = rmfield(t, 'TrialType');
end
if isfield(t, 'WordImageSequence')
    t = rmfield(t, 'WordImageSequence');
end

% reconvert the table into timing information
tsync = cat(1, t.SynchWithScanner_OffsetTime);
tons = [ ...
    cat(1, t.Instruct_OnsetTime), ...
    cat(1, t.Image_OnsetTime), ...
    cat(1, t.Rating_OnsetTime)] - repmat(tsync, 1, 3);
tdur = ones(numel(tsync), 1) * [2000, 6000, 3000];
ratp = cat(1, t.Rating_RT);

% subject and run ID
subid = t(1).Subject;
while subid > 999
    subid = floor(subid / 10);
end
runid = {t.Procedure_Block_}';

% create PRTs (trials should be in multiples of 20)
% - find unique letters
runs = uunion(runid, {});
for rc = 1:numel(runs)

    % new PRT
    prt = xff('new:prt');

    % find rows for this run
    rr = strcmp(runid, runs{rc});

    % get sub-data
    st = t(rr);
    stons = tons(rr, :);
    stdur = tdur(rr, :);
    sratp = ratp(rr);

    % get look, negative, food, alc rows
    look = strcmpi({st.Word}, 'look');
    nega = strcmpi({st.Word}, 'negative');
    food = strcmpi({st.Type}, 'food');
    alco = strcmpi({st.Type}, 'alc');

    % add instruction conditions
    prt.AddCond(conds{1}, [stons(look, 1), stons(look, 1) + stdur(look, 1)], condcols(1, :));
    prt.AddCond(conds{2}, [stons(nega, 1), stons(nega, 1) + stdur(nega, 1)], condcols(2, :));

    % add image conditions (2x2)
    prt.AddCond(conds{3}, [stons(food & look, 2), stons(food & look, 2) + stdur(food & look, 2)], condcols(3, :));
    prt.AddCond(conds{4}, [stons(food & nega, 2), stons(food & nega, 2) + stdur(food & nega, 2)], condcols(4, :));
    prt.AddCond(conds{5}, [stons(alco & look, 2), stons(alco & look, 2) + stdur(alco & look, 2)], condcols(5, :));
    prt.AddCond(conds{6}, [stons(alco & nega, 2), stons(alco & nega, 2) + stdur(alco & nega, 2)], condcols(6, :));

    % add rating conditions
    prt.AddCond(conds{7}, [stons(food, 3), stons(food, 3) + stdur(food, 3)], condcols(7, :));
    prt.AddCond(conds{8}, [stons(alco, 3), stons(alco, 3) + stdur(alco, 3)], condcols(8, :));

    % add RT parameter
    prt.Cond(7).Weights = ztrans(sratp(food));
    prt.Cond(8).Weights = ztrans(sratp(alco));

    % save as
    prt.SaveAs(sprintf('mot%3d_run%d.prt', subid, rc));

    % clear object
    prt.ClearObject;
end
