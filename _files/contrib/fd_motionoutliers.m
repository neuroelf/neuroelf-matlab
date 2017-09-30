% this script locates the PRTs, VTCs, and motion parameter files
% for the two longitudinal SOBC tasks, and then removes from those lists
% the files for which motion exceeds 10% of volumes (1.5mm displacement)

% basefolder
basefolder = '/Volumes/BeckPort/SOBC/longFunctional';

% NeuroElf library
n = neuroelf;

% use findfiles to identify prts, vtcs, and motion parameters; for both tasks
app_prts = n.findfiles([basefolder '/protocols/Appetitive/????*/Phase*'], 'P*.prt', '-d1');
app_vtcs = n.findfiles([basefolder '/preprocessed/????*/Phase*'],'*APP*.vtc', '-d1');
app_rps  = n.findfiles([basefolder '/preprocessed/????*/Phase*/func/run*APP'],'rp*.txt', '-d1');
ave_prts = n.findfiles([basefolder '/protocols/Aversive/????*/Phase*'], 'P*.prt', '-d1');
ave_vtcs = n.findfiles([basefolder '/preprocessed/????*/Phase*'],'*AVE*.vtc', '-d1');
ave_rps  = n.findfiles([basefolder '/preprocessed/????*/Phase*/func/run*AVER'],'rp*.txt', '-d1');
app_bad_vols = app_rps; 
ave_bad_vols = ave_rps;

% %% TEMPORARY FIX UNTIL PRTS AND VTCS MATCH!!
app_prts = app_vtcs;
ave_prts = ave_vtcs;
% %% TAKE OUT AFTER RECONCILING PRTS AND VTCS!!

% load rps and identify bad volumes
app_thresh = zeros(numel(app_rps), 1);
ave_thresh = zeros(numel(ave_rps), 1);
for r = 1:numel(app_rps)
    rp = load(app_rps{r});
    app_thresh(r) = size(rp, 1);
    rp(:, 4:6) = (180/pi) .*rp(:, 4:6);
    rpd = diff(rp);
    app_bad_vols{r} = 1 + find(any(abs(rpd(:, 1:3)) > 1.5, 2) | any(abs(rpd(:, 4:6)) > 2, 2));
end
for r = 1:numel(ave_rps)
    rp = load(ave_rps{r});
    ave_thresh(r) = size(rp, 1);
    rp(:, 4:6) = (180/pi) .*rp(:, 4:6);
    rpd = diff(rp);
    ave_bad_vols{r} = 1 + find(any(abs(rpd(:, 1:3)) > 1.5, 2) | any(abs(rpd(:, 4:6)) > 2, 2));
end

%Create output of # of bad volumes per run and which runs to cut, use cut
%off of 9 volumes for appetitive and 17 for aversive (10% of volumes)
app_rpnum = cellfun(@numel, app_bad_vols);
app_really_bad_runs = find(app_rpnum(:) > (0.1 .* app_thresh));
ave_rpnum = cellfun(@numel, ave_bad_vols);
ave_really_bad_runs = find(ave_rpnum(:) > (0.1 .* ave_thresh));

% display bad runs
disp(char(app_vtcs(app_really_bad_runs)));
disp(char(ave_vtcs(ave_really_bad_runs)));

% modify runtimevars for vtc files so as to not include bad volumes
for v = 1:numel(app_vtcs)
    
    % read VTC with data as transio (to save LOTS of time!)
    vtc = xff(app_vtcs{v}, 't');
    vtc.RunTimeVars.AutoSave = true;
    vtc.RunTimeVars.Discard = app_bad_vols{v}(:);
    vtc.SaveRunTimeVars;
    vtc.ClearObject;
end
for v = 1:numel(ave_vtcs)
    vtc = xff(ave_vtcs{v}, 't');
    vtc.RunTimeVars.AutoSave = true;
    vtc.RunTimeVars.Discard = ave_bad_vols{v}(:);
    vtc.SaveRunTimeVars;
    vtc.ClearObject;
end

% print out summary information about subjects
app_subs = unique(n.mfileparts(n.mfileparts(app_vtcs)));
ave_subs = unique(n.mfileparts(n.mfileparts(ave_vtcs)));
tot_subs = sort(union(app_subs, ave_subs));
[~, tot_subs] = n.mfileparts(tot_subs);

% print out summary information
for sc = 1:numel(tot_subs)
    
    % number of appetitive/aversive VTCs
    app_sub_runs = (~cellfun('isempty', regexpi(app_vtcs, tot_subs{sc})));
    app_runs_p1 = find(app_sub_runs & ~cellfun('isempty', regexpi(app_vtcs, 'Phase1')));
    app_runs_p2 = find(app_sub_runs & ~cellfun('isempty', regexpi(app_vtcs, 'Phase1')));
    app_badruns_p1 = intersect(app_runs_p1, app_really_bad_runs);
    app_badruns_p2 = intersect(app_runs_p2, app_really_bad_runs);
    ave_sub_runs = (~cellfun('isempty', regexpi(ave_vtcs, tot_subs{sc})));
    ave_runs_p1 = find(ave_sub_runs & ~cellfun('isempty', regexpi(ave_vtcs, 'Phase1')));
    ave_runs_p2 = find(ave_sub_runs & ~cellfun('isempty', regexpi(ave_vtcs, 'Phase1')));
    ave_badruns_p1 = intersect(ave_runs_p1, ave_really_bad_runs);
    ave_badruns_p2 = intersect(ave_runs_p2, ave_really_bad_runs);
    
    % print out information
    fprintf('%s: app - %d/%d (P1), %d/%d (P2)  // ave - %d/%d (P1), %d/%d (P2)\n', tot_subs{sc}, ...
        numel(app_badruns_p1), numel(app_runs_p1), ...
        numel(app_badruns_p2), numel(app_runs_p2), ...
        numel(ave_badruns_p1), numel(ave_runs_p1), ...
        numel(ave_badruns_p2), numel(ave_runs_p2))
end

% remove bad runs 
app_prts(app_really_bad_runs) = [];
app_vtcs(app_really_bad_runs) = [];
app_rps(app_really_bad_runs)  = [];
ave_prts(ave_really_bad_runs) = [];
ave_vtcs(ave_really_bad_runs) = [];
ave_rps(ave_really_bad_runs)  = [];

% create MDMs
% app_mdm = xff('new:mdm');
% app_mdm.XTC_RTC = [app_vtcs, app_prts];
% app_mdm.RunTimeVars.AutoSave = true;
% app_mdm.RunTimeVars.MotionParameters = app_rps;
% ave_mdm = xff('new:mdm');
% ave_mdm.XTC_RTC = [ave_vtcs, ave_prts];
% ave_mdm.RunTimeVars.AutoSave = true;
% ave_mdm.RunTimeVars.MotionParameters = ave_rps;
