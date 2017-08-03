
% load data
load /Volumes/mrp/Imaging/Jochen/group/RawData_allsubs_78runs_pp14

% load MDM
mdm = xff('/Volumes/mrp/Imaging/Jochen/group/MRP_16subs_OLS_SEG_redone.mdm');

% get list id subject IDs
sids = mdm.Subjects;
sall = mdm.Subjects(true);
snum = unique(data(:, 1));

% get unique IAPS numbers
iaps = [unique(data(data(:, 14) > 0, 16)); unique(data(data(:, 15) > 0, 16))];
iaps(iaps < 0) = [];

% iterate over subjects
for sc = 1:numel(sids)

    % get the portion of the data for this subject
    sdata = data(data(:, 1) == snum(sc), :);

    % get the number of runs
    rids = find(strcmp(sall, sids{sc}));
    rnum = unique(sdata(:, 4));
    if numel(rids) ~= numel(rnum)
        error('script:error', 'data mismatch for subject %s.', sids{sc});
    end

    % iterate over runs
    for rc = 1:numel(rids)

        % load PRT
        prt = xff(mdm.XTC_RTC{rids(rc), 2});

        % collapse across conditions
        prt.Collapse('.*_img_neg', 'img_neg');
        prt.Collapse('.*_img_neu', 'img_neu');

        % get condition names (for matching)
        prtc = prt.ConditionNames;
        cond = prt.Cond;
        cneg = cond(strcmp(prtc, 'img_neg'));
        cneu = cond(strcmp(prtc, 'img_neu'));

        % get portion of data pertaining this run
        rdata = sdata(sdata(:, 4) == rnum(rc), :);

        % reduce to neg and neu items
        dataneg = rdata(rdata(:, 14) > 0, :);
        dataneu = rdata(rdata(:, 15) > 0, :);

        % check
        if size(dataneg, 1) ~= size(cneg.OnOffsets, 1) || ...
            size(dataneu, 1) ~= size(cneu.OnOffsets, 1)
            error('script:error', 'data mismatch for subject %s, run %s.', ...
                sids{sc}, rc);
        end

        % create new protocol
        nprt = xff('new:prt');

        % copy main settings
        nprt.ResolutionOfTime = prt.ResolutionOfTime;

        % all IAPS conditions
        for ic = 1:30
            nprt.AddCond(sprintf('neg_%04d', iaps(ic)), zeros(0, 2), [255, 0, 0]);
        end
        for ic = 31:60
            nprt.AddCond(sprintf('neu_%04d', iaps(ic)), zeros(0, 2), [0, 255, 0]);
        end

        % then add conditions for temperatures
        nprt.AddCond('acc_tmp_neg', ...
            cond(strcmp(prtc, 'acc_tmp_neg')).OnOffsets, ...
            cond(strcmp(prtc, 'acc_tmp_neg')).Color);
        nprt.AddCond('rea_tmp_neg', ...
            cond(strcmp(prtc, 'acc_tmp_neg')).OnOffsets, ...
            cond(strcmp(prtc, 'acc_tmp_neg')).Color);
        nprt.AddCond('acc_tmp_neu', ...
            cond(strcmp(prtc, 'acc_tmp_neu')).OnOffsets, ...
            cond(strcmp(prtc, 'acc_tmp_neu')).Color);
        nprt.AddCond('rea_tmp_neu', ...
            cond(strcmp(prtc, 'acc_tmp_neu')).OnOffsets, ...
            cond(strcmp(prtc, 'acc_tmp_neu')).Color);

        % then fill conditions
        for ic = 1:size(dataneg, 1)
            iapc = find(iaps == dataneg(ic, 16));
            nprt.Cond(iapc).OnOffsets = cneg.OnOffsets(ic, :);
            nprt.Cond(iapc).NrOfOnOffsets = 1;
        end
        for ic = 1:size(dataneu, 1)
            iapc = find(iaps == dataneu(ic, 16));
            nprt.Cond(iapc).OnOffsets = cneu.OnOffsets(ic, :);
            nprt.Cond(iapc).NrOfOnOffsets = 1;
        end

        % then save under alternative filename
        nprt.SaveAs(strrep(prt.FilenameOnDisk, '.prt', '_negneu_iaps_st.prt'));

        % clear objects
        prt.ClearObject;
        nprt.ClearObject;
    end
end
