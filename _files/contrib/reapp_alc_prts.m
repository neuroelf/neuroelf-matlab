% reapp_alc_prts
%
% Create PRT files for the text (log) files from the REAPP_ALC
% project (ROCRPT_K23)
%
% The script automatically locates the text files (selected by
% a file pattern), and then attempts to create PRTs for the
% trial data within each of the files.
%
% Both subject and session ID will be extracted from the data
% file, and the information will be used to generate PRT file
% names.

% change into correct folder
cd('/Users/nasirnaqvi/Desktop/Projects by funding/Active funding/ROCRPT_K23/Data/ROC behavioral data');
filepattern = '*export*.txt';

% use neuroelf
n = neuroelf;

% conditions
condnames = { ...
    'LookCue', ...
    'NegCue', ...
    'LookFood', ...
    'LookAlc', ...
    'NegFood', ...
    'NegAlc', ...
    'Rating'};
condcols = [ ...
     64, 192,  64; ...
    192,  64,  64; ...
    192, 255,   0; ...
      0, 128, 255; ...
    255, 128,   0; ...
    255,   0,   0; ...
    128, 128, 192];

% locate files
logfiles = n.findfiles(pwd, filepattern, '-d1');

% iterate over files
for fc = 1:numel(logfiles)
    
    % try to read file
    try
        logfile = n.readeprimetextlog(logfiles{fc});
    catch ne_eo;
        warning(ne_eo.message);
        continue;
    end
    
    % subject and session
    subsess = sprintf('sub%04d_sess%d', logfile.Subject, logfile.Session);
    
    % get sync times
    synctimes = cat(1, logfile.Log.SynchWithScanner_OffsetTime);
    
    % runstarts
    runstarts = unique(synctimes, 'stable');
    numruns = numel(runstarts);
    
    % iterate over runs
    for rc = 1:numruns
        
        % generate PRT
        prt = xff('new:prt');
        
        % rows corresponding to run
        rs = runstarts(rc);
        runrows = find(synctimes == rs);
        runrows = runrows(:);
        
        % trial condition
        look = runrows(strcmpi({logfile.Log(runrows).Word}, 'look'));
        neg = runrows(strcmpi({logfile.Log(runrows).Word}, 'negative'));
        food = runrows(~cellfun('isempty', regexpi({logfile.Log(runrows).TrialType}, 'food')));
        alc = runrows(~cellfun('isempty', regexpi({logfile.Log(runrows).TrialType}, 'alc')));
        lookfood = intersect(look, food);
        lookalc = intersect(look, alc);
        negfood = intersect(neg, food);
        negalc = intersect(neg, alc);
        
        % add conditions
        prt.AddCond(condnames{1}, [cat(1, logfile.Log(look).Instruct_OnsetTime), cat(1, logfile.Log(look).Instruct_OffsetTime)] - rs, condcols(1, :));
        prt.AddCond(condnames{2}, [cat(1, logfile.Log(neg).Instruct_OnsetTime),  cat(1, logfile.Log(neg).Instruct_OffsetTime)]  - rs, condcols(2, :));
        prt.AddCond(condnames{3}, [cat(1, logfile.Log(lookfood).Image_OnsetTime), cat(1, logfile.Log(lookfood).Image_OffsetTime)] - rs, condcols(3, :));
        prt.AddCond(condnames{4}, [cat(1, logfile.Log(lookalc).Image_OnsetTime),  cat(1, logfile.Log(lookalc).Image_OffsetTime)]  - rs, condcols(4, :));
        prt.AddCond(condnames{5}, [cat(1, logfile.Log(negfood).Image_OnsetTime),  cat(1, logfile.Log(negfood).Image_OffsetTime)]  - rs, condcols(5, :));
        prt.AddCond(condnames{6}, [cat(1, logfile.Log(negalc).Image_OnsetTime),   cat(1, logfile.Log(negalc).Image_OffsetTime)]   - rs, condcols(6, :));
        prt.AddCond(condnames{7}, [cat(1, logfile.Log(runrows).Rating_OnsetTime), cat(1, logfile.Log(runrows).Rating_OffsetTime)] - rs, condcols(7, :));
        
        % save PRT
        prt.SaveAs(sprintf('%s_alc_reapp_run%d.prt', subsess, rc));
        fprintf('PRT saved: %s.\n', prt.FilenameOnDisk);
        prt.ClearObject;
        pause(0.001);
    end
end
