% this script renames VTC files that have the _RUN??_MNI.* pattern
% to contain the original run name from the func/ subfolder

% base folder
base = '/Volumes/BeckPort/SOBC/longFunctional/preprocessed';

% NeuroElf library
n = neuroelf;

% find subject folders
sfolders = n.findfiles([base, '/????*'], 'Phase*', '-d1D');

% iterate over subject (phase) folders
for fc = 1:numel(sfolders)
    
    % locate all files with _RUN??_MNI.* pattern
    runfiles = n.findfiles(sfolders{fc}, '*_RUN??_MNI.*', '-d1');
    runrenamed = runfiles;
    
    % locate all runs within the 'func/' subfolder
    runfolders = n.findfiles([sfolders{fc}, '/func'], 'run*', '-d1Dr');
    
    % make sure it's a multiple!
    if abs(mod(numel(runfiles) / numel(runfolders), 1)) > 1e-8
        warning('script:warning', 'Number of files vs. folders mismatch in %s.', sfolders{fc});
        continue;
    end
    
    % for each run file
    for rfc = 1:numel(runfiles)
        
        % which run is it from
        runnumber = regexpi(runfiles{rfc}, '_RUN(\d+)_MNI');
        runnumber = str2double(runfiles{rfc}(runnumber+4:runnumber+5));
        
        % replace
        runrenamed{rfc} = strrep(runfiles{rfc}, sprintf('_RUN%02d_MNI', runnumber), ...
            sprintf('_%s_MNI', runfolders{runnumber}));
    end
    
    % rename folders
    n.renamefile(runfiles, runrenamed);
end
