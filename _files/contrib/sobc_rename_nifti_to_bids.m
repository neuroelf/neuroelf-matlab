% switch to base folder
basefolder = 'F:\SOBC_unzipped';
cd(basefolder);

% find subjects (and strip trailing _E)
subs = dir('*_E');
subs = strrep({subs.name}, '_E', '');

% ~ (tilde for Italian keyboard copy-n-paste)

% arrange all the runs in the appropriate folders BIDS guidelines
for i = 1:numel(subs)
    
    % change into subject folder
    cd([basefolder filesep subs{i} '_E']);
    
    % find nifti folders
    nfolders = dir('nifti*');
    if isempty(nfolders)
        continue;
    end
    nfolders = {nfolders.name};
    
    % subject id
    subid = sprintf('sub-%s', subs{i});
    
    % iterate over nifti folders
    for sc = 1:numel(nfolders)
        
        % session folder
        sessid = sprintf('ses-t%d', sc);
        
        % continue if session folder exists
        if exist(sessid, 'dir') > 0
            continue;
        end
        mkdir(sessid);
        mkdir(sessid, 'anat');
        mkdir(sessid, 'func');
        mkdir(sessid, 'other');
        
        % locate MPRage
        mprage = dir([nfolders{sc} filesep '20*MPRage*.nii']);
        if isempty(mprage)
            warning('script:warning', 'No anatomical for %s/%s.', ...
                subid, sessid);
            continue;
        end
        system(sprintf('move "%s" "%s"', ...
            [basefolder filesep subs{i} '_E' filesep nfolders{sc} filesep mprage(1).name], ...
            [basefolder filesep subs{i} '_E' filesep sessid filesep 'anat' filesep subid '_' sessid '_T1w.nii'])); 
        
        % find BOLD runs
        bolds = dir([nfolders{sc} filesep '20*BOLDEMO*.nii']);
        if isempty(bolds)
            continue;
        end
        if mod(numel(bolds), 2) == 1
            warning('script:warning', 'Number of BOLD runs not even for %s/%s.', subid, sessid);
            continue;
        end
        bolds = {bolds.name};
        for rc = 2:2:numel(bolds)
            system(sprintf('move "%s" "%s"', ...
                [basefolder filesep subs{i} '_E' filesep nfolders{sc} filesep bolds{rc}], ...
                [basefolder filesep subs{i} '_E' filesep sessid filesep 'func' filesep subid '_' sessid '_' ...
                 'task-emogonogo_run-' num2str(round(rc/2)) '_bold.nii'])); 
        end
        
        % find all remaining files
        rfiles = dir([nfolders{sc} filesep '20*.*']);
        rfiles = {rfiles.name};
        for rc = 1:numel(rfiles)
            system(sprintf('move "%s" "%s"', ...
                [basefolder filesep subs{i} '_E' filesep nfolders{sc} filesep rfiles{rc}], ...
                [basefolder filesep subs{i} '_E' filesep  sessid filesep 'other' filesep subid '_' sessid '_' ...
                 'other-' rfiles{rc}(16:end)])); 
        end
        rfiles = dir([nfolders{sc} filesep '*20*.*']);
        rfiles = {rfiles.name};
        for rc = 1:numel(rfiles)
            system(sprintf('move "%s" "%s"', ...
                [basefolder filesep subs{i} '_E' filesep nfolders{sc} filesep rfiles{rc}], ...
                [basefolder filesep subs{i} '_E' filesep  sessid filesep 'other' filesep subid '_' sessid '_' ...
                 'other-' rfiles{rc}(16:end)])); 
        end
    end
end

