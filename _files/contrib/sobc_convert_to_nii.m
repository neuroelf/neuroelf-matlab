% base folder
basefolder = 'F:\SOBC_unzipped';
cd(basefolder);

% get task E folders
efolders = dir('*_E');
efolders = {efolders.name};

% iterate over subjects
for sc = 1:numel(efolders)
    
    % find dicom folders
    cd([basefolder filesep efolders{sc}]);
    dfolders = dir('dicom*');
    if isempty(dfolders)
        continue;
    end
    dfolders = {dfolders.name};
    
    % create synonymous nifti folders
    nfolders = strrep(dfolders, 'dicom', 'nifti');
    
    % iterate over dicom folders
    for dc = 1:numel(dfolders)
        
        % if nifti folder exists, do nothing
        if exist([pwd filesep nfolders{dc}], 'dir') > 0
            continue;
        end
        mkdir(nfolders{dc});
        
        % perform dcm2nii operation
        command = sprintf( ...
            '"C:\\Program Files (x86)\\MRIcron\\dcm2nii.exe" -4 y -g n -o "%s" "%s"', ...
            [basefolder filesep efolders{sc} filesep nfolders{dc}], ...
            [basefolder filesep efolders{sc} filesep dfolders{dc}]);
        fprintf('Processing %s\\%s...\n', efolders{sc}, dfolders{dc});
        [status, result] = system(command);
        if status > 0
            disp(result);
        end
    end
end

