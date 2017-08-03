function xnat_fmriquality(folder)

% load neuroelf library
n = neuroelf;

% look for scans in folder
scans = n.findfiles(folder, '*', '-d1Dr');

% for each scan folder
for sc = 1:numel(scans)

    % create functional file
    file = n.dcm2nii([folder '/' scans{sc}], folder);

    % run fmriquality
    n.fmriquality(file, struct('motcor', true, 'savefq', false, 'savefqres', [folder '/' scans{sc} '_qa']));

    % delete functional raw data and file
    system(sprintf('rm -rf "%s/%s"', folder, scans{sc}));
    if ~iscell(file)
        file = {file};
    end
    n.mdelete(file);
end

