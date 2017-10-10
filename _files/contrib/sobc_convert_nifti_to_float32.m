% this script reads in all nii files within a folder (and subfolders)
% and forces the datatype to float32 (instead of int)

% root folder
rootfolder = 'F:\SOBC_unzipped';

% locate subfolders
sdirs = dir(rootfolder);
sisdir = cat(1, sdirs.isdir);
sdirs(~sisdir) = [];
niifiles = dir([rootfolder filesep '*.nii']);

% add to files from subfolders
while ~isempty(sdirs)
    ssdirs = sdirs(1);
    ssdirs(:) = [];
    sniifiles = niifiles;
    sniifiles(:) = [];
    for sc = 1:numel(sdirs)
        if any(strcmp(sdirs(sc).name, {'.','..'}))
            continue;
        end
        snewdirs = dir([sdirs(sc).folder filesep sdirs(sc).name]);
        sisdirs = cat(1, snewdirs.isdir);
        snewdirs(~sisdirs) = [];
        ssdirs(end+1:end+numel(snewdirs)) = snewdirs;
        snewfiles = dir([sdirs(sc).folder filesep sdirs(sc).name filesep '*.nii']);
        if ~isempty(snewfiles)
            sniifiles(end+1:end+numel(snewfiles)) = snewfiles;
        end
    end
    sdirs = ssdirs;
    niifiles = [niifiles(:); sniifiles(:)];
end

% iterate over files
for fc = 1:numel(niifiles)
    
    % open file as SPM vol
    niifile = [niifiles(fc).folder filesep niifiles(fc).name];
    v = spm_vol(niifile);
    
    % check datatype
    if v(1).dt(1) == 16
        continue;
    end
    
    % information
    fprintf('Processing file %d/%d (%s)...\n', fc, numel(niifiles), niifile);
    
    % read in all the data (if this fails, the datatype is not supported!)
    try
        y = spm_read_vols(v);
    catch
        continue;
    end
    
    % get dimensions (product) and offset (to update headers)
    pdim = prod(v(1).dim(1:3));
    pi1 = v(1).pinfo(3);
    
    % write out a new file with correct settings
    for vc = 1:numel(v)
        v(vc).fname = [v(vc).fname(1:end-4) '_float32.nii'];
        v(vc).dt(1) = 16;
        v(vc).pinfo = [1; 0; pi1 + (vc - 1) * pdim * 4];
    end
    v = rmfield(v, 'private');
    for vc = 1:numel(v)
        vx = spm_write_vol(v(vc), y(:, :, :, vc));
    end
    
    % move file over
    movefile(vx.fname, niifile, 'f');
end
