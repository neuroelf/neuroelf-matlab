%% HEADER: parameters

% settings
maxvols = 2;        % only work on files with up to maxvols volumes
toslices = 66;      % cut outer slices and leave toslices middle slices
gzipfiles = true;   % gzip output files

% file search pattern (positive and negative)
filepattern = '*FieldMap*.nii.gz';
negpattern = '_to66';

% file name new extension (MUST NOT contain .gz !!)
filenewext = '_to66.nii';



%% CODE: initialization/checks

% using NeuroElf
n = neuroelf;

% change folder
searchpath = uigetdir(pwd, ...
    sprintf('Please select folder to cut ''%s'' files to %d slices...', ...
    filepattern, toslices));
if isempty(searchpath) || isequal(searchpath, 0)
    error('NYSPItools:error:noFolderChosen', 'No folder chosen.');
end
cd(searchpath);

% search for files
files = n.findfiles(pwd, filepattern);
files(~cellfun('isempty', regexpi(files, negpattern))) = [];
if isempty(files)
    error('NYSPItools:error:noFilesFound', 'No files found.');
end

% sanity check
doit = questdlg(sprintf('%d file(s) found. Continue?', numel(files)), ...
    'Do you want to continue?', 'Yes', 'No', 'Yes');
if ~strcmpi(doit, 'yes')
    error('NYSPItools:error:userAbort', 'User abort.');
end


%% CODE: work on files loop
nifti = [];

% iterate over files
for fc = 1:numel(files)

    % open file
    try
        nifti = xff(files{fc});

        % number of volumes too high
        if size(nifti.VoxelData, 4) > maxvols
            nifti.ClearObject;
            continue;
        end

        % get number of slices
        numslices = size(nifti.VoxelData, 3);

        % number of slices OK
        if numslices <= toslices
            nifti.ClearObject;
            continue;
        end

        % how many slices to cut off
        cutslices = numslices - toslices;
        topcut = ceil(cutslices / 2);
        botcut = cutslices - topcut;

        % get transformation value required to add
        transadd = [nifti.DataHist.NIftI1.AffineTransX(3), ...
            nifti.DataHist.NIftI1.AffineTransY(3), ...
            nifti.DataHist.NIftI1.AffineTransZ(3)];

        % compute additional vector
        addvect = botcut .* transadd;

        % add to values
        nifti.DataHist.NIftI1.QuatOffsetX = nifti.DataHist.NIftI1.QuatOffsetX + addvect(1);
        nifti.DataHist.NIftI1.QuatOffsetY = nifti.DataHist.NIftI1.QuatOffsetY + addvect(2);
        nifti.DataHist.NIftI1.QuatOffsetZ = nifti.DataHist.NIftI1.QuatOffsetZ + addvect(3);
        nifti.DataHist.NIftI1.AffineTransX(end) = nifti.DataHist.NIftI1.AffineTransX(end) + addvect(1);
        nifti.DataHist.NIftI1.AffineTransY(end) = nifti.DataHist.NIftI1.AffineTransY(end) + addvect(2);
        nifti.DataHist.NIftI1.AffineTransZ(end) = nifti.DataHist.NIftI1.AffineTransZ(end) + addvect(3);

        % cut slices
        nifti.VoxelData = nifti.VoxelData(:, :, (botcut + 1):(botcut + toslices), :);

        % patch header
        nifti.ImgDim.Dim(4) = toslices;

        % save
        [niftipath, niftifile] = fileparts(files{fc});
        niftifile = regexprep(niftifile, '\.nii$', '');
        newfile = [niftipath filesep niftifile filenewext];
        nifti.SaveAs(newfile);
        nifti.ClearObject;
        nifti = [];

        % gzip
        if gzipfiles
            system(sprintf('gzip -9 -f "%s"', newfile));
        end

    % handle errors
    catch eobj
        if numel(nifti) == 1 && isxff(nifti)
            nifti.ClearObject;
        end
        fprintf('%s\n', eobj.message);
    end
end
