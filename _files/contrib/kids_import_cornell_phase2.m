function kids_import_cornell

% use NeuroElf
n = neuroelf;

% hard-coded folders for now
rawdata = '/Volumes/kids_reapp/kids_grant/Imaging/raw_data';
subjects = '/tmp/subjects';

% select source folder
source = uigetdir(rawdata, 'Please select the raw data folder to import from...');
if isequal(source, 0) || ...
    isempty(source)
    return;
end

% select if Phase 1 or Phase 2
phase = input('Enter Phase 1 (1) or Phase 2 (2)):');

% the user selected a path (not canceled), get subject name
[subjectparent, subjectid] = fileparts(source);

% find all IMA files in the folder
importfiles = n.findfiles([source '/*1'], '*C*.IMA');

% find run numbers
% - first convert to char
filenames = char(importfiles);

% - then find the first character that changes
bjloc = findstr(filenames(1,:),'SI-BJ');
if isempty(bjloc)
    bjloc = findstr(filenames(1, :), 'HSNER');
end

% - and go back to the last dot before this
dotbefore = bjloc+5;

% - get the run numbers
runnumbers = filenames(:, dotbefore+1:dotbefore+4);

% - get first image in each run
firstinrun = [1; 1 + find(any(diff(double(runnumbers)) ~= 0, 2))];
numberofruns = numel(firstinrun);
firstinrun(end+1) = numel(importfiles) + 1;

% - get number of images in each run
imagesperrun = diff(firstinrun);

% create text for listbox
listboxtext = cell(numberofruns, 1);
xyresolution = zeros(numberofruns, 2);
for runcount = numberofruns:-1:1
    
    % read in first file
    dcm = xff(importfiles{firstinrun(runcount)}, 'dcm');
    
    % not a valid run
    if ~any(strcmpi(dcm.DataKeys, 'k_7fe0_0010'))
        firstinrun(runcount) = [];
        numberofruns = numberofruns - 1;
        imagesperrun(runcount) = [];
        listboxtext(runcount) = [];
        xyresolution(runcount, :) = [];
        dcm.ClearObject;
        continue;
    end
    
    % get resolution
    xyresolution(runcount, :) = [dcm.Value('Rows'), dcm.Value('Columns')];
    
    % create text
    listboxtext{runcount} = sprintf('Subject %s, Run %s (%d files, %dx%d resolution)', ...
        subjectid, runnumbers(firstinrun(runcount), :), imagesperrun(runcount), ...
        dcm.Value('Rows'), dcm.Value('Columns'));
    
    % clear object
    dcm.ClearObject;
end

% determine mosaic status
ismosaic = (xyresolution(:, 1) == xyresolution(:, 2)) & ...
    all(xyresolution(:, [1, 1, 1, 1]) ~= (ones(numberofruns, 1) * [64, 128, 256, 512]), 2);

% autodetect initial selection
importselection = imagesperrun;
for runcount = numberofruns:-1:2
    if importselection(runcount-1) == importselection(runcount)
        importselection(runcount-1) = 0;
    end
end
importselection = find(importselection ~= 0);

% show listbox
[importselection, selectok] = listdlg( ...
    'ListString',    listboxtext, ...
    'SelectionMode', 'multiple', ...
    'ListSize',      [640, 400], ...
    'InitialValue',  importselection, ...
    'Name',          'Select runs to import for this subject...');

% return if nothing needs to be done
if isequal(selectok, 0) || ...
    isempty(importselection)
    return;
end

% create target folder if necessary
target = [subjects '/' subjectid sprintf('/Phase%01d',phase)];

% if not exists, create
if exist(target, 'dir') == 0
    n.mkadir(target, '-p');
end

% import data
for runcount = importselection
    
    % for mosaic images
    if ismosaic(runcount)
    
        % load first file (again)
        dcm = xff(importfiles{firstinrun(runcount)}, 'dcm');

        % determine number of slices
        try
            imagedata = reshape(dcm.PixelData, xyresolution(runcount, :));
        catch
            warning('kids_reapp:ImportProblem', ...
                'Error reshaping data for run %d (%s).', runcount, ...
                importfiles{firstinrun(runcount)});
            dcm.ClearObject;
            continue;
        end
        
        % clear object
        dcm.ClearObject;
        
        % determine the base resolution
        if mod(xyresolution(runcount, 1), 64) == 0
            baseres = 64;
        elseif mod(xyresolution(runcount, 1), 96) == 0
            baseres = 96;
        elseif mod(xyresolution(runcount, 1), 80) == 0
            baseres = 80;
        else
            warning('kids_reapp:warning', 'Cannot determine mosaic resolution for run %d.', runcount);
            continue;
        end
        
        % unpack
        imagedata = n.unpackmosaic(imagedata, [1, 2], 3, xyresolution(runcount, 1) / baseres);
        nslices = size(imagedata, 3);
        xyres = [size(imagedata, 1), size(imagedata, 2)];
        
    % otherwise
    else
        
        % number of slices = number of images
        nslices = imagesperrun(runcount);
        
        % resolution directly as in image
        xyres = xyresolution(runcount, :);
    end

    % create target folder
    runtarget = sprintf('%s/run%02d_%04dimg', target, runcount, imagesperrun(runcount));
    if exist(runtarget, 'dir') == 0
        n.mkadir(runtarget,'-p');
    end
    runtargetfile = [runtarget '/vols.nii'];
    
    % create NII
    nii = n.dicom2nii(importfiles(firstinrun(runcount):firstinrun(runcount+1)-1), ...
        struct('mosaic', ismosaic(runcount), 'nslices', nslices, 'xyres', xyres));
    
    % re-set origin
    trf = nii.CoordinateFrame.Trf;
    offset = 0.5 .* (size(nii.VoxelData) - 1);
    origin = -(trf(1:3, 1:3) * n.lsqueeze(offset(1:3)));
    nii.DataHist.NIftI1.AffineTransX(end) = origin(1);
    nii.DataHist.NIftI1.AffineTransY(end) = origin(2);
    nii.DataHist.NIftI1.AffineTransZ(end) = origin(3);
    
    % flip data
    if ndims(nii.VoxelData) == 3
        nii.VoxelData = nii.VoxelData(end:-1:1, :, end:-1:1);
    else
        nii.VoxelData = nii.VoxelData(:, end:-1:1, :, :);
    end
    
    % store data
    nii.SaveAs(runtargetfile);
    
    % clear object
    nii.ClearObject;
end


vans = questdlg('Organize subject folder at this time?', 'Kids import - request', 'Yes');
if ischar(vans) && ...
    strcmpi(vans, 'yes') && ...
   ~isempty(importselection)

    %% Identify, rename and move aversive files
    source = target;
    importdirsavers = n.findfiles(source,'run*_0175img','dirs=1','depth=1');
    mkdir(source,'aversfunc');
    dest=[source '/aversfunc'];
    for i=1:numel(importdirsavers)
        n.renamefile(importdirsavers{i}, strrep(importdirsavers{i},'img','_avers_img'));
    end; 
    importdirsavers = n.findfiles(source,'run*avers*','dirs=1','depth=1');
    for i=1:numel(importdirsavers)
        movefile(importdirsavers{i},dest);
    end

    %% Identify, rename and move aversive files
    importdirsapp = n.findfiles(source,'run*_0098img','dirs=1','depth=1');
    mkdir(source,'appfunc');
    dest=[source '/appfunc'];
    for i=1:numel(importdirsapp)
        n.renamefile(importdirsapp{i},strrep(importdirsapp{i},'img','_app_img'));
    end; 
    importdirsapp = n.findfiles(source,'run*app*','dirs=1','depth=1');
    for i=1:numel(importdirsapp),
        movefile(importdirsapp{i},dest);
    end

    %% Make copies of spgr for appetitive and aversive preprocessing
    spgr_folder = n.findfiles(source, '*0160img','dirs=1','depth=1');
    spgr_file   = n.findfiles(spgr_folder{1}, 'v*.nii', 'depth=1');
    mkdir(source, 'spgr_app');
    copyfile(spgr_file{1}, sprintf('%s/spgr_app/vols.nii', source));
    mkdir(source,'spgr_avers');
    copyfile(spgr_file{1}, sprintf('%s/spgr_avers/vols.nii', source));
end
