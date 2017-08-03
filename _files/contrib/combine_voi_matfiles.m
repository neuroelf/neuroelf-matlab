% find subject folders (and remove files that match)
subjects = dir('*_SEMVOIS');
subisdir = cat(1, subjects.isdir);
subjects(~subisdir) = [];

% for each subject
for sc = 1:numel(subjects)

    % find mat files (for left and right)
    lvoi = dir([subjects(sc).name '/VOI_L*.mat']);
    rvoi = dir([subjects(sc).name '/VOI_R*.mat']);

    % load the mat files
    lvoitab = cell(1, numel(lvoi));
    for fc = 1:numel(lvoi)
        matcont = load([subjects(sc).name '/' lvoi(fc).name]);
        matfield = fieldnames(matcont);
        lvoitab{fc} = matcont.(matfield{1});
        lvoitab{fc} = lvoitab{fc}(:);
    end
    lvoitab = cat(2, lvoitab{:});

    % also for right hemisphere
    rvoitab = cell(1, numel(rvoi));
    for fc = 1:numel(rvoi)
        matcont = load([subjects(sc).name '/' rvoi(fc).name]);
        matfield = fieldnames(matcont);
        rvoitab{fc} = matcont.(matfield{1});
        rvoitab{fc} = rvoitab{fc}(:);
    end
    rvoitab = cat(2, rvoitab{:});

    % names
    lvoinames = {lvoi.name};
    rvoinames = {rvoi.name};

    % create output file names
    lfname = sprintf('%s_left.txt', strrep(subjects(sc).name, '_SEMVOIS', ''));
    rfname = sprintf('%s_right.txt', strrep(subjects(sc).name, '_SEMVOIS', ''));

    % write output files
    fid = fopen(lfname, 'w');
    fprintf(fid, '%s\n', sprintf('%s\t', lvoinames{:}));
    fprintf(fid, [repmat('%g\t', 1, size(lvoitab, 2)), '\n'], lvoitab');
    fclose(fid);
    fid = fopen(rfname, 'w');
    fprintf(fid, '%s\n', sprintf('%s\t', rvoinames{:}));
    fprintf(fid, [repmat('%g\t', 1, size(rvoitab, 2)), '\n'], rvoitab');
    fclose(fid);
end
