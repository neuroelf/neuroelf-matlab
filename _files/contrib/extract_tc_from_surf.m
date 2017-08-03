% FreeSurfer home folder (needs to be configured somehow later)
freesurfer_home = '/Applications/freesurfer';
subjects_dir = '/Data/reapptrain';
%subjects_dir = '/Users/yoshi/Desktop/reapptrain';

% folder containing the surfaces
anatfolder = '/Data/reapptrain/420/T2/t1';
anatfiles = findfiles(anatfolder, '+*.nii', 'relative=');
if numel(anatfiles) ~= 1
    error('Exactly one anatomical file required.');
end

% check that segmentation output files exist
if exist([anatfolder '/lh.white'], 'file') ~= 2 || ...
    exist([anatfolder '/rh.white'], 'file') ~= 2

    % find anatomical scans
    cd(anatfolder);

    % call freesurfer
    mkdir(subjects_dir, 'surfrecon');
    cmd = strcat( ...
        ['setenv FREESURFER_HOME ' freesurfer_home ';'], ...
         'setenv SUBJECTS_DIR /Users/yoshi/Desktop/reapptrain/;', ...
         'source $FREESURFER_HOME/SetUpFreeSurfer.csh;', ...
         ['recon-all -i ', anatfolder, '/', anatfiles{1}, ' -subjid surfrecon;'], ...
         'recon-all -all -subjid surfrecon');
    unix(cmd);

    % copy lh.white and rh.white from target folder
    copyfile([subjects_dir '/surf/lh.white'], anatfolder);
    copyfile([subjects_dir '/surf/rh.white'], anatfolder);

    % destroy recon folder
    cmd = ['rm -rf ' subjects_dir '/surfrecon'];
    unix(cmd);
end

% check that reduced meshes exist
if exist([anatfolder '/lh_spmspace_32kvert.mat'], 'file') ~= 2 || ...
    exist([anatfolder '/rh_spmspace_32kvert.mat'], 'file') ~= 2

    % load LH white matter surface (as FreeSurfer Binary Format = fsbf)
    lh = xff([anatfolder '/lh.white'], 'fsbf');
    rh = xff([anatfolder '/rh.white'], 'fsbf');
    lht = u8str2double(lh.REMAININGCONTENT, 10, 5);
    lht(:, 1:2) = [];
    rht = u8str2double(rh.REMAININGCONTENT, 10, 5);
    rht(:, 1:2) = [];

    % equal transformation!
    if ~isequaln(lht, rht)
        error('Transformation different for LH and RH mesh, invalid files.');
    end

    % get coordinates from surface and add column of 1's for 4x4 transformation
    lvc = lh.VertexCoordinate;
    rvc = rh.VertexCoordinate;

    % get required transformation parameters
    vm = find(all(lht == 256, 2));
    if all(lht(vm+1, :) == 1)
        lht = lht(vm+2:vm+5, :);
        mc = lht(4, :);

        % apply transformation
        lvc = lvc + ones(size(lvc, 1), 1) * mc;
        rvc = rvc + ones(size(rvc, 1), 1) * mc;
    else
        warning('Transformation matrix not found.');
    end

    % load normalization parameters (inverse, as we apply to coordinates!)
    np = load([anatfolder '/' anatfiles{1}(1:end-4) '_seg_inv_sn.mat']);

    % apply to coordinates
    lvc = applyspmsnc(lvc, np.Tr, np.VG(1).dim, inv(np.VG(1).mat), np.VF(1).mat * np.Affine);
    rvc = applyspmsnc(rvc, np.Tr, np.VG(1).dim, inv(np.VG(1).mat), np.VF(1).mat * np.Affine);

    % reduce vertices
    l = struct('vertices', lvc, 'faces', lh.TriangleVertex);
    lhred = reducepatch(l, 65000);
    save([anatfolder '/lh_spmspace_32kvert.mat'], 'lhred');
    l = struct('vertices', rvc, 'faces', rh.TriangleVertex);
    rhred = reducepatch(l, 65000);
    save([anatfolder '/rh_spmspace_32kvert.mat'], 'rhred');

    % clear objects
    lh.ClearObject;
    rh.ClearObject;

% otherwise load files
else
    load([anatfolder '/lh_spmspace_32kvert.mat']);
    load([anatfolder '/rh_spmspace_32kvert.mat']);
end

% load residuals (which come from normalized/warped but unsmoothed data!)
cd /Data/reapptrain/T2/func_07+RRS_T2_TASK_EPI_ME
res = xff(findfiles(pwd, 'ResI_*.hdr', 'depth=1'));

% apply coordinates transformation to normalized coordinates
vc = lhred.vertices;
vc(:, 4) = 1;
vci = (inv(res.CoordinateFrame.Trf) * vc')';

% create residuals (extracted) time course matrix
restc = zeros(size(vci, 1), size(res.VoxelData, 4));

% extract from volumes
for c = 1:size(res.VoxelData, 4)
    v = res.VoxelData(:, :, :, c);
    restc(:, c) = flexinterpn_method(v, vci(:, 1:3), 'cubic');
end

% apply butterworth (bandpass) filter
[b, a] = butter(5, [0.016, 0.18]);
restcf = filter(b, a, restc')';
