function extract_brain(file)

% neuroelf library
using(neuroelf, {'clustercoordsc', 'maxpos', 'renamefile', 'smoothdata3'});

% if result exists, exit
[fpath, fname, fext] = fileparts(file);
if isempty(fpath)
    fpath = pwd;
end
tfile = [fpath '/be/' fname fext];
if exist(tfile, 'file') > 0
    return;
end

% load job
load smooth_coreg_segment_brainextract_rawonly;

% adapt job
file = [file ',1'];
tfile = [fname '_be' fext];
matlabbatch{1}.spm.spatial.smooth.data = {file};
matlabbatch{2}.spm.spatial.coreg.estimate.other = {file};
matlabbatch{3}.spm.spatial.preproc.data = {file};

% initialize SPM
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% run jobs
cd(fpath);
spm_jobman('run', matlabbatch);

% mid-clean up
delete([fpath '/s' fname fext]);

% get data
modulated = xff([fpath '/m' fname fext]);
class1 = xff([fpath '/c1' fname fext]);
class2 = xff([fpath '/c2' fname fext]);
class3 = xff([fpath '/c3' fname fext]);

% create mask
class123 = (uint16(class1.VoxelData(:, :, :)) + uint16(class2.VoxelData(:, :, :)) + uint16(class3.VoxelData(:, :, :))) > uint16(127);

% close objects
class1.ClearObject;
class2.ClearObject;
class3.ClearObject;

% smooth and cluster mask
class123 = smoothdata3(double(class123), [2, 2, 2]) > 0.5;
[cs, class123] = clustercoordsc(class123, 2, 8);
csmax = maxpos(cs);
class123 = (class123 == csmax);
[cs, class123] = clustercoordsc(~class123, 2, 8);
csmax = maxpos(cs);
class123 = (class123 ~= csmax);

% load and mask data
modulated.LoadVoxelData;
modulated.VoxelData(~class123) = 0;

% save as brain extracted
modulated.SaveAs([fpath '/be/' fname fext]);
modulated.ClearObject;

% normalize brain
spm_write_sn([fpath '/be/' fname fext], [fpath '/' fname '_seg_sn.mat'], ...
    struct('interp', 3, 'wrap', [0, 0, 0], 'vox', [1, 1, 1], ...
    'bb', [-90, -126, -72; 90, 90, 102]));

% clean up...
delete([fpath '/c1' fname fext]);
delete([fpath '/c2' fname fext]);
delete([fpath '/c3' fname fext]);
delete([fpath '/m' fname fext]);
delete([fpath '/' fname '_seg_inv_sn.mat']);
delete([fpath '/' fname '_seg_sn.mat']);
renamefile([fpath '/be/w' fname fext], [fpath '/wbe/' fname fext]);
