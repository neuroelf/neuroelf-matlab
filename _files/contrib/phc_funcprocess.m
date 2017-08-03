function phc_funcprocess

% for now, some settings are fixed here, can be made optional later!
tr = 2;
ta = -1;
so = 'ai1';

% we need an anatomical scan (and from the name derive the files for
% normalization!)
[anatomical, anatpath] = uigetfile( ...
    {'*.img', 'Analyze 7.5 Image file (*.img)'; ...
     '*.nii', 'NIftI file (*.nii)'}, 'Please select the anatomical file for coregistration...', ...
     'MultiSelect', 'off');
if isequal(anatomical, 0)
    return;
end
anatomical = [anatpath anatomical ',1'];
[~, aname] = fileparts(anatomical);
segsn = [anatpath aname '_seg_sn.mat'];
seginvsn = [anatpath aname '_seg_inv_sn.mat'];

% overall script to process one run of functional data
funcfolder = uigetdir(pwd, 'Please select the folder containing the run to be processed...');
if isequal(funcfolder, 0)
    return;
end
if isempty(funcfolder)
    funcfolder = pwd;
end

% change into folder
cd(funcfolder);

% select files to be processed
[funcfiles, funcpath] = uigetfile( ...
    {'*.img', 'Analyze 7.5 Image files (*.img)'; ...
     '*.nii', 'NIftI files (*.nii)'}, 'Please select the files to be preprocessed...', ...
    'MultiSelect', 'on');
if isequal(funcfiles, 0)
    return;
end
if ischar(funcfiles)
    funcfiles = {funcfiles};
end

% resolve 4D NIftI file (and ask for skipped volumes)
if numel(funcfiles) == 1
    funcfiles{1} = [funcpath funcfiles{1}];
    if ~isempty(regexpi(funcfiles{1}, '\.nii$'))
        skip = inputdlg({'How many volumes to skip (discarded acquisitions)?'}, ...
            'Volumes to skip...', 1, {'4'}, struct('WindowStyle', 'modal'));
        if isempty(skip)
            skip = 0;
        else
            try
                skip = str2double(skip{1});
            catch e;
                disp(e.message);
                return;
            end
        end
        try
            funcfirst = funcfiles{1};
            nh = xff(funcfirst);
            nv = (1+skip):size(nh.VoxelData, 4);
            nh.ClearObject;
            funcfiles = cell(numel(nv), 1);
            for fc = 1:numel(nv)
                funcfiles{fc} = sprintf('%s,%d', funcfirst, nv(fc));
            end
        catch e;
            disp(e.message);
            return;
        end
    end
else
    for fc = 1:numel(funcfiles)
        funcfiles{fc} = [funcpath funcfiles{fc} ',1'];
    end
end

% request onsets (multiple conditions) file
[onsetfile, onsetpath] = uigetfile( ...
    {'*.mat', 'Onsets MAT file (*.mat)'}, 'Please select the onsets file...', ...
    'MultiSelect', 'off');
if isequal(onsetfile, 0)
    onsetfile = {};
else
    if isequal(onsetpath, 0) || ...
        isempty(onsetpath)
        onsetpath = [pwd filesep];
    end
    onsetfile = {[onsetpath onsetfile]};
end

% number of volumes for this run
nv = numel(funcfiles);

% get some settings
try
    nh = xff(funcfiles{1});
    nslices = size(nh.VoxelData, 3);
catch e;
    disp(e.message);
    return;
end

% auto-TA?
if ta <= 0
    ta = tr - (tr / nslices);
end

% slice acquisition order by string
if ischar(so)
    switch lower(so)
        case {'ai1'}
            so = [1:2:nslices, 2:2:nslices];
        otherwise
            error('phc:invalidsetting', 'Invalid slice-order token: %s.', so);
    end
end

% create matlab batch
matlabbatch = repmat({struct('spm', struct)}, 1, 4);

% slice-timing
matlabbatch{1}.spm.temporal.st = struct;
matlabbatch{1}.spm.temporal.st.scans = {funcfiles};
matlabbatch{1}.spm.temporal.st.nslices = nslices;
matlabbatch{1}.spm.temporal.st.tr = tr;
matlabbatch{1}.spm.temporal.st.ta = ta;
matlabbatch{1}.spm.temporal.st.so = so;
matlabbatch{1}.spm.temporal.st.refslice = 1;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

% realign (and write mean for later use)
funcfiles = filenameprep(funcfiles, 'a');
matlabbatch{2}.spm.spatial.realign.estwrite = struct;
matlabbatch{2}.spm.spatial.realign.estwrite.data = {funcfiles};
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions = ...
    struct('quality', 0.9, 'sep', 4, 'fwhm', 5, 'rtm', 1, 'interp', 2, ...
    'wrap', [0,0,0], 'weight', '');
matlabbatch{2}.spm.spatial.realign.estwrite.roptions = ...
    struct('which', [0, 1], 'interp', 4, 'wrap', [0,0,0], 'mask', 1, 'prefix', 'r');

% realignment parameters filename
[~, rpname] = fileparts(funcfiles{1});
rpname = [funcpath 'rp_' rpname '.txt'];

% co-register with anatomical (mean, then others)
matlabbatch{3}.spm.spatial.coreg.estimate = struct;
matlabbatch{3}.spm.spatial.coreg.estimate.ref = {anatomical};
matlabbatch{3}.spm.spatial.coreg.estimate.source = filenameprep(funcfiles(1), 'mean');
matlabbatch{3}.spm.spatial.coreg.estimate.other = funcfiles;
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions = ...
    struct('cost_fun', 'nmi', 'sep', [4,2], ...
    'tol', [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001], ...
    'fwhm', [7,7]);

% normalize
funcfiles = filenameprep(funcfiles, 'w');
matlabbatch{4}.spm.spatial.normalise.write = struct;
matlabbatch{4}.spm.spatial.normalise.write.subj = ...
    struct('matname', {{segsn}}, 'resample', {funcfiles});
matlabbatch{4}.spm.spatial.normalise.write.roptions = ...
    struct('preserve', 0, 'bb', [-78, -112, -50; 78, 76, 85], ...
    'vox', [3,3,3], 'interp', 1, 'wrap', [0,0,0], 'prefix', 'w');

% run job
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

% compute global signal (3-components)
% get three masks
[gm, wm, csf] = loadgwcmasks(anatomical);

% generate extract-array
xtc = zeros(nv, 3);

% iterate over functional files (volumes)
for vc = 1:nv

    % read volume
    v = spm_vol(funcfiles{vc});
    y = spm_read_vols(v);

    % extract gray matter voxels
    ys = y(gm);

    % create implicit mask with available data (not inf/nan/0)
    ysg = ~isinf(ys) & ~isnan(ys) & ys ~= 0;

    % average
    xtc(vc, 1) = sum(ys(ysg)) ./ sum(ysg(:));

    % repeat for WM and CSF
    ys = y(wm);
    ysg = ~isinf(ys) & ~isnan(ys) & ys ~= 0;
    xtc(vc, 2) = sum(ys(ysg)) ./ sum(ysg(:));
    ys = y(csf);
    ysg = ~isinf(ys) & ~isnan(ys) & ys ~= 0;
    xtc(vc, 3) = sum(ys(ysg)) ./ sum(ysg(:));
end

% detrend global signal components (will be added as a special term!)
xtc = detrend(xtc);

% create pseudo-derivative
xtcd = zeros(size(xtc));
xtcd(:, 1) = 10 .* (interp1(1:nv, xtc(:, 1), (1:nv)+0.05, 'spline') - interp1(1:nv, xtc(:, 1), (1:nv)-0.05, 'spline'));
xtcd(:, 2) = 10 .* (interp1(1:nv, xtc(:, 2), (1:nv)+0.05, 'spline') - interp1(1:nv, xtc(:, 2), (1:nv)-0.05, 'spline'));
xtcd(:, 3) = 10 .* (interp1(1:nv, xtc(:, 3), (1:nv)+0.05, 'spline') - interp1(1:nv, xtc(:, 3), (1:nv)-0.05, 'spline'));

% extend (detrended) global signal with derivatives
xtc(:, 4:6) = xtcd;

% load motion parameters
rrp = load(rpname);

% detrend motion parameters
rp = detrend(rrp);

% compute derivative (same as with xtc)
rpd = zeros(size(rp));
rpd(:, 1) = 10 .* (interp1(1:nv, rp(:, 1), (1:nv)+0.05, 'spline') - interp1(1:nv, rp(:, 1), (1:nv)-0.05, 'spline'));
rpd(:, 2) = 10 .* (interp1(1:nv, rp(:, 2), (1:nv)+0.05, 'spline') - interp1(1:nv, rp(:, 2), (1:nv)-0.05, 'spline'));
rpd(:, 3) = 10 .* (interp1(1:nv, rp(:, 3), (1:nv)+0.05, 'spline') - interp1(1:nv, rp(:, 3), (1:nv)-0.05, 'spline'));
rpd(:, 4) = 10 .* (interp1(1:nv, rp(:, 4), (1:nv)+0.05, 'spline') - interp1(1:nv, rp(:, 4), (1:nv)-0.05, 'spline'));
rpd(:, 5) = 10 .* (interp1(1:nv, rp(:, 5), (1:nv)+0.05, 'spline') - interp1(1:nv, rp(:, 5), (1:nv)-0.05, 'spline'));
rpd(:, 6) = 10 .* (interp1(1:nv, rp(:, 6), (1:nv)+0.05, 'spline') - interp1(1:nv, rp(:, 6), (1:nv)-0.05, 'spline'));
rp(:, 7:12) = rpd;

% combine global signal and motion parameters with trend parameter
gs = [xtc, rp, ztrans((1:nv)')];

% create new filename
gsname = filenameprep({rpname}, 'gs');
save(gsname{1}, 'gs', '-ascii');

% empty condition, regressors, factors structs
conds = struct('name', [], 'onset', [], 'duration', [], 'tmod', [], 'pmod', []);
conds(:) = [];
regresss = struct('name', [], 'val', []);
regresss(:) = [];
facts = struct('name', [], 'levels', []);
facts(:) = [];

% regress out global signal, motion parameters and their diffs -> will be
% residual
matlabbatch = {struct('spm', struct)};
matlabbatch{1}.spm.stats.fmri_spec = struct( ...
    'dir', {{funcpath}}, 'timing', struct('units', 'secs', 'RT', tr, 'fmri_t', 16, 'fmri_t0', 1), ...
    'sess', struct('scans', {funcfiles}, 'cond', conds, 'multi', {onsetfile}, ...
    'regress', regresss, 'multi_reg', {gsname}, 'hpf', Inf), ...
    'fact', facts, 'bases', struct('hrf', struct('derivs', [0, 0])), ...
    'volt', 1, 'global', 'None', 'mask', {{''}}, 'cvi', 'none');
spm_jobman('run', matlabbatch);
matlabbatch = {struct('spm', struct)};
matlabbatch{1}.spm.stats.fmri_est = struct( ...
    'spmmat', {{[funcpath 'SPM.mat']}}, 'method', struct('Classical', 1));
spm_jobman('run', matlabbatch);

% betas, mask, RPV, etc. not needed
files = {'^mask\..{3}$','^RPV\..{3}$', '^beta_.{4}\..{3}$', '^ess_.{4}\..{3}$'};
for fc = 1:numel(files)
    selectedfiles = spm_select('List', funcpath, files{fc});
    for sfc = 1:size(selectedfiles,1)
        spm_unlink(deblank(selectedfiles(sfc, :)));
    end
end

% load residual
residuals = xff([funcpath 'ResI*.hdr']);
residuals.LoadVoxelData;
vd = reshape(residuals.VoxelData, prod(residuals.ImgDim.Dim(2:4)), nv);

% filter (within mask)
[b, a] = butter(5, [0.016, 0.18]);
residmask = sum(abs(vd), 2) > 0;
vd(residmask, :) = filter(b, a, vd(residmask, :)')';

% and store back as 4D NIftI, then delete residual files
residuals.VoxelData = reshape(single(vd), size(residuals.VoxelData));
residuals.FileMagic = 'n+1';
residuals.NIIFileType = 2;
residuals.ImgDim.DataType = 16;
residuals.ImgDim.BitsPerPixel = 32;
funcfiles = filenameprep(funcfiles, 'f');
[~, filtfile] = fileparts(funcfiles{1});
residuals.SaveAs([funcpath filtfile '.nii']);
files = {'^ResI_.{4}\..{3}$'};
for fc = 1:numel(files)
    selectedfiles = spm_select('List', funcpath, files{fc});
    for sfc = 1:size(selectedfiles,1)
        spm_unlink(deblank(selectedfiles(sfc, :)));
    end
end

% estimate FD
drp = diff(rrp);
FD = sum(abs([drp(:, 1:3), 50 .* drp(:, 4:6)]), 2);

% estimate global signal for DVARS
% GM+WM mask
gmwm = gm | wm;
ngmwm = sum(gmwm(:));

% re-assess global signal
DVARS = zeros(nv, ngmwm);

% extract global signal
for vc = 1:size(DVARS, 1)
    v = residuals.VoxelData(:, :, :, vc);
    DVARS(vc, :) = v(gmwm);
end
DVARS = diff(psctrans(DVARS));

% compute DVARS
DVARS = sqrt(sum(DVARS .* DVARS, 2) ./ ngmwm);

% select "good" volumes
outliers = find(FD >= .3 | DVARS >= 3.5) + 1;

% add 2 timepoints before and after
outliers = unique([outliers(:) - 2; outliers(:) - 1; outliers(:); outliers(:) + 1; outliers(:) + 2]);
outliers(outliers < 1 | outliers > size(w.VoxelData, 4) - 4) = [];

% re-regress out nuisance (no more RP/temporal filtering, as that would be
% working against "discontinuities"

% extract from surface (already processed!!)

% and then the REAL FUN begins!!



%
% sub-function to alter filenames
function f = filenameprep(f, p)
fs = filesep;
for fc = 1:numel(f)
    [fp, fn, fe] = fileparts(f{fc});
    f{fc} = [fp fs p fn fe];
end


%
% sub-function to get 3 masks for gray matter, white matter, and CSF
function [gm, wm, csf] = loadgwcmasks(anatomical)

% construct filenames
gmfile = filenameprep({anatomical}, 'wrc1');
wmfile = filenameprep({anatomical}, 'wrc2');
csffile = {[fileparts(anatomical) '/wrminverse_c1_c2.img']};

% get three volumes representing warped (normalized) GM/WM/CSF estimates
v1 = spm_vol(gmfile{1});
v2 = spm_vol(wmfile{1});
v3 = spm_vol(csffile{1});

% read volumes into y1/y2/y3
y1 = spm_read_vols(v1);
y2 = spm_read_vols(v2);
y3 = spm_read_vols(v3);

% convert to binary (false/true) with .5 threshold
gm = y1 >= .5;
wm = y2 >= .5;
csf = y3 >= .5;
