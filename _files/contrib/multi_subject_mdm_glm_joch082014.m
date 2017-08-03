
% open Matlab (and make sure NeuroElf is on the path for the command line tools)

% change into the directory that contains the subject folders (s01_nifti, etc.)
subfolder = uigetdir(pwd, 'Please select the folder containing the subjects folders...');
if isequal(subfolder, 0) || ...
    isempty(subfolder) || ...
    exist(subfolder, 'dir') ~= 7
    error('No folder selected.');
end
cd(subfolder);

% locate the VTCs, PRTs, and realignment parameter (RP) files for the "image" task...
vtcs = findfiles([pwd '/s*nifti'], '*image*.vtc', 'depth=1');
prts = findfiles([pwd '/s*nifti/*protocol'], '*image*.prt', 'depth=1');
rps = findfiles([pwd '/s*nifti/*epibold_image'], 'rp*.txt', 'depth=1');

% make sure the numbers match
if numel(prts) ~= numel(vtcs) || ...
    numel(rps) ~= numel(vtcs)
    error('Numbers of files mismatch.');
end

% create MDM object
mdm = xff('new:mdm');

% store information
mdm.XTC_RTC = [vtcs, prts];
mdm.RunTimeVars.AutoSave = true;
mdm.RunTimeVars.MotionParameters = rps;

% set to fixed-effects
mdm.RFX_GLM = 0;
mdm.SeparatePredictors = 2; % or other value, depending on your preference!

% save MDM
mdm.SaveAs;

%To create SDMs (that can be used with BrainVoyager), use the following line of code (after all of the above!)
% run GLM, and *save* SDMs from PRTs
glm = mdm.ComputeGLM(struct('motpars', true, 'tfilter', 240, 'tfilttype', 'dct', 'savesdms', '.sdm'));


%and then alter the MDM to use the SDMs instead:
% replace the ".prt" part in the filenames with ".sdm"
mdm.XTC_RTC(:, 2) = strrep(prts, '.prt', '.sdm');

% save again (with a new name!)
mdm.SaveAs; % give it a new name!!!

