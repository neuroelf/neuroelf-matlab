function varargout = PT_glm_searchlight(part, outof)
%PT_glm_searchlight  Run Searchlight analysis within the PT GLM space
%   PT_glm_searchlight(PART, OUTOF) runs portion PART out of a total of
%   OUTOF parts, such that if PART is 12 and OUTOF is 100, it runs the
%   12th of 100 batches.
%
%   Results will be saved in a VMP that is GZIP-ed automatically.
%
%   This is particularly useful if access to a computational cluster is
%   available, in which case a possible script to submit jobs would be:
%
%   See also: gzip, neuroelf, @neuroelf/private/slsvmclassify,
%             xff, @xff/private/glm_Searchlight.

% preset output
varargout = cell(1, nargout);

% instantiate NeuroElf object as "n" (for searchlight function handle)
t = tic;
fprintf('Initializing NeuroElf (library)...\n');
n = neuroelf;

% load GLM
fprintf('Prepare to run searchlight batch %d of %d...\n', part, outof);
fprintf('Accessing GLM...\n');
glm = xff('GLM_st-psg_FFX.glm');

% create searchlight indices (8mm -> max. number of features/voxel := 81)
fprintf('Finding searchlight source indices...\n');
[midx, midxdist] = glm.Searchlight(struct('dist', 8, 'midxonly', true));

% load the mask (69,797 in-mask voxels)
fprintf('Loading mask...\n');
msk = xff('colin_searchlight.msk');

% subtract the "done" mask
fprintf('Loading and applying done mask...\n');
done = xff('done_searchlight.msk');
msk.Mask(done.Mask > 0) = 0;
done.ClearObject;

% run one VMP batch with options
% - condrep       -> replace the single-trial conditions appropriately
% - condsel       -> select the replaced conditions
% - dodebug       -> show debug messages...
% - mask          -> pass in the mask object
% - maxtime       -> drop out after 110 minutes
% - midx/midxdist -> pass in the indices generated from the default mask
% - nulld         -> run 20 null samples (we can combine those later!)
% - slfunc        -> use slsvmclassify, and simply pass in data and labels
% - subpart       -> specify which range of voxels to run
fprintf('Running searchlight...\n');
vmp = glm.Searchlight(struct( ...
    'condrep', {{'.*Res(Neut|Neg)_.*', 'Resilient'; '.*Sens(Neut|Neg)_.*', 'Sensitive'}}, ...
    'condsel', '^(Resilient|Sensitive)$', ...
    'debug', true, ...
    'mask', msk, ...
    'maxtime', 110 * 60, ...
    'midx', {midx}, 'midxdist', {midxdist}, ...
    'nulld', 20, ...
    'slfunc', {{@n.slsvmclassify, 'data', 'cond'}}, ...
    'subpart', [part, outof]));

% clear GLM/MSK
glm.ClearObject;
msk.ClearObject;

% save the VMP
fprintf('Saving VMP content...\n');
vmpc = getcont(vmp);
save(sprintf('PT_searchlight_Resilient_vs_Sensitive_part_%05dof%05d_vmp.mat', part, outof), 'vmpc');

% clear VMP
vmp.ClearObject;

% time elapsed
fprintf('Used %.0f:%.3f minutes.\n', floor(toc(t)/60), mod(toc(t), 60));
