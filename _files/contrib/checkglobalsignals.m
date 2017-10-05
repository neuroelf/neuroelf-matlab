% checks some global signal covariates to get a sense of which settings
% make the most sense given a dataset

% settings
psctrans = true; % PSC transform extracted voxels prior to averaging

% NeuroElf library
n = neuroelf;

% load MDM
mdm = xff('*.mdm', 'Please select an MDM to check global signals...');

% check MDM
if ~isxff(mdm, 'mdm')
    error('neuroelf:scriptError', 'No or invalid MDM file selected.');
end

% select global signal mask
msk = xff('*.msk', 'Please select a MSK to extract global signal from...');

% check MSK
if ~isxff(msk, 'msk')
    error('neuroelf:scriptError', 'No or invalid MSK file selected.');
end

% get mask voxels (in mm coordinates so they can be sampled later
mskvox = find(msk.Mask(:) > 0);
msktrf = msk.BoundingBox.QuatB2T;
[mskx, msky, mskz] = ind2sub(size(msk.Mask), mskvox);
mskxyz = [mskx, msky, mskz, ones(numel(mskx), 1)];
mskxyz = (msktrf * mskxyz')';
mskxyz(:, 4) = [];

% iterate over runs
vtcs = mdm.XTC_RTC(:, 1);
for rc = 1:numel(vtcs)
    
    % load VTC
    fprintf('Processing %s (%d/%d)...\n', vtcs{rc}, rc, numel(vtcs));
    vtc = xff(vtcs{rc});
    
    % same space as mask
    if isequal(msktrf, vtc.BoundingBox.QuatB2T)
        
        % simply access data
        globsigs = double(vtc.VTCData(:, mskvox));
        
    % different space as mask
    else
        
        % create necessary array
        globsigs = zeros(size(vtc.VTCData, 1), numel(mskvox));
        
        % iterate over time points
        for vc = 1:size(globsigs, 1)
            
            % sample data
            globsigs(vc, :) = vtc.SampleData3D(mskxyz, struct('mapvol', vc));
        end
    end
    
    % remove any invalid timecourses
    globsigs(:, any(isinf(globsigs) | isnan(globsigs))) = [];
    
    % PSC transform
    if psctrans
        globsigs = n.psctrans(globsigs);
    end
    
    % compute mean, and z-transform
    globsig = n.ztrans(mean(globsigs, 2));
    
    % compute square, first derivative, and squared derivative
    globsigsq = n.ztrans(globsig .* globsig);
    globsigd = n.ztrans([0; diff(globsig)]);
    globsigsqd = n.ztrans(globsigd .* globsigd);
    globsigd2 = n.ztrans([diff(globsig); 0]);
    globsigsqd2 = n.ztrans(globsigd2 .* globsigd2);
    
    % regress against VTC
    vmpgs = vtc.Correlate(globsig, struct('regnames', {{'globsig'}}));
    vmpgssq = vtc.Correlate(globsigsq, struct('regnames', {{'globsigsq'}}));
    vmpgsd = vtc.Correlate(globsigd, struct('regnames', {{'globsigd'}}));
    vmpgssqd = vtc.Correlate(globsigsqd, struct('regnames', {{'globsigsqd'}}));
    vmpgsd2 = vtc.Correlate(globsigd2, struct('regnames', {{'globsigd2'}}));
    vmpgssqd2 = vtc.Correlate(globsigsqd2, struct('regnames', {{'globsigsqd2'}}));
    
    % clear VTC from memory
    vtc.ClearObject;
    
    % either keep or add to existing VMPs
    if rc == 1
        tvmpgs = vmpgs;
        tvmpgssq = vmpgssq;
        tvmpgsd = vmpgsd;
        tvmpgssqd = vmpgssqd;
        tvmpgsd2 = vmpgsd2;
        tvmpgssqd2 = vmpgssqd2;
    else
        tvmpgs.Map(rc) = vmpgs.Map(1);
        tvmpgssq.Map(rc) = vmpgssq.Map(1);
        tvmpgsd.Map(rc) = vmpgsd.Map(1);
        tvmpgssqd.Map(rc) = vmpgssqd.Map(1);
        tvmpgsd2.Map(rc) = vmpgsd2.Map(1);
        tvmpgssqd2.Map(rc) = vmpgssqd2.Map(1);
        vmpgs.ClearObject;
        vmpgssq.ClearObject;
        vmpgsd.ClearObject;
        vmpgssqd.ClearObject;
        vmpgsd2.ClearObject;
        vmpgssqd2.ClearObject;
    end
end
