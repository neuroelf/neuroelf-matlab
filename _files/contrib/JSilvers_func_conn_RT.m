% Hey Jen, I am taking the liberty of replying in this M-file, which
% at least then means that all comments, etc. are in the same location;
% that way, you can easily see what I was getting at, and possibly make
% any necessary changes yourself...
%
% I started by copying the list of steps from the email (method 1) into
% the M-file as "high-level comments", and I am trying to comment
% additionally where needed.

% first, you would locate the VTC files
n = neuroelf;
VTC_FOLDER = pwd; % change as needed, e.g. '/some/folder/sub*/time*/'
VTC_MASK = '*.vtc'; % change as needed, e.g. {'*task1*.vtc', '*task2*.vtc'}
vtcs = n.findfiles(VTC_FOLDER, VTC_MASK, '-d1');

% specify a pattern that relates to how to get from the "mask" images
% needed to extract GM/WM/CSF signals
amaskpath = 'struct*/spgr*';
% this works by going from the VTC folder to a folder pattern that
% matches this relative path:
% VTC_FOLDER/structural/spgr176slices/

% temporal filtering cut-off (in seconds)
tfilt_cutoff = 100;

% load the func-connectivity (seed) VOI(s)
seed_voi = xff('*.voi');
if numel(seed_voi) ~= 1 || ~isxff(seed_voi, 'voi')
    error('No VOI selected.');
end

% target VMP filename (VTC extension replacement
vmpnamerep = '_seed_voi.vmp';

% iterate over VTCs
for vc = 1:numel(vtcs)

    % load VTC
    vtc = xff(vtcs);

    % %%% regress out nuisance of VTCs

    % locate the wc?x* files (masks) for the global signals
    wcs = n.findfiles([fileparts(vtc.FilenameOnDisk) '/' amaskpath], 'wc?x*.nii', '-d1');
    if isempty(wcs)
        warning('script:general:fileNotFound', 'wc?x*.nii files not found for %s.', ...
            vtc.FilenameOnDisk);
        vtc.ClearObject;
        continue;
    end
    
    % get VTC data from the coordinates for the three masks
    wcx = wcs;
    for wcc = 1:numel(wcs)

        % load wc?x file
        wc = xff(wcs{wcc});
        wc.LoadVoxelData;

        % get transformation matrix to real-space coordinates
        wctrf = wc.CoordinateFrame.Trf;

        % get voxel indices of "in-mask" voxels
        wcx{wcc} = find(wc.VoxelData(:) > 192); % these files are "byte"
        % and a threshold of 192 means 75% probability of this tissue type

        % convert to mm-coordinates
        [x, y, z] = ind2sub(size(wc.VoxelData), wcx{wcc});
        wcx{wcc} = (wctrf * [x(:)'; y(:)'; z(:)'; ones(1, numel(x))])';
        wcx{wcc}(:, 4) = [];

        % convert to VTC indices
        wcx{wcc} = n.bvcoordconv(wcx{wcc}, 'tal2bvx', vtc.BoundingBox);

        % get unique coordinates
        wcx{wcc}(isnan(wcx{wcc})) = [];
        wcx{wcc} = unique(wcx{wcc});

        % clear wc?x object
        wc.ClearObject;

        % access (and average) VTCData
        wcx{wcc} = n.ztrans(n.meannoinfnan(vtc.VTCData(:, wcx{wcc}), 2, true));
        if all(wcx{wcc} == wcx{wcc}(1))
            wcx{wcc} = [];
        end
    end

    % regress out nuisance and temporal filters
    vtc.Filter(struct('nuisreg', cat(2, wcx{:}), 'temp', true, ...
        'tempdct', tfilt_cutoff));
    
    % %%% extract the ROIs from the VTCs
    seed_vtc = vtc.VOITimeCourse(seed_voi);

    % %%% correlate the filtered time-series with the ROI time courses
    seed_vmp = vtc.Correlate(seed_vtc, struct('regnames', {seed_voi.VOINames}));
    % add 'trobust', true to the list if you wish robust regression here!

    % %%% transform to Z values
    for mc = 1:numel(seed_vmp.Map)
        seed_vmp.Map(mc).VMPData(isinf(seed_vmp.Map(mc).VMPData) | isnan(seed_vmp.Map(mc).VMPData)) = 0;
        seed_vmp.Map(mc).VMPData = single(n.fisherr2z(double(seed_vmp.Map(mc).VMPData)));
    end

    % save
    seed_vmp.SaveAs(strrep(vtc.FilenameOnDisk, '.vtc', vmpnamerep));

    % clear both
    seed_vmp.ClearObject;
    vtc.ClearObject;
end

% at this stage, my suggestion would be to then load all VMPs
vmps = n.findfiles(VTC_FOLDER, strrep(VTC_MASK, '.vtc', vmpnamerep), '-d1');

% collect maps
for vc = 1:numel(vmps)
    vmp = xff(vmps{vc});
    vmp.Map = vmp.Map(:);
    if vc == 1
        tvmp = vmp;
    else
        tvmp.Map(end+1:end+numel(vmp.Map)) = vmp.Map;
        vmp.ClearObject;
    end
end

% browse
tvmp.Browse;

% and then you can perform the t-test (incl. robust regression) with the
% compute-formula function...
%
% Hope this makes some sense...?
% /jochen
