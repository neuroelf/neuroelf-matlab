function [gmap, fmaps] = generategsmask(mdm, opts)

% argument check
if nargin < 1 || ...
    isempty(mdm) || ...
   (~ischar(mdm) && ...
    (numel(mdm) ~= 1 || ...
     ~isxff(mdm, 'mdm')))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'cutoff') || ...
   ~isa(opts.cutoff, 'double') || ...
    numel(opts.cutoff) ~= 1 || ...
    isinf(opts.cutoff) || ...
    isnan(opts.cutoff) || ...
    opts.cutoff <= eps
    opts.cutoff = eps;
else
    opts.cutoff = min(1e-3, opts.cutoff);
end
if ~isfield(opts, 'thresh') || ...
   ~isa(opts.thresh, 'double') || ...
    numel(opts.thresh) ~= 1 || ...
    isinf(opts.thresh) || ...
    isnan(opts.thresh) || ...
    opts.thresh <= eps
    opts.thresh = 1e-6;
else
    opts.thresh = min(1e-3, opts.thresh);
end

% load MDM
mdms = {};
if ischar(mdm)
    try
        mdms{1} = xff(mdm(:)', 't');
        mdm = mdms{1};
        if ~isxff(mdm, 'mdm')
            error( ...
                'neuroelf:BadArgument', ...
                'Not an MDM filename.' ...
            );
        end
    catch ne_eo;
        clearxffobjects(mdms);
        rethrow(ne_eo);
    end
end

% VTC data?
if ~strcmpi(mdm.TypeOfFunctionalData, 'vtc')
    clearxffobjects(mdms);
    error( ...
        'neuroelf:BadObject', ...
        'MDM must be VTC-based.' ...
    );
end
[mdmpath, mdmfile] = fileparts(mdm.FilenameOnDisk);

% big try/catch (to deal with errors)
msko = {[]};
try

    % load (create) SDMs
    csdm = {[], []};
    sdms = {};
    vtcs = repmat({[]}, 1, 24);
    sdms = mdm.SDMs(opts);

    % get subjects list
    subids = mdm.Subjects;
    runsubid = mdm.Subjects(true);

    % generate fmaps container
    fmaps = cell(numel(subids), 4);

    % for each subject
    for sc = 1:numel(subids)

        % get associated runs
        runidx = find(strcmpi(runsubid, subids{sc}));

        % load VTCs (as transio)
        nvol = 1;
        tvol = zeros(numel(runidx), 2);
        vtcsz = zeros(numel(runidx), 4);
        for rc = 1:numel(runidx)
            tvol(rc, 1) = nvol;
            vtcs{rc} = xff(mdm.XTC_RTC{runidx(rc), 1}, 't');
            vtcsz(rc, :) = size(vtcs{rc}.VTCData);
            if isempty(msko{1})
                msk = xff('new:msk');
                msk.Resolution = vtcs{rc}.Resolution;
                msk.XStart = vtcs{rc}.XStart;
                msk.XEnd = vtcs{rc}.XEnd;
                msk.YStart = vtcs{rc}.YStart;
                msk.YEnd = vtcs{rc}.YEnd;
                msk.ZStart = vtcs{rc}.ZStart;
                msk.ZEnd = vtcs{rc}.ZEnd;
                gmsk = zeros(vtcsz(1, 2:4));
                msk.Mask = uint8(gmsk);
                msko{1} = msk;
            end
            nvol = nvol + vtcsz(rc, 1);
            tvol(rc, 2) = nvol - 1;
        end
        if any(any(diff(vtcsz(:, 2:4))))
            error( ...
                'neuroelf:BadObject', ...
                'VTCs must match in size.' ...
            );
        end
        nvol = nvol - 1;

        % combine SDMs
        csdm{1} = sdms{runidx(1)}.CopyObject;
        for rc = 2:numel(runidx)
            csdm{2} = csdm{1}.Concatenate(sdms{runidx(rc)});
            csdm{1}.ClearObject;
            csdm{1} = [];
            csdm = csdm(1, [2, 1]);
        end
        fcp = csdm{1}.FirstConfoundPredictor;

        % create fmap
        fmap = zeros(vtcsz(1, 2:4));

        % data array
        slicedata = zeros([nvol, vtcsz(1, 2:3)]);

        % process along slices
        for zc = 1:vtcsz(1, 4)

            % then runs
            for rc = 1:numel(runidx)
                slicedata(tvol(rc, 1):tvol(rc, 2), :, :) = ...
                    vtcs{rc}.VTCData(:, :, :, zc);
            end

            % add to mask
            gmsk(:, :, zc) = gmsk(:, :, zc) + squeeze(meannoinfnan(slicedata, 1));

            % then compute statistic
            [fmap(:, :, zc), df1, df2] = modelcomp(csdm{1}.SDMMatrix,  ...
                csdm{1}.SDMMatrix(:, fcp:end), slicedata, 1);
        end

        % clear component VTCs and combined SDM
        clearxffobjects(vtcs(1:rc));
        vtcs(1:rc) = {[]};
        csdm{1}.ClearObject;
        csdm{1} = [];

        % remove values without merit
        fmap(isinf(fmap) | isnan(fmap)) = 0;

        % limit data
        cutoff = sdist('finv', opts.cutoff, df1, df2, true);
        fmap = limitrangec(fmap, 0, cutoff, 0);
        pmap = -log(sdist('fcdf', fmap, df1, df2, true));

        % store results
        fmaps(sc, :) = {single(fmap), df1, df2, pmap};
    end

% deal with errors
catch ne_eo;

    % make sure objects are deleted
    clearxffobjects(vtcs);
    clearxffobjects(csdm);
    clearxffobjects(sdms);
    clearxffobjects(mdms);
    clearxffobjects(msko);

    % then bail out
    rethrow(ne_eo);
end

% make sure objects are deleted
clearxffobjects(vtcs);
clearxffobjects(csdm);
clearxffobjects(sdms);
clearxffobjects(mdms);

% compute group map
df1s = cat(1, fmaps{:, 2});
df2s = cat(1, fmaps{:, 3});
dftotal = sum(df2s);
gmap = zeros(size(fmaps{1, 4}));
try

    % weigted p-map
    for sc = 1:size(fmaps, 1)
        gmap = gmap + (df2s(sc) / dftotal) .* fmaps{sc, 4};
    end

    % back-transform
    df1 = harmmean(df1s);
    df2 = harmmean(df2s);
    gmap = sdist('finv', exp(-gmap), df1, df2, true);
    gmap(isinf(gmap) | isnan(gmap) | gmap < 0) = 0;

    % group threshold
    gthresh = sdist('finv', opts.thresh, df1, df2, true);
    gmapmsk = (gmap > gthresh);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    warning( ...
        'neuroelf:SizeMismatch', ...
        'Maps mismatch in size across subjects.' ...
    );
end

% group mask
gmsk = uint8(gmsk > mean(gmsk(gmsk > 0)));

% sample colin to mask space
try
    colin = {[]};
    msks = {gmsk, [], [], []};
    mskn = {'', 'GM', 'WM', 'CSF'};
    colin{1} = neuroelf_file('c', 'colin_brain_ICBMnorm.vmr');
    colin{1}.SmoothData3D(msk.Resolution / sqrt(3), [Inf; -127; 1; 127] * ones(1, 3));
    bbx = msk.BoundingBox;
    cmsk = colin{1}.SampleBVBox(bbx, 1);
    colin{1}.ClearObject;
    colin{1} = [];

    % threshold to GM, WM and CSF
    msks{2} = uint8(cmsk > 120 & cmsk < 150);
    msks{3} = uint8(cmsk > 170);
    msks{4} = uint8(cmsk > 20 & cmsk < 80);

    % smooth the CSF mask and get rid of voxels on the outside
    msks{4} = (smoothdata3(double(msks{4}), 2.5) >= 0.4);
    msks{4}(msks{2} > 0 | msks{3} > 0) = 0;

    % cluster
    [csfs, csfv] = clustercoordsc(msks{4} > 0);
    msks{4} = uint8(csfv == maxpos(csfs));
catch ne_eo;
    warning( ...
        'neuroelf:BadObject', ...
        'Error sampling Colin brain: %s.', ...
        ne_eo.message ...
    );
    msks(2:end) = [];
end
clearxffobjects(colin);

% write out masks
msk.Mask = msks{1};
msk.SaveAs(sprintf('%s/%s_brain.msk', mdmpath, mdmfile));
msk.Mask(gmapmsk) = 0;
msk.SaveAs(sprintf('%s/%s_brain_noGS.msk', mdmpath, mdmfile));
for mc = 2:numel(msks)
    if ~isempty(msks{mc})
        msk.Mask = msks{mc};
        msk.SaveAs(sprintf('%s/%s_colin_%s.msk', mdmpath, mdmfile, mskn{mc}));
        msk.Mask(gmapmsk) = 0;
        msk.SaveAs(sprintf('%s/%s_colin_%s_noGS.msk', mdmpath, mdmfile, mskn{mc}));
    end
end
msk.ClearObject;
