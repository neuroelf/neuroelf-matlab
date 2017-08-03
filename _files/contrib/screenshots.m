% get global access to NeuroElf GUI config (ne_gcfg)
global ne_gcfg;

% list of maps to cluster and show effects of
maps = [1, 2, 11, 12, 21, 22, 31, 32];
dmaps = 51:58;
% get handle to xff class
x = xff;

% which file to use as underlay for display
try
    vmr = x('colin_brain_ICBMnorm.vmr');
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    vmr = neuroelf_file('c', 'colin_brain_ICBMnorm.vmr');
end
vmr.Browse;

% load GLM
glm = xff('*.glm');
glm.Browse;

% open beta plot window
neuroelf_gui('glmplotbetas');
drawnow;

% get tag and tagstruct
tstr = glm.Handles.PlotFig.Tag;
tstr = tstr(1:8);
tstruct = glm.Handles.PlotFig.TagStruct;

% set contrasts
tstruct.(['LB_' tstr '_Cons']).Value = (1:3)';
neuroelf_gui('glmplotbetasgui', tstr, 'Contrasts');

% and hide options
tstruct.(['CB_' tstr '_OVis']).Value = 0;
neuroelf_gui('glmplotbetasrsz', tstr, 'OptVis');

% load VMP
vmp = xff('38subs_fullruns_prt_sh250_oldconf_singletrial_MSKcolinbrain_ztrans_agecorr.vmp');
vmp.Browse;

% open montage window
neuroelf_gui('vismontage');
drawnow;

% options for montage images
miopt = struct( ...
    'atrans',    false, ...
    'atranscol', [0, 0, 0], ...
    'blx',       [1, 1], ...
    'brds',      0, ...
    'drc',       1, ...
    'drs',       256, ...
    'filename',  './X.png', ...
    'flp',       false, ...
    'flx',       false, ...
    'frame',     [128, 128, 128; -127.99, -127.99, -127.9999], ...
    'hFig',      xfigure('Wnd_NeuroElf_vismontage'), ...
    'imeth',     'cubic', ...
    'imetha',    'cubic', ...
    'join',      true, ...
    'ppv',       4, ...
    'showinfig', false, ...
    'slvar',     vmr, ...
    'stalp',     1, ...
    'stthr',     [0.1, 1.0], ...
    'stvar',     vmp, ...
    'stvix',     1, ...
    'sws',       false, ...
    'tpvol',     1);

% get (Matlab UI) handle of all beta plot windows
cfigs = fieldnames(ne_gcfg.cc);
bpfig = cfigs(~cellfun('isempty', regexpi(cfigs, '^BP')));
bpglm = bpfig;
bpglmf = bpglm;
bpglmvb = bpglm;
for fc = 1:numel(bpfig)
    bpglm{fc} = ne_gcfg.cc.(bpfig{fc}).Config.glm;
    bpfig{fc} = ne_gcfg.cc.(bpfig{fc}).SatelliteMLH;
    [glmpath, bpglmf{fc}] = fileparts(bpglm{fc}.FilenameOnDisk);
end
bpfig = cat(1, bpfig{:});

% set figure size(s)
set(bpfig, 'Position', [100, 100, 1024, 768]);

% iterate over maps
for mc = 1:numel(maps)

    % get map number and name
    mn = maps(mc);
    mnm = makelabel(vmp.Map(mn).Name);

    % enable clustering thresholds and two-tailed test
    vmp.Map(mn).ClusterSize = 57;
    vmp.Map(mn).EnableClusterCheck = 1;
    vmp.Map(mn).ShowPositiveNegativeFlag = 3;

    % set as current map
    neuroelf_gui('setcstatmap', mn);

    % enforce correct threshold
    neuroelf_gui('setstatthrpval', 0.01);

    % create cluster table
    neuroelf_gui('clustertable');

    % continue if no clusters
    if isempty(ne_gcfg.voi.VOI)
        continue;
    end

    % set map in options for montage image
    drawnow;
    miopt.stvix = dmaps(mc);
    miopt.stthr = [vmp.Map(mn).LowerThreshold, vmp.Map(mn).UpperThreshold];

    % restrict clusters to 8mm spheres
    neuroelf_gui('setcluster', 'set', 1:numel(ne_gcfg.voi.VOI));
    neuroelf_gui('limitclusters', {8, 's'});

    % get VOI object handle
    voi = ne_gcfg.voi;

    % extract betas
    for fc = 1:numel(bpglm)
        bpglmvb{fc} = bpglm{fc}.VOIBetas(voi);
    end

    % switch to display map
    neuroelf_gui('setcstatmap', dmaps(mc));
    drawnow;

    % iterate over clusters
    for cc = 2:2:numel(voi.VOI);

        % get coordinate and text
        svoi = voi.VOI(cc);
        svoic = svoi.Voxels(1, :);
        svoin = regexprep(svoi.Name, '^.*_(.H_.*)$', '$1');
        svoin = regexprep(svoin, '_sph.*$', '');
        svoin = regexprep(svoin, '_+$', '');

        % build cluster name
        clname = sprintf('%s%s_%d_%d_%d', mnm, svoin, svoic(1), svoic(2), svoic(3));

        % set cluster
        neuroelf_gui('setcluster', 'set', cc);

        % draw and update
        pause(0.01);
        drawnow;
        pause(0.01);

        % handle barplots
        for fc = 1:numel(bpfig)

            % get frame(s) from figure(s)
            bplot = getframe(bpfig(fc));

            % save as PNGs
            imwrite(bplot.cdata, sprintf('%s_%s_plot.png', bpglmf{fc}, clname));

            % save extract data as well
            vb = bpglmvb{fc}(:, [2, 4, 6], cc);
            save(sprintf('%s_%s_extract.txt', bpglmf{fc}, clname), 'vb', '-ascii');
        end

        % create three panel views
        miopt.drc = 1; % SAG
        miopt.filename = sprintf('%s_%s_SAG.png', bpglmf{fc}, clname);
        miopt.frame = [128, 128, 128; svoic(1), -127.99, -127.9999];
        neuroelf_gui('vismontage_create_ex', miopt);
        miopt.drc = 2; % COR
        miopt.filename = sprintf('%s_%s_COR.png', bpglmf{fc}, clname);
        miopt.frame = [128, 128, 128; -127.99, svoic(2), -127.9999];
        neuroelf_gui('vismontage_create_ex', miopt);
        miopt.drc = 3; % TRA
        miopt.filename = sprintf('%s_%s_TRA.png', bpglmf{fc}, clname);
        miopt.frame = [128, 128, 128; -127.99, -127.99, svoic(3)];
        neuroelf_gui('vismontage_create_ex', miopt);

        % show progress
        disp(sprintf('%d/%d, %d/%d, %s-%s', mc, numel(maps), ...
            0.5 .* cc, 0.5 .* numel(voi.VOI), bpglmf{fc}, clname));
        pause(0.01);
    end
end
