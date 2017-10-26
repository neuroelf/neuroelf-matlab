function neuroelf_makefiles(opts)
% neuroelf_makefiles  - create additional files for NeuroElf's GUI
%
% FORMAT:       neuroelf_makefiles([opts])
%
% Input fields:
%
%       opts        optional settings
%        .all       create all files (false)
%        .bmaskvmr  also create brain-masked versions of VMRs (true)
%        .bvsph     create BV-compatible SPH files of colin (true)
%        .darttmp   create DARTEL-based templates (false)
%        .dhsrf     dual-hemisphere surfaces (false)
%        .findold   find files from an older installation and reuse (true)
%        .flatmaps  create flatmap SRF files (false)
%        .headmesh  create head meshes from VMRs (true)
%        .hiresseg  create 0.5mm segmentation VMRs of colin brain (false)
%        .hiressph  create 0.5mm segmentation-based SPH files (true)
%        .hiressrf  create 0.5mm resolution-based SRF files (true)
%        .hiressrfb create backprojection VMR files from hires SRFs (false)
%        .hiresvmr  create 0.5mm resolution VMR of colin brain (false)
%        .icbmimg   also create ICBM-norm version of colin.img (true)
%        .icbmsph   also create ICBM-norm versions of SPH files (true)
%        .icbmsrf   also create ICBM-norm versions of SRF files (true)
%        .icbmtal   also create ICBM-norm version of talairach.nii (true)
%        .icbmvmr   also create ICBM-norm versions of VMR files (true)
%        .joined40k also create hemisphere-joined version of 40k srfs (true)
%        .recreate  delete all files and start anew (false)
%        .rendskull create file useful for rendering glass-skulls (false)
%        .shenfiles create Shen 268-parcellation atlas files (true)
%        .subcort   also create SRF files for subcortical structures (false)
%        .wmaskvmr  also create (anti-)WM-masked brain VMRs (false)
%
% No output fields.

% Version:  v1.1
% Build:    16060819
% Date:     Jun-08 2016, 7:23 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% argument
if nargin > 0 && ischar(opts) && strcmpi(opts(:)', 'all')
    opts = struct('all', true);
elseif nargin > 0 && ischar(opts) && strcmpi(opts(:)', 'ask')
    dfa = {'  no', '  yes'};
    mfa = inputdlg({ ...
        'Create brain-masked versions of VMRs:', ...
        'Create BV-compatible SPH files of Colin brain:', ...
        'Create dual-hemisphere surfaces:', ...
        'Create flatmap surface files:', ...
        'Create head-mesh surface files:', ...
        'Create 0.5mm segmentation VMRs of Colin brain:', ...
        'Create 0.5mm segmentation-based SPH surface files:', ...
        'Create 0.5mm resolution-VMR-based surfaces:', ...
        'Create backprojected VMR files from hires surfaces:', ...
        'Create 0.5mm resolution VMR of Colin brain:', ...
        'Create ICBM-normalized version of files:', ...
        'Create file useful for rendering glass-skulls:', ...
        'Create joined surfaces from 40k-vertices hemispheres:', ...
        'Create surfaces of subcortical structures:', ...
        'Create Shen et al. 268-parcel Atlas files:', ...
        'Create (anti-) white-matter-masked brain VMRs:', ...
        'Create DARTEL-based templates:'}, ...
        'NeuroElf - makefiles', 1, ...
        dfa(1, [2, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1]));
    if ~iscell(mfa) || numel(mfa) ~= 16
        return;
    end
    mfa = ddeblank(mfa);
    if any(cellfun('isempty', mfa))
        return;
    end
    opts = struct;
    opts.bmaskvmr = (lower(mfa{1}(1)) == 'y');
    opts.bvsph    = (lower(mfa{2}(1)) == 'y');
    opts.dhsrf    = (lower(mfa{3}(1)) == 'y');
    opts.flatmaps = (lower(mfa{4}(1)) == 'y');
    opts.headmesh = (lower(mfa{5}(1)) == 'y');
    opts.hiresseg = (lower(mfa{6}(1)) == 'y');
    opts.hiressph = (lower(mfa{7}(1)) == 'y');
    opts.hiressrf = (lower(mfa{8}(1)) == 'y');
    opts.hiressrfb = (lower(mfa{9}(1)) == 'y');
    opts.hiresvmr = (lower(mfa{10}(1)) == 'y');
    opts.icbmimg  = (lower(mfa{11}(1)) == 'y');
    opts.icbmsph  = (lower(mfa{11}(1)) == 'y');
    opts.icbmsrf  = (lower(mfa{11}(1)) == 'y');
    opts.icbmtal  = (lower(mfa{11}(1)) == 'y');
    opts.icbmvmr  = (lower(mfa{11}(1)) == 'y');
    opts.rendskull = (lower(mfa{12}(1)) == 'y');
    opts.joined40k = (lower(mfa{12}(1)) == 'y');
    opts.subcort  = (lower(mfa{14}(1)) == 'y');
    opts.shenfiles = (lower(mfa{15}(1)) == 'y');
    opts.wmaskvmr = (lower(mfa{16}(1)) == 'y');
    opts.darttmp = (lower(mfa{17}(1)) == 'y');
end
if nargin < 1 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if isfield(opts, 'all') && islogical(opts.all) && numel(opts.all) == 1 && opts.all
    opts.bmaskvmr = true;
    opts.bvsph = true;
    opts.darttmp = true;
    opts.dhsrf = true;
    opts.flatmaps = true;
    opts.headmesh = true;
    opts.hiresseg = true;
    opts.hiressph = true;
    opts.hiressrf = true;
    opts.hiressrfb = true;
    opts.hiresvmr = true;
    opts.icbmimg = true;
    opts.icbmsph = true;
    opts.icbmsrf = true;
    opts.icbmtal = true;
    opts.icbmvmr = true;
    opts.joined40k = true;
    opts.rendskull = true;
    opts.shenfiles = true;
    opts.subcort = true;
    opts.wmaskvmr = true;
end
if ~isfield(opts, 'bmaskvmr') || ~islogical(opts.bmaskvmr) || numel(opts.bmaskvmr) ~= 1
    opts.bmaskvmr = true;
end
if ~isfield(opts, 'bvsph') || ~islogical(opts.bvsph) || numel(opts.bvsph) ~= 1
    opts.bvsph = true;
end
if ~isfield(opts, 'darttmp') || ~islogical(opts.darttmp) || numel(opts.darttmp) ~= 1
    opts.darttmp = false;
end
if ~isfield(opts, 'dhsrf') || ~islogical(opts.dhsrf) || numel(opts.dhsrf) ~= 1
    opts.dhsrf = false;
end
if ~isfield(opts, 'findold') || ~islogical(opts.findold) || numel(opts.findold) ~= 1
    opts.findold = true;
end
if ~isfield(opts, 'flatmaps') || ~islogical(opts.flatmaps) || numel(opts.flatmaps) ~= 1
    opts.flatmaps = false;
end
if ~isfield(opts, 'headmesh') || ~islogical(opts.headmesh) || numel(opts.headmesh) ~= 1
    opts.headmesh = true;
end
if ~isfield(opts, 'hiresseg') || ~islogical(opts.hiresseg) || numel(opts.hiresseg) ~= 1
    opts.hiresseg = false;
end
if ~isfield(opts, 'hiressph') || ~islogical(opts.hiressph) || numel(opts.hiressph) ~= 1
    opts.hiressph = true;
end
if ~isfield(opts, 'hiressrf') || ~islogical(opts.hiressrf) || numel(opts.hiressrf) ~= 1
    opts.hiressrf = true;
end
if ~isfield(opts, 'hiressrfb') || ~islogical(opts.hiressrfb) || numel(opts.hiressrfb) ~= 1
    opts.hiressrfb = false;
end
if ~isfield(opts, 'hiresvmr') || ~islogical(opts.hiresvmr) || numel(opts.hiresvmr) ~= 1
    opts.hiresvmr = false;
end
if ~isfield(opts, 'icbmimg') || ~islogical(opts.icbmimg) || numel(opts.icbmimg) ~= 1
    opts.icbmimg = true;
end
if ~isfield(opts, 'icbmsph') || ~islogical(opts.icbmsph) || numel(opts.icbmsph) ~= 1
    opts.icbmsph = true;
end
if ~isfield(opts, 'icbmsrf') || ~islogical(opts.icbmsrf) || numel(opts.icbmsrf) ~= 1
    opts.icbmsrf = true;
end
if ~isfield(opts, 'icbmtal') || ~islogical(opts.icbmtal') || numel(opts.icbmtal) ~= 1
    opts.icbmtal = true;
end
if ~isfield(opts, 'icbmvmr') || ~islogical(opts.icbmvmr) || numel(opts.icbmvmr) ~= 1
    opts.icbmvmr = true;
end
if ~isfield(opts, 'joined40k') || ~islogical(opts.joined40k) || numel(opts.joined40k) ~= 1
    opts.joined40k = true;
end
if ~isfield(opts, 'recreate') || ~islogical(opts.recreate) || numel(opts.recreate) ~= 1
    opts.recreate = false;
end
if ~isfield(opts, 'rendskull') || ~islogical(opts.rendskull) || numel(opts.rendskull) ~= 1
    opts.rendskull = false;
end
if ~isfield(opts, 'shenfiles') || ~islogical(opts.shenfiles) || numel(opts.shenfiles) ~= 1
    opts.shenfiles = true;
end
if ~isfield(opts, 'subcort') || ~islogical(opts.subcort) || numel(opts.subcort) ~= 1
    opts.subcort = false;
end
if ~isfield(opts, 'wmaskvmr') || ~islogical(opts.wmaskvmr) || numel(opts.wmaskvmr) ~= 1
    opts.wmaskvmr = false;
end

% base path
cpath = [neuroelf_path('colin') filesep];
ppath = [neuroelf_path('shen') filesep];
spath = [neuroelf_path('spm') filesep];
tpath = [neuroelf_path('tal') filesep];

% locate an older installation's files
if opts.findold

    % locate older installation folders
    nefullp = neuroelf_path;
    [neparent, necurrent] = fileparts(nefullp);
    oldneis = findfiles(neparent, 'NeuroElf*', '-d1D');
    oldneis(~cellfun('isempty', regexpi(oldneis, necurrent))) = [];

    % folders and file patterns
    fspec = { ...
        ['_files' filesep 'colin'], {'colin.v*', 'colin_*.*'}, '-d1'; ...
        ['_files' filesep 'neurosynth' filesep 'rawdata'], 'neurosynth*.*', '-d1'; ...
        ['_files' filesep 'neurosynth' filesep 'terms'], '*.nii.gz', '-d1'; ...
        ['_files' filesep 'shenparcel'], {'shen_parcel*.*'}, '-d1'; ...
        ['_files' filesep 'tal'], {'talairach_ICBM*.*'}, '-d1'; ...
        ['_files' filesep 'spm'], {'spm8_dartel_template.mat'}, '-d1'; ...
        ['_files' filesep 'spm'], {'dart*'}, '-d1D'};

    % iterate over those folders
    for ofc = 1:numel(oldneis)

        % indicate what's going in
        idisp(sprintf('Moving files from older NeuroElf (%s) into current one...', oldneis{ofc}));

        % iterate over findfiles spec
        for fsc = 1:size(fspec, 1)

            % locate files in colin folder
            cfiles = findfiles([oldneis{ofc} filesep fspec{fsc, 1}], fspec{fsc, 2:3});

            % iterate over files
            for cfc = 1:numel(cfiles)

                % if file exists, don't move
                [sfolder, sfile, sfext] = fileparts(cfiles{cfc});
                netfile = [nefullp filesep fspec{fsc, 1} filesep sfile sfext];
                if exist(netfile, 'file') == 2
                    continue;
                end

                % otherwise move
                movefile(cfiles{cfc}, netfile);
            end
        end
    end
end

% flatmaps not available
if opts.flatmaps && exist([cpath 'flatmap.mat'], 'file') ~= 2
    warning('neuroelf:makefile:fileNotAvailable', ...
        'The required flatmap.mat file is not available.');
    opts.flatmaps = false;
end

% flatmaps require hiressph and hiressrf
flat = 0;
if opts.flatmaps
    try
        flat = load([cpath 'flatmap.mat']);
        flat = flat.flat;
        ucc = flat.ucc;
        lhcrd = flat.lh;
        lhcrh = flat.lhhi;
        lhcrm = flat.lhmid;
        rhcrd = flat.rh;
        rhcrh = flat.rhhi;
        rhcrm = flat.rhmid;
        flat = 0;
        opts.hiressph = true;
        opts.hiressrf = true;
    catch ne_eo;
        warning('neuroelf:makefiles:fileContentInvalid', ...
            'Invalid flatmap.mat file content (%s).', ne_eo.message);
        opts.flatmaps = false;
        flat = 0;
    end
end

% dual-hemisphere require hiressrf preparation
if opts.dhsrf
    opts.hiresseg = true;
    opts.hiressrf = true;
end

% rendskull requires wmaskvmr
if opts.rendskull
    opts.wmaskvmr = true;
end

% hiressrf not selected
if ~opts.hiressrf

    % then hiressrfb is not available!
    opts.hiressrfb = false;

% otherwise, if selected
elseif opts.hiressrfb

    % create options
    btvopt = struct('fillmode', 'linear', 'nfrom', -0.5, 'nstep', 0.5, 'nto', 0.5, 'res', 0.5, 'triovsmp', 3);
end

% load required file
load([cpath 'sph.mat']);

% recreate
if opts.recreate

    % find all files
    cf = findfiles(cpath, 'colin_*.*');

    % delete each in turn
    for fc = 1:numel(cf)
        try
            delete(cf{fc});
        catch ne_eo;
            error('neuroelf:makefiles:writeError', ...
                'Error deleting files: %s', ne_eo.message);
        end
    end
end

% pre-set objects (for later checks)
chdr = 0;
colin = 0;
dtmp = 0;
hemi = 0;
hcrv = 0;
icbm = 0;
jhem = {[], [], [], [], [], []};
reco = 0;
recoi = 0;
recos = {};
res = 0;
sph = {[], [], []};
skull = 0;
vmr = 0;
pvmr = 0;

% load colin, but remove V16 data
try
    colin = neuroelf_file('c', 'colin.vmr');
    colin.VMRData16 = [];
catch ne_eo;

    % nothing based on colin.vmr
    opts.bmaskvmr = false;
    opts.icbmvmr = false;
    opts.rendskull = false;
    opts.wmaskvmr = false;
    warning('neuroelf:makefiles:internalError', ...
        'Error loading colin.vmr file from NeuroElf installation: %s.', ne_eo.message);
end

% for white-matter mask we need brain mask!
if opts.wmaskvmr
    opts.bmaskvmr = true;
end

% load segmentation-based warping parameters
try
    if opts.icbmimg || opts.icbmvmr
        sn = load([cpath 'seg_sn.mat']);
    end
    if opts.icbmsph || opts.icbmsrf
        sni = load([cpath 'seg_inv_sn.mat']);
    end
catch ne_eo;

    % we can't do any of this then!
    opts.icbmimg = false;
    opts.icbmsph = false;
    opts.icbmsrf = false;
    opts.icbmvmr = false;
    warning('neuroelf:makefiles:internalError', ...
        'Error loading spatial-normalization warping parameters: %s.', ne_eo.message);
end

% for each thing we do, use error handling!
try

    % check if masked version exists
    if opts.bmaskvmr && exist([cpath 'colin_brain.vmr'], 'file') ~= 2

        % information
        idisp('Creating brain-masked version of colin.vmr -> colin_brain.vmr ...');

        % load mask
        cbm = load([cpath 'brainmask.mat']);
        cbm = cbm.cbm;

        % get data for sampling
        vd = max(uint8(1), colin.VMRData);

        % mask brain
        colin.VMRData = uint8(cbm) .* vd;

        % compute fringes coordinates
        cbmx = dilate3d(cbm);
        fringe = find(cbmx(:) & ~cbm(:));
        [frx, fry, frz] = ind2sub(size(vd), fringe);

        % those will be smoothed versions of the masked data
        colin.VMRData(fringe) = ...
            flexinterpn(colin.VMRData, [frx, fry, frz], smoothkern(1.5), 1, 0);

        % then compute inner fringe coordinates
        cbmx = erode3d(cbm);
        fringe = find(cbm(:) & ~cbmx(:));
        [frx, fry, frz] = ind2sub(size(vd), fringe);

        % those are slightly smoothed versions of the current data
        colin.VMRData(fringe) = ...
            flexinterpn(colin.VMRData, [frx, fry, frz], smoothkern(1), 1, 0);

        % save masked version
        colin.SaveAs([cpath 'colin_brain.vmr']);

    % still required for white-matter mask
    elseif opts.wmaskvmr && exist([cpath 'colin_brain_WMmasked.vmr'], 'file') ~= 2

        % load colin
        if isxff(colin, true)
            colin.ClearObject;
            colin = 0;
        end
        colin = xff([cpath 'colin_brain.vmr']);
    end

    % white matter mask
    if opts.wmaskvmr && exist([cpath 'colin_brain_WMmasked.vmr'], 'file') ~= 2

        % information
        idisp('Creating WM-masked version of colin_brain.vmr -> colin_brain_WMmasked.vmr ...');

        % get colin data clustered
        [cols, colv] = clustercoordsc(colin.VMRData(:, :, :) >= 165, 1, 100);
        colv = (colv == maxpos(cols));
        [colv, colvo, colvs] = minarray(colv, 0, 2, 1);
        colv = (smoothdata3(single(colv), [2, 2, 2]) > 0.5);

        % 2 passes of erode and smooth
        for spc = 1:2
            colv = (smoothdata3(single(erode3d(colv)), [2, 2, 2]) > 0.25);
        end

        % set to 0 in data
        colin.VMRData = colin.VMRData(:, :, :);
        colin.VMRData(colvo(1):colvo(1)+colvs(1)-1, ...
            colvo(2):colvo(2)+colvs(2)-1, colvo(3):colvo(3)+colvs(3)-1) = ...
            colin.VMRData(colvo(1):colvo(1)+colvs(1)-1, ...
            colvo(2):colvo(2)+colvs(2)-1, colvo(3):colvo(3)+colvs(3)-1) .* ...
            uint8(~colv);

        % set ventricles to zero also!
        colin.VMRData(smoothdata3(single(colin.VMRData > 0 & colin.VMRData < 60), ...
            [2.5, 2.5, 2.5]) > 0.5) = 0;

        % and establish a minimum value of 30
        colin.VMRData(colin.VMRData > 0 & colin.VMRData < 30) = 30;

        % smooth 0/1 values with a 3mm kernel
        colin.SmoothData3D(3, [Inf, Inf, Inf; -80, -112, -76; 1, 1, 1; 80, 80, 88], ...
            struct('range', [0, 1]));

        % save as
        colin.SaveAs([cpath 'colin_brain_WMmasked.vmr']);

    % still required
    elseif opts.rendskull && exist([cpath 'colin_brain_rendskull.vmr'], 'file') ~= 2

        % load colin
        if isxff(colin, true)
            colin.ClearObject;
            colin = 0;
        end
        colin = xff([cpath 'colin_brain_WMmasked.vmr']);
    end

    % create rendskull VMR
    if opts.rendskull && exist([cpath 'colin_brain_rendskull.vmr'], 'file') ~= 2

        % information
        idisp('Creating render-skull version of colin_brain_WMmasked.vmr -> colin_brain_rendskull.vmr ...');

        % create copy
        skull = colin.CopyObject;

        % load original colin
        vmr = xff([cpath 'colin.vmr']);

        % smooth initial content (lower values)
        skull.SmoothData3D(1.5, [Inf; -127; 1; 127] * ones(1, 3), struct('range', [1, 55]));

        % mix data
        skull.VMRData(:, :, :) = max(uint8(round(0.67 .* double(skull.VMRData))), ...
            uint8(round(0.05 .* (double(vmr.VMRData) + 8))));

        % clear original VMR
        vmr.ClearObject;
        vmr = 0;

        % smooth ever so slightly (everything)
        skull.SmoothData3D(1.25, [Inf; -127; 1; 127] * ones(1, 3));

        % update and set scaling window
        skull.SetScalingWindow([0, 600], true);
        skull.RunTimeVars.ScalingWindow = [3, 450];

        % save and clear
        skull.SaveAs([cpath 'colin_brain_rendskull.vmr']);
        skull.SaveRunTimeVars;
        skull.ClearObject;
        skull = 0;
    end

    % clear colin
    if isxff(colin, true)
        colin.ClearObject;
    end
    colin = 0;

catch ne_eo;
    warning('neuroelf:makefiles:internalError', ...
        'Error creating colin_brain.vmr: %s.', ne_eo.message);
end

try

    % colin normalization as Analyze file (colin_ICBMnorm.img)
    if opts.icbmimg && exist([cpath 'colin_ICBMnorm.img'], 'file') ~= 2

        % information
        idisp('Creating ICBM-norm version of colin.img -> colin_ICBMnorm.img ...');

        % load header/image
        chdr = xff([cpath 'colin.hdr']);
        chdr.LoadVoxelData;

        % normalize within bounding box
        bb = [ones(1, 4) * chdr.CoordinateFrame.Trf'; ...
            [chdr.ImgDim.Dim(2:4), 1] * chdr.CoordinateFrame.Trf'];
        nd = uint8(limitrangec(flexinterpn_method(chdr.VoxelData, ...
            applyspmsnc([Inf, Inf, Inf; bb(1, 1:3); 1, 1, 1; bb(2, 1:3)], ...
            sn.Tr, sn.VG(1).dim, inv(sn.VG(1).mat), sn.Affine), 'lanczos3'), 0, 254));

        % put back into object and save
        chdr.VoxelData = reshape(nd, chdr.ImgDim.Dim(2:4));
        chdr.SaveAs([cpath 'colin_ICBMnorm.hdr']);
        chdr.ClearObject;
        chdr = 0;
    end
catch ne_eo;
    warning('neuroelf:makefiles:internalError', ...
        'Error creating colin_ICBMnorm.img: %s.', ne_eo.message);
end

% creating BV-compatible colin SPH files
try

    % SPH files needed
    if opts.bvsph

        % left hemisphere SPH
        if exist([cpath 'colin_LH_SPH.srf'], 'file') ~= 2

            % information
            idisp('Creating colin_LH_SPH.srf ...');

            % create hemisphere SPH (for triangles and neighbors)
            hemi = spheresrf(100, 6);

            % set extended neighbors to 0, and make other settings
            hemi.ExtendedNeighbors = 0;
            hemi.NrOfVertices = 40962;
            hemi.NrOfTriangles = 81920;
            hemi.MeshCenter = [128, 128, 128];
            hemi.ConvexRGBA = colin_sph.ConvexRGBA;
            hemi.ConcaveRGBA = colin_sph.ConcaveRGBA;
            hemi.NrOfTriangleStrips = 0;
            hemi.TriangleStripSequence = zeros(0, 1);
            hemi.AutoLinkedSRF = '';

            % set coordinates
            hemi.VertexCoordinate = (1 / 256) .* double(colin_sph.LH.VertexCoordinate);

            % zero out normals (reserve space), then compute
            hemi.VertexNormal = zeros(40962, 3);
            hemi.RecalcNormals;

            % and create curvature map (for coloring)
            hcrv = hemi.CurvatureMap;
            hemi.VertexColor = [double(hcrv.Map.SMPData >= 0.01), zeros(40962, 3)];
            hcrv.ClearObject;
            hcrv = 0;

            % then save
            hemi.SaveAs([cpath 'colin_LH_SPH.srf']);

        % if still requested
        elseif opts.icbmsph

            % load
            hemi = xff([cpath 'colin_LH_SPH.srf']);
        end

        % use later?
        if opts.joined40k
            jhem{1} = hemi.CopyObject;
        end

        % normalized left hemisphere
        if opts.icbmsph && exist([cpath 'colin_LH_SPH_ICBMnorm.srf'], 'file') ~= 2

            % information
            idisp('Creating colin_LH_SPH_ICBMnorm.srf ...');

            % apply normalization to coordinates (with INV file -> forward)
            nc = applyspmsnc(128 - hemi.VertexCoordinate(:, [3, 1, 2]), ...
                sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

            % set coordinates, recalc normals
            hemi.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
            hemi.RecalcNormals;

            % then re-create curvature information
            hcrv = hemi.CurvatureMap;
            hemi.VertexColor = [double(hcrv.Map.SMPData >= 0.01), zeros(40962, 3)];
            hcrv.ClearObject;
            hcrv = 0;

            % and save file
            hemi.SaveAs([cpath 'colin_LH_SPH_ICBMnorm.srf']);
            if opts.joined40k
                jhem{4} = hemi.CopyObject;
            end
        elseif opts.icbmsph && opts.joined40k
            jhem{4} = xff([cpath 'colin_LH_SPH_ICBMnorm.srf']);
        end
        if isxff(hemi, true)
            hemi.ClearObject;
            hemi = 0;
        end

        % right hemisphere
        if exist([cpath 'colin_RH_SPH.srf'], 'file') ~= 2

            % information
            idisp('Creating colin_RH_SPH.srf ...');

            % create hemisphere SPH (for triangles and neighbors)
            if ~isxff(hemi, 'srf')
                hemi = spheresrf(100, 6);
            end

            % set extended neighbors to 0, and make other settings
            hemi.ExtendedNeighbors = 0;
            hemi.NrOfVertices = 40962;
            hemi.NrOfTriangles = 81920;
            hemi.MeshCenter = [128, 128, 128];
            hemi.ConvexRGBA = colin_sph.ConvexRGBA;
            hemi.ConcaveRGBA = colin_sph.ConcaveRGBA;
            hemi.NrOfTriangleStrips = 0;
            hemi.TriangleStripSequence = zeros(0, 1);
            hemi.AutoLinkedSRF = '';

            % set coordinates
            hemi.VertexCoordinate = (1 / 256) .* double(colin_sph.RH.VertexCoordinate);

            % zero out normals (reserve space), then compute
            hemi.VertexNormal = zeros(40962, 3);
            hemi.RecalcNormals;

            % and create curvature map (coloring)
            hcrv = hemi.CurvatureMap;
            hemi.VertexColor = [double(hcrv.Map.SMPData >= 0.01), zeros(40962, 3)];
            hcrv.ClearObject;
            hcrv = 0;

            % then save
            hemi.SaveAs([cpath 'colin_RH_SPH.srf']);

        % if still requested
        elseif opts.icbmsph

            % load
            hemi = xff([cpath 'colin_RH_SPH.srf']);
        end

        % use later?
        if opts.joined40k
            jhem{2} = hemi.CopyObject;
        end

        % normalized right hemisphere
        if opts.icbmsph && ...
            exist([cpath 'colin_RH_SPH_ICBMnorm.srf'], 'file') ~= 2

            % information
            idisp('Creating colin_RH_SPH_ICBMnorm.srf ...');

            % normalization to coordinates (with INV file -> forward)
            nc = applyspmsnc(128 - hemi.VertexCoordinate(:, [3, 1, 2]), ...
                sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

            % set coordinates, recalc normals
            hemi.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
            hemi.RecalcNormals;

            % then re-create curvature information
            hcrv = hemi.CurvatureMap;
            hemi.VertexColor = [double(hcrv.Map.SMPData >= 0.01), zeros(40962, 3)];
            hcrv.ClearObject;
            hcrv = 0;

            % and save file
            hemi.SaveAs([cpath 'colin_RH_SPH_ICBMnorm.srf']);
            if opts.joined40k
                jhem{5} = hemi.CopyObject;
            end
        elseif opts.icbmsph && opts.joined40k
            jhem{5} = xff([cpath 'colin_RH_SPH_ICBMnorm.srf']);
        end

        % joined
        if opts.joined40k && exist([cpath 'colin_BH_SPH_joined.srf'], 'file') ~= 2
            jhem{3} = jhem{1}.Combine(jhem{2});
            jhem{3}.SaveAs([cpath 'colin_BH_SPH_joined.srf']);
        end
        if opts.icbmsph && opts.joined40k && exist([cpath 'colin_BH_SPH_joined_ICBMnorm.srf'], 'file') ~= 2
            jhem{6} = jhem{4}.Combine(jhem{5});
            jhem{6}.SaveAs([cpath 'colin_BH_SPH_joined_ICBMnorm.srf']);
        end
        if opts.joined40k
            clearxffobjects(jhem);
            jhem = {[], [], [], [], [], []};
        end
    end
catch ne_eo;
    warning( ...
        'neuroelf:xffError', ...
        'Error creating colin_xH_SPH.srf file/s: %s.', ...
        ne_eo.message ...
    );
end

try

    % high-resolution VMR
    if opts.hiresseg || ...
        opts.hiressph || ...
        opts.hiressrf

        % information
        idisp('Preparing hi-resolution segmentation VMR (for various files)...');

        % load file
        seg_05mm = struct;
        load([cpath 'seg_05mm.mat']);

        % create VMR
        vmr = xff('new:vmr');

        % settings
        vmr.FramingCube = 512;
        vmr.PosInfoVerified = 1;
        vmr.Slice1CenterX = -127.75;
        vmr.SliceNCenterX = 127.75;
        vmr.NRows = 512;
        vmr.NCols = 512;
        vmr.SliceThickness = 0.5;
        vmr.VoxResX = 0.5;
        vmr.VoxResY = 0.5;
        vmr.VoxResZ = 0.5;
        vmr.VoxResInTalairach = 1;
        vmr.VoxResVerified = 1;

        % fill with LH information
        vmr.VMRData = uint8(240) .* uint8(seg_05mm.lh.vox);

        % offsets
        vmr.OffsetX = double(seg_05mm.lh.off(1));
        vmr.OffsetY = double(seg_05mm.lh.off(2));
        vmr.OffsetZ = double(seg_05mm.lh.off(3));

        % VMR file?
        if opts.hiresseg && ...
            exist([cpath 'colin_hires_LH_seg.vmr'], 'file') ~= 2

            % information
            idisp('Saving colin_hires_LH_seg.vmr ...');
            vmr.SaveAs([cpath 'colin_hires_LH_seg.vmr']);
        end

        % SRF file?
        if (opts.hiressph || ...
            opts.hiressrf) && ...
            exist([cpath 'colin_hires_LH_seg_RECOSM.srf'], 'file') ~= 2

            % information
            idisp('Creating colin_hires_LH_seg_RECOSM.srf ...');

            % reconstruct
            reco = vmr.DBReco;
            reco.ConvexRGBA = colin_sph.ConvexRGBA;
            reco.ConcaveRGBA = colin_sph.ConcaveRGBA;

            % smooth (a bit)
            reco.Smooth(100, 0.1, struct('show', false));

            % save
            if opts.hiressrf

                % first, curvature map (coloring)
                hcrv = reco.CurvatureMap;
                reco.VertexColor = [double(hcrv.Map.SMPData >= 0.002), ...
                    zeros(reco.NrOfVertices, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % then save
                reco.SaveAs([cpath 'colin_hires_LH_seg_RECOSM.srf']);

                % smooth, inflate, and save again
                recovc = reco.VertexCoordinate;
                recovn = reco.VertexNormal;
                reco.Smooth(500, 0.5, struct('show', false, 'distwsq', true));
                reco.Inflate(125, 0.5, struct('show', false));
                reco.SaveAs([cpath 'colin_hires_LH_seg_RECOSM_sm500_05_inf125_05.srf']);
                reco.Smooth(500, 0.5, struct('show', false, 'distwsq', true));
                reco.Inflate(125, 0.5, struct('show', false));
                reco.SaveAs([cpath 'colin_hires_LH_seg_RECOSM_sm1000_05_inf250_05.srf']);

                % then restore unsmoothed coordinates
                reco.VertexCoordinate = recovc;
                reco.VertexNormal = recovn;
            end

            % ICBM-norm version
            if opts.icbmsrf && ...
                exist([cpath 'colin_hires_LH_seg_RECOSM_ICBMnorm.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_hires_LH_seg_RECOSM_ICBMnorm.srf ...');

                % make copy (so as not to overwrite coordinates!)
                recoi = reco.CopyObject;

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - reco.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                recoi.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                recoi.RecalcNormals;

                % then re-create curvature information
                hcrv = recoi.CurvatureMap;
                recoi.VertexColor = [double(hcrv.Map.SMPData >= 0.002), ...
                    zeros(recoi.NrOfVertices, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % and save file
                recoi.AutoLinkedSRF = '';
                recoi.SaveAs([cpath 'colin_hires_LH_seg_RECOSM_ICBMnorm.srf']);

                % smooth, inflate, and save again
                recoi.Smooth(500, 0.5, struct('show', false, 'distwsq', true));
                recoi.Inflate(125, 0.5, struct('show', false));
                recoi.SaveAs([cpath 'colin_hires_LH_seg_RECOSM_ICBMnorm_sm500_05_inf125_05.srf']);
                recoi.Smooth(500, 0.5, struct('show', false, 'distwsq', true));
                recoi.Inflate(125, 0.5, struct('show', false));
                recoi.SaveAs([cpath 'colin_hires_LH_seg_RECOSM_ICBMnorm_sm1000_05_inf250_05.srf']);

                % clear again
                recoi.ClearObject;
                recoi = 0;
            end

        % still required
        elseif (opts.hiressph && ...
           (exist([cpath 'colin_hires_LH_SPH160k.srf'], 'file') ~= 2 || ...
            exist([cpath 'colin_midres_LH_SPH40k.srf'], 'file') ~= 2 || ...
            exist([cpath 'colin_lowres_LH_SPH10k.srf'], 'file') ~= 2)) || ...
           (opts.flatmaps && ...
            exist([cpath 'colin_hires_LH_seg_FLAT.srf'], 'file') ~= 2) || ...
           (opts.hiressrfb && ...
            exist([cpath 'colin_hires_LH_segback.vmr'], 'file') ~= 2)

            % load
            reco = xff([cpath 'colin_hires_LH_seg_RECOSM.srf']);
        end

        % hires SRF backprojection VMR
        if opts.hiressrfb && ...
            exist([cpath 'colin_hires_LH_segback.vmr'], 'file') ~= 2

            % information
            idisp('Creating colin_hires_LH_segback+hull.vmr');

            % create back projection
            reco.BackToVMR(btvopt);

            % create new hiresvmr
            hvmr = newhiresvmr;
            hvmr.VMRData = reco.VertexVMRData(77:512, 53:404, 77:440);
            reco.VertexVMRData(:) = [];
            hvmr.OffsetX = 76;
            hvmr.OffsetY = 52;
            hvmr.OffsetZ= 76;
            hvmr.SmoothData3D(1.25);
            hvmr.SetScalingWindow;
            hvmr.RunTimeVars.ScalingWindow = [33, 150];
            hvmr.RunTimeVars.ScalingWindowLim = [33, 300];
            hvmr.SaveAs([cpath 'colin_hires_LH_segback_hull.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.VMRData(floodfill3(hvmr.VMRData <= 50, 219, 169, 233)) = 150;
            hvmr.SmoothData3D(1.25, [], struct('range', [0, 125]));
            hvmr.SetScalingWindow;
            hvmr.RunTimeVars.ScalingWindow = [33, 100];
            hvmr.SaveAs([cpath 'colin_hires_LH_segback.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.ClearObject;
        end

        % hires SRF backprojection VMR
        if opts.hiressrfb && ...
            opts.icbmsrf && ...
            exist([cpath 'colin_hires_LH_segback_ICBMnorm.vmr'], 'file') ~= 2 && ...
            exist([cpath 'colin_hires_LH_seg_RECOSM_ICBMnorm.srf'], 'file') == 2

            % information
            idisp('Creating colin_hires_LH_segback_ICBMnorm+hull.vmr');

            % create back projection
            recoi = xff([cpath 'colin_hires_LH_seg_RECOSM_ICBMnorm.srf']);
            recoi.BackToVMR(btvopt);

            % create new hiresvmr
            hvmr = newhiresvmr;
            hvmr.VMRData = recoi.VertexVMRData(77:512, 53:404, 77:440);
            recoi.ClearObject;
            recoi = 0;
            hvmr.OffsetX = 76;
            hvmr.OffsetY = 52;
            hvmr.OffsetZ= 76;
            hvmr.SmoothData3D(1.25);
            hvmr.SetScalingWindow;
            hvmr.RunTimeVars.ScalingWindow = [33, 150];
            hvmr.RunTimeVars.ScalingWindowLim = [33, 300];
            hvmr.SaveAs([cpath 'colin_hires_LH_segback_ICBMnorm_hull.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.VMRData(floodfill3(hvmr.VMRData <= 47, 219, 169, 233)) = 150;
            hvmr.SmoothData3D(1.25, [], struct('range', [0, 125]));
            hvmr.SetScalingWindow;
            hvmr.RunTimeVars.ScalingWindow = [33, 100];
            hvmr.SaveAs([cpath 'colin_hires_LH_segback_ICBMnorm.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.ClearObject;
        end

        % flatmaps
        if opts.flatmaps && ...
            exist([cpath 'colin_hires_LH_seg_FLAT.srf'], 'file') ~= 2

            % information
            idisp('Creating colin_hires_LH_seg_FLAT.srf ...');

            % copy object
            flat = reco.CopyObject;

            % fill color first
            flat.VertexColor(:) = 0;
            flat.VertexColor(:, 1) = ...
                ucc(1 + floor((1 / 2 ^ 28) .* double(lhcrd)));

            % then remove color from coordinates information
            lhcrd = mod(lhcrd, uint32(2 ^ 28));

            % and fill coordinates
            flat.VertexCoordinate(:, 1) = (1 / 64) .* ...
                floor((1 / 2 ^ 14) .* double(lhcrd));
            flat.VertexCoordinate(:, 2) = (1 / 64) .* ...
                double(mod(lhcrd, uint32(2 ^ 14)));
            flat.VertexCoordinate(:, 3) = 128;

            % recalc normals
            flat.RecalcNormals;

            % save and clear
            flat.SaveAs([cpath 'colin_hires_LH_seg_FLAT.srf']);
            flat.ClearObject;
            flat = 0;
        end

        % hi-res SPH version
        if opts.hiressph

            % re-sampled surface -> 160k LH
            if exist([cpath 'colin_hires_LH_SPH160k.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_hires_LH_SPH160k.srf ...');

                % generate sph's
                sph{3} = spheresrf(140, 7);

                % create TSM
                tsm = xff('new:tsm');
                nv = sph{3}.NrOfVertices;
                tsm.NrOfTargetVertices = nv;
                tsm.NrOfSourceVertices = reco.NrOfVertices;
                tsm.NrOfSourceTriangles = reco.NrOfTriangles;

                % get data
                tsmd = double(seg_05mm.lh.tsm(1:nv));
                teld = floor(tsmd ./ (2 ^ 22));
                tsm.SourceTriangleOfTarget = mod(tsmd, 2 ^ 22);
                tsm.TriangleEdgeLengths = ...
                    0.05 .* ([floor(teld ./ 32), mod(teld, 32)] - 10);

                % re-sample surface (160k)
                res = reco.ApplyTSM(tsm, sph{3});
                tsm.ClearObject;

                % smooth (a bit)
                res.Smooth(15, 0.02, struct('show', false));

                % curvature map (coloring)
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.005), zeros(163842, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % save
                res.SaveAs([cpath 'colin_hires_LH_SPH160k.srf']);

            % load anyway
            elseif (opts.icbmsph && ...
                exist([cpath 'colin_hires_LH_SPH160k_ICBMnorm.srf'], 'file') ~= 2) || ...
               (opts.flatmaps && ...
                exist([cpath 'colin_hires_LH_SPH160k_FLAT.srf'], 'file') ~= 2)

                % load
                res = xff([cpath 'colin_hires_LH_SPH160k.srf']);
            end

            % flatmaps
            if opts.flatmaps && ...
                exist([cpath 'colin_hires_LH_SPH160k_FLAT.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_hires_LH_SPH160k_FLAT.srf ...');

                % copy object
                flat = res.CopyObject;

                % fill color first
                flat.VertexColor(:) = 0;
                flat.VertexColor(:, 1) = ...
                    ucc(1 + floor((1 / 2 ^ 28) .* double(lhcrh)));

                % then remove color from coordinates information
                lhcrh = mod(lhcrh, uint32(2 ^ 28));

                % and fill coordinates
                flat.VertexCoordinate(:, 1) = (1 / 64) .* ...
                    floor((1 / 2 ^ 14) .* double(lhcrh));
                flat.VertexCoordinate(:, 2) = (1 / 64) .* ...
                    double(mod(lhcrh, uint32(2 ^ 14)));
                flat.VertexCoordinate(:, 3) = 128;

                % recalc normals
                flat.RecalcNormals;

                % save and clear
                flat.SaveAs([cpath 'colin_hires_LH_SPH160k_FLAT.srf']);
                flat.ClearObject;
                flat = 0;
            end

            % ICBM-normalized version
            if opts.icbmsph && ...
                exist([cpath 'colin_hires_LH_SPH160k_ICBMnorm.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_hires_LH_SPH160k_ICBMnorm.srf ...');

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - res.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                res.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                res.RecalcNormals;

                % then re-create curvature information
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.005), zeros(163842, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % and save file
                res.SaveAs([cpath 'colin_hires_LH_SPH160k_ICBMnorm.srf']);
            end

            % clear object
            clearxffobjects({res});
            res = 0;

            % re-sampled surface -> 40k LH
            if exist([cpath 'colin_midres_LH_SPH40k.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_midres_LH_SPH40k.srf ...');

                % generate sph's
                sph{2} = spheresrf(140, 6);

                % create TSM
                tsm = xff('new:tsm');
                nv = sph{2}.NrOfVertices;
                tsm.NrOfTargetVertices = nv;
                tsm.NrOfSourceVertices = reco.NrOfVertices;
                tsm.NrOfSourceTriangles = reco.NrOfTriangles;
                tsmd = double(seg_05mm.lh.tsm(1:nv));
                teld = floor(tsmd ./ (2 ^ 22));
                tsm.SourceTriangleOfTarget = mod(tsmd, 2 ^ 22);
                tsm.TriangleEdgeLengths = ...
                    0.05 .* ([floor(teld ./ 32), mod(teld, 32)] - 10);

                % re-sample surface (40k)
                res = reco.ApplyTSM(tsm, sph{2});
                tsm.ClearObject;

                % smooth (a bit)
                res.Smooth(20, 0.015, struct('show', false));

                % curvature map (coloring)
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.01), zeros(40962, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % save
                res.SaveAs([cpath 'colin_midres_LH_SPH40k.srf']);

            % load anyway
            elseif (opts.icbmsph && ...
                exist([cpath 'colin_midres_LH_SPH40k_ICBMnorm.srf'], 'file') ~= 2) || ...
               (opts.flatmaps && ...
                exist([cpath 'colin_midres_LH_SPH40k_FLAT.srf'], 'file') ~= 2)

                % load
                res = xff([cpath 'colin_midres_LH_SPH40k.srf']);
            end

            % flatmaps
            if opts.flatmaps && ...
                exist([cpath 'colin_midres_LH_SPH40k_FLAT.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_midres_LH_SPH40k_FLAT.srf ...');

                % copy object
                flat = res.CopyObject;

                % fill color first
                flat.VertexColor(:) = 0;
                flat.VertexColor(:, 1) = ...
                    ucc(1 + floor((1 / 2 ^ 28) .* double(lhcrm)));

                % then remove color from coordinates information
                lhcrm = mod(lhcrm, uint32(2 ^ 28));

                % and fill coordinates
                flat.VertexCoordinate(:, 1) = (1 / 64) .* ...
                    floor((1 / 2 ^ 14) .* double(lhcrm));
                flat.VertexCoordinate(:, 2) = (1 / 64) .* ...
                    double(mod(lhcrm, uint32(2 ^ 14)));
                flat.VertexCoordinate(:, 3) = 128;

                % recalc normals
                flat.RecalcNormals;

                % save and clear
                flat.SaveAs([cpath 'colin_midres_LH_SPH40k_FLAT.srf']);
                flat.ClearObject;
                flat = 0;
            end

            % ICBM-normalized version
            if opts.icbmsph && ...
                exist([cpath 'colin_midres_LH_SPH40k_ICBMnorm.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_midres_LH_SPH40k_ICBMnorm.srf ...');

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - res.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                res.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                res.RecalcNormals;

                % then re-create curvature information
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.01), zeros(40962, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % and save file
                res.SaveAs([cpath 'colin_midres_LH_SPH40k_ICBMnorm.srf']);
            end

            % clear object
            clearxffobjects({res});
            res = 0;

            % re-sampled surface -> 10k LH
            if exist([cpath 'colin_lowres_LH_SPH10k.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_lowres_LH_SPH10k.srf ...');

                % generate sph's
                sph{1} = spheresrf(140, 5);

                % create TSM
                tsm = xff('new:tsm');
                nv = sph{1}.NrOfVertices;
                tsm.NrOfTargetVertices = nv;
                tsm.NrOfSourceVertices = reco.NrOfVertices;
                tsm.NrOfSourceTriangles = reco.NrOfTriangles;
                tsmd = double(seg_05mm.lh.tsm(1:nv));
                teld = floor(tsmd ./ (2 ^ 22));
                tsm.SourceTriangleOfTarget = mod(tsmd, 2 ^ 22);
                tsm.TriangleEdgeLengths = ...
                    0.05 .* ([floor(teld ./ 32), mod(teld, 32)] - 10);

                % re-sample surface (10k)
                res = reco.ApplyTSM(tsm, sph{1});
                tsm.ClearObject;

                % smooth (a bit)
                res.Smooth(25, 0.01, struct('show', false));

                % curvature map (coloring)
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.02), zeros(10242, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % save
                res.SaveAs([cpath 'colin_lowres_LH_SPH10k.srf']);

            % load anyway
            elseif opts.icbmsph && ...
                exist([cpath 'colin_lowres_LH_SPH10k_ICBMnorm.srf'], 'file') ~= 2

                % load
                res = xff([cpath 'colin_lowres_LH_SPH10k.srf']);
            end

            % ICBM-normalized version
            if opts.icbmsph && ...
                exist([cpath 'colin_lowres_LH_SPH10k_ICBMnorm.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_lowres_LH_SPH10k_ICBMnorm.srf ...');

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - res.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                res.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                res.RecalcNormals;

                % then re-create curvature information
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.02), zeros(10242, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % and save file
                res.SaveAs([cpath 'colin_lowres_LH_SPH10k_ICBMnorm.srf']);
            end

            % clear object
            clearxffobjects({res});
            res = 0;

        end

        % clear object
        clearxffobjects({reco});
        reco = 0;

        % fill with RH information
        vmr.VMRData = uint8(240) .* uint8(seg_05mm.rh.vox);

        % offsets
        vmr.OffsetX = double(seg_05mm.rh.off(1));
        vmr.OffsetY = double(seg_05mm.rh.off(2));
        vmr.OffsetZ = double(seg_05mm.rh.off(3));

        % VMR file?
        if opts.hiresseg && ...
            exist([cpath 'colin_hires_RH_seg.vmr'], 'file') ~= 2

            % information
            idisp('Saving colin_hires_RH_seg.vmr ...');
            vmr.SaveAs([cpath 'colin_hires_RH_seg.vmr']);
        end

        % SRF file?
        if (opts.hiressph || ...
            opts.hiressrf) && ...
            exist([cpath 'colin_hires_RH_seg_RECOSM.srf'], 'file') ~= 2

            % information
            idisp('Creating colin_hires_RH_seg_RECOSM.srf ...');

            % reconstruct
            reco = vmr.DBReco;
            reco.ConvexRGBA = colin_sph.ConvexRGBA;
            reco.ConcaveRGBA = colin_sph.ConcaveRGBA;

            % smooth (a bit)
            reco.Smooth(100, 0.1, struct('show', false));

            % save
            if opts.hiressrf

                % first, curvature map (coloring)
                hcrv = reco.CurvatureMap;
                reco.VertexColor = [double(hcrv.Map.SMPData >= 0.002), ...
                    zeros(reco.NrOfVertices, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % then save
                reco.SaveAs([cpath 'colin_hires_RH_seg_RECOSM.srf']);

                % smooth, inflate, and save again
                recovc = reco.VertexCoordinate;
                recovn = reco.VertexNormal;
                reco.Smooth(500, 0.5, struct('show', false, 'distwsq', true));
                reco.Inflate(125, 0.5, struct('show', false));
                reco.SaveAs([cpath 'colin_hires_RH_seg_RECOSM_sm500_05_inf125_05.srf']);
                reco.Smooth(500, 0.5, struct('show', false, 'distwsq', true));
                reco.Inflate(125, 0.5, struct('show', false));
                reco.SaveAs([cpath 'colin_hires_RH_seg_RECOSM_sm1000_05_inf250_05.srf']);

                % then restore unsmoothed coordinates
                reco.VertexCoordinate = recovc;
                reco.VertexNormal = recovn;
            end

            % ICBM-norm version
            if opts.icbmsrf && ...
                exist([cpath 'colin_hires_RH_seg_RECOSM_ICBMnorm.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_hires_RH_seg_RECOSM_ICBMnorm.srf ...');

                % copy object (so as not to overwrite coordinates)
                recoi = reco.CopyObject;

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - reco.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                recoi.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                recoi.RecalcNormals;

                % then re-create curvature information
                hcrv = recoi.CurvatureMap;
                recoi.VertexColor = [double(hcrv.Map.SMPData >= 0.002), ...
                    zeros(recoi.NrOfVertices, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % and save file
                recoi.AutoLinkedSRF = '';
                recoi.SaveAs([cpath 'colin_hires_RH_seg_RECOSM_ICBMnorm.srf']);

                % smooth, inflate, and save again
                recoi.Smooth(500, 0.5, struct('show', false, 'distwsq', true));
                recoi.Inflate(125, 0.5, struct('show', false));
                recoi.SaveAs([cpath 'colin_hires_RH_seg_RECOSM_ICBMnorm_sm500_05_inf125_05.srf']);
                recoi.Smooth(500, 0.5, struct('show', false, 'distwsq', true));
                recoi.Inflate(125, 0.5, struct('show', false));
                recoi.SaveAs([cpath 'colin_hires_RH_seg_RECOSM_ICBMnorm_sm1000_05_inf250_05.srf']);

                % clear object
                recoi.ClearObject;
                recoi = 0;
            end

        % still required
        elseif (opts.hiressph && ...
           (exist([cpath 'colin_hires_RH_SPH160k.srf'], 'file') ~= 2 || ...
            exist([cpath 'colin_midres_RH_SPH40k.srf'], 'file') ~= 2 || ...
            exist([cpath 'colin_lowres_RH_SPH10k.srf'], 'file') ~= 2)) || ...
           (opts.flatmaps && ...
            exist([cpath 'colin_hires_RH_seg_FLAT.srf'], 'file') ~= 2) || ...
           (opts.hiressrfb && ...
            exist([cpath 'colin_hires_RH_segback.vmr'], 'file') ~= 2)

            % load
            reco = xff([cpath 'colin_hires_RH_seg_RECOSM.srf']);
        end

        % hires SRF backprojection VMR
        if opts.hiressrfb && ...
            exist([cpath 'colin_hires_RH_segback.vmr'], 'file') ~= 2

            % information
            idisp('Creating colin_hires_RH_segback+hull.vmr');

            % create back projection
            reco.BackToVMR(btvopt);

            % create new hiresvmr
            hvmr = newhiresvmr;
            hvmr.VMRData = reco.VertexVMRData(77:512, 53:404, 77:440);
            reco.VertexVMRData(:) = [];
            hvmr.OffsetX = 76;
            hvmr.OffsetY = 52;
            hvmr.OffsetZ= 76;
            hvmr.SmoothData3D(1.25);
            hvmr.SetScalingWindow;
            hvmr.RunTimeVars.ScalingWindow = [33, 150];
            hvmr.RunTimeVars.ScalingWindowLim = [33, 300];
            hvmr.SaveAs([cpath 'colin_hires_RH_segback_hull.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.VMRData(floodfill3(hvmr.VMRData <= 50, 219, 181, 131)) = 150;
            hvmr.SmoothData3D(1.25, [], struct('range', [0, 125]));
            hvmr.SetScalingWindow;
            hvmr.RunTimeVars.ScalingWindow = [33, 100];
            hvmr.SaveAs([cpath 'colin_hires_RH_segback.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.ClearObject;
        end

        % combined object
        if opts.hiressrfb && ...
            exist([cpath 'colin_hires_LH_segback.vmr'], 'file') == 2 && ...
            exist([cpath 'colin_hires_RH_segback.vmr'], 'file') == 2 && ...
            exist([cpath 'colin_hires_BH_segback.vmr'], 'file') ~= 2

            % information
            idisp('Creating colin_hires_BH_segback.vmr');

            % load both and combine
            hvmr = xff([cpath 'colin_hires_LH_segback.vmr']);
            hvm2 = xff([cpath 'colin_hires_RH_segback.vmr']);
            hvmr.VMRData = max(hvmr.VMRData, hvm2.VMRData);
            hvm2.ClearObject;
            hvmr.SaveAs([cpath 'colin_hires_BH_segback.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.ClearObject;
        end

        % hires SRF backprojection VMR
        if opts.hiressrfb && ...
            opts.icbmsrf && ...
            exist([cpath 'colin_hires_RH_segback_ICBMnorm.vmr'], 'file') ~= 2 && ...
            exist([cpath 'colin_hires_RH_seg_RECOSM_ICBMnorm.srf'], 'file') == 2

            % information
            idisp('Creating colin_hires_RH_segback_ICBMnorm+hull.vmr');

            % create back projection
            recoi = xff([cpath 'colin_hires_RH_seg_RECOSM_ICBMnorm.srf']);
            recoi.BackToVMR(btvopt);

            % create new hiresvmr
            hvmr = newhiresvmr;
            hvmr.VMRData = recoi.VertexVMRData(77:512, 53:404, 77:440);
            recoi.ClearObject;
            hvmr.OffsetX = 76;
            hvmr.OffsetY = 52;
            hvmr.OffsetZ= 76;
            hvmr.SmoothData3D(1.25);
            hvmr.SetScalingWindow;
            hvmr.RunTimeVars.ScalingWindow = [33, 150];
            hvmr.RunTimeVars.ScalingWindowLim = [33, 300];
            hvmr.SaveAs([cpath 'colin_hires_RH_segback_ICBMnorm_hull.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.VMRData(floodfill3(hvmr.VMRData <= 47, 219, 181, 131)) = 150;
            hvmr.SmoothData3D(1.25, [], struct('range', [0, 125]));
            hvmr.SetScalingWindow;
            hvmr.RunTimeVars.ScalingWindow = [33, 100];
            hvmr.SaveAs([cpath 'colin_hires_RH_segback_ICBMnorm.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.ClearObject;
        end

        % combined object
        if opts.hiressrfb && ...
            exist([cpath 'colin_hires_LH_segback_ICBMnorm.vmr'], 'file') == 2 && ...
            exist([cpath 'colin_hires_RH_segback_ICBMnorm.vmr'], 'file') == 2 && ...
            exist([cpath 'colin_hires_BH_segback_ICBMnorm.vmr'], 'file') ~= 2

            % information
            idisp('Creating colin_hires_BH_segback_ICBMnorm.vmr');

            % load both and combine
            hvmr = xff([cpath 'colin_hires_LH_segback_ICBMnorm.vmr']);
            hvm2 = xff([cpath 'colin_hires_RH_segback_ICBMnorm.vmr']);
            hvmr.VMRData = max(hvmr.VMRData, hvm2.VMRData);
            hvm2.ClearObject;
            hvmr.SaveAs([cpath 'colin_hires_BH_segback_ICBMnorm.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.ClearObject;
        end

        % flatmaps
        if opts.flatmaps && ...
            exist([cpath 'colin_hires_RH_seg_FLAT.srf'], 'file') ~= 2

            % information
            idisp('Creating colin_hires_RH_seg_FLAT.srf ...');

            % copy object
            flat = reco.CopyObject;

            % fill color first
            flat.VertexColor(:) = 0;
            flat.VertexColor(:, 1) = ...
                ucc(1 + floor((1 / 2 ^ 28) .* double(rhcrd)));

            % then remove color from coordinates information
            rhcrd = mod(rhcrd, uint32(2 ^ 28));

            % and fill coordinates
            flat.VertexCoordinate(:, 1) = (1 / 64) .* ...
                floor((1 / 2 ^ 14) .* double(rhcrd));
            flat.VertexCoordinate(:, 2) = (1 / 64) .* ...
                double(mod(rhcrd, uint32(2 ^ 14)));
            flat.VertexCoordinate(:, 3) = 128;

            % recalc normals
            flat.RecalcNormals;

            % save and clear
            flat.SaveAs([cpath 'colin_hires_RH_seg_FLAT.srf']);
            flat.ClearObject;
            flat = 0;
        end

        % hi-res SPH version
        if opts.hiressph

            % re-sampled surface -> 160k RH
            if exist([cpath 'colin_hires_RH_SPH160k.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_hires_RH_SPH160k.srf ...');

                % generate sph's
                if isempty(sph{3})
                    sph{3} = spheresrf(140, 7);
                end

                % create TSM
                tsm = xff('new:tsm');
                nv = sph{3}.NrOfVertices;
                tsm.NrOfTargetVertices = nv;
                tsm.NrOfSourceVertices = reco.NrOfVertices;
                tsm.NrOfSourceTriangles = reco.NrOfTriangles;
                tsmd = double(seg_05mm.rh.tsm(1:nv));
                teld = floor(tsmd ./ (2 ^ 22));
                tsm.SourceTriangleOfTarget = mod(tsmd, 2 ^ 22);
                tsm.TriangleEdgeLengths = ...
                    0.05 .* ([floor(teld ./ 32), mod(teld, 32)] - 10);

                % re-sample surface (160k)
                res = reco.ApplyTSM(tsm, sph{3});
                tsm.ClearObject;

                % smooth (a bit)
                res.Smooth(15, 0.02, struct('show', false));

                % curvature map (coloring)
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.005), zeros(163842, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % save
                res.SaveAs([cpath 'colin_hires_RH_SPH160k.srf']);

            % load anyway
            elseif (opts.icbmsph && ...
                exist([cpath 'colin_hires_RH_SPH160k_ICBMnorm.srf'], 'file') ~= 2) || ...
               (opts.flatmaps && ...
                exist([cpath 'colin_hires_RH_SPH160k_FLAT.srf'], 'file') ~= 2)

                % load
                res = xff([cpath 'colin_hires_RH_SPH160k.srf']);
            end

            % flatmaps
            if opts.flatmaps && ...
                exist([cpath 'colin_hires_RH_SPH160k_FLAT.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_hires_RH_SPH160k_FLAT.srf ...');

                % copy object
                flat = res.CopyObject;

                % fill color first
                flat.VertexColor(:) = 0;
                flat.VertexColor(:, 1) = ...
                    ucc(1 + floor((1 / 2 ^ 28) .* double(rhcrh)));

                % then remove color from coordinates information
                rhcrh = mod(rhcrh, uint32(2 ^ 28));

                % and fill coordinates
                flat.VertexCoordinate(:, 1) = (1 / 64) .* ...
                    floor((1 / 2 ^ 14) .* double(rhcrh));
                flat.VertexCoordinate(:, 2) = (1 / 64) .* ...
                    double(mod(rhcrh, uint32(2 ^ 14)));
                flat.VertexCoordinate(:, 3) = 128;

                % recalc normals
                flat.RecalcNormals;

                % save and clear
                flat.SaveAs([cpath 'colin_hires_RH_SPH160k_FLAT.srf']);
                flat.ClearObject;
                flat = 0;
            end

            % ICBM-normalized version
            if opts.icbmsph && ...
                exist([cpath 'colin_hires_RH_SPH160k_ICBMnorm.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_hires_RH_SPH160k_ICBMnorm.srf ...');

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - res.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                res.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                res.RecalcNormals;

                % then re-create curvature information
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.005), zeros(163842, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % and save file
                res.SaveAs([cpath 'colin_hires_RH_SPH160k_ICBMnorm.srf']);
            end

            % clear object
            clearxffobjects({res});
            res = 0;

            % re-sampled surface -> 40k RH
            if exist([cpath 'colin_midres_RH_SPH40k.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_midres_RH_SPH40k.srf ...');

                % generate sph's
                if isempty(sph{2})
                    sph{2} = spheresrf(140, 6);
                end

                % create TSM
                tsm = xff('new:tsm');
                nv = sph{2}.NrOfVertices;
                tsm.NrOfTargetVertices = nv;
                tsm.NrOfSourceVertices = reco.NrOfVertices;
                tsm.NrOfSourceTriangles = reco.NrOfTriangles;
                tsmd = double(seg_05mm.rh.tsm(1:nv));
                teld = floor(tsmd ./ (2 ^ 22));
                tsm.SourceTriangleOfTarget = mod(tsmd, 2 ^ 22);
                tsm.TriangleEdgeLengths = ...
                    0.05 .* ([floor(teld ./ 32), mod(teld, 32)] - 10);

                % re-sample surface (40k)
                res = reco.ApplyTSM(tsm, sph{2});
                tsm.ClearObject;

                % smooth (a bit)
                res.Smooth(20, 0.015, struct('show', false));

                % curvature map (coloring)
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.01), zeros(40962, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % save
                res.SaveAs([cpath 'colin_midres_RH_SPH40k.srf']);

            % load anyway
            elseif (opts.icbmsph && ...
                exist([cpath 'colin_midres_RH_SPH40k_ICBMnorm.srf'], 'file') ~= 2) || ...
               (opts.flatmaps && ...
                exist([cpath 'colin_midres_RH_SPH40k_FLAT.srf'], 'file') ~= 2)

                % load
                res = xff([cpath 'colin_midres_RH_SPH40k.srf']);
            end

            % flatmaps
            if opts.flatmaps && ...
                exist([cpath 'colin_midres_RH_SPH40k_FLAT.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_midres_RH_SPH40k_FLAT.srf ...');

                % copy object
                flat = res.CopyObject;

                % fill color first
                flat.VertexColor(:) = 0;
                flat.VertexColor(:, 1) = ...
                    ucc(1 + floor((1 / 2 ^ 28) .* double(rhcrm)));

                % then remove color from coordinates information
                rhcrm = mod(rhcrm, uint32(2 ^ 28));

                % and fill coordinates
                flat.VertexCoordinate(:, 1) = (1 / 64) .* ...
                    floor((1 / 2 ^ 14) .* double(rhcrm));
                flat.VertexCoordinate(:, 2) = (1 / 64) .* ...
                    double(mod(rhcrm, uint32(2 ^ 14)));
                flat.VertexCoordinate(:, 3) = 128;

                % recalc normals
                flat.RecalcNormals;

                % save and clear
                flat.SaveAs([cpath 'colin_midres_RH_SPH40k_FLAT.srf']);
                flat.ClearObject;
                flat = 0;
            end

            % ICBM-normalized version
            if opts.icbmsph && ...
                exist([cpath 'colin_midres_RH_SPH40k_ICBMnorm.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_midres_RH_SPH40k_ICBMnorm.srf ...');

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - res.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                res.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                res.RecalcNormals;

                % then re-create curvature information
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.01), zeros(40962, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % and save file
                res.SaveAs([cpath 'colin_midres_RH_SPH40k_ICBMnorm.srf']);
            end

            % clear object
            clearxffobjects({res});
            res = 0;

            % re-sampled surface -> 10k RH
            if exist([cpath 'colin_lowres_RH_SPH10k.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_lowres_RH_SPH10k.srf ...');

                % generate sph's
                if isempty(sph{1})
                    sph{1} = spheresrf(140, 5);
                end

                % create TSM
                tsm = xff('new:tsm');
                nv = sph{1}.NrOfVertices;
                tsm.NrOfTargetVertices = nv;
                tsm.NrOfSourceVertices = reco.NrOfVertices;
                tsm.NrOfSourceTriangles = reco.NrOfTriangles;
                tsmd = double(seg_05mm.rh.tsm(1:nv));
                teld = floor(tsmd ./ (2 ^ 22));
                tsm.SourceTriangleOfTarget = mod(tsmd, 2 ^ 22);
                tsm.TriangleEdgeLengths = ...
                    0.05 .* ([floor(teld ./ 32), mod(teld, 32)] - 10);

                % re-sample surface (10k)
                res = reco.ApplyTSM(tsm, sph{1});
                tsm.ClearObject;

                % smooth (a bit)
                res.Smooth(25, 0.01, struct('show', false));

                % curvature map (coloring)
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.02), zeros(10242, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % save
                res.SaveAs([cpath 'colin_lowres_RH_SPH10k.srf']);

            % load anyway
            elseif opts.icbmsph && ...
                exist([cpath 'colin_lowres_RH_SPH10k_ICBMnorm.srf'], 'file') ~= 2

                % load
                res = xff([cpath 'colin_lowres_RH_SPH10k.srf']);
            end

            % ICBM-normalized version
            if opts.icbmsph && ...
                exist([cpath 'colin_lowres_RH_SPH10k_ICBMnorm.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_lowres_RH_SPH10k_ICBMnorm.srf ...');

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - res.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                res.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                res.RecalcNormals;

                % then re-create curvature information
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.02), zeros(10242, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % and save file
                res.SaveAs([cpath 'colin_lowres_RH_SPH10k_ICBMnorm.srf']);
            end

            % clear object
            clearxffobjects({res});
            res = 0;

        end

        % BH files
        if opts.dhsrf && ...
            exist([cpath 'colin_hires_BH_seg.vmr'], 'file') ~= 2

            % information
            idisp('Creating colin_hires_BH_seg.vmr ...');

            % load LH-VMR + RH-VMR
            if isxff(vmr, true)
                vmr.ClearObject;
            end
            vmr = xff([cpath 'colin_hires_LH_seg.vmr']);
            pvmr = xff([cpath 'colin_hires_RH_seg.vmr']);

            % reframe parts
            vmr.Reframe([0, 0, 0; 511, 511, 511]);
            pvmr.Reframe([0, 0, 0; 511, 511, 511]);

            % put together
            vmr.VMRData = max(vmr.VMRData, pvmr.VMRData);

            % clear second object
            pvmr.ClearObject;
            pvmr = 0;

            % add missing slice
            bhsiz = size(seg_05mm.bh.vox);
            bhoff = double(seg_05mm.bh.off);
            mslice = vmr.VMRData(bhoff(1)+1:bhoff(1)+bhsiz(1), bhoff(2)+1:bhoff(2)+bhsiz(2), 257);
            mslice(seg_05mm.bh.vox) = 240;
            vmr.VMRData(bhoff(1)+1:bhoff(1)+bhsiz(1), bhoff(2)+1:bhoff(2)+bhsiz(2), 257) = mslice;

            % reframe to smaller size again
            vmr.MinBox;

            % then add other voxels
            vmr.VMRData(seg_05mm.bh.xvx) = 240;

            % save
            if opts.hiresseg
                vmr.SaveAs([cpath 'colin_hires_BH_seg.vmr']);
            end

        % other wise
        elseif opts.dhsrf

            % load VMR
            if isxff(vmr, true)
                vmr.ClearObject;
            end
            vmr = xff([cpath 'colin_hires_BH_seg.vmr']);
        end

        % SRF file?
        if opts.dhsrf && ...
            exist([cpath 'colin_hires_BH_seg_RECOSM.srf'], 'file') ~= 2

            % information
            idisp('Creating colin_hires_BH_seg_RECOSM.srf ...');

            % reconstruct
            reco = vmr.DBReco;
            reco.ConvexRGBA = colin_sph.ConvexRGBA;
            reco.ConcaveRGBA = colin_sph.ConcaveRGBA;

            % smooth (a bit)
            reco.Smooth(100, 0.1, struct('show', false));

            % save
            if opts.hiressrf

                % first, curvature map (coloring)
                hcrv = reco.CurvatureMap;
                reco.VertexColor = [double(hcrv.Map.SMPData >= 0.002), ...
                    zeros(reco.NrOfVertices, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % then save
                reco.SaveAs([cpath 'colin_hires_BH_seg_RECOSM.srf']);

                % smooth, inflate, and save again
                recovc = reco.VertexCoordinate;
                recovn = reco.VertexNormal;
                reco.Smooth(500, 0.5, struct('show', false, 'distwsq', true));
                reco.Inflate(125, 0.5, struct('show', false));
                reco.SaveAs([cpath 'colin_hires_BH_seg_RECOSM_sm500_05_inf125_05.srf']);
                reco.Smooth(500, 0.5, struct('show', false, 'distwsq', true));
                reco.Inflate(125, 0.5, struct('show', false));
                reco.SaveAs([cpath 'colin_hires_BH_seg_RECOSM_sm1000_05_inf250_05.srf']);

                % then restore unsmoothed coordinates
                reco.VertexCoordinate = recovc;
                reco.VertexNormal = recovn;
            end

            % ICBM-norm version
            if opts.icbmsrf  && ...
                opts.dhsrf && ...
                exist([cpath 'colin_hires_BH_seg_RECOSM_ICBMnorm.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_hires_BH_seg_RECOSM_ICBMnorm.srf ...');

                % make copy (so as not to overwrite coordinates!)
                recoi = reco.CopyObject;

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - reco.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                recoi.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                recoi.RecalcNormals;

                % then re-create curvature information
                hcrv = recoi.CurvatureMap;
                recoi.VertexColor = [double(hcrv.Map.SMPData >= 0.002), ...
                    zeros(recoi.NrOfVertices, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % and save file
                recoi.AutoLinkedSRF = '';
                recoi.SaveAs([cpath 'colin_hires_BH_seg_RECOSM_ICBMnorm.srf']);
                recoi.Smooth(500, 0.5, struct('show', false, 'distwsq', true));
                recoi.Inflate(125, 0.5, struct('show', false));
                recoi.SaveAs([cpath 'colin_hires_BH_seg_RECOSM_ICBMnorm_sm500_05_inf125_05.srf']);
                recoi.Smooth(500, 0.5, struct('show', false, 'distwsq', true));
                recoi.Inflate(125, 0.5, struct('show', false));
                recoi.SaveAs([cpath 'colin_hires_BH_seg_RECOSM_ICBMnorm_sm1000_05_inf250_05.srf']);

                % clear again
                recoi.ClearObject;
                recoi = 0;
            end

        % still required
        elseif opts.dhsrf

            % load
            reco = xff([cpath 'colin_hires_BH_seg_RECOSM.srf']);
        end

        % hires SRF backprojection VMR
        if opts.hiressrfb && ...
            opts.dhsrf && ...
            exist([cpath 'colin_hires_BH_segback.vmr'], 'file') ~= 2

            % information
            idisp('Creating colin_hires_BH_segback+hull.vmr');

            % create back projection
            reco.BackToVMR(btvopt);

            % create new hiresvmr
            hvmr = newhiresvmr;
            hvmr.VMRData = reco.VertexVMRData(77:512, 53:404, 77:440);
            reco.VertexVMRData(:) = [];
            hvmr.OffsetX = 76;
            hvmr.OffsetY = 52;
            hvmr.OffsetZ= 76;
            hvmr.SmoothData3D(1.25);
            hvmr.SetScalingWindow;
            hvmr.RunTimeVars.ScalingWindow = [33, 150];
            hvmr.RunTimeVars.ScalingWindowLim = [33, 300];
            hvmr.SaveAs([cpath 'colin_hires_BH_segback_hull.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.VMRData(floodfill3(hvmr.VMRData <= 50, 219, 169, 233)) = 150;
            hvmr.SmoothData3D(1.25, [], struct('range', [0, 125]));
            hvmr.SetScalingWindow;
            hvmr.RunTimeVars.ScalingWindow = [33, 100];
            hvmr.SaveAs([cpath 'colin_hires_BH_segback.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.ClearObject;
        end

        % hires SRF backprojection VMR
        if opts.hiressrfb && ...
            opts.icbmsrf && ...
            opts.dhsrf && ...
            exist([cpath 'colin_hires_BH_segback_ICBMnorm.vmr'], 'file') ~= 2 && ...
            exist([cpath 'colin_hires_BH_seg_RECOSM_ICBMnorm.srf'], 'file') == 2

            % information
            idisp('Creating colin_hires_BH_segback_ICBMnorm+hull.vmr');

            % create back projection
            recoi = xff([cpath 'colin_hires_BH_seg_RECOSM_ICBMnorm.srf']);
            recoi.BackToVMR(btvopt);

            % create new hiresvmr
            hvmr = newhiresvmr;
            hvmr.VMRData = recoi.VertexVMRData(77:512, 53:404, 77:440);
            recoi.ClearObject;
            recoi = 0;
            hvmr.OffsetX = 76;
            hvmr.OffsetY = 52;
            hvmr.OffsetZ= 76;
            hvmr.SmoothData3D(1.25);
            hvmr.SetScalingWindow;
            hvmr.RunTimeVars.ScalingWindow = [33, 150];
            hvmr.RunTimeVars.ScalingWindowLim = [33, 300];
            hvmr.SaveAs([cpath 'colin_hires_BH_segback_ICBMnorm_hull.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.VMRData(floodfill3(hvmr.VMRData <= 50, 219, 169, 233)) = 150;
            hvmr.SmoothData3D(1.25, [], struct('range', [0, 125]));
            hvmr.SetScalingWindow;
            hvmr.RunTimeVars.ScalingWindow = [33, 100];
            hvmr.SaveAs([cpath 'colin_hires_BH_segback_ICBMnorm.vmr']);
            hvmr.SaveRunTimeVars;
            hvmr.ClearObject;
        end

        % hi-res SPH version
        if opts.hiressph && ...
            opts.dhsrf

            % re-sampled surface -> 160k BH
            if exist([cpath 'colin_hires_BH_SPH160k.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_hires_BH_SPH160k.srf ...');

                % generate sph's
                if isempty(sph{3})
                    sph{3} = spheresrf(140, 7);
                end

                % create TSM
                tsm = xff('new:tsm');
                nv = sph{3}.NrOfVertices;
                tsm.NrOfTargetVertices = nv;
                tsm.NrOfSourceVertices = reco.NrOfVertices;
                tsm.NrOfSourceTriangles = reco.NrOfTriangles;

                % get data
                tsmd = double(seg_05mm.bh.tsm(1:nv));
                teld = floor(tsmd ./ (2 ^ 22));
                tsm.SourceTriangleOfTarget = mod(tsmd, 2 ^ 22);
                tsm.TriangleEdgeLengths = ...
                    0.05 .* ([floor(teld ./ 32), mod(teld, 32)] - 10);

                % re-sample surface (160k)
                res = reco.ApplyTSM(tsm, sph{3});
                tsm.ClearObject;

                % smooth (a bit)
                res.Smooth(15, 0.02, struct('show', false));

                % curvature map (coloring)
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.005), zeros(163842, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % save
                res.SaveAs([cpath 'colin_hires_BH_SPH160k.srf']);

            % load anyway
            elseif opts.icbmsph && ...
                opts.dhsrf && ...
                exist([cpath 'colin_hires_BH_SPH160k_ICBMnorm.srf'], 'file') ~= 2

                % load
                res = xff([cpath 'colin_hires_BH_SPH160k.srf']);
            end

            % ICBM-normalized version
            if opts.icbmsph && ...
                opts.dhsrf && ...
                exist([cpath 'colin_hires_BH_SPH160k_ICBMnorm.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_hires_BH_SPH160k_ICBMnorm.srf ...');

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - res.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                res.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                res.RecalcNormals;

                % then re-create curvature information
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.005), zeros(163842, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % and save file
                res.SaveAs([cpath 'colin_hires_BH_SPH160k_ICBMnorm.srf']);
            end

            % clear object
            clearxffobjects({res});
            res = 0;

            % re-sampled surface -> 40k BH
            if opts.dhsrf && ...
                exist([cpath 'colin_midres_BH_SPH40k.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_midres_BH_SPH40k.srf ...');

                % generate sph's
                if isempty(sph{2})
                    sph{2} = spheresrf(140, 6);
                end

                % create TSM
                tsm = xff('new:tsm');
                nv = sph{2}.NrOfVertices;
                tsm.NrOfTargetVertices = nv;
                tsm.NrOfSourceVertices = reco.NrOfVertices;
                tsm.NrOfSourceTriangles = reco.NrOfTriangles;
                tsmd = double(seg_05mm.bh.tsm(1:nv));
                teld = floor(tsmd ./ (2 ^ 22));
                tsm.SourceTriangleOfTarget = mod(tsmd, 2 ^ 22);
                tsm.TriangleEdgeLengths = ...
                    0.05 .* ([floor(teld ./ 32), mod(teld, 32)] - 10);

                % re-sample surface (40k)
                res = reco.ApplyTSM(tsm, sph{2});
                tsm.ClearObject;

                % smooth (a bit)
                res.Smooth(20, 0.015, struct('show', false));

                % curvature map (coloring)
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.01), zeros(40962, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % save
                res.SaveAs([cpath 'colin_midres_BH_SPH40k.srf']);

            % load anyway
            elseif opts.icbmsph && ...
                opts.dhsrf && ...
                exist([cpath 'colin_midres_BH_SPH40k_ICBMnorm.srf'], 'file') ~= 2

                % load
                res = xff([cpath 'colin_midres_BH_SPH40k.srf']);
            end

            % ICBM-normalized version
            if opts.icbmsph && ...
                opts.dhsrf && ...
                exist([cpath 'colin_midres_BH_SPH40k_ICBMnorm.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_midres_BH_SPH40k_ICBMnorm.srf ...');

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - res.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                res.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                res.RecalcNormals;

                % then re-create curvature information
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.01), zeros(40962, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % and save file
                res.SaveAs([cpath 'colin_midres_BH_SPH40k_ICBMnorm.srf']);
            end

            % clear object
            clearxffobjects({res});
            res = 0;

            % re-sampled surface -> 10k BH
            if opts.dhsrf && ...
                exist([cpath 'colin_lowres_BH_SPH10k.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_lowres_BH_SPH10k.srf ...');

                % generate sph's
                if isempty(sph{1})
                    sph{1} = spheresrf(140, 5);
                end

                % create TSM
                tsm = xff('new:tsm');
                nv = sph{1}.NrOfVertices;
                tsm.NrOfTargetVertices = nv;
                tsm.NrOfSourceVertices = reco.NrOfVertices;
                tsm.NrOfSourceTriangles = reco.NrOfTriangles;
                tsmd = double(seg_05mm.bh.tsm(1:nv));
                teld = floor(tsmd ./ (2 ^ 22));
                tsm.SourceTriangleOfTarget = mod(tsmd, 2 ^ 22);
                tsm.TriangleEdgeLengths = ...
                    0.05 .* ([floor(teld ./ 32), mod(teld, 32)] - 10);

                % re-sample surface (10k)
                res = reco.ApplyTSM(tsm, sph{1});
                tsm.ClearObject;

                % smooth (a bit)
                res.Smooth(25, 0.01, struct('show', false));

                % curvature map (coloring)
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.02), zeros(10242, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % save
                res.SaveAs([cpath 'colin_lowres_BH_SPH10k.srf']);

            % load anyway
            elseif opts.icbmsph && ...
                opts.dhsrf && ...
                exist([cpath 'colin_lowres_BH_SPH10k_ICBMnorm.srf'], 'file') ~= 2

                % load
                res = xff([cpath 'colin_lowres_BH_SPH10k.srf']);
            end

            % ICBM-normalized version
            if opts.icbmsph && ...
                opts.dhsrf && ...
                exist([cpath 'colin_lowres_BH_SPH10k_ICBMnorm.srf'], 'file') ~= 2

                % information
                idisp('Creating colin_lowres_BH_SPH10k_ICBMnorm.srf ...');

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - res.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                res.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                res.RecalcNormals;

                % then re-create curvature information
                hcrv = res.CurvatureMap;
                res.VertexColor = [double(hcrv.Map.SMPData >= 0.02), zeros(10242, 3)];
                hcrv.ClearObject;
                hcrv = 0;

                % and save file
                res.SaveAs([cpath 'colin_lowres_BH_SPH10k_ICBMnorm.srf']);
            end

            % clear object
            clearxffobjects({res});
            res = 0;
        end

        % clear object
        clearxffobjects({reco});
        reco = 0;
    end
catch ne_eo;
    warning( ...
        'neuroelf:InternalError', ...
        'Error creating hi-res segmentation based file: %s.', ...
        ne_eo.message ...
    );
end

try

    % ICBM-norm version VMR
    if opts.icbmvmr && ...
        exist([cpath 'colin_ICBMnorm.vmr'], 'file') ~= 2

        % information
        idisp('Creating colin_ICBMnorm.vmr ...');

        % reload colin
        colin = neuroelf_file('c', 'colin.vmr');
        colin.VMRData16 = [];

        % create transformation argument for VMR space
        itrf = [0, -1, 0, 129; 0, 0, -1, 129; -1, 0, 0, 129; 0, 0, 0, 1];

        % resample VMRData
        nd = uint8(limitrangec(flexinterpn_method(colin.VMRData(:, :, :), ...
            applyspmsnc([Inf; 128; -1; -127] * ones(1, 3), ...
            sn.Tr, sn.VG(1).dim, inv(sn.VG(1).mat), itrf*sn.VF(1).mat*sn.Affine), ...
            'lanczos3'), 0, 225));
        icbm = colin.CopyObject;
        icbm.VMRData = permute(reshape(nd, [256, 256, 256]), [2, 3, 1]);
        icbm.SaveAs([cpath 'colin_ICBMnorm.vmr']);

    % stil required
    elseif opts.bmaskvmr && ...
        opts.icbmvmr

        % load normalized VMR
        icbm = xff([cpath 'colin_ICBMnorm.vmr']);
    end

    % colin brain in ICBM normalization
    if opts.bmaskvmr && ...
        opts.icbmvmr && ...
        exist([cpath 'colin_brain_ICBMnorm.vmr'], 'file') ~= 2

        % information
        idisp('Creating colin_brain_ICBMnorm.vmr ...');

        % load mask
        cbm = load([cpath 'brainmask.mat']);
        cbm = cbm.cbm;

        % smooth mask
        cbm = smoothdata3(single(cbm), [2, 2, 2]);

        % apply normalization
        itrf = [0, -1, 0, 129; 0, 0, -1, 129; -1, 0, 0, 129; 0, 0, 0, 1];
        cbm = permute(reshape(flexinterpn_method(cbm, ...
            applyspmsnc([Inf; 128; -1; -127] * ones(1, 3), ...
            sn.Tr, sn.VG(1).dim, inv(sn.VG(1).mat), itrf*sn.VF(1).mat*sn.Affine), ...
            'lanczos3') >= 0.5, [256, 256, 256]), [2, 3, 1]);

        % get data for sampling
        vd = max(uint8(1), icbm.VMRData(:, :, :));

        % mask brain
        icbm.VMRData = uint8(cbm) .* vd;

        % compute fringes coordinates
        cbmx = dilate3d(cbm);
        fringe = find(cbmx(:) & ~cbm(:));
        [frx, fry, frz] = ind2sub(size(vd), fringe);

        % those will be smoothed versions of the masked data
        icbm.VMRData(fringe) = ...
            flexinterpn(icbm.VMRData, [frx, fry, frz], smoothkern(1.5), 1, 0);

        % then compute inner fringe coordinates
        cbmx = erode3d(cbm);
        fringe = find(cbm(:) & ~cbmx(:));
        [frx, fry, frz] = ind2sub(size(vd), fringe);

        % those are slightly smoothed versions of the current data
        icbm.VMRData(fringe) = ...
            flexinterpn(icbm.VMRData, [frx, fry, frz], smoothkern(1), 1, 0);

        % save masked version
        icbm.SaveAs([cpath 'colin_brain_ICBMnorm.vmr']);

    % still required for WM masked version
    elseif opts.wmaskvmr && ...
        opts.icbmvmr && ...
        exist([cpath 'colin_brain_WMmasked_ICBMnorm.vmr'], 'file') ~= 2

        % load ICBMnorm brain
        icbm = xff([cpath 'colin_brain_ICBMnorm.vmr']);
    end

    % WM-masked version of ICBM brain
    if opts.wmaskvmr && ...
        opts.icbmvmr && ...
        exist([cpath 'colin_brain_WMmasked_ICBMnorm.vmr'], 'file') ~= 2

        % information
        idisp('Creating WM-masked version of colin_brain_ICBMnorm.vmr -> colin_brain_WMmasked_ICBMnorm.vmr ...');

        % get colin data clustered
        [cols, colv] = clustercoordsc(icbm.VMRData(:, :, :) >= 165, 1, 100);
        colv = (colv == maxpos(cols));
        [colv, colvo, colvs] = minarray(colv, 0, 2, 1);
        colv = (smoothdata3(single(colv), [2, 2, 2]) > 0.5);

        % 2 passes of erode and smooth
        for spc = 1:2
            colv = (smoothdata3(single(erode3d(colv)), [2, 2, 2]) > 0.25);
        end

        % set to 0 in data
        icbm.VMRData = icbm.VMRData(:, :, :);
        icbm.VMRData(colvo(1):colvo(1)+colvs(1)-1, ...
            colvo(2):colvo(2)+colvs(2)-1, colvo(3):colvo(3)+colvs(3)-1) = ...
            icbm.VMRData(colvo(1):colvo(1)+colvs(1)-1, ...
            colvo(2):colvo(2)+colvs(2)-1, colvo(3):colvo(3)+colvs(3)-1) .* ...
            uint8(~colv);

        % set ventricles to zero also!
        icbm.VMRData(smoothdata3(single(icbm.VMRData > 0 & icbm.VMRData < 60), ...
            [2.5, 2.5, 2.5]) > 0.5) = 0;

        % and establish a minimum value of 30
        icbm.VMRData(icbm.VMRData > 0 & icbm.VMRData < 30) = 30;

        % smooth 0/1 values with a 3mm kernel
        icbm.SmoothData3D(3, [Inf, Inf, Inf; -80, -112, -76; 1, 1, 1; 80, 80, 88], ...
            struct('range', [0, 1]));

        % save as
        icbm.SaveAs([cpath 'colin_brain_WMmasked_ICBMnorm.vmr']);

    % still required
    elseif opts.rendskull && ...
        opts.icbmvmr && ...
        exist([cpath 'colin_brain_rendskull_ICBMnorm.vmr'], 'file') ~= 2

        % load colin
        icbm = xff([cpath 'colin_brain_WMmasked_ICBMnorm.vmr']);
    end

    % create rendskull VMR
    if opts.rendskull && ...
        opts.icbmvmr && ...
        exist([cpath 'colin_brain_rendskull_ICBMnorm.vmr'], 'file') ~= 2

        % information
        idisp('Creating render-skull version of colin_brain_WMmasked_ICBMnorm.vmr -> colin_brain_rendskull_ICBMnorm.vmr ...');

        % create copy
        skull = icbm.CopyObject;

        % load original colin
        vmr = xff([cpath 'colin_ICBMnorm.vmr']);

        % smooth initial content (lower values)
        skull.SmoothData3D(1.5, [Inf; -127; 1; 127] * ones(1, 3), struct('range', [1, 55]));

        % mix data
        skull.VMRData(:, :, :) = max( ...
            uint8(round(0.67 .* double(skull.VMRData))), ...
            uint8(round(0.05 .* (double(vmr.VMRData) + 8))));

        % clear original VMR
        vmr.ClearObject;
        vmr = 0;

        % smooth ever so slightly (everything)
        skull.SmoothData3D(1.25, [Inf; -127; 1; 127] * ones(1, 3));

        % set/update scaling window
        skull.SetScalingWindow([0, 600], true);
        skull.RunTimeVars.ScalingWindow = [3, 450];

        % save and clear
        skull.SaveAs([cpath 'colin_brain_rendskull_ICBMnorm.vmr']);
        skull.SaveRunTimeVars;
        skull.ClearObject;
        skull = 0;
    end

    % clear colin
    if isxff(icbm, true)
        icbm.ClearObject;
    end
    icbm = 0;

catch ne_eo;
    warning( ...
        'neuroelf:xffError', ...
        'Error creating colin_brain.vmr: %s.', ...
        ne_eo.message ...
    );
end

try

    % ICBM-version of talairach.nii
    t = {[], []};
    if opts.icbmtal && ...
        exist([tpath, 'talairach_ICBMnorm.nii'], 'file') ~= 2

        % information
        idisp('Creating talairach_ICBMnorm.nii ...');

        % open talairach nii
        t{1} = xff([tpath 'talairach.nii']);

        % create new HDR/NII object
        t{2} = xff('new:hdr');
        t{2}.FileMagic = 'n+1';
        t{2}.NIIFileType = 2;
        t{2}.ImgDim.Dim = [3, 181, 217, 181, 1, 1, 1, 1];
        t{2}.ImgDim.DataType = 4;
        t{2}.ImgDim.BitsPerPixel = 16;
        t{2}.ImgDim.PixSpacing(:) = 1;
        t{2}.DataHist.Description = 'Talairach ICBM-transformed';
        t{2}.DataHist.NIftI1.QFormCode = 1;
        t{2}.DataHist.NIftI1.SFormCode = 1;
        t{2}.DataHist.NIftI1.QuaternionB = 0;
        t{2}.DataHist.NIftI1.QuaternionC = 0;
        t{2}.DataHist.NIftI1.QuaternionD = 0;
        t{2}.DataHist.NIftI1.QuatOffsetX = -90;
        t{2}.DataHist.NIftI1.QuatOffsetY = -126;
        t{2}.DataHist.NIftI1.QuatOffsetZ = -72;
        t{2}.DataHist.NIftI1.AffineTransX = [1, 0, 0,  -90];
        t{2}.DataHist.NIftI1.AffineTransY = [0, 1, 0, -126];
        t{2}.DataHist.NIftI1.AffineTransZ = [0, 0, 1,  -72];

        % sample data
        t{2}.VoxelData = int16(t{1}.SampleTalBox(struct( ...
            'BBox', [-90, -126, -72; 90, 90, 108], 'ResXYZ', [1, 1, 1]), ...
            1, 'nearest', [], load([tpath 'talairach_seg_sn.mat'])));

        % transform labels
        tallabels = splittocell( ...
            strrep(char(t{1}.IntermedData(13:end)), '.', ','), ...
            char(10), 1);
        tallabels(1) = [];
        tallabels(end) = [];
        t{2}.RunTimeVars.AutoSave = true;
        t{2}.RunTimeVars.TALLabels = tallabels;

        % save data
        t{2}.SaveAs([tpath 'talairach_ICBMnorm.nii']);

        % clear data
        clearxffobjects(t);
    end
catch ne_eo;
    warning('neuroelf:xffError', ...
        'Error creating talairach_ICBMnorm.nii: %s.', ...
        ne_eo.message ...
    );
    clearxffobjects(t);
end

try

    % hi-res VMR
    if opts.hiresvmr && ...
        exist([cpath 'colin_hires.vmr'], 'file') ~= 2

        % clear all objects!
        clear xff;

        % information
        idisp('Creating colin_hires.vmr ...');

        % create hi-res version (e.g. for re-segmentation)
        if isxff(colin, true)
            colin.ClearObject;
            colin = 0;
        end
        colin = importvmrfromanalyze( ...
            [cpath 'colin.hdr'], 'lanczos8', [0.005, 0.999], 0.5, ...
            [], [38, 26, 38; 256, 202, 220]);

        % save as
        colin.SaveAs([cpath 'colin_hires.vmr']);
    end

    % hires ICBM
    if opts.hiresvmr && ...
        opts.icbmvmr && ...
        exist([cpath 'colin_hires_ICBMnorm.vmr'], 'file') ~= 2

        % clear old objects
        clear xff;

        % information
        idisp('Creating colin_hires_ICBMnorm.vmr ...');

        % create hi-res version (e.g. for re-segmentation)
        if isxff(colin, true)
            colin.ClearObject;
            colin = 0;
        end
        colin = importvmrfromanalyze( ...
            [cpath 'colin_ICBMnorm.hdr'], 'lanczos8', [0.005, 0.999], 0.5, ...
            [], [38, 26, 38; 256, 202, 220]);

        % save as
        colin.SaveAs([cpath 'colin_hires_ICBMnorm.vmr']);
    end

    % brain-masked versions
    if opts.bmaskvmr

        % hi-res VMR
        if opts.hiresvmr && ...
            exist([cpath 'colin_hires_brain.vmr'], 'file') ~= 2

            % clear all objects!
            clear xff;

            % information
            idisp('Creating colin_hires_brain.vmr ...');

            % load required objects
            if isxff(colin, true)
                colin.ClearObject;
                colin = 0;
            end
            colin = xff([cpath 'colin_hires.vmr']);
            colindata = colin.VMRData(:, :, :);
            colinmsk = xff([cpath 'colin_brain.vmr']);
            msk = uint8(colinmsk.VMRData(:, :, :) > 10);
            colinmsk.ClearObject;

            % get part of mask we need
            msk = msk(39:256, 27:202, 39:220);

            % apply slice by slice
            smk = smoothkern(16, 0.01);
            for sc = 1:size(colindata, 1)
                msksl = flexinterpn(msk, [Inf, Inf, Inf; ...
                    0.5 * sc, 1, 1; 1, 0.5, 0.5; 0.5 * (sc + 4), 176.5, 182.5]);
                msksl = flexinterpn(msksl, [Inf, Inf, Inf; ...
                    2, 1, 1; 1, 1, 1; 2, size(msksl, 2), size(msksl, 3)], smk, 4);
                colindata(sc, :, :) = ...
                    uint8(double(colindata(sc, :, :)) .* msksl);
            end

            % save as
            colin.VMRData = colindata;
            colin.VMRData16 = [];
            colin.SaveAs([cpath 'colin_hires_brain.vmr']);
        end

        % hires ICBM
        if opts.hiresvmr && ...
            opts.icbmvmr && ...
            exist([cpath 'colin_hires_brain_ICBMnorm.vmr'], 'file') ~= 2

            % clear old objects
            clear xff;

            % information
            idisp('Creating colin_hires_brain_ICBMnorm.vmr ...');

            % load required objects
            colin = xff([cpath 'colin_hires_ICBMnorm.vmr']);
            colindata = colin.VMRData(:, :, :);
            colinmsk = xff([cpath 'colin_brain_ICBMnorm.vmr']);
            msk = uint8(colinmsk.VMRData(:, :, :) > 10);
            colinmsk.ClearObject;

            % get part of mask we need
            msk = msk(39:256, 27:202, 39:220);

            % apply slice by slice
            smk = smoothkern(16, 0.01);
            for sc = 1:size(colindata, 1)
                msksl = flexinterpn(msk, [Inf, Inf, Inf; ...
                    0.5 * sc, 1, 1; 1, 0.5, 0.5; 0.5 * (sc + 4), 176.5, 182.5]);
                msksl = flexinterpn(msksl, [Inf, Inf, Inf; ...
                    2, 1, 1; 1, 1, 1; 2, size(msksl, 2), size(msksl, 3)], smk, 4);
                colindata(sc, :, :) = ...
                    uint8(double(colindata(sc, :, :)) .* msksl);
            end

            % save as
            colin.VMRData = colindata;
            colin.VMRData16 = [];
            colin.SaveAs([cpath 'colin_hires_brain_ICBMnorm.vmr']);
        end
    end
catch ne_eo;
    warning( ...
        'neuroelf:InternalError', ...
        'Error creating hi-res VMR version: %s.', ...
        ne_eo.message ...
    );
end

% sub-cortical structures
try
    if opts.subcort && ...
        exist([cpath 'colin_subcort.srf'], 'file') ~= 2

        % display information
        idisp('Creating subcortical surfaces ...');

        % load required file
        chdr = neuroelf_file('c', 'subcortical.nii');
        chv = chdr.VoxelData(:, :, :);
        chdr.ClearObject;

        % for each of the surfaces
        scodes = [10:13, 16:18, 26, 49:54, 58];
        snames = { ...
            'lThal', 'lCaud', 'lPut', 'lGPal', 'BStem', 'lHipp', 'lAmyg', 'lNAcc', ...
            'rThal', 'rCaud', 'rPut', 'rGPal', 'rHipp', 'rAmyg', 'rNAcc'};
        recos = cell(numel(scodes), 2);
        for sc = 1:numel(scodes)

            % get points, triangles, neighbors
            [p, t, n] = mesh_reconstruct(smoothdata3( ...
                double(chv == scodes(sc) | chv == (scodes(sc) + 100)), ...
                [2, 2, 2]) >= 0.5);

            % add offset to p
            p = 128 + (0.5 .* p) + ones(size(p, 1), 1) * [-25, -25.5, -40];

            % create new surface
            recos{sc} = xff('new:srf');

            % fill with data
            recos{sc}.NrOfVertices = size(p, 1);
            recos{sc}.NrOfTriangles = size(t, 1);
            recos{sc}.VertexCoordinate = p;
            recos{sc}.ConvexRGBA = [208, 208, 208, 255] ./ 255;
            recos{sc}.ConcaveRGBA = [176, 176, 176, 255] ./ 255;
            recos{sc}.VertexColor = zeros(size(p, 1), 4);
            recos{sc}.Neighbors = n;
            recos{sc}.TriangleVertex = t;
            recos{sc}.RecalcNormals;
            recos{sc}.Smooth(100, 0.1, struct('show', false));
            recos{sc}.SaveAs([cpath 'colin_subcort_' snames{sc} '.srf']);

            % combined surface
            if sc == 1
                reco = recos{1};
            else
                recoc = reco.Combine(recos{sc});
                reco.ClearObject;
                reco = recoc;
                if sc == numel(scodes)
                    reco.SaveAs([cpath 'colin_subcort_combined.srf']);
                    reco.ClearObject;
                    reco = 0;
                end
            end

            % ICBM-normalized version
            if opts.icbmsph

                % copy object
                recos{sc, 2} = recos{sc}.CopyObject;

                % apply normalization to coordinates (with INV file -> forward)
                nc = applyspmsnc(128 - recos{sc}.VertexCoordinate(:, [3, 1, 2]), ...
                    sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

                % set coordinates, recalc normals
                recos{sc, 2}.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
                recos{sc, 2}.RecalcNormals;

                % and save file
                recos{sc, 2}.SaveAs([cpath 'colin_subcort_' snames{sc} '_ICBMnorm.srf']);

                % combined surface
                if sc == 1
                    recoi = recos{1, 2};
                else
                    recoic = recoi.Combine(recos{sc, 2});
                    recoi.ClearObject;
                    recoi = recoic;
                    if sc == numel(scodes)
                        recoi.SaveAs([cpath 'colin_subcort_combined_ICBMnorm.srf']);
                        recoi.ClearObject;
                        recoi = 0;
                    end
                end
            end
        end

        % create actually combined surface
        [p, t, n] = mesh_reconstruct(smoothdata3( ...
            double(chv > 0), [2, 2, 2]) >= 0.5);

        % add offset to p
        p = 128 + (0.5 .* p) + ones(size(p, 1), 1) * [-25, -25.5, -40];

        % create new surface
        reco = xff('new:srf');

        % fill with data
        reco.NrOfVertices = size(p, 1);
        reco.NrOfTriangles = size(t, 1);
        reco.VertexCoordinate = p;
        reco.ConvexRGBA = [208, 208, 208, 255] ./ 255;
        reco.ConcaveRGBA = [176, 176, 176, 255] ./ 255;
        reco.VertexColor = zeros(size(p, 1), 4);
        reco.Neighbors = n;
        reco.TriangleVertex = t;
        reco.RecalcNormals;
        reco.Smooth(100, 0.1, struct('show', false));
        reco.SaveAs([cpath 'colin_subcort.srf']);

        % ICBM-normalized version
        if opts.icbmsph

            % apply normalization to coordinates (with INV file -> forward)
            nc = applyspmsnc(128 - reco.VertexCoordinate(:, [3, 1, 2]), ...
                sni.Tr, sni.VG(1).dim, inv(sni.VG(1).mat), sni.VF(1).mat * sni.Affine);

            % set coordinates, recalc normals
            reco.VertexCoordinate = 128 - nc(:, [2, 3, 1]);
            reco.RecalcNormals;

            % and save file
            reco.SaveAs([cpath 'colin_subcort_ICBMnorm.srf']);
        end

        % clear
        reco.ClearObject;
        reco = 0;
    end
catch ne_eo;
    warning( ...
        'neuroelf:InternalError', ...
        'Error creating sub-cortical surfaces: %s.', ...
        ne_eo.message ...
    );
end

% Shen et al. Atlas files
try
    if opts.shenfiles && (exist([ppath 'shen_parcels_BH.srf'], 'file') ~= 2 || ...
        exist([ppath 'shen_parcels_LH.srf'], 'file') ~= 2 || ...
        exist([ppath 'shen_parcels_RH.srf'], 'file') ~= 2 || ...
        exist([ppath 'shen_parcels.voi'], 'file') ~= 2)
        neuroelf_makeshenfiles;
    end
catch ne_eo;
    warning('neuroelf:makeFiles:internalError', ...
        'Error creating Shen et al. Atlas files: %s.', ne_eo.message);
end

% DARTEL-based templates
hmc = false;
try

    % requested and required MAT file
    if opts.darttmp && ...
        exist([spath 'spm8_dartel_template.mat'], 'file') ~= 2

        % download template file
        hmc = true;
        neuroelf_gui('httpget', 'http://neuroelf.net/spm8_dartel_template.mat', ...
            [spath 'spm8_dartel_template.mat'], struct('steps', [276480, Inf]));
    end

    % requested and required template files
    if opts.darttmp && ...
       (exist([spath 'dartel/Template1mm'], 'dir') ~= 7 || ...
        exist([spath 'dartel/Template1mm/Template1mm_7.nii'], 'file') ~= 2)
        if exist([spath 'dartel/Template1mm'], 'dir') ~= 7
            mkadir([spath 'dartel/Template1mm'], '-p');
        end
        if exist([spath 'dartel/Template1mm'], 'dir') ~= 7
            error( ...
                'neuroelf:MkadirFailed', ...
                'Error creating DARTEL Template folder.' ...
            );
        end
        h0c = load([spath 'spm8_dartel_template.mat']);
        if ~isstruct(h0c) || ...
           ~isfield(h0c, 'h0c') || ...
           ~isstruct(h0c.h0c)
            error( ...
                'neuroelf:BadFileContent', ...
                'Bad DARTEL Template file content.' ...
            );
        end
        h0c = h0c.h0c;
        dtmp = xff('new:hdr');
        dtmp.VoxelData = single(0);
        dtmp.VoxelData(181, 217, 175, 2) = 0;
        imgdata = dtmp.VoxelData(:, :, :, 1);
        dtmp.Endian = h0c.Endian;
        dtmp.FileMagic = h0c.FileMagic;
        dtmp.NIIFileType = h0c.NIIFileType;
        dtmp.HdrKey = h0c.HdrKey;
        dtmp.ImgDim = h0c.ImgDim;
        dtmp.DataHist = h0c.DataHist;
        dtmp.IntermedData = h0c.IntermedData;
        doff = h0c.RunTimeVars.mskoff;
        smsk = h0c.RunTimeVars.mskcont;
        smsz = size(smsk);
        dmsk = false(size(imgdata));
        dmsk(doff(1):doff(1)+smsz(1)-1, doff(2):doff(2)+smsz(2)-1, doff(3):doff(3)+smsz(3)-1) = smsk;
        for tc = 1:8
            imgdata(dmsk) = (1 / 65535) .* double(h0c.VoxelData(tc, :)');
            dtmp.VoxelData(:, :, :, 1) = imgdata;
            imgdata(dmsk) = (1 / 65535) .* double(h0c.VoxelData(tc + 8, :)');
            dtmp.VoxelData(:, :, :, 2) = imgdata;
            dtmp.SaveAs(sprintf('%sdartel/Template1mm/Template1mm_%d.nii', spath, tc - 1));
        end
        dtmp.ClearObject;
        dtmp = 0;
        mni = h0c.RunTimeVars.WarpToMNI.mni;
        save(sprintf('%sdartel/Template1mm/Template1mm_7_2mni.mat', spath), 'mni', '-v6');
    end
catch ne_eo;
    warning( ...
        'neuroelf:InternalError', ...
        'Error creating DARTEL templates: %s.', ...
        ne_eo.message ...
    );
end

% head meshes
try

    % create head mesh
    if opts.headmesh && ...
        exist([cpath 'colin_head.srf'], 'file') ~= 2 && ...
        exist([cpath 'colin.vmr'], 'file') == 2

        % display information
        idisp('Creating colin_head.srf...');

        % load VMR
        if isxff(colin, true)
            colin.ClearObject;
            colin = 0;
        end
        colin = xff([cpath 'colin.vmr']);

        % use the GUI
        neuroelf_gui;
        hmc = true;
        colin.Browse;
        neuroelf_gui('srf_tools', 'findintensity', {'77', '0.75', '85', 'in'});

        % then save resulting surface
        srf = neuroelf_gui('srf_save', [cpath 'colin_head.srf']);
        srf.ClearObject;
        colin.ClearObject;
        colin = 0;
    end

    % create ICBMnorm head mesh
    if opts.headmesh && ...
        exist([cpath 'colin_head_ICBMnorm.srf'], 'file') ~= 2 && ...
        exist([cpath 'colin_ICBMnorm.vmr'], 'file') == 2

        % display information
        idisp('Creating colin_head_ICBMnorm.srf...');

        % load VMR
        if isxff(colin, true)
            colin.ClearObject;
            colin = 0;
        end
        colin = xff([cpath 'colin_ICBMnorm.vmr']);

        % use the GUI
        neuroelf_gui;
        hmc = true;
        colin.Browse;
        neuroelf_gui('srf_tools', 'findintensity', {'77', '0.75', '85', 'in'});

        % then save resulting surface
        srf = neuroelf_gui('srf_save', [cpath 'colin_head_ICBMnorm.srf']);
        srf.ClearObject;
        colin.ClearObject;
        colin = 0;
    end

    % create dura mesh
    if opts.headmesh && ...
        exist([cpath 'colin_dura.srf'], 'file') ~= 2 && ...
        exist([cpath 'colin_brain.vmr'], 'file') == 2

        % display information
        idisp('Creating colin_dura.srf...');

        % load VMR
        if isxff(colin, true)
            colin.ClearObject;
            colin = 0;
        end
        colin = xff([cpath 'colin_brain.vmr']);

        % use the GUI
        neuroelf_gui;
        hmc = true;
        colin.Browse;
        neuroelf_gui('srf_tools', 'findintensity', {'88', '0.75', '75', 'in'});

        % then save resulting surface
        srf = neuroelf_gui('srf_save', [cpath 'colin_dura.srf']);
        srf.ClearObject;
        colin.ClearObject;
        colin = 0;
    end

    % create ICBMnorm dura mesh
    if opts.headmesh && ...
        exist([cpath 'colin_dura_ICBMnorm.srf'], 'file') ~= 2 && ...
        exist([cpath 'colin_brain_ICBMnorm.vmr'], 'file') == 2

        % display information
        idisp('Creating colin_dura_ICBMnorm.srf...');

        % load VMR
        if isxff(colin, true)
            colin.ClearObject;
            colin = 0;
        end
        colin = xff([cpath 'colin_brain_ICBMnorm.vmr']);

        % use the GUI
        neuroelf_gui;
        hmc = true;
        colin.Browse;
        neuroelf_gui('srf_tools', 'findintensity', {'77', '0.75', '75', 'in'});

        % then save resulting surface
        srf = neuroelf_gui('srf_save', [cpath 'colin_dura_ICBMnorm.srf']);
        srf.ClearObject;
        colin.ClearObject;
        colin = 0;
    end
    
catch ne_eo;
    warning( ...
        'neuroelf:InternalError', ...
        'Error creating head meshes: %s.', ...
        ne_eo.message ...
    );
end

% clear objects
clearxffobjects(sph);
clearxffobjects( ...
    {chdr, colin, dtmp, flat, hemi, hcrv, icbm, pvmr, reco, recoi, res, skull, vmr});
clearxffobjects(jhem);
clearxffobjects(recos(:)');

% NeuroElf GUI was used (for head-meshes, etc.)
if hmc
    neuroelf_gui('closewindow', 'yes');
end
clear xff;

% sub-function
function idisp(t)
disp(t);
drawnow;
