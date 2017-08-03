function outfile = nyspi_recon_mux_gpuArray(pfile, outfile, sense, fermi, opts)
% nyspi_recon_mux_gpuArray  - reconstruct multiband sequence to NII file
%
% FORMAT:       outfile = nyspi_recon_mux_gpuArray(pfile, outfile, sense, fermi)
%
% Input fields:
%
%       pfile       P-file name (GE scanner raw data file)
%       outfile     output NII file (create filename if empty)
%       fermi       apply circular Fermi filter (default: false)
%       sense       use SENSE (default: false -> use GRAPPA)
%       opts        optional settings
%        .debugout  if 1x1 true, print some information to console (false)
%        .parpool   if 1x1 true, attempt using parpool (false)
%        .scaling   if 1x1 true, apply data scaling (receiver gains, false)
%
% Output fields:
%
%       outfile     filename of output file
%
% Main program for slice-multiplexing EPI reconstruction.
%
% (c) Kangrong Zhu,     Stanford University                      Aug 2012
% (c) Jochen Weber,     New York State Psychiatric Institute     Dec 2014

% contributions: all main algorithm developed by Kangrong Zhu, GPUarray and
% optimizations by Jochen Weber

% main input argument
if nargin < 1 || ~ischar(pfile) || isempty(pfile) || exist(pfile(:)', 'file') ~= 2
    error('nyspi:mux_recon:argumentError', 'No pfile specified!');
end

% options
if ~exist('opts', 'var') || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'debugout') || ~islogical(opts.debugout) || numel(opts.debugout) ~= 1
    opts.debugout = false;
end
if ~isfield(opts, 'parpool') || ~islogical(opts.parpool) || numel(opts.parpool) ~= 1
    opts.parpool = false;
end

% setting for gpuArray or not
nr_usegpu = false;
try
    g = gpuDevice;
    if g.AvailableMemory >= 8e9
        disp('Using GPU device.');
        pause(0.01);
        nr_usegpu = true;
    end
catch nr_eo;
    if opts.debugout
        fprintf('Error while attempting to access a GPU: %s.\n', nr_eo.message);
    end
end

% main output arguments
if nargin < 2 || ~ischar(outfile) || isempty(outfile) || numel(outfile) < 5 || size(outfile, 1) > 1
    outfile = '';
else
    outfile = outfile(:)';
end

% other inputs
if ~exist('sense', 'var')
    sense = false;
end
if ~exist('fermi', 'var')
    fermi = false;
end

% We do lots of ffts on arrays of the same size, so it's worth letting fftw
% measure the optimal algorithm the first time we run one.
fftw('planner', 'measure');

% construct required filenames
[pfile, ref, vrgf, refp] = get_pfile_name(pfile);

% read all parameters of the mux epi scan and for the reconstruction into a structure
p = mux_epi_params(pfile, [], [], [], sense);
[pgoodvols, pisgood] = full_epi_vols(pfile, p);
if ~pisgood
    p.nt_to_recon = pgoodvols - (p.pfile_header.npasses - p.nt_to_recon);
    if p.nt_to_recon > 0 && opts.debugout
        fprintf('PFile %s data incomplete; truncating to %d volumes!\n', pfile, p.nt_to_recon);
    else
        error('PFile %s MUX REF data incomplete; cannot reconstruct!\n', pfile);
    end
end

% construct default output filename?
if isempty(outfile)
    outfile = sprintf('E%d_S%d_mux%d_arc%d.nii', ...
        p.pfile_header.fullcont.exam.ex_no, p.pfile_header.fullcont.series.se_no, ...
        p.mux, p.pfile_header.fullcont.rdb_hdr.ileaves);
    if opts.debugout
        fprintf('Writing output to %s.\n', outfile);
    end
end

% additional files must exist
if ~exist(ref, 'file') && strcmp(p.ref_for_default_ecc, 'ref.dat')
    error('Cannot find the ref.dat file!');
end
if ~exist(vrgf, 'file') && p.vrgf
    error('Cannot find the vrgf.dat file!');
end
if ~exist(refp, 'file') && (strcmp(p.ref_for_default_ecc, 'ref pfile') || ...
    (p.md_ecc && strcmp(p.ref_for_md_ecc, 'ref pfile')))
    error('Cannot find the reference scan pfile!');
end

% calibration points and output size
cps = p.mux * p.num_mux_cycle;
outsize = [p.nx_pres, p.ny_pres, p.num_unmuxed_slices, p.nt_to_recon + p.num_mux_cycle];

% write output file
if ~generate_nifti(outfile, outsize, p)
    error('Cannot generate output file');
end

% scaling (added between Aug'14 and Oct'14!)
if ~isfield(opts, 'scaling') || ~islogical(opts.scaling) || numel(opts.scaling) ~= 1 || ~opts.scaling
    imscaling = 1;
else
    imscaling = 0.1 * 0.5 ^ ((p.r2 - 15) + 0.5 * (p.r1 - 14));
    if opts.debugout
        fprintf('Applying a receiver-gain scaling factor of %.5g\n', imscaling);
    end
end

% get number of threads
warning('off','MATLAB:maxNumCompThreads:Deprecated');
nt = maxNumCompThreads;

% with GPU
if nr_usegpu
    
    % check memory again, and, if requested, try opening parpool
    g = gpuDevice;
    try
        usepp = false;
        gmem = g.AvailableMemory;
        if opts.parpool && gmem >= 1.1e10
            pp = parpool(floor(gmem / 5.45e9));
            usepp = true;
        else
            pp = [];
        end
    catch
        if opts.debugout
            fprintf('parpool failed.\n');
        end
        pp = [];
    end

    % with A LOT of memory (and parpool option)
    if usepp
        recstep = 8;
        ppsize = floor(g.AvailableMemory / 5.45e9);
        ppnt = max(1, nt / ppsize);

        % par-for
        parfor sl = 1:numel(p.slices_to_recon)

            % give workers a slight temporal discrepancy to avoid
            % GPU out-of-memory errors
            pause((labindex - 1) * 12);

            % set fermfilt and gker to empty
            fermfilt = [];
            gker = [];
            warning('off','MATLAB:maxNumCompThreads:Deprecated');
            maxNumCompThreads(ppnt);

            % iterate over volumes (actual time points)
            for tpoint = [1, 2:recstep:p.nt_to_recon]

                % read calibration data
                if tpoint == 1
                    icps = 0;
                    rcps = cps;
                    mpoints = 0;
                    epoint = tpoint;
                else
                    icps = cps;
                    rcps = 0;
                    mpoints = p.num_mux_cycle;
                    epoint = min(p.nt_to_recon, tpoint + recstep - 1);
                end

                % load the EPI time series data
                [dat, p_recon] = epi_load_tseries(pfile, ref, vrgf, refp, ...
                    sl, p.sl_acq_order(sl), [], tpoint + icps, epoint + icps, rcps, p, nr_usegpu);
                % p_recon will be used as the parameter structure in the reconstruction
                % several fields may have been changed in function epi_load_tseries
                % parameters may change later if external calibration or coil compression
                % is used; can't modify values in p because it is reused by multiple slices

                % get data size
                szdat = [size(dat{1}), numel(dat)];

                % data whitening
                % after the whitening transformation, the coil noise covariance matrix
                % will become identity and can be ignored in further data processing
                % (SOS will directly be SNR optimal coil combination and the sample
                % covariance matrix in SENSE recon will become identity)
                if ~isempty(p.coil_noise_std)
                    datsz = get_dat_sz(dat, p);
                    % dim: [Coil, Kx-Ky-Echo-Slice-Time]
                    for dc = 1:numel(dat)
                        dat{dc} = reshape(permute(dat{dc}, [5,1,2,3,4]), ...
                            [datsz.c, datsz.x*datsz.y*datsz.ec*datsz.sl]);

                        % if covariance matrix M of the random vector X is positive definite,
                        % then M^(-1/2)X has covariance matrix I. Reference:
                        % Matt A. Bernstein et al, MRM 1994, 32:330-334, Eq[1]
                        dat{dc} = diag(1./p.coil_noise_std) * dat{dc};

                        % dim: [Kx, Ky, Echo, Slice, Coil, Time]
                        dat{dc} = permute(reshape(dat{dc}, ...
                            [datsz.c, datsz.x, datsz.y, datsz.ec, datsz.sl]), [2,3,4,5,1]);
                    end
                end

                % Fermi filter
                if fermi
                    if isempty(fermfilt)
                        fermfilt = gpuArray(gen_fermi_filter([size(dat{1},1), size(dat{1},2)], ...
                        p_recon.fermi_radius, p_recon.fermi_width, 'circ'));
                    end
                    for dc = 1:numel(dat)
                        dat{dc} = dat{dc} .* repmat(fermfilt, [1, 1, szdat(3:5)]);
                    end
                end

                % GRAPPA
                [im, gker] = mux_epi_process_data_grappa(dat, p_recon, gker, nr_usegpu);

                % -- Fix phase-encode direction for pepolar scans
                if bitand(p.dacq_ctrl, 4) % 3rd bit indicates odd-echo phase-flip (i.e., "pepolar")
                    im = flip(im, p.PE_DIM);
                end

                % final reshape and scaling
                im = imscaling .* reshape(abs(im), [p.nx_pres, p.ny_pres, p.mux, size(im, p.T_DIM)]);

                % -- Fix slice ordering.
                if p.descend_acq                                                     % If the multiplexed slices are acquired in descending order. Need this because the unmuxed slices are always in ascending order (i.e., from bottom to top) for the current DFT encoding scheme in the PSD.
                    sl_loc = (p.num_slices*(p.mux-1)+sl) : -p.num_slices : 1;
                else
                    sl_loc = sl : p.num_slices : p.num_unmuxed_slices;
                end

                % then write into NII
                nfid = fopen(outfile, 'r+', 'l');
                if isa(im, 'gpuArray')
                    im = gather(im);
                end
                write_into_nifti(nfid, outsize, im(end:-1:1, :, :, :), sl_loc, tpoint + mpoints);
                fclose(nfid);
            end
        end
        delete(pp);

    % GPU with less memory
    else
        
        % regular-for
        recstep = 16;
        for sl = 1:numel(p.slices_to_recon)

            % set gker to empty and reset GPU
            fermfilt = [];
            gker = [];

            % iterate over volumes (actual time points)
            for tpoint = [1, 2:recstep:p.nt_to_recon]

                % read calibration data
                if tpoint == 1
                    icps = 0;
                    rcps = cps;
                    mpoints = 0;
                    epoint = tpoint;
                    if opts.debugout
                        fprintf('Working on slice %d, volume 1 (+kernel, etc.)\n', sl);
                    end
                else
                    icps = cps;
                    rcps = 0;
                    mpoints = p.num_mux_cycle;
                    epoint = min(p.nt_to_recon, tpoint + recstep - 1);
                    if opts.debugout
                        fprintf('Working on slice %d, volumes %d through %d\n', sl, tpoint, epoint);
                    end
                end

                % load the EPI time series data
                [dat, p_recon] = epi_load_tseries(pfile, ref, vrgf, refp, ...
                    sl, p.sl_acq_order(sl), [], tpoint + icps, epoint + icps, rcps, p, nr_usegpu);
                % p_recon will be used as the parameter structure in the reconstruction
                % several fields may have been changed in function epi_load_tseries
                % parameters may change later if external calibration or coil compression
                % is used; can't modify values in p because it is reused by multiple slices

                % get data size
                szdat = [size(dat{1}), numel(dat)];

                % data whitening
                % after the whitening transformation, the coil noise covariance matrix
                % will become identity and can be ignored in further data processing
                % (SOS will directly be SNR optimal coil combination and the sample
                % covariance matrix in SENSE recon will become identity)
                if ~isempty(p.coil_noise_std)
                    datsz = get_dat_sz(dat, p);
                    % dim: [Coil, Kx-Ky-Echo-Slice-Time]
                    for dc = 1:numel(dat)
                        dat{dc} = reshape(permute(dat{dc}, [5,1,2,3,4]), ...
                            [datsz.c, datsz.x*datsz.y*datsz.ec*datsz.sl]);

                        % if covariance matrix M of the random vector X is positive definite,
                        % then M^(-1/2)X has covariance matrix I. Reference:
                        % Matt A. Bernstein et al, MRM 1994, 32:330-334, Eq[1]
                        dat{dc} = diag(1./p.coil_noise_std) * dat{dc};

                        % dim: [Kx, Ky, Echo, Slice, Coil, Time]
                        dat{dc} = permute(reshape(dat{dc}, ...
                            [datsz.c, datsz.x, datsz.y, datsz.ec, datsz.sl]), [2,3,4,5,1]);
                    end
                end

                % Fermi filter
                if fermi
                    if isempty(fermfilt)
                        fermfilt = gpuArray(gen_fermi_filter([size(dat{1},1), size(dat{1},2)], ...
                        p_recon.fermi_radius, p_recon.fermi_width, 'circ'));
                    end
                    for dc = 1:numel(dat)
                        dat{dc} = dat{dc} .* repmat(fermfilt, [1, 1, szdat(3:5)]);
                    end
                end

                % GRAPPA
                [im, gker] = mux_epi_process_data_grappa(dat, p_recon, gker, nr_usegpu);

                % -- Fix phase-encode direction for pepolar scans
                if bitand(p.dacq_ctrl, 4) % 3rd bit indicates odd-echo phase-flip (i.e., "pepolar")
                    im = flip(im, p.PE_DIM);
                end

                % final reshape and scaling
                im = imscaling .* reshape(abs(im), [p.nx_pres, p.ny_pres, p.mux, size(im, p.T_DIM)]);

                % -- Fix slice ordering.
                if p.descend_acq                                                     % If the multiplexed slices are acquired in descending order. Need this because the unmuxed slices are always in ascending order (i.e., from bottom to top) for the current DFT encoding scheme in the PSD.
                    sl_loc = (p.num_slices*(p.mux-1)+sl) : -p.num_slices : 1;
                else
                    sl_loc = sl : p.num_slices : p.num_unmuxed_slices;
                end

                % then write into NII
                nfid = fopen(outfile, 'r+', 'l');
                if isa(im, 'gpuArray')
                    im = gather(im);
                end
                write_into_nifti(nfid, outsize, im(end:-1:1, :, :, :), sl_loc, tpoint + mpoints);
                fclose(nfid);
            end
        end
    end

% no GPU
else
    recstep = 8;
    for sl = 1:numel(p.slices_to_recon)

        % set gker to empty and reset GPU
        % reset(gpuDevice(1));
        fermfilt = [];
        gker = [];

        % iterate over volumes (actual time points)
        for tpoint = [1, 2:recstep:min(40, p.nt_to_recon)]

            % read calibration data
            if tpoint == 1
                icps = 0;
                rcps = cps;
                mpoints = 0;
                epoint = tpoint;
                if opts.debugout
                    fprintf('Working on slice %d, volume 1 (+kernel, etc.)\n', sl);
                end
            else
                icps = cps;
                rcps = 0;
                mpoints = p.num_mux_cycle;
                epoint = min(p.nt_to_recon, tpoint + recstep - 1);
                if opts.debugout
                    fprintf('Working on slice %d, volumes %d through %d\n', sl, tpoint, epoint);
                end
            end

            % load the EPI time series data
            [dat, p_recon] = epi_load_tseries(pfile, ref, vrgf, refp, ...
                sl, p.sl_acq_order(sl), [], tpoint + icps, epoint + icps, rcps, p, nr_usegpu);
            % p_recon will be used as the parameter structure in the reconstruction
            % several fields may have been changed in function epi_load_tseries
            % parameters may change later if external calibration or coil compression
            % is used; can't modify values in p because it is reused by multiple slices

            % get data size
            szdat = [size(dat{1}), numel(dat)];

            % data whitening
            % after the whitening transformation, the coil noise covariance matrix
            % will become identity and can be ignored in further data processing
            % (SOS will directly be SNR optimal coil combination and the sample
            % covariance matrix in SENSE recon will become identity)
            if ~isempty(p.coil_noise_std)
                datsz = get_dat_sz(dat, p);
                % dim: [Coil, Kx-Ky-Echo-Slice-Time]
                for dc = 1:numel(dat)
                    dat{dc} = reshape(permute(dat{dc}, [5,1,2,3,4]), ...
                        [datsz.c, datsz.x*datsz.y*datsz.ec*datsz.sl]);

                    % if covariance matrix M of the random vector X is positive definite,
                    % then M^(-1/2)X has covariance matrix I. Reference:
                    % Matt A. Bernstein et al, MRM 1994, 32:330-334, Eq[1]
                    dat{dc} = diag(1./p.coil_noise_std) * dat{dc};

                    % dim: [Kx, Ky, Echo, Slice, Coil, Time]
                    dat{dc} = permute(reshape(dat{dc}, ...
                        [datsz.c, datsz.x, datsz.y, datsz.ec, datsz.sl]), [2,3,4,5,1]);
                end
            end

            % Fermi filter
            if fermi
                if isempty(fermfilt)
                    fermfilt = gen_fermi_filter([size(dat{1},1), size(dat{1},2)], ...
                    p_recon.fermi_radius, p_recon.fermi_width, 'circ');
                end
                for dc = 1:numel(dat)
                    dat{dc} = dat{dc} .* repmat(fermfilt, [1, 1, szdat(3:5)]);
                end
            end

            % GRAPPA
            [im, gker] = mux_epi_process_data_grappa(dat, p_recon, gker, nr_usegpu);

            % -- Fix phase-encode direction for pepolar scans
            if bitand(p.dacq_ctrl, 4) % 3rd bit indicates odd-echo phase-flip (i.e., "pepolar")
                im = flip(im, p.PE_DIM);
            end

            % final reshape and scaling
            im = imscaling .* reshape(abs(im), [p.nx_pres, p.ny_pres, p.mux, size(im, p.T_DIM)]);

            % -- Fix slice ordering.
            if p.descend_acq                                                     % If the multiplexed slices are acquired in descending order. Need this because the unmuxed slices are always in ascending order (i.e., from bottom to top) for the current DFT encoding scheme in the PSD.
                sl_loc = (p.num_slices*(p.mux-1)+sl) : -p.num_slices : 1;
            else
                sl_loc = sl : p.num_slices : p.num_unmuxed_slices;
            end

            % then write into NII
            nfid = fopen(outfile, 'r+', 'l');
            write_into_nifti(nfid, outsize, im(end:-1:1, :, :, :), sl_loc, tpoint + mpoints);
            fclose(nfid);
        end
    end
end


% sub functions (three empty lines before, from file of first function)



function [pfile, ref, vrgf, refp, pfile_name] = get_pfile_name(d)
% function [pfile, ref, vrgf, refp, pfile_name] = get_pfile_name(d)
%
% Get the full paths to the pfile and the associated ref.dat, vrgf.dat and reference scan pfile.
%
% Input
%   d          - Filename of the pfile, or the directory containing a pfile.
%
% Outputs
%   pfile      - Filename of the p-file.
%   ref        - Filename of the ref.dat file associated with the p-file.
%   vrgf       - Filename of the vrgf.dat file associated with the p-file.
%   refp       - Filename of reference scan pfile associated with the p-file.
%   pfile_name - Name of pfile, PXXXXX.7.
%
% (c) Kangrong Zhu      Stanford University     Oct 2013

% P-file
pfile = d;
[pfile_dir, name, ext] = fileparts(pfile);
if isempty(pfile_dir)
    pfile_dir = pwd;
end
pfile_name = [name ext];

% ref.dat file
ref = get_pfile_associated_file(pfile_dir, pfile_name, 'ref.dat');

% vrgf.dat file
vrgf = get_pfile_associated_file(pfile_dir, pfile_name, 'vrgf.dat');

% Reference scan pfile
refp = get_pfile_associated_file(pfile_dir, pfile_name, 'refscan.7');

function file = get_pfile_associated_file(pfile_dir, pfile_name, associated_file_type)

file = fullfile(pfile_dir, [pfile_name '_' associated_file_type]);
if ~exist(file, 'file')
    file = fullfile(pfile_dir, ['_' pfile_name '_' associated_file_type]);
end
if ~exist(file, 'file')
    file = fullfile(pfile_dir, associated_file_type);
end
if ~exist(file, 'file')
    files_in_dir = dir(fullfile(pfile_dir, ['*' associated_file_type]));
    if ~isempty(files_in_dir)
        file = files_in_dir(1).name;
        file = fullfile(pfile_dir, file);
    end
end
if ~exist(file, 'file')
    file = '';
end



function p = mux_epi_params(pfile, slices, nt_to_recon, n_vcoils, use_sense)
%
% p = mux_epi_params(pfile[, slices, nt_to_recon, n_vcoils, use_sense])
%
% Generates all parameters needed for reconstructing slice-multiplexed EPI data.
%
% Inputs:
%   pfile       - The p-file (need information in p-file header).
%   slices      - Indices of the (muxed) slices to reconstruct. Default: all slices.
%   nt_to_recon - Number of time points to reconstruct, excluding the first
%                 few mux phase cycling time points. Default: all time points.
%   n_vcoils    - Number of virtual coils for coil compression. Set to [] for no coil compression.
%   use_sense   - True: use SENSE; False: use GRAPPA. Default: False.
%                 This flag will be ignored for MICA-type acquisition (always uses SENSE).
%
% Output:
%   p           - A structure with all parameters.
%
% (c) Bob Dougherty <bobd@stanford.edu> Stanford University     September 2012
% Modified by Kangrong Zhu              Stanford University     September 2012

% -- Load full p-file header.
hdr = read_MR_headers(pfile);

% Dimension of the k-/image-space data
p = struct;
p.FE_DIM = 1;        % Dimension for frequency encoding
p.PE_DIM = 2;        % Dimension for phase encoding
p.EC_DIM = 3;        % Dimension for echoes
p.SL_DIM = 4;        % Dimension for slices
p.C_DIM  = 5;        % Dimension for coils
p.T_DIM  = 6;        % Dimension for time points
p.KZ_DIM = 7;        % Dimension for Kz, after reformulating the slice-multiplexed data into a 3D data set

% For the 'repmat' function
p.KEEP_ORIG_SZ = 1;  % Keep the original size in a dimension when repeating a matrix.

% For ky traversal direction
p.TOP_DOWN   = 0;
p.CENTER_OUT = 1;
p.BOTTOM_UP  = 2;

% Number of coefficients per ky line for EPI eddy current correction (i.e. phase correction in x-ky space)
p.NUM_COE_PER_KY_LINE = 2;

% For virtual coil concept
p.ONLY_ACTUAL_COILS        = 0;
p.ACTUAL_AND_VIRTUAL_COILS = 1;
p.ONLY_VIRTUAL_COILS       = 2;

% For Homodyne partial-k recon
RHTYPFRACTNEX = 16;  % The bit indicating fractional nex in rhtype(int16). GE software manual for rhtype: bit4 (16) Set for fractional processing.

% -- Set parameters using p-file header.
% Save some header fields into a structure to avoid loading them again when loading the k-space data using load_raw_tseries.m
p.pfile_header = struct;
p.pfile_header.fullcont   = hdr;
p.pfile_header.version    = hdr.rdb_hdr.rdbm_rev;
p.pfile_header.npasses    = hdr.rdb_hdr.npasses;
p.pfile_header.nslices    = hdr.rdb_hdr.nslices;
p.pfile_header.nframes    = hdr.rdb_hdr.nframes;
p.pfile_header.nechoes    = hdr.rdb_hdr.nechoes;
p.pfile_header.hnover     = hdr.rdb_hdr.hnover;
p.pfile_header.frsize     = hdr.rdb_hdr.frame_size;
p.pfile_header.ncoils     = hdr.rdb_hdr.dab(2)-hdr.rdb_hdr.dab(1)+1;
p.pfile_header.ptsize     = hdr.rdb_hdr.point_size;
p.pfile_header.rawhdrsize = hdr.rdb_hdr.off_data;
p.pfile_header.rawsize    = hdr.rdb_hdr.raw_pass_size;

% Scan params
p.num_echoes      = hdr.rdb_hdr.nechoes;                     % Total number of echoes
p.num_slices      = hdr.image.slquant;                       % Total number of multiplexed slices
p.num_coils       = hdr.rdb_hdr.dab(2)-hdr.rdb_hdr.dab(1)+1; % Total number of receiving coils
p.num_passes      = hdr.rdb_hdr.npasses;                     % Total number of time points in the mux epi scan, including the first few mux phase cycling time points
p.inplane_R       = hdr.rdb_hdr.ileaves;                     % In-plane reduction factor (e.g. ARC 2 corresponds to inplane_R = 2)
p.partial_ky      = (hdr.rdb_hdr.hnover ~= 0);               % True: used partial ky acquisition
p.dacq_ctrl       = hdr.rdb_hdr.dacq_ctrl;                   % Bits contianing information about the DACQ read-out
p.kissoff_views   = hdr.rdb_hdr.kissoff_views;               % Not sure what this is, but seems important for RDS PE-line interleaving
p.vrgf            = hdr.rdb_hdr.vrgf;                        % True: ramp sampling on
p.nx_pres         = hdr.rdb_hdr.vrgfxres;                    % Prescribed size in FE. Pass an empty matrix to use the default value. Default: the same as the size in PE of the raw k-space data
p.kydir           = hdr.rdb_hdr.kydir;                       % Ky traversal direction. Value equal to TOP_DOWN, CENTER_OUT or BOTTOM_UP.
p.ny_pres         = hdr.image.dim_Y;                         % Prescribed size in PE
p.slthick         = hdr.image.slthick;                       % Slice thickness, in mm
p.slspacing       = hdr.image.scanspacing;                   % Slice spacing (gap between two adjacent prescribed slices), in mm.
p.start_loc       = hdr.series.start_loc;                    % The z location of the 1st prescribed slices, in mm
p.mux_excited     = hdr.rdb_hdr.user6;                       % Number of simultaneously excited slices
p.num_mux_cycle   = hdr.rdb_hdr.user7;                       % Number of complete slice phase cycling at the beginning of the slice-multiplexed scan
p.sldist          = hdr.rdb_hdr.user8;                       % Distance between two adjacent simultaneously excited slices, in mm
p.use_gzblips     = hdr.rdb_hdr.user13;                      % True: used Gz blips for slice encoding; False: no Gz blips, direct aliasing of all simultaneous slices
p.cap_random_flag = hdr.rdb_hdr.user17;                      % Type of Gz blip phase encoding. 0: CAIPI; 1: Bit-reversed MICA; 2: Random MICA.
p.swappf          = ~hdr.image.swappf;                       % Whether to swap the phase and frequency encoding directions in the reconed images. True: swap; False: don't swap.
p.flip_angle      = hdr.image.mr_flip;                       % Flip angle, in degrees
p.slice_shuffling = hdr.rdb_hdr.user30;                      % True: slice shuffling during acquisition; false: normal slice order
p.extra_tr        = hdr.rdb_hdr.user31;                      % Wait time after the last slice in each pass
p.pseq            = hdr.rdb_hdr.user32;                      % Type of pulse sequence
p.r1              = hdr.rdb_hdr.ps_aps_r1;                   % Receive gain R1
p.r2              = hdr.rdb_hdr.ps_aps_r2;                   % Receive gain R2

% Fermi filter params
p.fermi_radius   = hdr.rdb_hdr.fermi_radius;
p.fermi_width    = hdr.rdb_hdr.fermi_width;
p.fermi_ecc      = hdr.rdb_hdr.fermi_ecc;

switch p.num_slices
    case hdr.rdb_hdr.nslices/hdr.rdb_hdr.npasses
        p.acq_order = 'interleaved';                         % Interleaved acquisition, hdr.data_acq_tab.pass_number is always 0.
        p.sl_acq_order = hdr.data_acq_tab.slice_in_pass(1 : p.num_slices); % The acquisition order of the prescribed slices
    case hdr.rdb_hdr.nslices/hdr.rdb_hdr.reps                % For sequential acquisition, reps = opfphases is the number of time points.
        p.acq_order = 'sequential';                          % Sequential acquisition, hdr.data_acq_tab.slice_in_pass is always 1.
        p.sl_acq_order = hdr.data_acq_tab.pass_number(1 : p.num_slices) + 1;
        p.pfile_header.npasses = hdr.rdb_hdr.reps;
        p.num_passes = hdr.rdb_hdr.reps;
    otherwise
        error('Inconsistent slice information in header.');
end

p.internal_cal = hdr.rdb_hdr.user5 && (p.num_mux_cycle > 0); % True: acquired mux phase cycling data as internal calibration; False: need external calibration from a separate calibration scan.

if hdr.series.start_loc > hdr.series.end_loc                 % Header indicates the 1st slice has a higher location in the slice select direction than the last slice.
    p.descend_acq = true;                                    % True: acquired slices in descending order, i.e. from higher frequency to lower frequency in the slice select direction; False: acquired in ascending order.
else
    p.descend_acq = false;
end

if hdr.series.start_ras=='S' || hdr.series.start_ras=='I'
    p.scan_orient = 'axial';
elseif hdr.series.start_ras=='A' || hdr.series.start_ras=='P'
    p.scan_orient = 'coronal';
else
    p.scan_orient = 'sagittal';
end

% RF attributes
p.multiband_TBW = hdr.rdb_hdr.user9;                         % Time-bandwidth product of the multiband RF pulse
p.multiband_pw  = hdr.rdb_hdr.user10/1000.0;                 % Duration of the multiband RF pulse

% Flag indicating the raw data order. When we use the RDS to save a pfile, the data are ordered differently than a normal p-file.
% (Note: for rdsdata saved before 2013.10.23, rhuser16=1 indicates rds data order.)
p.rds_data_order = ( (hdr.rdb_hdr.recon>=9000) || (hdr.rdb_hdr.user16==1) );

% Is diffusion scan or not
p.isdifscan = (hdr.rdb_hdr.numdifdirs > 1);                  % True: The scan was a diffusion scan acquired with the muxarcepi2 sequence; False: The scan was acquired with the muxarcepi sequence.

% Params for calibration time points
p.dftz_type = 'inverse';                                     % Type of DFTz encoding conducted on the simultaneous slices. 'inverse': inverse DFT; 'forward': forward DFT.
p.cap_blip_start_cal = hdr.rdb_hdr.user22;                   % Starting index of the kz blips for the calibration time points. 0~(abs(cap_fov_shift_cal)-1) corrrespond to -kzmax~kzmax.
if hdr.rdb_hdr.user25 == 0                                   % Deal with older version sequence
    if p.isdifscan
        % This is set to the # of T2 images for diffusion scans, so default to the mux factor if this is a dwi scan.
        p.cap_fov_shift_cal = p.mux_excited;
    else
        p.cap_fov_shift_cal = hdr.rdb_hdr.user23;                % Older version sequence used rhuser23 to save cap_fov_shift_cal
    end
else
    p.cap_fov_shift_cal = hdr.rdb_hdr.user25;                % CAIPI FOV shift for the calibration time points. A positive integer.
end
if p.cap_fov_shift_cal == 0
    p.cap_fov_shift_cal = p.mux_excited;
end
if (p.mux_excited == 1) && (p.cap_fov_shift_cal ~= 1)
    p.cap_fov_shift_cal = 1;
end
p.mux_encoded = p.cap_fov_shift_cal;
if p.mux_encoded < p.mux_excited
    error('mux_encoded cannot be smaller than mux_excited.');
end
p.recon_cal_use_mux_extra = 0;                               % True: when mux_encoded>mux_excited, use extra encoded slices in recon calibration; False: only use the desired mux_excited slices in recon calibration and ignore the extra slices even if there were any.
if ~(p.use_gzblips && (p.mux_encoded > p.mux_excited))       % Use Gz blips && Encode more slices than excited in calibration time points
    p.recon_cal_use_mux_extra = 0;
end
if p.recon_cal_use_mux_extra
    p.mux = p.mux_encoded;                                   % p.mux is the mux factor in recon
else
    p.mux = p.mux_excited;
end
p.cap_fov_shift_cal = (-1)^strcmp(p.dftz_type, 'inverse') * p.cap_fov_shift_cal; % The type of DFTz encoding for calibration data will be passed as the sign of p.cap_fov_shift_cal.

% Type of acquisition
p.caipi = (p.use_gzblips && p.cap_random_flag == 0);         % CAIPI
p.mica_br = (p.use_gzblips && p.cap_random_flag == 1);       % Bit-reversed MICA
p.mica_rand = (p.use_gzblips && p.cap_random_flag == 2);     % Random MICA
p.mica_perturbed_caipi = (p.use_gzblips && p.cap_random_flag == 3); % MICA, whose kz encoding is CAIPI kz encoding wih random kz perturbation added
p.mica_poisson = (p.use_gzblips && p.cap_random_flag == 4);  % MICA with Poisson-disc-like undersampling pattern on the ky-kz plane

% Params needed only for CAIPI, MICA with randomly perturbed CAIPI scheme, and MICA with Poisson-disc ky-kz pattern
if p.caipi || p.mica_perturbed_caipi || p.mica_poisson
    p.cap_blip_start = hdr.rdb_hdr.user14;                   % Starting index of the kz blips for the accelerated time points. 0~(abs(cap_fov_shift)-1) correspond to -kzmax~kzmax.
    p.cap_blip_inc   = hdr.rdb_hdr.user15;                   % Increment of the kz blip index for adjacent acquired ky lines
    p.cap_fov_shift  = hdr.rdb_hdr.user21;                   % CAIPI FOV shift for the accelerated time points. A positive integer.
    if p.cap_fov_shift == 0
        p.cap_fov_shift = p.mux_excited;
    end
    p.cap_fov_shift = (-1)^strcmp(p.dftz_type, 'inverse') * p.cap_fov_shift; % The shifting in the PE direction of successive slices is FOVy/p.cap_fov_shift. The type of DFTz encoding for accelerated data will be passed as the sign of p.cap_fov_shift.
end

% Params needed only for MICA
if p.mica_br || p.mica_rand || p.mica_perturbed_caipi
    p.cap_seed_shift = hdr.rdb_hdr.user18;                   % Seed shift for subsequent images in subsequent time points
end
if p.mica_perturbed_caipi
    p.cap_kz_rand_pert_range = hdr.rdb_hdr.user26;
end

% Params for eddy current correction
p.md_ecc = 0;                                                % True: use matrix-decoding eddy current correction; False: use the default correction (either a single-slice or a slice-averaged correction).
p.ref_for_default_ecc = 'ref.dat';                           % Which data to use for the default eddy current correction. 'ref.dat': use ref.dat file; 'ecc data': use ECC data collected by the mux sequence; 'ref pfile': use reference scan pfile.
if p.md_ecc
    p.ref_for_md_ecc = 'ref pfile';                          % Which data to use for the matrix-decoding eddy current correction. 'ecc data': use ECC data collected by the mux sequence; 'ref pfile': use reference scan pfile.
end
p.cap_get_ecc = hdr.rdb_hdr.user19;
if p.cap_get_ecc
    p.cap_get_ecc_z = hdr.rdb_hdr.user20;
    p.pcslice       = hdr.rdb_hdr.pcspacial;                 % Slice averaging when calculating x-ky phase correction coefficients. 0: no averaging across slice; >=1 && <=nslics: use one of the slices' coefficients for all slices; -1: average across slices.
    if p.recon_cal_use_mux_extra
        p.ref_for_default_ecc = 'ref.dat';                   % TODO: Enable use of ecc data
        if p.md_ecc
            p.md_ecc = false;                                % TODO: Enable matrix-decoding ghosting correction
        end
    else
        p.ref_for_default_ecc = 'ref.dat';
        if p.md_ecc
            p.ref_for_default_ecc = 'ecc data';
            p.ref_for_md_ecc = 'ecc data';
        end
    end
end
if ~strcmp(p.ref_for_default_ecc, 'ref.dat') || p.md_ecc
    if p.use_gzblips && p.cap_get_ecc && p.cap_get_ecc_z     % Use Gz blips && Measured eddy current effects with Gz blips on
        p.smooth_pha_coe_acc = 0;                            % True: smooth the x-ky phase correction coefficients along the ky direction for the accelerated data; False: no smoothing along ky.
    else
        p.smooth_pha_coe_acc = 1;
    end
end

p.output_ecc_data = (p.cap_get_ecc && (strcmp(p.ref_for_default_ecc, 'ecc_data') || (p.md_ecc && strcmp(p.ref_for_md_ecc, 'ecc data')))); % True: output ECC data collected by the mux sequence along with the reconstructed images; False: don't output the ECC data.

% Receiver coil noise standard deviation
if isfield(hdr, 'psc')
    p.coil_noise_std = hdr.psc.rec_std(1:p.num_coils);       % Receiver coil noise standard deviation. Dim: [1, p.num_coils]
    p.coil_noise_std = 1./ ((1./p.coil_noise_std) ./ sqrt(sum((1./p.coil_noise_std).^2))); % Normalize so that sum_k(1/sigma_k.^2)=1, where sigma_k is the noise standard deviation in the k-th coil. Reference: Matt A. Bernstein et al, MRM 1994, 32:330-334, Eq[1]
else
    p.coil_noise_std = [];                                   % Pass empty matrix to not do data whitening
end

% Params for Homodyne partial-k recon
if p.partial_ky
    rhtype = hdr.rdb_hdr.data_collect_type;
    use_homodyne = ( bitand(uint16(RHTYPFRACTNEX), rhtype) == RHTYPFRACTNEX ); % True: use Homodyne; False: use zero-filling
    if use_homodyne
        p.homodyne_niter = 4;                                % Number of iterations for homodyne partial-k recon
    else
        p.homodyne_niter = 0;
    end
    p.homodyne_ntran = hdr.rdb_hdr.ntran;                    % Transition width for homodyne partial-k recon
end

% -- Specify Y chopping
p.ychop = true;                                              % True: used y-chopping in the acquisition; False: no y-chopping. TODO: Get ychop from the pfile header (Whether to use ychop was not set by bit 0 of rhtype)

% -- Specify whether to remove the encoding added to each individual slice for debugging
p.decode_each_slice = false;

% -- Specify params for coil compression
if isempty(n_vcoils) || (n_vcoils == 0) || (n_vcoils >= p.num_coils)
    p.coil_compress = false;                                 % True: conduct coil compression before parallel imaging recon
else
    p.coil_compress = true;
    p.gcc_slwin     = 5;                                     % Odd number of window size in space for computing the coil compression matrices in Geometric decomposition Coil Compression (GCC)
    if n_vcoils >= p.mux * p.inplane_R
        p.num_vcoils = n_vcoils;                             % Number of virtual coils for coil compression
    else
        p.num_vcoils = p.mux * p.inplane_R;                  % To conduct parallel imaging, must have num_vcoils >= Total_acceleration_factor.
    end
end

% -- Specify params for reconstruction
p.frames             = [];                                   % Always pass an empty matrix to reconstruct all frames.
p.coils              = [];                                   % Always pass an empty matrix to reconstruct all coils.
p.echoes             = [];                                   % Always pass an empty matrix to reconstruct all echoes.
p.num_unmuxed_slices = p.num_slices * p.mux;                 % Total # of unmuxed slices.
p.cal_dat_tpoints    = p.mux*(p.num_mux_cycle-1)+1 : p.mux*p.num_mux_cycle; % Time points corresponding to the calibration data, i.e. the last group of mux phase cycling

% -- Specify reconstruction algorithm to use
if use_sense
    p.mux_recon_method = 'sense';                            % 'sense', '1Dgrappa' or 'slice-grappa'
else
    p.mux_recon_method = '1Dgrappa';
end
if p.mica_br || p.mica_rand || p.mica_perturbed_caipi || p.md_ecc % Always use SENSE if MICA OR if using matrix decoding eddy current correction
    p.mux_recon_method = 'sense';
end

% -- Specify params for SENSE reconstruction
if strcmp(p.mux_recon_method, 'sense')
    p.crop_smap          = 0;                                % True: crop the sensitivity maps to make SENSE recon work better.
    p.sense_comb_cal_dat = false;                            % How to combine the coil images of the fully sampled calibration data. True: use SENSE; False: use SOS (square-root-of-sum-of-squares).
    p.sense_lambda       = [];                               % The lambda for Tikhonov regularization in the SENSE recon. Pass an empty matrix to use the default of norm(encode_mtx' * encode_mtx, 'fro') / size(encode_mtx,2) * p.sense_lambda_default_ratio.
    p.sense_lambda_default_ratio = 0.02;                     % Default regularization coefficient for the SENSE recon is norm(encode_mtx' * encode_mtx, 'fro') / size(encode_mtx,2) * p.sense_lambda_default_ratio.
    p.sense_slwin        = 1;                                % Sliding window length along x for SENSE recon

    % Algorithm to calculate sensitivity maps
    p.smap_type = 'espirit';                                 % Algorithm to calculate the sensitivity maps. 'espirit' or 'coil_over_sos'.

    % Params for ESPIRiT sensitivity map calculations
    if strcmp(p.smap_type, 'espirit')
        p.smap_acssz = struct('x', {64}, 'y', {64});         % K-space size to use to calculate the sensitivity maps. Maximum size possible will be used if empty matrices are passed or if the specified size is larger than the maximum size available.
        p.espirit_nmap = 1;                                  % Number of sensitivity map sets to use in the recon
        p.espirit_ksize = [6, 6];                            % Kernel size to calculate the ESPIRiT sensitivity maps. Dim: [Kx, Ky] (The sensitivity maps will be calculated for each individual slice)

        p.espirit_use_c_code = 0;                            % True: Use Martin Uecker's compiled C code for ESPIRiT
        if ispc
            p.espirit_use_c_code = 0;                        % No compiled C code for Windows
        end
        if p.espirit_use_c_code
            p.espirit_c_code_dir = '~/Dropbox/ESPIRiT_matlab_and_c_code/recon_v0.1.10';% Directory to the compiled C code for ESPIRiT
            if ~exist(p.espirit_c_code_dir, 'dir')
                p.espirit_use_c_code = 0;
            end
        end
        if ~p.espirit_use_c_code                             % Use Michael Lustig's matlab code for ESPIRiT
            p.espirit_eigThresh_1 = 0.02;                    % Threshold to define the nullspace from the 1st SVD in the k-space
            setenv('TOOLBOX_PATH', '');                      % Set environmental variable used by compiled C code to be empty
            p.espirit_tmp_dir = '';                          % Set variable used by compiled C code to be empty
        else                                                 % Use Martin Uecker's compiled C code for ESPIRiT
            p.espirit_eigThresh_1 = 0.001;
            addpath([p.espirit_c_code_dir '/matlab']);       % Add directory to the matlab interface functions
            setenv('TOOLBOX_PATH', p.espirit_c_code_dir);    % Set environmental variable to work around a Matlab bug
            setenv('PATH', strcat(getenv('TOOLBOX_PATH'), ':', getenv('PATH'))); % Set environmental variable to work around a Matlab bug
            if ismac                                         % On Mac
                setenv('DYLD_LIBRARY_PATH', '');
            else                                             % On unix
                setenv('LD_LIBRARY_PATH', '');
            end
            p.espirit_tmp_dir = 'tmp_data/';                 % Temporary directory to store ESPIRiT intermediate results
        end

        if (p.espirit_nmap > 1) || p.espirit_use_c_code
            p.crop_map = 1;
        end
        if p.crop_smap
            p.espirit_eigThresh_2 = 0.95;                    % Threshold on the eigenvalues for cropping the sensitivity maps
        else
            p.espirit_eigThresh_2 = 0;
        end
    end

    % Params for coil-over-sos sensitivity map calculations
    if strcmp(p.smap_type, 'coil_over_sos')
        p.smap_acssz = struct('x', {[]}, 'y', {[]});         % Always use all data available, better performance than using only the central part of k-space
    end
end

% -- Specify params for GRAPPA reconstruction
% The kernel size represents the number of acquired points used in each
% dimension to interpolate one missing point. The ACS size represents the
% number of points used in each dimension for the ACS area.
if strcmp(p.mux_recon_method, '1Dgrappa') || strcmp(p.mux_recon_method, 'slice-grappa')
    p.grappa_domain = 'imSpace';                             % 'imSpace' or 'kSpace', specifies the domain in which the GRAPPA recon will be carried out.

    p.inplane_kersz = struct( 'x',{7},  'y',{4}  );          % Kernel size for solving inplane acceleration. Pass an empty matrix for a field to use the default value. Default: [7, 4] for [x, y].
    p.inplane_acssz = struct( 'x',{[]}, 'y',{[]} );          % ACS size for solving inplane acceleration. Pass an empty matrix for a field to use the default value. Default: use all available calibration k-space data.

    switch p.mux_recon_method
        case '1Dgrappa'
            p.mux_kersz = struct( 'x',{7},  'y',{4}  );      % Kernel size for solving slice multiplexing. Pass an empty matrix for a field to use the default value. Default: [7, 4] for [x, y].
            p.mux_acssz = struct( 'x',{32}, 'y',{p.ny_pres * (p.mux-1)} ); % ACS size for solving slice multiplexing. Pass an empty matrix for a field to use the default value. Default: [32, p.ny_pres * (p.mux-1)] for [x, y].

            p.zpad_1Dgrappa = 'minZpad';                     % How to zero pad in 1D GRAPPA recon. 'minZpad': pad minimum size zero matrices between adjacent slices; 'normalZpad': pad zero matrices of the same size as the original image matrix.
            if p.caipi && (abs(p.cap_fov_shift) ~= p.mux) && strcmp(p.zpad_1Dgrappa, 'minZpad');
                p.zpad_1Dgrappa = 'normalZpad';              % Minimum zero padding does not work for CAIPI if abs(p.cap_fov_shift) ~= p.mux
            end
        case 'slice-grappa'
            p.mux_kersz = struct( 'x',{7}, 'y',{7} );        % Kernel size for slice-GRAPPA. Pass an empty matrix for a field to use the default value. Default: [7, 7] for [x, y]
            p.mux_acssz = struct( 'x',{[]}, 'y',{[]} );      % ACS size for slice-GRAPPA. Pass an empty matrix for a field to use the default value. Default: Use all available k-space data.
            p.slice_grappa_odd_even_fit = 1;                 % True: Fit and apply different kernels to the odd and even ky lines in slice-GRAPPA.
    end
end

% -- Specify whether to add virtual coil concept in recon
% The virtual concept(Reference: Martin Blaimer et al. MRM 2009; 61:93-102)
% is different from coil compression: It synthesizes an additional set of
% data to better utilize the encoding power of the phase(including background phase) in the sensitivity profiles.
% Encoding scheme for virtual coils: Take conjugate of measured data <=>
% Take conjugate of sensitivity, use negative frequencies, and take conjugate of eddy current effect phase term.
p.add_vcc = p.ONLY_ACTUAL_COILS;                             % p.ONLY_ACTUAL_COILS: Use data in only actual coils; p.ACTUAL_AND_VIRTUAL_COILS: Use data in both actual and virtual coils; p.ONLY_VIRTUAL_COILS: Use data in only virtual coils.
if p.add_vcc
    switch p.mux_recon_method
        case 'sense'
            switch p.smap_type
                case 'coil_over_sos'
                    p.smap_acssz = struct('x', {[]}, 'y', {[]}); % Must use all data available, must capture rapid phase change for virtual coil concept to work well
                case 'espirit'
                    p.add_vcc = p.ONLY_ACTUAL_COILS;
            end
        case '1Dgrappa'
            p.add_vcc = p.ONLY_ACTUAL_COILS;
    end
end

% -- Specify whether to return complex (coil) image(s) or coil-combined magnitude image
p.return_comb_mag_im = true;                                 % True: return the coil-combined magnitude image; False: if using GRAPPA, return complex coil images, if using SENSE, return the coil-combined complex image.
if p.partial_ky
    p.return_comb_mag_im = true;                             % If partial ky acquisition was used, reconstructed images are real images, can only return the coil-combined magnitude image.
end

% -- Reconstruction params input from the user
p.debug = false;                                             % True: calculate intermediate results and print out messages.

if isempty(slices)
    slices = 1 : p.num_slices;
end
[~, nslices_inrange] = checkrange(slices, 1, p.num_slices);
if nslices_inrange < length(slices)
    error('Acquired %d slices. Slices to recon contain out of range slice indices.', p.num_slices);
end
p.slices_to_recon = slices;                                  % Indices of the (muxed) slices to reconstruct

p = mux_epi_params_set_tpoints(p, nt_to_recon);



function [v, fileok] = full_epi_vols(fname, params)

% number of expected volumes
v = params.pfile_header.npasses;
fileok = true;

% actual filesize (use dir output)
pfsize = dir(fname);
pfsize = pfsize.bytes;

% compute expected number of bytes per full volume (pass) as...
% number of coils (different acquisition signals) *
% number of slices *
% number of frames (+ hangover) / inplane_R (acceleration) *
% framesize (factor) *
% number of echoes *
% 2 (real+imag) * number of bytes per datapoint (bytes per value)
passbytes = params.pfile_header.ncoils * params.num_slices * ...
    ((params.pfile_header.nframes + params.pfile_header.hnover) / params.inplane_R) * ...
    params.pfile_header.frsize * params.pfile_header.nechoes * 2 * params.pfile_header.ptsize;
 
% compute expected file size (from number of full volumes + mux_ref
% volumes, together with raw header size)
expected_num_bytes = params.pfile_header.npasses * passbytes + ...
    params.mux * (params.inplane_R - 1) * params.num_mux_cycle * passbytes + ...
    params.pfile_header.rawhdrsize;

% if filesize doesn't match
if expected_num_bytes ~= pfsize
    
    % recalculate as the smaller of known volumes and the difference
    % between filesize - (header + mux_ref data)
    v = min(v, floor((pfsize - ...
        (params.mux * (params.inplane_R - 1) * params.num_mux_cycle * passbytes + ...
        params.pfile_header.rawhdrsize)) / passbytes));
    if expected_num_bytes > pfsize
        fileok = false;
    end
end



function nifti_ok = generate_nifti(outfile, outsize, p)

% assume it didn't work
nifti_ok = false;

% prepare output
iprm = p.pfile_header.fullcont.image;
outres = [iprm.dfov / iprm.dim_X, iprm.dfov_rect / iprm.dim_Y, ...
    iprm.slthick + iprm.scanspacing, iprm.tr / 1e6];
sltime = outres(4) / p.num_slices;
desc = sprintf('te=%.2f;ti=%.2f;fa=%.2f;ec=%.4f;acq=[%d,%d];mt=%g;rp=%g;', ...
    iprm.te / 1000, iprm.ti / 1000, iprm.mr_flip, iprm.effechospace / 1e6, ...
    p.pfile_header.fullcont.rdb_hdr.rc_xres, p.pfile_header.fullcont.rdb_hdr.rc_yres, ...
    iprm.offsetfreq, p.pfile_header.fullcont.rdb_hdr.ileaves);
tlhc = [iprm.tlhc_R; iprm.tlhc_A; iprm.tlhc_S];
trhc = [iprm.trhc_R; iprm.trhc_A; iprm.trhc_S];
brhc = [iprm.brhc_R; iprm.brhc_A; iprm.brhc_S];
slnorm = [-iprm.norm_R; -iprm.norm_A; iprm.norm_S];
numbnd = p.pfile_header.fullcont.rdb_hdr.user6;
bndspc = p.pfile_header.fullcont.rdb_hdr.user8;
lrdiff = tlhc - trhc;
sidiff = trhc - brhc;
rowcos = lrdiff ./ sqrt(sum(lrdiff .* lrdiff));
colcos = sidiff ./ sqrt(sum(sidiff .* sidiff));
imgpos = tlhc - (0.5 * (numbnd - 1) * bndspc) .* slnorm - ...
    (rowcos + colcos) .* (0.5 .* outres(1:3)');
slnorm = [-iprm.norm_R; -iprm.norm_A; -iprm.norm_S];
m = generate_4by4matrix(imgpos, rowcos, colcos, slnorm, outres(1:3));
q = compute_quaternion(m);

% open file
nfid = fopen(outfile, 'w+', 'l');
if nfid < 1
    return;
end
fseek(nfid, 0, -1);

% determine some header parameters
s11 = mod(outsize(1), 256);
s12 = floor(outsize(1) / 256);
s21 = mod(outsize(2), 256);
s22 = floor(outsize(2) / 256);
s31 = mod(outsize(3), 256);
s32 = floor(outsize(3) / 256);
s41 = mod(outsize(4), 256);
s42 = floor(outsize(4) / 256);
desc = uint8(desc(:));
desc(80) = 0;
if numel(desc) > 80
    desc(81:end) = [];
end
desc(104) = 0;

% write header
fwrite(nfid, uint8( ...
    [ 92,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, ...
       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, ...
       0,   0,   0,   0,   0,   0, 114,  57,   4,   0, s11, s12, s21, s22, s31, s32, ...
     s41, s42,   1,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0, ...
       0,   0,   0,   0,   0,   0,  16,   0,  32,   0,   0,   0]), 'uint8');
fwrite(nfid, single([1, outres, 1, 1, 1, 352, 1, 0]), 'single');
fwrite(nfid, uint8([outsize(3) - 1,   0,   1,  10]), 'uint8');
fwrite(nfid, single([1000, 1, sltime, 0, 0, 0]), 'single');
fwrite(nfid, desc, 'uint8');
fwrite(nfid, uint8([1, 0, 1, 0]), 'uint8');
fwrite(nfid, single([q(1), q(2), q(3), ...
    m(1, 4), m(2, 4), m(3, 4), m(1, :), m(2, :), m(3, :), 0, 0, 0, 0]), 'single');
fwrite(nfid, uint8([110, 43, 49, 0, 0, 0, 0, 0]), 'uint8');

% extend file as necessary 
e = forceseek(nfid, 352 + 4 * prod(outsize), -1);
if e < 0
    fclose(nfid);
    return;
end
fclose(nfid);
nifti_ok = true;


function write_ok = write_into_nifti(nfid, outsize, d, sl_loc, tpoint)
write_ok = false;
slsize = prod(outsize(1:2));
try
    for tc = 1:size(d, 4)
        for sc = 1:numel(sl_loc)
            fseek(nfid, 352 + 4 * ((tpoint - 1) * prod(outsize(1:3)) + ...
                (sl_loc(sc) - 1) * slsize), -1);
            fwrite(nfid, reshape(single(d(:, :, sc, tc)), slsize, 1), 'single');
        end
        tpoint = tpoint + 1;
    end
catch errobj
    warning(errobj.message);
    return;
end
write_ok = true;


function m = generate_4by4matrix(t, r1, r2, r3, s)

% translation
t44 = eye(4);
t44(1:3, 4) = t(:);

% rotation
r44 = [ ...
   -r1(1), -r2(1), -r3(1), 0; ...
   -r1(2), -r2(2), -r3(2), 0; ...
   -r1(3), -r2(3), -r3(3), 0; ...
     0, 0, 0, 1];

% scaling
s44 = diag([s, 1]);

% complete matrix
m = t44 * r44 * s44;

function Q = compute_quaternion(M)
% Convert from rotation matrix to quaternion form
% See: http://skal.planet-d.net/demo/matrixfaq.htm

% computation
R = M(1:3, 1:3);
vx = sqrt(sum(M(1:3, 1:3) .^ 2));
vx(vx == 0) = 1;
R = R * diag(1 ./ vx);
[U, S, V] = svd(R);
R = U * V';
if det(R) <= 0
    R = R * diag([1, 1, -1]);
end
d = diag(R);
t = sum(d) + 1;
if t > 0.5
    s = sqrt(t) * 2;
    Q = [(R(3, 2) - R(2, 3)) / s, ...
         (R(1, 3) - R(3, 1)) / s, ...
         (R(2, 1) - R(1, 2)) / s, ...
         0.25 * s]';
else
    t = find(d == max(d));
    t = t(1);
    switch (t)
        case {1}
            s = 2 * sqrt(1 + R(1, 1) - R(2, 2) - R(3, 3));
            Q = [0.25 * s, ...
                 (R(1, 2) + R(2, 1)) / s, ...
                 (R(3, 1) + R(1, 3)) / s, ...
                 (R(3, 2) - R(2, 3)) / s]';
        case {2}
            s = 2 * sqrt(1 + R(2, 2) - R(1, 1) - R(3, 3));
            Q = [(R(1, 2) + R(2, 1)) / s, ...
                 0.25 * s, ...
                 (R(2, 3) + R(3, 2)) / s, ...
                 (R(1, 3) - R(3, 1)) / s]';
        case {3}
            s = 2 * sqrt(1 + R(3, 3) - R(1, 1) - R(2, 2));
            Q = [(R(3, 1) + R(1, 3)) / s, ...
                 (R(2, 3) + R(3, 2)) / s, ...
                 0.25 * s, ...
                 (R(2, 1) - R(1, 2)) / s]';
    end
end
if Q(4) < 0
    Q = -Q;
end


function e = forceseek(fid, p, o)
e = fseek(fid, p, o);
if e < 0 && o == -1
    fseek(fid, 0, 1);
    eop = p - ftell(fid);
    bck = zeros(8192, 1);
    while eop >= 65536
        fwrite(fid, bck, 'double');
        eop = eop - 65536;
    end
    if eop > 0
        fwrite(fid, uint8(zeros(eop, 1)), 'uint8');
    end
    if ftell(fid) ~= p
        e = -1;
    else
        e = 0;
    end
end


function [dat, gker] = mux_epi_process_data_grappa(dat, p, gker, nr_usegpu)
%
% function [dat, gker] = mux_epi_process_data_grappa(dat, p, gker, nr_usegpu)
%
% Process slice-multiplexed EPI data. Use GRAPPA for reconstruction.
%
% Inputs
%   dat - MUX EPI data, already phase and ramp-sample corrected. The first
%         p.mux_encoded*p.num_mux_cycle time points are mux phase
%         cycling time points(totally p.num_mux_cycle mux phase cycles).
%         The rest time points are accelerated MUX EPI data. Dim: [Kx, Ky,
%         Echo, Slice, Coil, Time].
%   p   - Parameter structure. See mux_epi_params.m for details. The following
%         fields are used in this function: inplane_R, mux, partial_ky,
%         ny_pres, num_mux_cycle, debug, inplane_kersz, inplane_acssz,
%         mux_kersz, mux_acssz, PE_DIM, C_DIM, T_DIM, caipi, cal_dat_tpoints,
%         cap_fov_shift_cal, cap_fov_shift, grappa_domain, return_comb_mag_im.
%   gker- kernel used for GRAPPA reconstruction (if empty, compute)
%   nr_usegpu- use GPU flag
%
% Output
%   dat - Reconstructed images. Dim: [FE, PE, Echo, Slice*SimultaneousSlice, Coil(=1 if
%         p.return_comb_mag_im == true; =nc if p.return_comb_mag_im == false), Time].
%
% Calibration data for GRAPPA
%   (A) Inplane
%         time points in dat to recon                                         -- time points in dat used as calibration
%         1 : abs(cap_fov_shift_cal)*num_mux_cycling(Not accelerated inplane) -- NONE
%         abs(cap_fov_shift_cal)*num_mux_cycling+1 : end                      -- abs(cap_fov_shift_cal)*(num_mux_cycling-1)+1 (Zblip Off) or synthesized from abs(cap_fov_shift_cal)*(num_mux_cycling-1)+1 : abs(cap_fov_shift_cal)*num_mux_cycling (Zblip On)
%   (B) Through-plane
%         time points in dat to recon                                         -- time points in dat used as calibration
%         1 : abs(cap_fov_shift_cal)*num_mux_cycling                          -- NONE
%         abs(cap_fov_shift_cal)*num_mux_cycling+1 : end                      -- abs(cap_fov_shift_cal)*(num_mux_cycling-1)+1 : abs(cap_fov_shift_cal)*num_mux_cycling
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

% required size
sz_dat = get_dat_sz(dat, p);

% -- For through-plane accelerated scans
if p.mux >1
    % - Default parameters for through-plane GRAPPA recon
    if isempty(p.mux_kersz.x)
        p.mux_kersz.x = 7;
    end
    if isempty(p.mux_kersz.y)
        p.mux_kersz.y = 4;
    end
    if isempty(p.mux_acssz.x)
        p.mux_acssz.x = 32;
    end
    if isempty(p.mux_acssz.y)
        p.mux_acssz.y = p.ny_pres * (p.mux-1);
    end
end

% -- If partial ky acquisition was used, zero-pad k-space to the prescribed size.
if p.partial_ky
    if nr_usegpu
        for dc = 1:numel(dat)
            dat{dc} = cat(p.PE_DIM, dat{dc}, gpuArray.zeros(sz_dat.x, ...
                p.ny_pres-sz_dat.y, sz_dat.ec, sz_dat.sl, sz_dat.c));
        end
    else
        for dc = 1:numel(dat)
            dat{dc} = cat(p.PE_DIM, dat{dc}, zeros(sz_dat.x, ...
                p.ny_pres-sz_dat.y, sz_dat.ec, sz_dat.sl, sz_dat.c));
        end
    end

    % # of acquired ky lines in partial ky acquisition
    ny_partk = sz_dat.y;
    sz_dat.y = p.ny_pres;
end

% -- Through-plane GRAPPA recon.
if p.mux > 1
    
    % Dim:  [FE, PE, Echo, Slice, Coil, TemporalPhase]
    [dat, gker] = mux_recon_1Dgrappa(dat, p, gker, nr_usegpu);
end

% -- Change to image space if p.mux=1 or if through-plane GRAPPA returned kspace data
if strcmp(p.grappa_domain, 'kSpace') || (p.mux == 1)
    dat = ifft2c(dat, nr_usegpu);
end

% -- If partial k-space acquisition was used, recover missing k-space data.
if p.partial_ky
    % Data input to and output from homodyne are all image space data
    dat = homodyne(dat, ny_partk, p.homodyne_ntran, p.homodyne_niter, 'xy', 'xy', nr_usegpu);
end

% -- Final images
if p.return_comb_mag_im
    dat = sos(dat, p.C_DIM);
end



function [reconed, gker] = mux_recon_1Dgrappa(ksp, p, gker, nr_usegpu)
%
% function [reconed, genkern] = mux_recon_1Dgrappa(ksp, p, genkern, nr_usegpu)
% 
% Reconstruct slice-multiplexed data using 1D GRAPPA.
% The last group of mux phase cycling is used as the reference images.
%
% Inputs
%   ksp     - Slice-multiplexed k-space data. The first 'p.mux_encoded
%             *p.num_mux_cycle' time points are the mux phase cycling time points.
%             Dim: [Kx, Ky, Echo, Slice(with slice multiplexing), Coil, Time].
%   p       - Parameter structure. See mux_epi_params.m for details.
%             The following fields are used in this function: mux_kersz,
%             mux_acssz, mux, num_mux_cycle, PE_DIM, T_DIM, KEEP_ORIG_SZ,
%             caipi, cap_fov_shift, grappa_domain, zpad_1Dgrappa, debug.
%   genkern - kernel required for GRAPPA reconstruction (generated if empty)
%   nr_usegpu- use GPU flag
%
% Outputs
%   reconed - Reconctructed data. This is k-space data if p.grappa_domain is
%             'kSpace', and is image space data if p.grappa_domain is 'imSpace'.
%             Dim: [FE, PE, Echo, Slice(solved for slice multiplexing), Coil, Time].
%   sz_out  - Size of 'reconed'. A structure array with fields: x, y, ec, sl, c, t, kz.
%
% (c) Kangrong Zhu,     Stanford University     July 2012

% -- Variables deducted from inputs.
if p.caipi
    if (abs(p.cap_fov_shift) ~= p.mux) && strcmp(p.zpad_1Dgrappa, 'minZpad');
        p.zpad_1Dgrappa = 'normalZpad';
    end
    if strcmp(p.zpad_1Dgrappa, 'minZpad')   % Use minimum zero padding, only zero pad a total of 1*datsz.y pixels in the PE direction.
        mux_R = p.mux + 1;                  % Reduction factor for the 1D GRAPPA recon.
    else                                    % Use an easy way for zero padding, zero pad (p.mux-1)*datsz.y pixels in the PE direction.
        mux_R = 1 + (p.mux-1)*2;
    end
else
    mux_R = p.mux;
end
datsz = get_dat_sz(ksp, p);
if nargin < 3 || ...
   ~iscell(gker) || ...
   ~isequal(size(gker), [datsz.sl, datsz.ec])
    gker = cell(datsz.sl, datsz.ec);
    kgiven = false;
    idxplus = p.mux * p.num_mux_cycle;
    nt_solve = datsz.t - idxplus; % # of temporal phases in 'ksp' to solve for.
else
    kgiven = true;
    idxplus = 0;
    nt_solve = datsz.t;
end

if p.caipi
    if strcmp(p.zpad_1Dgrappa, 'minZpad')   % Minimum zero padding
        ny_zpad = floor(datsz.y/p.mux);                       % # of y lines in the zero matrices for zero padding. Zero padding using matrices can only account for shifts with interger number of pixels.
        ny_zpad_residual_total = datsz.y - p.mux*ny_zpad;     % # of y lines lacking after the zero matrix padding.
        ny_zpad_residual_each = ny_zpad_residual_total/p.mux; % The difference between the fractional number and the floored integer. This will be used to correct shifts with fractional number of pixels.
        slice_shift_amt = sign(p.cap_fov_shift) * (0:p.mux-1) * ny_zpad_residual_each/(datsz.y+ny_zpad); % The amount of shift to apply to each slice, as a fraction of the FOV in y. This will be used to correct shifts with fractional number of pixels.
    else                                    % Easy way for zero padding
        ny_zpad = datsz.y;                                    % # of y lines in the zero matrices for zero padding.
        slice_indices = - floor(p.mux/2) : 1 : (ceil(p.mux/2)-1);       % Slice indices, must center around index 0. [-1,0,1] for mux=3, [-2,-1,0,1] for mux=4. ([-1,0,1] and [2,0,1] are the same when abs(cap_fov_shift)==p.mux, but different when abs(cap_fov_shift)~=p.mux).
        slice_indices = ifftshift(slice_indices);                       % Make the first slice the 0-th indexed slice.
        slice_shift_amt_raw = mod(slice_indices, abs(p.cap_fov_shift)); % [0,1,2] for p.mux=3 and abs(p.cap_fov_shift)=3; [0,1,1] for p.mux=3 and abs(p.cap_fov_shift)=2; [0,1,2,3] for p.mux=4 and abs(p.cap_fov_shift)=4; [0,1,1,2] for p.mux=4 and abs(p.cap_fov_shift)=3.
        slice_shift_amt = slice_shift_amt_raw /(p.cap_fov_shift*2);     % The amount of shift to apply to each slice, as a fraction of the FOV in y. This will be used to synthesize the shifting for the reference image.
    end
end
if nt_solve < 0
    error('The number of temporal phases to solve after mux phase cycling can''t be negative.');
end

% -- The temporal phases for mux phase cycling.
if ~kgiven
    ref_dat = cat(6, ksp{1:p.mux_encoded*p.num_mux_cycle});
    ref_dat = reshape(ref_dat, [datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, p.mux_encoded, p.num_mux_cycle]);
    ref_dat = mux_dftz(ref_dat, p.T_DIM, p.cap_fov_shift_cal, p.mux, 'decode'); % Dim: [Kx, Ky, Echo, Slice, Coil, SimultaneousSlice, MultiplexingPhaseCycle]
    if (p.pseq == 2) && mod(p.mux, 2)==0 && (~p.caipi || (p.caipi && ~p.cal_gzblips))         % Adjust the slice ordering for the non-caipi case, if the number of simultaneous bands is even. Somehow this is needed for the muxarcepi GRE sequence, but not for the muxarcepi2 sequence. (changed between Aug'14 and Oct'14)
        ref_dat = circshift(ref_dat, [0, 0, 0, 0, 0, -1, 0]);
    end
end

% -- GRAPPA reconstruction.
reconed = complex(zeros(datsz.x, datsz.y, datsz.ec, datsz.sl, p.mux, datsz.c, nt_solve));

if nt_solve > 0
    % Raw reference images.
    if ~kgiven
        refic_raw = ifft2c(ref_dat(:, :, :, :, :, :, p.num_mux_cycle), nr_usegpu); % Dim: [FE, PE, Echo, Slice, Coil, SimultaneousSlice, MultiplexingPhaseCycle(=1)]
        if nr_usegpu
            fft(complex(gpuArray.randn(4, 1), gpuArray.randn(4, 1)));
        end
        if p.caipi
            if nr_usegpu
                if p.cap_fov_shift > 0
                    refic_raw = cat(p.PE_DIM, gpuArray.zeros(datsz.x, ny_zpad, datsz.ec, datsz.sl, datsz.c, p.mux), refic_raw);
                else % p.cap_fov_shift <0
                    refic_raw = cat(p.PE_DIM, refic_raw, gpuArray.zeros(datsz.x, ny_zpad, datsz.ec, datsz.sl, datsz.c, p.mux));
                end
            else
                if p.cap_fov_shift > 0
                    refic_raw = cat(p.PE_DIM, zeros(datsz.x, ny_zpad, datsz.ec, datsz.sl, datsz.c, p.mux), refic_raw);
                else % p.cap_fov_shift <0
                    refic_raw = cat(p.PE_DIM, refic_raw, zeros(datsz.x, ny_zpad, datsz.ec, datsz.sl, datsz.c, p.mux));
                end
            end
            refic_raw = mux_shift(refic_raw, p.T_DIM, slice_shift_amt, 'im', nr_usegpu); % Note the negative y direction in matlab is the positive y direction in the physical image space. If (CAIPI && Minimum zero padding), this step is the 1st step to correct the shifts with fractional number of pixels.
        end
    end
   
    % The position for the ACS data.
    acspos = grappa_acspos( struct('x',{datsz.x},'y',{datsz.y*mux_R}), p.mux_acssz,...
        struct('x',{[1,datsz.x]}, 'y',{[1,datsz.y*mux_R]}), false);
    
    % Recon.
    for slice = 1 : datsz.sl
        for echo = 1 : datsz.ec
            if ~kgiven
                refic = refic_raw( :, :, echo, slice, :, :); % Dim: [FE, PE, Echo, Slice, Coil, SimultaneousSlice]
                refic = permute(refic, [1,2,6,5,3,4]);       % Dim: [FE, PE, SimultaneousSlice, Coil]
                if p.caipi && strcmp(p.zpad_1Dgrappa, 'minZpad') && p.cap_fov_shift > 0
                    refic = flip(refic, 3);
                end
                refic = reshape( refic, [datsz.x, size(refic,2)*size(refic,3), datsz.c]); % Dim: CAIPI && (Minimum zero padding): [nx, ny*mux+(ny-ny_zpad_residual_total), nc]; CAIPI && (easy way for zero padding): [nx, ny*mux*2, nc]; Non-CAIPI: [nx, ny*mux, nc].
                if p.caipi
                    if strcmp(p.zpad_1Dgrappa, 'minZpad')    % Minimum zero padding
                        if p.cap_fov_shift > 0
                            refic = cat(p.PE_DIM, zeros(datsz.x, ny_zpad_residual_total, datsz.c), refic); % 2nd step to correct the shifts with fractional number of pixels.
                        else % p.cap_fov_shift < 0
                            refic = cat(p.PE_DIM, refic, zeros(datsz.x, ny_zpad_residual_total, datsz.c));
                        end
                    else                                     % Use an easy way for zero padding
                        if p.cap_fov_shift > 0
                            refic = refic(:, datsz.y+1:end, :);
                        else % p.cap_fov_shift < 0
                            refic = cat(p.PE_DIM, refic(:, 1:datsz.y, :), refic(:, datsz.y*2+1:end, :));
                        end
                    end
                end
                refic = fft2c(refic, nr_usegpu);
                refic = refic(acspos.x, acspos.y, :);
                if nr_usegpu
                    fft(complex(gpuArray.randn(4, 1), gpuArray.randn(4, 1)));
                    ker = gpuArray(grappa_kernel(refic, [p.mux_kersz.x, p.mux_kersz.y], ...  
                        mux_R, p.grappa_domain, [datsz.x, datsz.y*mux_R], [], nr_usegpu));
                    fft(complex(gpuArray.randn(4, 1), gpuArray.randn(4, 1)));
                else
                    ker = grappa_kernel(refic, [p.mux_kersz.x, p.mux_kersz.y], ...  
                        mux_R, p.grappa_domain, [datsz.x, datsz.y*mux_R], [], nr_usegpu);
                end
                gker{slice, echo} = ker;
            else
                ker = gker{slice, echo};
            end
            
            % Solve the multiplexed slices for each temporal frame.
            for tphase = 1 : nt_solve
                if nr_usegpu
                    us_dat = gpuArray.zeros(datsz.x, datsz.y*mux_R, datsz.c);
                else
                    us_dat = zeros(datsz.x, datsz.y*mux_R, datsz.c);
                end
                us_dat(:, 1:mux_R:end, :) = reshape( ksp{idxplus + tphase}(:, :, echo, slice, :), ...
                    [datsz.x, datsz.y, datsz.c]);        % Dim: [FE, PE, Coil]
                
                if mod(mux_R, 2) == 0                    % For even 'mux_R'(must be a non-caipi case or (CAIPI && Minimum zero padding) case), need to adjust the aliasing pattern to match the concatenated reference image.
                    us_dat(:, 1 : mux_R*2 : end, :) = -us_dat(:, 1 : mux_R*2 : end, :); % Dim: [FE, PE, Coil]
                end
            
                method = p.grappa_domain;
%                 if ~exist('method', 'var') || isempty(method)
%                     method = 'imSpace';
%                 end

%                 if ~(strcmp(method, 'kSpace') || strcmp(method, 'imSpace'))
%                     error('''method'' must be either ''kSpace'' or ''imSpace''.');
%                 end

                % -- Data size
                [n_fe, n_pe, n_coil] = size(us_dat);

                % -- GRAPPA reconstruction
                if nr_usegpu
                    rec_dat = complex(gpuArray.zeros(n_fe, n_pe, n_coil));
                else
                    rec_dat = complex(zeros(n_fe, n_pe, n_coil));
                end

                if strcmp(method, 'kSpace')
                    for tgt_coil = 1:n_coil
                        for src_coil = 1:n_coil
                            rec_dat(:, :, tgt_coil) = rec_dat(:, :, tgt_coil) + conv2( ...
                                us_dat(:, :, src_coil), ker(:, :, src_coil, tgt_coil), 'same');
                        end
                    end
                end

                if strcmp(method, 'imSpace')
                    if (size(ker, 1) ~= n_fe) || (size(ker,2) ~= n_pe)
                        error('The input kernel size is not consistent with the undersampled k-space data size.');
                    end

                    us_dat = ifft2c(us_dat, nr_usegpu);                                    % undersampled coil images. The variable name 'us_dat' is kept unchanged to reduce memory use.
                    if nr_usegpu
                        fft(complex(gpuArray.randn(4, 1), gpuArray.randn(4, 1)));
                    end
                    for coil = 1 : n_coil
                        rec_dat(:, :, coil) = sum(us_dat .* ker(:, :, :, coil), 3); % reconstructed coil images
                    end
                end  
             
                if p.caipi
                    if strcmp(p.grappa_domain, 'kSpace')
                        rec_dat = ifft2c(rec_dat, nr_usegpu);                    % Reconstructed coil images
                        if nr_usegpu
                            fft(complex(gpuArray.randn(4, 1), gpuArray.randn(4, 1)));
                        end
                    end
                    
                    if strcmp(p.zpad_1Dgrappa, 'minZpad')% Minimum zero padding
                        if p.cap_fov_shift >0
                            rec_dat = rec_dat(:, ny_zpad_residual_total+1 : end, :); % Undo the effects of the 2nd step to correct the shifts with fractional number of pixels.
                        else % p.cap_fov_shift <0
                            rec_dat = rec_dat(:, 1 : end-ny_zpad_residual_total, :);
                        end
                        rec_dat = reshape(rec_dat, [datsz.x, datsz.y + ny_zpad, p.mux, datsz.c]);
                        if p.cap_fov_shift > 0
                            rec_dat = flip(rec_dat, 3);
                        end
                    else                                 % Use an easy way for zero padding
                        if p.cap_fov_shift > 0
                            rec_dat = cat(p.PE_DIM, zeros(datsz.x, datsz.y, datsz.c), rec_dat);
                        else % p.cap_fov_shift < 0
                            rec_dat = cat(p.PE_DIM, cat(p.PE_DIM, rec_dat(:, 1:datsz.y, :), zeros(datsz.x, datsz.y, datsz.c)), rec_dat(:, datsz.y+1:end, :));
                        end
                        rec_dat = reshape(rec_dat, [datsz.x, 2*datsz.y, p.mux, datsz.c]);
                    end
                   
                    rec_dat = fftshift(mux_shift(rec_dat, 3, -slice_shift_amt, 'im', nr_usegpu), 3);           % If minimum zero padding, this undos the effects of the 1st step to correct the shifts with fractional number of pixels.
                    if nr_usegpu
                        fft(complex(gpuArray.randn(4, 1), gpuArray.randn(4, 1)));
                        if p.cap_fov_shift > 0
                            reconed(:, :, echo, slice, :, :, tphase) = ...
                                gather(rec_dat(:, ny_zpad+1:end, :, :)); % reconed Dim: [FE, PE, Echo, Slice, Mux, Coil, TemporalPhase]
                        else % p.cap_fov_shift < 0
                            reconed(:, :, echo, slice, :, :, tphase) = ...
                                gather(rec_dat(:, 1:datsz.y, :, :));
                        end
                        if strcmp(p.grappa_domain, 'kSpace')
                            reconed(:, :, echo, slice, :, :, tphase) = ...
                                gather(fft2c(reconed(:, :, echo, slice, :, :, tphase), nr_usegpu));
                        end
                    else
                        if p.cap_fov_shift > 0
                            reconed(:, :, echo, slice, :, :, tphase) = ...
                                rec_dat(:, ny_zpad+1:end, :, :); % reconed Dim: [FE, PE, Echo, Slice, Mux, Coil, TemporalPhase]
                        else % p.cap_fov_shift < 0
                            reconed(:, :, echo, slice, :, :, tphase) = ...
                                rec_dat(:, 1:datsz.y, :, :);
                        end
                        if strcmp(p.grappa_domain, 'kSpace')
                            reconed(:, :, echo, slice, :, :, tphase) = ...
                                fft2c(reconed(:, :, echo, slice, :, :, tphase), nr_usegpu);
                        end
                    end
                else     % Non-CAIPI
                    if nr_usegpu
                        if strcmp(p.grappa_domain, 'imSpace')
                            reconed(:, :, echo, slice, :, :, tphase) = ...                              % reconed Dim: [FE, PE, Echo, Slice, Mux, Coil, TemporalPhase]
                                gather(fftshift(reshape(rec_dat, [datsz.x, datsz.y, p.mux, datsz.c]), 3));
                        else % p.grappa_domain is 'kSpace'
                            reconed(:, :, echo, slice, :, :, tphase) = ...
                                gather(fft2c( fftshift(reshape(ifft2c(rec_dat, nr_usegpu), [datsz.x, datsz.y, p.mux, datsz.c]), 3), nr_usegpu));
                        end
                    else
                        if strcmp(p.grappa_domain, 'imSpace')
                            reconed(:, :, echo, slice, :, :, tphase) = ...                              % reconed Dim: [FE, PE, Echo, Slice, Mux, Coil, TemporalPhase]
                                fftshift(reshape(rec_dat, [datsz.x, datsz.y, p.mux, datsz.c]), 3);
                        else % p.grappa_domain is 'kSpace'
                            reconed(:, :, echo, slice, :, :, tphase) = ...
                                fft2c(fftshift(reshape(ifft2c(rec_dat), [datsz.x, datsz.y, p.mux, datsz.c]), 3), nr_usegpu);
                        end
                    end
                end
                if nr_usegpu
                    fft(complex(gpuArray.randn(4, 1), gpuArray.randn(4, 1)));
                end
            end
        end
    end
end

reconed = reconed .* sqrt(p.mux/mux_R);  % Adjust the scaling in ifft2c to account for the upsampling. The scaling in ifft2c should be consistent with the scaling in fft2c.
if ~kgiven
    ref_dat = permute(fftshift(ref_dat,6), [1, 2, 3, 4, 6, 5, 7]);
    if strcmp(p.grappa_domain, 'imSpace')
        ref_dat = ifft2c(ref_dat, nr_usegpu);
    end
    if nr_usegpu
        fft(complex(gpuArray.randn(4, 1), gpuArray.randn(4, 1)));
        reconed = cat(7, gather(ref_dat), reconed);
    else
        reconed = cat(7, ref_dat, reconed);
    end
    reconed = reshape( reconed, [datsz.x, datsz.y, datsz.ec, datsz.sl*p.mux, ...
        datsz.c, size(reconed, 7)]); % Dim: [FE, PE, Echo, Slice, Coil, TemporalPhase]
else
    reconed = reshape( reconed, [datsz.x, datsz.y, datsz.ec, datsz.sl*p.mux, ...
        datsz.c, nt_solve]); % Dim: [FE, PE, Echo, Slice, Coil, TemporalPhase]
end



function dat = homodyne(dat, ny_part, ntran, niter, in_dom, out_dom, nr_usegpu)
%
% function dat = homodyne(dat, ny_part, [ntran=2], [niter=4], [in_dom='KxKy'], [out_dom=in_dom], nr_usegpu)
%
% Homodyne detection reconstruction of partially acquired k-space.
%
% Inputs
%   dat     - Partially acquired k-space, with full matrix size. The
%             acquired k-space data is FOLLOWED by zeros in the PE direction.
%             Dim: [FE, PE, Echo, Slice, Coil, TemporalPhase].
%   ny_part - # of acquired ky lines in the partial k-space acquisition.
%   ntran   - Number of points in the transition band of the filters.
%   niter   - Number of iterations.
%             0: Use zero-filling.
%             1,2,3...: Use Homodyne. Iterative when niter > 1.
%   in_dom  - A string specifying the domain of the input 'dat'. 'xy', 'xKy', 'KxY' or 'KxKy'(uppercase or lowercase doesn't matter).
%   out_dom - A string specifying the domain of the output 'dat'.
%
% Output
%   dat     - Reconstructed k-space data, having the same size and in the same domain as the input 'dat'.
%
% From GE company. Based on Homodyne.cpp.
% Modified by Robert Dougherty and Kangrong Zhu, Stanford University, July 2014

% persistent gpuArray
persistent gpuad;
if isempty(gpuad)
    gpuad = {-1, [], [], []};
end

if ~exist('niter', 'var') || isempty(niter)
    niter = 4;
end
if ~exist('ntran', 'var') || isempty(ntran)
    ntran = 2;
end
if ~exist('in_dom', 'var') || isempty(in_dom)
    in_dom = 'KxKy';
end
if ~exist('out_dom', 'var') || isempty(out_dom)
    out_dom = in_dom;
end
in_dom = lower(in_dom);
out_dom = lower(out_dom);

FE_DIM = 1;
PE_DIM = 2;

% Data size
datsz = size(dat);
ny_pres = datsz(PE_DIM);
if length(datsz) < 7
    datsz(end+1 : 7) = 1;        % To make sure datsz(7) exists when it is called below
end

% Transform input data into (x/Kx, Ky)
if any(strcmp(in_dom, {'kxy', 'xy'}))
    if gpuad{1} ~= size(dat, PE_DIM)
        if nr_usegpu
            gpuad = {size(dat, PE_DIM), gpuArray(sqrt(1 / size(dat, PE_DIM))), ...
                fftshift(gpuArray((1:ny_pres)')), ifftshift(gpuArray((1:ny_pres)'))};
        else
            gpuad = {size(dat, PE_DIM), sqrt(1 / size(dat, PE_DIM)), ...
                fftshift((1:ny_pres)'), ifftshift((1:ny_pres)')};
        end
    end
    dat = dat(:, gpuad{4}, :, :, :, :, :);
    dat = fft(dat, [], PE_DIM);
    dat = dat(:, gpuad{3}, :, :, :, :, :);
    dat = gpuad{2} .* dat;
end

% Keep a copy of the input data in (x, Ky)
dat_input = dat;                                   % If in_dom is 'xy' or 'xky', dat_input is in (x, Ky); If in_dom is 'kxy' or 'kxky', dat_input is in (Kx, Ky)
if strcmp(in_dom, 'kxy') || strcmp(in_dom, 'kxky')
    dat_input = ifftc(dat_input, FE_DIM, nr_usegpu);          % dat_input is in (x, Ky)
end

% Generate filters
[lowPassFilter, highPassFilter, mergeWindow] = generate_filters(ntran, ny_part, ny_pres, nr_usegpu);

lowPassFilter(ny_part+1 : ny_pres) = 0.0;
highPassFilter(ny_part+1 : ny_pres) = 0.0;

lowPassFilter = repmat( lowPassFilter,[size(dat,1) 1 size(dat,3) size(dat,4) size(dat,5)]);
highPassFilter = repmat(highPassFilter,[size(dat,1) 1 size(dat,3) size(dat,4) size(dat,5)]);

% Calculate high and low pass images in (x/Kx, Ky)
if nr_usegpu
    lowPass = gpuArray.zeros(datsz);
else
    lowPass = zeros(datsz);
end
for i67 = 1 : prod(datsz(6:7))
    % Get low pass views
    lowPass(:, :, :, :, :, i67)  = dat(:, :, :, :, :, i67) .* lowPassFilter;
    % Get high pass views
    dat(:, :, :, :, :, i67) = dat(:, :, :, :, :, i67) .* highPassFilter;
end

% Inverse FT both high and low pass images into (x, y)
if strcmp(in_dom, 'xy') || strcmp(in_dom, 'xky')        % dat in (x, Ky)
    lowPass = ifftc(lowPass, PE_DIM, nr_usegpu);
    dat = ifftc(dat, PE_DIM, nr_usegpu);
elseif strcmp(in_dom, 'kxy') || strcmp(in_dom, 'kxky') % dat in (Kx, Ky)
    lowPass = ifft2c(lowPass, nr_usegpu);
    dat = ifft2c(dat, nr_usegpu);
end

% Apply Phase Correction to high pass images - use real part only
scalar = 0.000001;
phaseCorrection = abs(lowPass);
scalar = lowPass + scalar;
lpImagePhase = scalar ./ phaseCorrection;
phaseCorrection = phaseCorrection ./ scalar;
dat = dat .* phaseCorrection;
dat = real(dat); % real part

i1mergeWindow = repmat(mergeWindow, [datsz(1), 1, datsz(3:7)]);
for iter = 1 : niter
    % Reinsert phase
    dat = dat .* lpImagePhase;
    
    % Transform high-pass image back to (x, Ky)
    dat = sqrt(1 / size(dat, PE_DIM)) .* fftshift(fft(ifftshift(dat, PE_DIM), [], PE_DIM), PE_DIM);
    
    % Here the output HP image is as follows:
    % Orig*Merge + HP*(1-Merge) = (Orig-Hp)*Merge + HP;
    dat = (dat_input - dat) .* i1mergeWindow + dat;
    
    % Transform high-pass image back to (x,y)
    dat = ifftc(dat, PE_DIM, nr_usegpu);
    
    % Correct phase
    dat = real(dat .* phaseCorrection);
end

% Load output
switch out_dom
    case {'kxy', 'xky'}
        dat = sqrt(1 / size(dat, PE_DIM)) .* fftshift(fft(ifftshift(dat, PE_DIM), [], PE_DIM), PE_DIM);
    case 'kxky'
        dat = fft2c(dat, nr_usegpu);
end



function [lowPassFilter, highPassFilter, mergeWindow] = generate_filters(ntran, nacq, npres, nr_usegpu)
%
% function [lowPassFilter, highPassFilter, mergeWindow] = generate_filters(ntran, nacq, npres, nr_usegpu)
%
% Inputs
%   ntran          - Number of points in transition band.
%   nacq           - Number of partial-k acquired lines.
%   npres          - Number of total prescribed lines.
%
% Outputs
%   lowPassFilter  - Low-pass filter.
%   highPassFilter - High-pass filter.
%   mergeWindow    - Window for data merging.
%
% Generate filters for Homodyne detection reconstruction.
%

symmetricSamples = nacq - npres/2;
transitionWidth = ntran;
mergeWindowSize = npres;
filterSize = nacq;
dataCenter = filterSize - symmetricSamples;
dataOffset = mergeWindowSize/2 - dataCenter;

windowCenter2 = filterSize - (2*transitionWidth);
windowCenter1 = dataCenter - (windowCenter2 - dataCenter);

if nr_usegpu
    exponential1 = gpuArray(generate_exponential(filterSize, windowCenter1, transitionWidth));
    exponential2 = gpuArray(generate_exponential(filterSize, windowCenter2, transitionWidth));
else
    exponential1 = generate_exponential(filterSize, windowCenter1, transitionWidth);
    exponential2 = generate_exponential(filterSize, windowCenter2, transitionWidth);
end

lowPassFilter = exponential2 - exponential1;
highPassFilter = exponential2 + exponential1;

if nr_usegpu
    mergeWindow = gpuArray.zeros(1, mergeWindowSize);
else
    mergeWindow = zeros(1, mergeWindowSize);
end
mergeWindow((dataOffset : dataOffset + length(exponential2) -1) + 1) = exponential2;


function vector = generate_exponential(vecSize,winWidth,transWidth)
%
% function vector = generate_exponential(vecSize,winWidth,transWidth)
%
% Generate exponentials for Homodyne detection reconstruction.
%

lowThreshold = -87.0;
highThreshold = 87.0;

vector = 1:vecSize;
vector = vector - winWidth;
vector = vector ./ transWidth;

vector(vector>highThreshold) = highThreshold;
vector(vector<lowThreshold) = lowThreshold;

vector = exp(vector) + 1.0;
vector =  1.0./vector;



function dat_out = mux_shift(dat, shift_dim, shift_amt, dat_type, nr_usegpu)
%
% function dat_out = mux_shift(dat, shift_dim, shift_amt [, dat_type], nr_usegpu)
%
% Shift each slice in the specified dimension by a certain amount in the PE direction.
%
% Inputs : 
%   dat       - Input k-space or image-space data. Dim: [FE, PE, ...].
%   shift_dim - The dimension for the slice shifting.
%   shift_amt - A vector specifying the amount of shift which will be
%               applied to each slice. The m-th (m=0,1,2...) slice in the 
%               specified 'shift_dim' dimension will be shifted by 
%               shift_amt(m)*FOVy in the image space in the PE direction.
%   dat_type  - 'ksp'(Default): 'dat' is k-space data.
%               'im'          : 'dat' is image-space data.
%
% Output: 
%   dat_out   - K-space or image-space data after slice shifting, having
%                the same size and data type as the input 'dat'.
%
% (c) Kangrong Zhu,     Stanford University,    Oct 2012

if ~exist('dat_type', 'var') || isempty(dat_type)
    dat_type = 'ksp';
end

if strcmp(dat_type, 'im')
    dat = fft2c(dat, nr_usegpu);
else if ~strcmp(dat_type, 'ksp')
        error('The input ''dat_type'' must be either ''ksp'' or ''im''.');
    end
end

PE_DIM = 2;                   % The dimension for phase encoding in the data matrix
KEEP_ORIG_SZ = 1;             % Keep original matrix size when using 'repmat'

sz = size(dat);
ny = sz(PE_DIM);              % Number of ky lines.
nim = sz(shift_dim);          % Number of images to shift.
ndim = length(sz);            % Number of dimensions in 'dat'.

shift_amt = (shift_amt(:)).'; % Make shift_amt a row vector

% Basic phases to apply.
if nr_usegpu
    pha = gpuArray(2*pi * (1:ny).'); % Applying a phase of -n*pha in the PE direction in k-space will shift the image by +n*FOVy in the image space.
    pha = repmat(pha, [KEEP_ORIG_SZ, nim]) .* repmat(shift_amt, [ny, KEEP_ORIG_SZ]);
else
    pha = repmat(2*pi * (1:ny).', [KEEP_ORIG_SZ, nim]) .* repmat(shift_amt, [ny, KEEP_ORIG_SZ]);
end

% Reshape and repeat to the same size as the input 'dat'.
reshape_sz = ones(1, ndim);
reshape_sz(PE_DIM) = ny;
reshape_sz(shift_dim) = nim;

rep_sz = sz;
rep_sz(PE_DIM) = KEEP_ORIG_SZ;
rep_sz(shift_dim) = KEEP_ORIG_SZ;

pha = repmat( reshape(pha, reshape_sz), rep_sz);

% Shift the images
if nr_usegpu
    dat_out = dat .* exp(gpuArray(1i) .* pha);
else
    dat_out = dat .* exp(1i .* pha);
end

if strcmp(dat_type, 'im')
    dat_out = ifft2c(dat_out, nr_usegpu);
    if nr_usegpu
        fft(complex(gpuArray.randn(4, 1), gpuArray.randn(4, 1)));
    end
end



function x = fftc(x, dim, nr_usegpu)
if iscell(x)
    if nr_usegpu
        for cc = 1:numel(x)
            x{cc} = sqrt(gpuArray(1 / size(x{cc}, dim))) .* fftshift(fft(ifftshift(x{cc}, dim), [], dim), dim);
        end
    else
        for cc = 1:numel(x)
            x{cc} = sqrt(1 / size(x{cc}, dim)) .* fftshift(fft(ifftshift(x{cc}, dim), [], dim), dim);
        end
    end
else
    if nr_usegpu
        x = sqrt(gpuArray(1 / size(x, dim))) .* fftshift(fft(ifftshift(x, dim), [], dim), dim);
    else
        x = sqrt(1 / size(x, dim)) .* fftshift(fft(ifftshift(x, dim), [], dim), dim);
    end
end

function x = fft2c(x, nr_usegpu)
if iscell(x)
    for cc = 1:numel(x)
        sz = size(x{cc});
        if nr_usegpu
            x{cc} = sqrt(gpuArray(1 / prod(sz(1:2)))) .* fftshift(fft2(ifftshift(x{cc})));
        else
            x{cc} = sqrt(1 / prod(sz(1:2))) .* fftshift(fft2(ifftshift(x{cc})));
        end
    end
else
    sz = size(x);
    if nr_usegpu
        x = sqrt(gpuArray(1 / prod(sz(1:2)))) .* fftshift(fft2(ifftshift(x)));
    else
        x = sqrt(1 / prod(sz(1:2))) .* fftshift(fft2(ifftshift(x)));
    end
end

function x = ifftc(x, dim, nr_usegpu)
if iscell(x)
    if nr_usegpu
        for cc = 1:numel(x)
            x{cc} = sqrt(gpuArray(size(x{cc}, dim))) .* fftshift(ifft(ifftshift(x{cc}, dim), [], dim), dim);
        end
    else
        for cc = 1:numel(x)
            x{cc} = sqrt(size(x{cc}, dim)) .* fftshift(ifft(ifftshift(x{cc}, dim), [], dim), dim);
        end
    end
else
    if nr_usegpu
        x = sqrt(gpuArray(size(x, dim))) .* fftshift(ifft(ifftshift(x, dim), [], dim), dim);
    else
        x = sqrt(size(x, dim)) .* fftshift(ifft(ifftshift(x, dim), [], dim), dim);
    end
end

function x = ifft2c(x, nr_usegpu)
if iscell(x)
    for cc = 1:numel(x)
        sz = size(x{cc});
        if nr_usegpu
            x{cc} = sqrt(gpuArray(prod(sz(1:2)))) .* ifftshift(ifft2(fftshift(x{cc})));
        else
            x{cc} = sqrt(prod(sz(1:2))) .* ifftshift(ifft2(fftshift(x{cc})));
        end
    end
else
    sz = size(x);
    if nr_usegpu
        x = sqrt(gpuArray(prod(sz(1:2)))) .* ifftshift(ifft2(fftshift(x)));
    else
        x = sqrt(prod(sz(1:2))) .* ifftshift(ifft2(fftshift(x)));
    end
end



function ker = grappa_kernel(acs_dat, kersz, R, method, kspsz, kersz_min, nr_usegpu)
%
% function ker = grappa_kernel(acs_dat, kersz, R [, method, kspsz, kersz_min], nr_usegpu)
%
% Calculate reconstruction kernel for 1D GRAPPA.
%
% Inputs:
%   acs_dat   - AutoCalibrating Signal in k-space. Dim: [FE, PE, n_coil].
%   kersz     - KERnel SiZe, representing # of acquired points used in each
%               dimension to interpolate one missing point. Dim: [FE, PE].
%   R         - Reduction factor in the phase encoding direction.
%   method    - 'imSpace'(Default): the reconstruction will be carried out 
%                                   by multiplication in the image space.
%               'kSpace': the reconstruction will be carried out by 
%                         convolution in the k-space.
%   kspsz     - The full size of the k-space data, used only when 'method' 
%               is 'imSpace'. Dim: [n_fe(full), n_pe(full), n_coil].
%   kersz_min - If there aren't enough ACS points for the specified 'kersz', 
%               'kersz' will be reduced automatically. 'kersz_min' is the 
%               minimum kernel size allowed when 'kersz' is reduced. 
%               Default: [3, 2], dim: [FE, PE].
% Output:
%   ker       - The interpolation KERnel.
%               Dim: [FE(source), PE(source), n_coil(source), n_coil(target)].
%               If 'method' is 'kSpace', size of 'ker' is [kersz(1),
%               R*ceil(kersz(2)/2)*2-1, n_coil, n_coil].
%               If 'method' is 'imSpace', size of 'ker' is [kspsz(1),
%               kspsz(2), n_coil, n_coil].
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

n_coil = size(acs_dat,3);

if ~exist('kersz_min', 'var') || isempty(kersz_min)
    kersz_min = [3, 2];
end

if strcmp(method, 'imSpace') && (~exist('kspsz', 'var') || isempty(kspsz))
    error('Must specify the full k-space size if the reconstruction will be conducted in the image space.');
end

acs_sz = [size(acs_dat,1), size(acs_dat,2)]; % Dim: [FE, PE]
[n_instances, acs_blk_sz] = calc_instances(kersz, acs_sz, R);

if ( sum(acs_sz < acs_blk_sz) ) > 0
    fprintf('Not enough ACS lines in %d of the FE and PE dimensions.\n',...
        sum(acs_sz < acs_blk_sz) );
    error('Not enough ACS lines.');
end

% -- Reduce the kernel size if there aren't enough ACS points for the specified kernel size.
if n_instances < prod(kersz)*n_coil
    fprintf(['   The specified kernel size is [%d, %d] in [FE, PE]. ' ...
        'Not enough ACS points for the specified kernel size.\n'], kersz(1), kersz(2) );
    fprintf('   Reducing the kernel size...\n');
    
    for dim = 1 : 2                 % First reduce kersz in FE, then in PE
        while (n_instances < prod(kersz)*n_coil) && (kersz(dim) > kersz_min(dim))
            kersz(dim) = kersz(dim) - 1;
            [n_instances, acs_blk_sz] = calc_instances(kersz, acs_sz, R);
        end
    end
    
    if n_instances < prod(kersz)*n_coil
        fprintf('   The minimum kernel size allowed is [%d, %d] in [FE, PE].\n', kersz_min(1), kersz_min(2));
        fprintf('   Total number of instances available in the ACS area to fit the minimum kernel is %d.\n', n_instances);
        fprintf('   Total number of coefficients to fit for the minimum kernel is prod(kersz)*n_coil = %d.\n', prod(kersz)*n_coil);
        error('Unable to estimate the minimum kersz size allowed with the current ACS area size. Please increase the ACS area.');
    else
        fprintf('   The reduced kernel size is [%d, %d] in [FE, PE].\n', kersz(1), kersz(2));
    end
end

% -- Set up a matrix with all of the calibration data.
if nr_usegpu
    acs_mtx = cell(n_coil, 1);
    % fprintf('the total coils are %d and siz(acs_mtx)=%d\n', n_coil, numel(acs_mtx));

    for idx = 1 : n_coil
         acs_mtx{idx} = gather(im2colgpu(acs_dat(:,:,idx), acs_blk_sz, 'sliding'));        % Dim: [acs_blk_sz(1)*acs_blk_sz(2), n_instances, n_coil]
    end

    acs_mtx = reshape(cat(3, acs_mtx{:}), [acs_blk_sz(1), acs_blk_sz(2), n_instances, n_coil]); % Dim: [acs_blk_sz(1)(i.e. kersz(1)), acs_blk_sz(2), n_instances, n_coil]
    
    % -- Source matrix.
    src_mtx = gpuArray(acs_mtx(:, 1:R:end, :, :));                                             % Dim: [kersz(1), kersz(2), n_instances, n_coil]
else
    acs_mtx = im2col(acs_dat(:,:,1), acs_blk_sz, 'sliding');
    acs_mtx(1, 1, n_coil) = 0;
    for idx = 2 : n_coil
         acs_mtx(:, :, idx) = im2col(acs_dat(:,:,idx), acs_blk_sz, 'sliding');        % Dim: [acs_blk_sz(1)*acs_blk_sz(2), n_instances, n_coil]
    end

    acs_mtx = reshape(acs_mtx, [acs_blk_sz(1), acs_blk_sz(2), n_instances, n_coil]); % Dim: [acs_blk_sz(1)(i.e. kersz(1)), acs_blk_sz(2), n_instances, n_coil]
    
    % -- Source matrix.
    src_mtx = acs_mtx(:, 1:R:end, :, :);                                             % Dim: [kersz(1), kersz(2), n_instances, n_coil]
end
src_mtx = reshape(permute(src_mtx,[3,1,2,4]), [n_instances, kersz(1)*kersz(2)*n_coil]); % Dim: [n_instances, kersz(1)*kersz(2)*n_coil]

% -- Target matrix.
tgt_fe_pos = floor(kersz(1)/2)+1;
tgt_pe_pos = R*floor((kersz(2)-1)/2)+2 : R*floor((kersz(2)-1)/2)+R;
if nr_usegpu
    tgt_mtx = gpuArray(acs_mtx(tgt_fe_pos, tgt_pe_pos, :, :));                            % Dim: [1, R-1, n_instances, n_coil]
else
    tgt_mtx = acs_mtx(tgt_fe_pos, tgt_pe_pos, :, :);                            % Dim: [1, R-1, n_instances, n_coil]
end
tgt_mtx = reshape(permute(tgt_mtx,[3,2,4,1]), [n_instances, (R-1)*n_coil]); % Dim: [n_instances, (R-1)*n_coil]

acs_mtx = [];

% -- Solve the kernel.
% -- Raw kernel
AtA = src_mtx' * src_mtx; 
lambda = norm(AtA, 'fro') / size(src_mtx,2) * 0.02;            % Regularization parameter, from Miki Lustig's code
if nr_usegpu
    ker = gather(inv(AtA + eye(size(AtA))*lambda) * src_mtx' * tgt_mtx);   % Dim: [kersz(1)*kersz(2)*n_coil, (R-1)*n_coil]
else
    ker = inv(AtA + eye(size(AtA))*lambda) * src_mtx' * tgt_mtx;   % Dim: [kersz(1)*kersz(2)*n_coil, (R-1)*n_coil]
end

acs_dat = [];
AtA = [];
lambda = [];
src_mtx = [];
tgt_mtx = [];

ker = reshape(ker, [kersz(1), kersz(2), n_coil, R-1, n_coil]); % Dim: [kersz(1)(source), kersz(2)(source), n_coil(source), (R-1)(target), n_coil(target)]

% -- Autocorrelation kernel
conv_ker = zeros(kersz(1), R*ceil(kersz(2)/2)*2-1, n_coil, n_coil); % Dim: [FE(src), PE(src), n_coil(src), n_coil(tgt)]
for tgt_coil = 1 : n_coil
    for tgt_pe = 1 : R-1
        conv_ker( :, R-tgt_pe:R:R*kersz(2)-1, :, tgt_coil ) = ker(:, :, :, tgt_pe, tgt_coil);
    end
    conv_ker( floor(kersz(1)/2)+1, R*floor((kersz(2)-1)/2)+R, ...
        tgt_coil, tgt_coil) = 1;
end

% -- Convolution kernel
conv_ker = flip(flip(conv_ker, 1), 2); % Dim: [FE(src), PE(src), n_coil(src), n_coil(tgt)]

if strcmp(method, 'kSpace')
    ker = conv_ker;                              % Dim: [FE(src), PE(src), n_coil(src), n_coil(tgt)]
end

% -- Image space multiplication kernel
if strcmp(method, 'imSpace')
    ker = zpad(sqrt(kspsz(1)*kspsz(2))*conv_ker, [kspsz(1), kspsz(2), n_coil, n_coil]);
    conv_ker = [];

%    ker = ifft2c(ker, nr_usegpu);                           % Dim: [FE(src), PE(src), n_coil(src), n_coil(tgt)]
    sz = size(ker);
%    n_2d = prod( sz(3 : end) ); % # of 2D arrays to perform 2D iFFT on.

%    for idx = 1 : n_2d
    ker = sqrt(prod(sz(1:2))) * ifftshift(ifft2(fftshift(ker)));
%    end
end



function [n_instances, acs_blk_sz] = calc_instances(kersz, acs_sz, R)
%
% Calculates the total number of instances for kernel fitting in GRAPPA.
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

acs_blk_sz = [kersz(1), R*(kersz(2)-1) + 1];   % Size of one block for kernel fitting, dim: [FE, PE]
n_instances_indv = acs_sz - acs_blk_sz + 1;    % # of instances for kernel fitting in [FE,PE]
n_instances = prod(n_instances_indv);          % Total # of instances for kernel fitting



function [dat, p] = epi_load_tseries(pfile, ref, vrgf, refp, slices, slices_pfile, slices_ref, tpointstart, tpointend, cali_points, p, nr_usegpu)
%
% function [dat, p] = epi_load_tseries(pfile, ref, vrgf, refp, slices, slices_pfile, slices_ref, p, nr_usegpu)
%
% Load and correct the data for an EPI time series.
%
% Inputs:
%   pfile        - Path to pfile.
%   ref          - Path to 'ref.dat'.
%   vrgf         - Path to 'vrgf.dat'.
%   refp         - Path to reference scan pfile.
%   slices       - Desired slices(Normally ordered indices) to load. Default: all slices.
%   slices_pfile - Slice indices in the p-file corresponding to 'slices'.
%                  'slices_pfile' and 'slices' are different for interleaved
%                  slice acquisition. Default: the same as 'slices'.
%   slices_ref   - Slice indices(Normally ordered indices, not pfile indices)
%                  in the reference scan corresponding to 'slices' in the actual scan.
%                  Only needed if the eddy current correction algorithm will use reference scan pfile.
%   tpointstart  - The index of the start volume to be reconstructed
%   tpointend    - The index of the end volume to be reconstructed, if set as default value 0,
%                  then all the volumes inside the pfile will be reconstructed
%   cali_points  - Indicate either external calibration (-1) or internal calibration (>0), and if internally calibrated, 
%                  the number of calibration volumes
%   p            - Parameter structure. See mux_epi_params.m for details.
%                  The following fields are used in this function: ychop,
%                  frames, echoes, coils, num_slices, num_coils, nx_pres,
%                  tpoints_to_load, pfile_header(used in function load_raw_tseries).
%
% Outputs:
%   dat          - K-space data. Dim: [FE, PE, Echo, Slice, Coil, Time].
%   p            - The output parameter structure, having the same fields as the input 'p'.
%                  The following fields might be changed inside this function if 
%                  the scan was aborted midway and the p-file was only
%                  partially collected: num_mux_cycle, num_passes, 
%                  nt_to_recon, tpoints_to_load.
%                  The following field will be added if p.md_ecc==true:
%                  pha_coe(EPI x-ky phase correction coefficients, will be
%                  passed to mux_epi_process_data_sense.m. Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Echo,
%                  Slice, SimultaneousSlice Z(indices [ifftshift(-floor(mux/2):1:ceil(mux/2)-1])),
%                  Coil(=1 if p.coil_compress==true; =nc if p.coil_compress==false)]
%                  The following fields will be changed if ((p.mux_encoded
%                  > p.mux_excited) && ~p.recon_cal_use_mux_extra):
%                  cap_fov_shift_cal, mux_encoded.
%   dat_ecc      - if p.cap_get_ecc==true, this is the reference k-space data for the
%                  eddy current correction which is collected by the mux sequence.
%                  Dim: [Kx, Ky, Echo, Slice*SimultaneousSlice, Coil(=nc), Time(=1)]
%                  if p.cap_get_ecc==false, this is [].
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

% -- Load k-space data
tstart = tpointstart;
tend = tpointend;
if tstart < 1
    tstart = 1;
end
if cali_points > 0
    tstart = tstart + cali_points;
    tend = tend + cali_points;
end
if(tend > numel(p.tpoints_to_load))
     tend = numel(p.tpoints_to_load);
end

% Dim: [FE, PE, Echo, Slice, Coil, Time]
[dat, p] = load_raw_tseries(pfile, p, slices_pfile, tstart:tend);

if cali_points > 0
    dat = cat(1, load_raw_tseries(pfile, p, slices_pfile, 1:cali_points), dat);
end

datsz = get_dat_sz(dat, p);

if p.vrgf && (datsz.x == p.nx_pres)
    error('Header indicates ramp sampling on, but data are not.');
end

% -- Account for y chopping
if p.ychop
    for dc = 1:numel(dat)
        dat{dc}(:, 2:2:end, :, :, :) = -dat{dc}(:, 2:2:end, :, :, :);
    end
end

% to GPU?
if nr_usegpu
    for dc = 1:numel(dat)
        dat{dc} = gpuArray(dat{dc});
    end
end

% ghost correction
switch p.ref_for_default_ecc
    case 'ref pfile'
        dat = default_ecc(dat, p, refp, slices_ref, nr_usegpu);                % Use reference scan pfile in default correction
    case 'ref.dat'
        dat = default_ecc(dat, p, ref, slices, nr_usegpu);                     % Use ref.dat file in default correction
end

% -- Ramp sampling On
if p.vrgf
    
    % -- Ramp sampling correction.
    ramp_flt = rawload_vrgf(datsz.x, p.nx_pres, vrgf);                  % Load the filter contained in the 'vrgf.dat' file.
    
    if size(ramp_flt, 1) ~= p.nx_pres;
        error('Wrong size of ramp filter.');
    end
    dat = epi_vrgf_correct(dat, ramp_flt, nr_usegpu);
    datsz.x = p.nx_pres;
    
    % -- Nyquist ghost correction: matrix decoding correction for ramp on data
    if p.md_ecc
        [dat, p.pha_coe] = md_ecc(dat, p, md_ref, pha_coe_default, ramp_flt, slices_ref);
    end
end



% function [dat, pha_coe_acc] = md_ecc(dat, p, ref, pha_coe_1stpass, ramp_flt, slices_ref)
% %
% % Matrix-decoding eddy current correction. Corrects fully sampled calibration data
% % and returns the x-ky phase correction coefficients for accelerated data.
% %
% % Inputs
% %   dat             - K-space data, including both mux phase cycling and accelerated data.
% %                     Dim: [Kx(=nx), Ky(=ny), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=mux*num_mux_cycle+num_accelerated_tpoints)].
% %   p               - Parameter structure. See mux_epi_params.m for details.
% %   ref             - Reference scan k-space data matrix(Dim: [Kx(=nx), Ky(=ny), Echo(=1), Slice(=nsl*mux), Coil(=nc)]),
% %                     or ref.dat file name('*.dat'), or reference scan pfile name('*refscan.7').
% %   pha_coe_1stpass - First-pass phase correction coefficients. See rawload_ref_data.m for details.
% %   ramp_flt        - First-pass ramp sampling correction filter. See rawload_ref_data.m for details.
% %   slices_ref      - Slice indices(Normally order indices, not pfile indices) in reference scan which correspond to the input 'dat'.
% %                     Only needed if the type of 'ref' is reference scan pfile name.
% %
% % Outputs
% %   dat             - K-space data with the calibration data corrected by the matrix decoding correction.
% %                     Dim: [Kx(=nx), Ky(=ny), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=mux*num_mux_cycle+num_accelerated_tpoints)].
% %   pha_coe_acc     - EPI x-ky phase correction coefficients for the accelerated data. Dim: [2(0thOrderPhase, 1stOrderPhase), Ky(=ny), Echo(=nec),
% %                     Slice(=nsl), SimultaneousSlice Z(=mux, indices [ifftshift(-floor(mux/2):1:ceil(mux/2)-1])),
% %                     Coil(=1 if p.coil_compress==true; =nc Otherwise)].
% %
% % (c) Kangrong Zhu,     Stanford University     Nov 2013
% 
% MUX_DIM = 5;
% 
% datsz = get_dat_sz(dat, p);
% ref_type = get_ref_type(ref);
% 
% if (~exist('slices_ref', 'var') || isempty(slices_ref)) && strcmp(ref_type, 'ref_pfile')
%     error('No slice indices for reference scan pfile');
% end
% 
% if strcmp(ref_type, 'ksp_data')
%     if p.inplane_R > 1                            % Multi-shot
%         ny_pershot = size(ref, p.PE_DIM) / p.inplane_R;
%         for shot = 2 : p.inplane_R
%             ref(:, shot:ny_pershot:end, :, :, :) = ref(:, 1:ny_pershot:end, :, :, :); % Repeat because process_refdat.m averages over shots
%         end
%     end
% end
% 
% % Load eddy current correction coefficients
% pccoil_cal = 0;
% do_quad_final_cal = true;
% if p.coil_compress
%     pccoil_acc = -1;                              % Average coefficients across coils to get one set of coefficients for all coils if coil compression will be used
% else
%     pccoil_acc = 0;                               % No averaging across coil
% end
% do_quad_final_acc = p.smooth_pha_coe_acc;
% debug = false;
% switch ref_type
%     case 'ref_pfile'
%         if p.debug
%             fprintf('   Using reference scan file in matrix-decoding ghosting correction.\n');
%         end
%         [pha_coe_cal, descend_acq_ref] = rawload_ref_pfile(p.frames, slices_ref, p.coils, ref, pha_coe_1stpass, ramp_flt, pccoil_cal, do_quad_final_cal, debug); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Slice(num_muxed_slice*mux), Coil(=nc)]
%         pha_coe_acc = rawload_ref_pfile(p.frames, slices_ref, p.coils, ref, pha_coe_1stpass, ramp_flt, pccoil_acc, do_quad_final_acc, debug); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Slice(num_muxed_slice*mux), Coil(=1 if p.coil_compress==true; =nc Otherwise)]
%     case 'ksp_data'
%         if p.debug
%             fprintf('   Using ecc data collected by mux sequence in matrix-decoding ghosting correction.\n');
%         end
%         pha_coe_cal = rawload_ref_data(ref, p, pha_coe_1stpass, ramp_flt, pccoil_cal, do_quad_final_cal, debug); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Slice(num_muxed_slice*mux), Coil(=nc)]
%         pha_coe_acc = rawload_ref_data(ref, p, pha_coe_1stpass, ramp_flt, pccoil_acc, do_quad_final_acc, debug); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Slice(num_muxed_slice*mux), Coil(=1 if p.coil_compress==true; =nc Otherwise]
% end
% 
% % Eddy current correction coefficients for accelerated data
% pha_coe_acc = reshape(pha_coe_acc, [p.NUM_COE_PER_KY_LINE, size(pha_coe_acc, 2), 1, size(pha_coe_acc,3)/p.mux, p.mux, size(pha_coe_acc,4)]); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Echo, Slice, SimultaneousSlice Z(indices [-floor(mux/2):1:ceil(mux/2)-1] for ascending acquisition; [ceil(mux/2)-1:-1:-floor(mux/2)] for descending acquisition), Coil(=1 if p.coil_compress==true; =nc Otherwise)]
% if strcmp(ref_type, 'refp_file') && descend_acq_ref
%     pha_coe_acc = flipdim(pha_coe_acc, MUX_DIM);  % Make the slice indices from negative to positive, for mux=3, the slice ordering is [-1, 0, 1].
% end
% pha_coe_acc = ifftshift(pha_coe_acc, MUX_DIM);    % Make the 0-th indexed slice the first slice, for mux=3, the slice ordering is [0, 1, -1]. pha_coe_acc Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Echo, Slice, SimultaneousSlice Z(indices [ifftshift(-floor(mux/2):1:ceil(mux/2)-1])), Coil(=1 if p.coil_compress==true; =nc Otherwise)]
% 
% % Eddy current correction coefficients for calibration data
% pha_coe_cal = reshape(pha_coe_cal, [p.NUM_COE_PER_KY_LINE, size(pha_coe_cal, 2), 1, size(pha_coe_cal,3)/p.mux, p.mux, size(pha_coe_cal,4)]); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Echo, Slice, SimultaneousSlice Z(indices [-floor(mux/2):1:ceil(mux/2)-1]), Coil(=nc)]
% if strcmp(ref_type, 'refp_file') && descend_acq_ref
%     pha_coe_cal = flipdim(pha_coe_cal, MUX_DIM);
% end
% pha_coe_cal = ifftshift(pha_coe_cal, MUX_DIM);    % For mux=3, the slice ordering is [0, 1, -1]. pha_coe_cal Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Echo, Slice, SimultaneousSlice Z(indices [ifftshift(-floor(mux/2):1:ceil(mux/2)-1])), Coil(=nc)]
% 
% % Correct fully sampled calibration data
% eddy = get_pha_flt( - pha_coe_cal, datsz.x);      % Dim: [X, Ky, Echo, Slice, SimultaneousSlice Z, Coil]
% cal_samp_ky = 1 : 1 : datsz.y;                    % All ky lines are sampled
% 
% cal_dat = dat(:, :, :, :, :, p.cal_dat_tpoints);
% cal_dat_nt = length(p.cal_dat_tpoints)/p.mux;
% cal_dat = reshape(cal_dat, [datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, p.mux, cal_dat_nt]); % Calibration data for recon. Dim: [Kx, Ky, Echo, Slice, Coil, Kz(=p.mux), Time]
% cal_dat = permute(cal_dat, [1,2,3,4,5,7,6]);      % Dim: [Kx, Ky, Echo, Slice, Coil, Time, Kz(=p.mux)]
% 
% [cal_dat, ranks] = mux_epi_md_ecc_full_kz(cal_dat, cal_samp_ky, eddy, p);
% 
% cal_dat = reshape(permute(cal_dat, [1,2,3,4,5,7,6]), [datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, p.mux*cal_dat_nt]);
% dat(:, :, :, :, :, p.cal_dat_tpoints) = cal_dat;


function [dat, pha_coe] = default_ecc(dat, p, ref, slices, nr_usegpu)
%
% Default(single-slice or slice-averaged) eddy current correction on both calibration and accelerated data.
%
% Inputs
%   dat     - K-space data to correct. Dim: [Kx(=nx), Ky(=ny), Echo(=nec), Slice(=nsl), Coil(=nc), Time].
%   p       - Parameter structure. See mux_epi_params.m for details.
%   ref     - Reference scan k-space data matrix(Dim: [Kx(=nx), Ky(=ny), Echo(=1), Slice(=nsl*mux), Coil(=nc)]),
%             or ref.dat file name('*.dat'), or reference scan pfile name('*refscan.7').
%   slices  - Slice indices(Normally order indices, not pfile indices) the
%             input 'dat' corresponds to. Only needed if the type of 'ref'
%             is ref.dat file name or reference scan pfile name.
%
% Outputs
%   dat     - Corrected k-space data. Size same as the input 'dat'.
%   pha_coe - EPI x-ky phase correction coefficients. Dim: [2(0th order, 1st order), ny, nsl, nc].
%
% (c) Kangrong Zhu,     Stanford University     Nov 2013

ref_type = get_ref_type(ref);

if strcmp(ref_type, 'ref.dat_file') || strcmp(ref_type, 'ksp_data')
    datsz = get_dat_sz(dat, p);
end
if (~exist('slices', 'var') || isempty(slices)) && (strcmp(ref_type, 'ref.dat_file') || strcmp(ref_type, 'ref_pfile'))
    error('No slice indices');
end

if strcmp(ref_type, 'ref.dat_file') % Load coefficients from ref.dat file
    if p.debug
        fprintf('   Using ref.dat file in default ghosting correction.\n');
    end
    pha_coe = rawload_ref(datsz.y, p.num_slices, p.num_coils, p.frames, slices, p.coils, ref);
else                                % Calculate coefficients from reference scan k-space data or from reference scan pfile
    pha_coe_1stpass = [];           % []: No first-pass phase correction
    ramp_flt = [];                  % []: No first-pass ramp sampling correction
    pccoil = 0;                     % 0: No averaging across coil
    do_quad_final = true;           % True: Conduct smoothing along ky on the least squares fitted coefficients
    debug = false;                  % True: Display correction coefficient maps
    
    switch ref_type
        case 'ref_pfile'
            if p.debug
                fprintf('   Using reference scan pfile in default ghosting correction.\n');
            end
            pha_coe = rawload_ref_pfile(p.frames, slices, p.coils, ref, pha_coe_1stpass, ramp_flt, pccoil, do_quad_final, debug);
        case 'ksp_data'
            if p.debug
                fprintf('   Using ecc data collected by mux sequence in default ghosting correction.\n');
            end
            ref = ref(:, :, :, datsz.sl*floor(p.mux/2)+1 : datsz.sl*(floor(p.mux/2)+1), :); % Reference scan data for the middle band. Dim: [Kx, Ky, Echo(=1), Slice(=nsl), Coil(=nc)].
% % %             ref = ref(:, :, :, 1 : datsz.sl*(floor(p.mux/2)), :); % Reference scan data for the 1st band. Dim: [Kx, Ky, Echo(=1), Slice(=nsl), Coil(=nc)].
            pha_coe = rawload_ref_data(ref, p, pha_coe_1stpass, ramp_flt, pccoil, do_quad_final, debug);
    end
end

% dat = epi_pha_correct(dat, pha_coe, p);
datsz = get_dat_sz(dat, p);
if nr_usegpu
    pha_flt = gpuArray(get_pha_flt(pha_coe, datsz.x));                              % Dim: [X, Ky, Slice, Coil]
    datfac = sqrt(gpuArray(size(dat{1}, p.FE_DIM)));
else
    pha_flt = get_pha_flt(pha_coe, datsz.x);                              % Dim: [X, Ky, Slice, Coil]
    datfac = sqrt(size(dat{1}, p.FE_DIM));
end
pha_flt = reshape(pha_flt, [datsz.x, datsz.y, 1, datsz.sl, datsz.c]); % Dim: [X, Ky, Echo(=1), Slice, Coil]
for dc = 1:numel(dat)
    dat{dc} = ifftshift(dat{dc}, p.FE_DIM);
    dat{dc} = ifft(dat{dc}, [], p.FE_DIM);
    dat{dc} = datfac .* fftshift(dat{dc}, p.FE_DIM);
    dat{dc} = bsxfun(@times, dat{dc}, pha_flt);
    dat{dc} = fftc(dat{dc}, p.FE_DIM, nr_usegpu);
end



function ref_type = get_ref_type(ref)
%
% Determines the the data type for eddy current correction.
%
% Input
%   ref      - Reference scan k-space data, or '*.dat' file name, or '*refscan.7' file name.
%
% Output
%   ref_type - What the input 'ref' contains.
%              'ksp_data'     - Reference scan k-space data.
%              'ref.dat_file' - Name of ref.dat file.
%              'ref_pfile'    - Name of reference scan pfile.
%
% (c) Kangrong Zhu,     Stanford University     Nov 2013

if isnumeric(ref)
    ref_type = 'ksp_data';
else
    if strcmp(ref(end-3:end), '.dat')
        ref_type = 'ref.dat_file';
    else
        if strcmp(ref(end-1:end), '.7')
            ref_type = 'ref_pfile';
        else
            error('ref must be kspace data or *.dat file or *.7 file.');
        end
    end
end



% function [ksp, pha_flt] = epi_pha_correct(ksp, pha_coe, p)
%
% % Correct the 0-th and 1st order phases in the x-ky space for EPI data.
% % Correspondingly, the constant phases of and the misalignment between the
% % odd and even echoes in the k-space are corrected.
% %
% % Inputs
% %   ksp     - Uncorrected k-space data. Dim: [Kx, Ky, Echo, Slice, Coil, Time]
% %   pha_coe - Coefficients for x-ky phase correction. See get_pha_flt.m for details.
% %             Dim: [2(0th order, 1st order), other dimensions(e.g. ny, nsl, nc].
% %   p       - Parameter structure. See mux_epi_params.m for details. The following
% %             fields are used in this function: FE_DIM, PE_DIM, EC_DIM, SL_DIM, C_DIM, T_DIM.
% %
% % Outputs
% %   ksp     - Corrected k-space data (note that the correction is done in-place).
% %   pha_flt - X-ky phase filter. xKy_corrected = xKy_uncorrected .* pha_flt.
% %             Dim: [X, Ky, Echo(=1), Slice, Coil].
% %
% % (c) Kangrong Zhu,     Stanford University     July 2012
% 
% datsz = get_dat_sz(ksp, p);
% if nr_usegpu
%     pha_flt = gpuArray(get_pha_flt(pha_coe, datsz.x));                              % Dim: [X, Ky, Slice, Coil]
% else
%     pha_flt = get_pha_flt(pha_coe, datsz.x);                              % Dim: [X, Ky, Slice, Coil]
% end
% pha_flt = reshape(pha_flt, [datsz.x, datsz.y, 1, datsz.sl, datsz.c]); % Dim: [X, Ky, Echo(=1), Slice, Coil]
% pha_flt = repmat(pha_flt, [1, 1, datsz.ec, 1, 1, datsz.t]);
% ksp = fftc(ifftc(ksp, p.FE_DIM) .* pha_flt, p.FE_DIM, nr_usegpu); % The same coefficients are used for every echo and every time point.



function [data, params] = load_raw_tseries(fname,params,slices,passes)
%
% Modified from rawloadX.m, to load k-space data of a time series conntained
% in a p-file.
%
% function [data, params] = load_raw_tseries(fname,params,slices,passes);
%
%	Function reads selected (or all) data from a P-file,
%	using given lists of slices and/or passes (i.e., timepoints).
%	Reads 12.x, 14.x and 20.x EPIC files.
%
%	INPUT:
%		fname   = path and file name.
%       params  = The structure with all parameters of the scan and
%                 for the reconstruction. Please refer to mux_epi_params.m
%                 for the definition of each field in this structure.
%                 The 'pfile_header' field in 'params' is used in
%                 this function. If 'params.pfile_header' is empty, the
%                 header information will be read using function 'rawheadX.m'.
%		slices  = slice numbers to read (1...).
%       passes  = pass numbers to read (1...). For the fMRI scan,
%                 this must be continuous from 1 to the last pass
%                 number to read. For an external calibration scan,
%                 this should correspond to the last group of slice
%                 phase cycling.

%
%	OUTPUT:
%		data = data array (up to 6 dimensions, including all
%			specified frames, echoes, slices, coils, passes.
%               params = The output parameter structure, having the same
%                       fields as the input 'params'. The following fields might
%                       have been changed inside this function if the scan was aborted
%                       midway and the p-file was only partially collected: num_mux_cycle,
%                       num_passes, nt_to_recon, tpoints_to_load.
%
%
%	NOTES:
%	1) If ranges are beyond those of the header, then they are
%		limited by ranges in the header.
%	2) If not all expected frames are read, then the data are returned
%		as a 2D array, with all frames read.
%
%
%	B.Hargreaves -- April 2008.
%       Modified by Kangrong Zhu,   Stanford University       July 2012

% -- Header information
if isempty(params.pfile_header)
    header = rawheadX(fname);
else
    header = params.pfile_header;
end

% -- Set defaults, if not provided.
if (~exist('slices','var') || isempty(slices))
    slices = 1:(header.nslices/header.npasses);
end;
if (~exist('passes','var') || isempty(passes))
    passes = 1:header.npasses;
end;

rds_data_order = params.rds_data_order;
inplane_R = params.inplane_R;

% -- Open file.
d = dir(fname);
file_bytes = d.bytes;
nframes = header.nframes+header.hnover;
nechos = header.nechoes;
ncoils  = header.ncoils;
num_acquired_slices = params.num_slices;

slices = rtcheckrange(slices,0,num_acquired_slices,'Slices');
passes = rtcheckrange(passes,0,header.npasses,'Passes');

if (passes(1)==0) && (passes(end) ~= length(passes)-1) % Assuming num_mux_cycle is always larger than 1, so the p-file must be from an fMRI scan, not from an external calibration scan, when passes(1)==1
    error('The pass numbers to read must be continuous from 1 to the last pass number to read.');
end

% --- Allocate array for data, based on passed arguments or defaults.
data = cell(length(passes), 1);

ptsize = header.ptsize;			% Sample size (bytes)
rawhdrsize = header.rawhdrsize;

aborted_scan = false;                   % true: the scan was aborted midway and the k-space data was collected only partially.

dattypes = {'int16','int32'};

if(rds_data_order)
    % It took a year, but I think I finally cracked the flip even/odd thing! Time will tell as the scan onslaught rages on...
    % The important rule seems to be: if the number of acquired echos is even, flip even frames. If it's odd, flip odd frames.
    % And note that the number of acquired echos is nframes/inplane_R.
    if(mod(nframes/inplane_R, 2))
        flip_even_frames = false;
    else
        flip_even_frames = true;
    end
    fip = fopen(fname,'r','l');
    if fip == -1
        error('File %s not found\n',fn);
    end
    % The following hack is no longer needed-- the rdsclient ensures that
    % the header size is set correctly.
    % if params.pfile_header.version<20.0065
    %     hdr_size = 149800;
    % else
    %     hdr_size = 149808;
    % end
    % if(rawhdrsize ~= hdr_size) %file reports 149788, ese22 was 149800, file size suggests 149808 for ese23
    %     fprintf('  Fixing raw header size.\n');
    %     rawhdrsize = hdr_size;
    % end
    
    % Data are arranged (from fastest to slowest moving):
    % coil - echo - slice - pass (timeframe)
    % note that in here a "frame" refers to a line of kspace (freq-encode line).
    framesize = header.frsize;
    viewsize = ncoils * nechos * framesize;
    viewbytes = viewsize * 2 * ptsize;
    % Size of one slice (bytes)
    slicebytes = viewbytes * (header.nframes+header.hnover)/inplane_R;
    passbytes = slicebytes * num_acquired_slices;      % Size of one pass (bytes)
    start_pass = 1;
end
% we might have decided that this isn't an rds file after counting bytes...
if(rds_data_order)
    passes = passes - 1;
    slices = slices - 1;
    for passindex = start_pass:length(passes)
        if(passes(passindex) < params.mux*(inplane_R-1)*params.num_mux_cycle)
            % These are calibration scans.
            acquired_pass = passes(passindex) + passes(passindex) * (inplane_R-1);
            for shot = 0:inplane_R-1
                for sliceindex = 1:length(slices)
                    fseek(fip, rawhdrsize + (acquired_pass + shot)*passbytes + slices(sliceindex)*slicebytes, -1);
                    [d,nread] = fread(fip, slicebytes/ptsize, dattypes{ptsize/2});
                    if (nread == slicebytes/ptsize)
                        d = d(1:2:end) + 1i*d(2:2:end);
%                         for frame = 0:inplane_R:nframes-1
%                             for echo = 0:nechos-1
%                                 for coil = 0:ncoils-1
%                                     startind = coil*framesize + echo*ncoils*framesize + frame/inplane_R*viewsize;
%                                     endind = startind + framesize;
%                                     even = mod(floor(frame/inplane_R),2)==0;
%                                     if((even && flip_even_frames) || (~even && ~flip_even_frames))
%                                         data(1:end, nframes-frame-shot, echo+1, sliceindex, coil+1, passindex) = d(endind:-1:startind+1);
%                                     else
%                                         data(1:end, nframes-frame-shot, echo+1, sliceindex, coil+1, passindex) = d(startind+1:endind);
%                                     end
%                                 end
%                             end
%                         end
                        ddata = reshape(d, [framesize, ncoils, nframes]);
                        ddata = permute(ddata, [1, 3, 2]);
                        ddata = reshape(ddata, [framesize, nframes, nechos, 1, ncoils]);
                        if flip_even_frames
                            ddata(:, 1:2:end, :, :, :) = ddata(end:-1:1, 1:2:end, :, :, :);
                        else
                            ddata(:, 2:2:end, :, :, :) = ddata(end:-1:1, 2:2:end, :, :, :);
                        end
                        data{passindex} = ddata(:, end:-1:1, :, :, :);
                    else
                        aborted_scan = true;
                        acquired_num_passes = passindex - 1;    % Assuming the data until the previous pass has been collected correctly.
                        break;                                  % Once we know the p-file was only collected partially, we don't need to read in any more data.
                    end;
                end;  % --Slice loop.
                if aborted_scan
                    break;
                end
            end % -- Shot loop
            if aborted_scan
                break;
            end
        else
            % Done with calibration scans-- load the normal scans here.
            acquired_pass = passes(passindex) + params.mux*(inplane_R-1)*params.num_mux_cycle;
            for sliceindex = 1:length(slices)
                % Seek to the correct slice.
                fseek(fip, rawhdrsize + acquired_pass*passbytes + slices(sliceindex)*slicebytes, -1);
                % Read the entire slice (all coils, all echos, all frames)
                [d,nread] = fread(fip, slicebytes/ptsize, dattypes{ptsize/2});
                if (nread == slicebytes/ptsize)
                    % convert to complex
                    d = d(1:2:end) + 1i*d(2:2:end);
                    % We only get the actual acquired data. We'll fill the acquired lines and leave the rest zero-padded.
%                     for frame = 0:inplane_R:nframes-1
%                         for echo = 0:nechos-1
%                             for coil = 0:ncoils-1
%                                 % EPI acquires data in a raster pattern, so we need to flip every other FE line.
%                                 startind = coil*framesize + echo*ncoils*framesize + frame/inplane_R*viewsize;
%                                 endind = startind + framesize;
%                                 even = mod(floor(frame/inplane_R),2)==0;
%                                 if((even && flip_even_frames) || (~even && ~flip_even_frames))
%                                     data(:, nframes-frame-(inplane_R-1), echo+1, sliceindex, coil+1, passindex) = d(endind:-1:startind+1);
%                                 else
%                                     data(:, nframes-frame-(inplane_R-1), echo+1, sliceindex, coil+1, passindex) = d(startind+1:endind);
%                                 end
%                             end
%                         end
%                     end
                    ddata = reshape(d, [framesize, ncoils, nframes]);
                    ddata = permute(ddata, [1, 3, 2]);
                    ddata = reshape(ddata, [framesize, nframes, nechos, 1, ncoils]);
                    if flip_even_frames
                        ddata(:, 1:2:end, :, :, :) = ddata(end:-1:1, 1:2:end, :, :, :);
                    else
                        ddata(:, 2:2:end, :, :, :) = ddata(end:-1:1, 2:2:end, :, :, :);
                    end
                    data{passindex} = ddata(:, end:-1:1, :, :, :);
                else
                    aborted_scan = true;
                    acquired_num_passes = passindex - 1;    % Assuming the data untill the previous pass has been collected correctly.
                    break;                                  % Once we know the p-file was only collected partially, we don't need to read in any more data.
                end;
            end;  % --Slice loop.
            if aborted_scan
                break;
            end
        end % --if(passes...
    end % --Pass loop
    % -- No need to skip remaining passes, as no more data to read!
    fclose(fip);	% -- Close file.
    if params.debug
        fprintf('  RDS-format data loaded.\n');
    end
else
    % TODO: rewrite this code here. It can be simplified. Also, we should load the data into a more compute-friendly
    % array order so that when we start the calculations we aren't thrashing all over RAM.
    [data, params] = rawloadX_tseries(fname,[],[],slices,[],passes,params);
end


% -- Deal with partially acquired p-files.
if aborted_scan
    if params.debug
        fprintf('  The scan was aborted midway and the p-file was only collected partially.\n');
    end

    if ~exist('params', 'var')
        error('The parameter structure doesn''t exist.');
    end

    if acquired_num_passes < 1          % p-file contains no data.
        error('The p-file doesn''t contain any data.');
    end

    if acquired_num_passes < params.mux % Data in p-file is less than one mux phase cycling.
        error(['The number of temporal frames contained in the p-file is %d', ...
            'mux is %d. Can''t conduct any reconstruction in this case.'], acquired_num_passes, params.mux);
    end

    if acquired_num_passes <= params.mux*params.num_mux_cycle; % Data in p-file are all mux phase cycling frames.
        params.num_mux_cycle = floor(acquired_num_passes/params.mux);
        params.num_passes = params.mux*params.num_mux_cycle;
        params.nt_to_recon = 0;
        params.tpoints_to_load = 1 : params.num_passes;
        data = data(params.tpoints_to_load);
        if params.debug
            fprintf(['  Number of frames contained in the p-file is %d, mux is %d, num_mux_cycle is %d.', ...
                'Only %d calibration frames will be reconstructed.\n'], ...
                acquired_num_passes, params.mux, params.num_mux_cycle, params.num_passes);
        end
    end

    if acquired_num_passes > params.mux*params.num_mux_cycle; % Some accelerated frames are contained in the p-file.
        params.num_passes = acquired_num_passes;
        params.nt_to_recon = acquired_num_passes - params.mux*params.num_mux_cycle;
        params.tpoints_to_load = 1 : min(params.num_passes, numel(data));
        data = data(params.tpoints_to_load);
        if params.debug
            fprintf('  Only the acquired %d frames will be loaded and reconstructed.\n', params.num_passes);
        end
    end
end



function arrout = rtcheckrange(arr,amin,amax,atype)
%%%% -- Check values in array are within range, and remove ones that
	%%%are not!
f = find(arr >= amin);
arrout = arr(f);
f = find(arrout <= amax);
arrout = arrout(f);

if (length(arrout) < length(arr))
  arr = arr(:);
  f1 = find(arr < amin);
  f2 = find(arr > amax);
  tt = sprintf('%d %s out of range: ',length(f1)+length(f2),atype); disp(tt);
  disp([arr(f1); arr(f2)].');
end;
%%%% --------



% Copyright (C) 2004 Josep Mones i Teixidor
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.
%
% -*- texinfo -*-
% @deftypefn {Function File} {@var{B} = } im2col (@var{A}, [@var{m},@var{n}], @var{block_type})
% @deftypefnx {Function File} {@var{B} = } im2col (@var{A}, [@var{m},@var{n}])
% @deftypefnx {Function File} {@var{B} = } im2col (@var{A}, 'indexed', ...)
% Rearranges image blocks into columns.
%
% @code{B=im2col(A, [m, n], blocktype)} rearranges blocks in @var{A}
% into columns in a way that's determined by @var{block_type}, which
% can take the following values:
%
% @table @code
% @item distinct
% Rearranges each distinct @var{m}-by-@var{n} block in image @var{A}
% into a column of @var{B}. Blocks are scanned from left to right and
% the up to bottom in @var{A}, and columns are added to @var{B} from
% left to right. If @var{A}'s size is not multiple @var{m}-by-@var{n}
% it is padded.
% @item sliding
% Rearranges any @var{m}-by-@var{n} sliding block of @var{A} in a
% column of @var{B}, without any padding, so only sliding blocks which
% can be built using a full @var{m}-by-@var{n} neighbourhood are taken.
% In consequence, @var{B} has @var{m}*@var{n} rows and
% (@var{mm}-@var{m}+1)*(@var{nn}-@var{n}+1) columns (where @var{mm} and
% @var{nn} are the size of @var{A}).
%
% This case is thought to be used applying operations on columns of
% @var{B} (for instance using sum(:)), so that result is a
% 1-by-(@var{mm}-@var{m}+1)*(@var{nn}-@var{n}+1) vector, that is what
% the complementary function @code{col2im} expects.
% @end table
%
% @code{B=im2col(A,[m,n])} takes @code{distinct} as a default value for
% @var{block_type}.
%
% @code{B=im2col(A,'indexed',...)} will treat @var{A} as an indexed
% image, so it will pad using 1 if @var{A} is double. All other cases
% (incluing indexed matrices with uint8 and uint16 types and
% non-indexed images) will use 0 as padding value.
%
% Any padding needed in 'distinct' processing will be added at right
% and bottom edges of the image.
%
% @seealso{col2im}
% @end deftypefn
%
% Author:  Josep Mones i Teixidor <jmones@puntbarra.com>

function B = im2col(A, varargin)

% check [m,n]
m=varargin{1}(1);
n=varargin{1}(2);

indices = (1:m)' * ones(1, n) + ones(m, 1) * (0:size(A,1):(size(A,1) * n - 1));
access = (0:(size(A, 1)-m))' * ones(1, size(A, 2) + 1 - n) + ones(size(A, 1) + 1 -m, 1) * (0:size(A, 1):(size(A, 1) * (size(A, 2) - n)));
indices = reshape(indices, numel(indices), 1) * ones(1, numel(access)) + ...
    ones(numel(indices), 1) * reshape(access, 1, numel(access));
B = reshape(A(reshape(indices, numel(indices), 1)), size(indices));



function B = im2colgpu(A, varargin)

% check [m,n]
m=varargin{1}(1);
n=varargin{1}(2);

indices = gpuArray((1:m)') * gpuArray.ones(1, n) + gpuArray.ones(m, 1) * gpuArray(0:size(A,1):(size(A,1) * n - 1));
access = gpuArray((0:(size(A, 1)-m))') * gpuArray.ones(1, size(A, 2) + 1 - n) + gpuArray.ones(size(A, 1) + 1 -m, 1) * gpuArray(0:size(A, 1):(size(A, 1) * (size(A, 2) - n)));
indices = reshape(indices, numel(indices), 1) * gpuArray.ones(1, numel(access)) + ...
    gpuArray.ones(numel(indices), 1) * reshape(access, 1, numel(access));
B = reshape(A(reshape(indices, numel(indices), 1)), size(indices));



function dat = epi_vrgf_correct(dat, ramp_flt, nr_usegpu)
% function dat_out = epi_vrgf_correct(dat, ramp_flt)
%
% Correct ramp sampling in EPI data. Correction will be conducted image by image.
%
% Inputs
%   dat      - K-space data with ramp sampling. Dim: [Kx(=nxi), Ky, ...]
%   ramp_flt - Correction filter. When not exist or empty, no correction will be conducted. Dim: [nxo, nxi]
%
% Output
%   dat_out  - Corrected k-space data. Dim: [Kx(=nxo), Ky, ...]
%
% (c) Kangrong Zhu,     Stanford University     Sep 2013

if ~exist('ramp_flt', 'var') || isempty(ramp_flt)
    return;
end
if nr_usegpu
    for dc = 1:numel(dat)
        dat{dc} = pagefun(@mtimes, ramp_flt, dat{dc});
    end
else
    datsz = size(dat{1});
    cdat = zeros([size(ramp_flt, 1), datsz(2:end)]);
    for dc = 1:numel(dat)
        for pc = 1:prod(datsz(3:end))
            cdat(:, :, pc) = ramp_flt * dat{dc}(:, :, pc);
        end
        dat{dc} = cdat;
    end
end



function header = read_MR_headers(filename)

hdr_szs = [24, 157276, 4096, 16384, 16384, 98304, ...
    2052, 2052, 2048, 1500, 1960, 2560, 2448, 7488];
rdbm_rev = hdr_szs(1);
formatID = ['pfile', num2str(rdbm_rev)];
endianID = 'ieee-le';
header_lengths = hdr_szs(3:numel(hdr_szs))';
header_list = 1:length(header_lengths);
header = struct;
header.endian = endianID;
header.format = formatID;

fid = fopen( filename, 'r', header.endian );
offset = 0;

% Step through all the headers present in the file, reading them in order
for i = 1:length(header_list)
    fseek( fid, offset, 'bof');
    % fprintf('  Header: %3d  Offset: %6d\n', header_list(i), offset);
    switch header_list(i)
        case 1 % rdb_hdr
            header.rdb_hdr = read_rdb_hdr( fid, rdbm_rev );
        case  2 %  2 per_pass (not implemented)
        case  3 %  3 unlock_raw (not implemented)
        case  4 %  4 data_acq_tab
            %                    header.data_acq_tab = read_data_acq_tab( fid, header.rdb_hdr.nslices);
            tmp_struct.remove_me = 1;
            for slicenum = 1:header.rdb_hdr.nslices/header.rdb_hdr.npasses % RFD: was "header.rdb_hdr.nslices", but we don't need to read the slices for all time points
                my_struct = read_data_acq_tab(fid, rdbm_rev);
                fnames = fieldnames(my_struct);
                %Assume 2D array at most
                for fdnum = 1:length(fnames)
                    xnum = size(my_struct.(fnames{fdnum}), 1);
                    ynum = size(my_struct.(fnames{fdnum}), 2);
                    for xvar = 1:xnum
                        for yvar = 1:ynum
                            tmp_struct = setfield(tmp_struct,fnames{fdnum},{slicenum,xvar,yvar},getfield(my_struct,fnames{fdnum},{xvar,yvar}));
                        end
                    end
                end
            end
            tmp_struct = rmfield(tmp_struct,'remove_me');
            header.data_acq_tab = tmp_struct;
            clear xnum ynum xvar yvar slicenum tmp_struct my_struct

        case  5 %  5 nex_tab (not implemented)
        case  6 %  6 nex_abort_tab (not implemented)
        case  7 %  7 tool (not implemented)
        case  8 %  8 prescan 
            header.psc  = read_psc_header(fid, rdbm_rev); 
        case  9 %  9 exam
            header.exam  = read_exam_header( fid, rdbm_rev );
        case 10 % 10 series
            header.series  = read_series_header( fid, rdbm_rev );
        case 11 % 11 image
            header.image  = read_image_header( fid, rdbm_rev );
        %case 12 % 12 grad_data
         %   header.grad_data = read_grad_header( fid, rdbm_rev );
    end
    offset = offset + header_lengths(header_list(i));
end
header.total_length = offset;

fclose(fid);



function im = sos(ic, dim)
%
% im = sos(ic [, dim])
%
% Square root of sum of squares(SOS) reconstruction of multicoil images.
% 
% Inputs:
%   ic  - Single coil images. 
%   dim - The dimension for coil in matrix 'ic'. Default: the last dimension of 'ic'.
%
% Output:
%   im  - The SOS combined images.
%
% (c) Kangrong Zhu,     Stanford University     2011

if ~exist('dim','var')
    dim = length( size(ic) );
end

im = sqrt( sum( ic.*conj(ic), dim) );



function a = read_series_header( my_file, rdbm_rev )
%read_series_header - Read GE series header
%
%  a = read_series_header( my_file, rdbm_rev );
%    my_file - string indicating file name to read
%    rdbm_rev - raw header (RDBM) revision number
%    a - structure with header values
%

% Copyright (c) 2012 by General Electric Company. All rights reserved.


% RDBM revision 24.000
if rdbm_rev == 24.000 
  for id = 1 : 32
    a.double_padding(id) = fread(my_file, 1, 'float64');
  end
  a.se_pds_a = fread(my_file, 1, 'float32');
  a.se_pds_c = fread(my_file, 1, 'float32');
  a.se_pds_u = fread(my_file, 1, 'float32');
  a.lmhor = fread(my_file, 1, 'float32');
  a.start_loc = fread(my_file, 1, 'float32');
  a.end_loc = fread(my_file, 1, 'float32');
  a.echo1_alpha = fread(my_file, 1, 'float32');
  a.echo1_beta = fread(my_file, 1, 'float32');
  a.echo2_alpha = fread(my_file, 1, 'float32');
  a.echo2_beta = fread(my_file, 1, 'float32');
  a.echo3_alpha = fread(my_file, 1, 'float32');
  a.echo3_beta = fread(my_file, 1, 'float32');
  a.echo4_alpha = fread(my_file, 1, 'float32');
  a.echo4_beta = fread(my_file, 1, 'float32');
  a.echo5_alpha = fread(my_file, 1, 'float32');
  a.echo5_beta = fread(my_file, 1, 'float32');
  a.echo6_alpha = fread(my_file, 1, 'float32');
  a.echo6_beta = fread(my_file, 1, 'float32');
  a.echo7_alpha = fread(my_file, 1, 'float32');
  a.echo7_beta = fread(my_file, 1, 'float32');
  a.echo8_alpha = fread(my_file, 1, 'float32');
  a.echo8_beta = fread(my_file, 1, 'float32');
  a.landmark = fread(my_file, 1, 'float32');
  a.tablePosition = fread(my_file, 1, 'float32');
  a.pure_lambda = fread(my_file, 1, 'float32');
  a.pure_tuning_factor_surface = fread(my_file, 1, 'float32');
  a.pure_tuning_factor_body = fread(my_file, 1, 'float32');
  a.pure_derived_cal_fraction = fread(my_file, 1, 'float32');
  a.pure_derived_cal_reapodization = fread(my_file, 1, 'float32');
  for id = 1 : 25
    a.float_padding(id) = fread(my_file, 1, 'float32');
  end
  a.se_complete = fread(my_file, 1, 'int32');
  a.se_numarch = fread(my_file, 1, 'int32');
  a.se_imagect = fread(my_file, 1, 'int32');
  a.se_numimages = fread(my_file, 1, 'int32');
  a.se_delta_cnt = fread(my_file, 1, 'int32');
  a.se_numunimg = fread(my_file, 1, 'int32');
  a.se_toarchcnt = fread(my_file, 1, 'int32');
  for id = 1 : 33
    a.int_padding1(id) = fread(my_file, 1, 'int32');
  end
  a.se_datetime = fread(my_file, 1, 'int32');
  a.se_actual_dt = fread(my_file, 1, 'int32');
  a.position = fread(my_file, 1, 'int32');
  a.entry = fread(my_file, 1, 'int32');
  a.se_lndmrkcnt = fread(my_file, 1, 'int32');
  a.se_lastmod = fread(my_file, 1, 'int32');
  a.ExpType = fread(my_file, 1, 'int32');
  a.TrRest = fread(my_file, 1, 'int32');
  a.TrActive = fread(my_file, 1, 'int32');
  a.DumAcq = fread(my_file, 1, 'int32');
  a.ExptTimePts = fread(my_file, 1, 'int32');
  a.cal_pass_set_vector = fread(my_file, 1, 'int32');
  a.cal_nex_vector = fread(my_file, 1, 'int32');
  a.cal_weight_vector = fread(my_file, 1, 'int32');
  a.pure_filtering_mode = fread(my_file, 1, 'int32');
  for id = 1 : 29
    a.int_padding2(id) = fread(my_file, 1, 'int32');
  end
  a.se_exno = fread(my_file, 1, 'uint16');
  a.echo1_window = fread(my_file, 1, 'uint16');
  a.echo2_window = fread(my_file, 1, 'uint16');
  a.echo3_window = fread(my_file, 1, 'uint16');
  a.echo4_window = fread(my_file, 1, 'uint16');
  a.echo5_window = fread(my_file, 1, 'uint16');
  a.echo6_window = fread(my_file, 1, 'uint16');
  a.echo7_window = fread(my_file, 1, 'uint16');
  a.echo8_window = fread(my_file, 1, 'uint16');
  a.echo8_level = fread(my_file, 1, 'int16');
  a.echo7_level = fread(my_file, 1, 'int16');
  a.echo6_level = fread(my_file, 1, 'int16');
  a.echo5_level = fread(my_file, 1, 'int16');
  a.echo4_level = fread(my_file, 1, 'int16');
  a.echo3_level = fread(my_file, 1, 'int16');
  a.echo2_level = fread(my_file, 1, 'int16');
  a.echo1_level = fread(my_file, 1, 'int16');
  a.se_no = fread(my_file, 1, 'int16');
  a.se_typ = fread(my_file, 1, 'int16');
  a.se_source = fread(my_file, 1, 'int16');
  a.se_plane = fread(my_file, 1, 'int16');
  a.scan_type = fread(my_file, 1, 'int16');
  a.se_uniq = fread(my_file, 1, 'int16');
  a.se_contrast = fread(my_file, 1, 'int16');
  a.se_pseq = fread(my_file, 1, 'int16');
  a.se_sortorder = fread(my_file, 1, 'int16');
  a.se_nacq = fread(my_file, 1, 'int16');
  a.xbasest = fread(my_file, 1, 'int16');
  a.xbaseend = fread(my_file, 1, 'int16');
  a.xenhst = fread(my_file, 1, 'int16');
  a.xenhend = fread(my_file, 1, 'int16');
  a.table_entry = fread(my_file, 1, 'int16');
  a.SwingAngle = fread(my_file, 1, 'int16');
  a.LateralOffset = fread(my_file, 1, 'int16');
  a.GradientCoil = fread(my_file, 1, 'int16');
  a.se_subtype = fread(my_file, 1, 'int16');
  a.BWRT = fread(my_file, 1, 'int16');
  a.assetcal_serno = fread(my_file, 1, 'int16');
  a.assetcal_scnno = fread(my_file, 1, 'int16');
  a.content_qualifn = fread(my_file, 1, 'int16');
  a.purecal_serno = fread(my_file, 1, 'int16');
  a.purecal_scnno = fread(my_file, 1, 'int16');
  a.ideal = fread(my_file, 1, 'int16');
  a.verify_corners = fread(my_file, 1, 'int16');
  a.asset_cal_type = fread(my_file, 1, 'int16');
  a.pure_compatible = fread(my_file, 1, 'int16');
  a.purecal_type = fread(my_file, 1, 'int16');
  for id = 1 : 29
    a.short_padding(id) = fread(my_file, 1, 'int16');
  end
  for id = 1 : 2
    a.se_verscre(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 2
    a.se_verscur(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 4
    a.se_suid(id) = fread(my_file, 1, 'char');
  end
  a.se_alloc_key = freadc(my_file, 13);
  a.se_diskid = fread(my_file, 1, 'char');
  a.se_desc = freadc(my_file, 65);
  a.pr_sysid = freadc(my_file, 9);
  a.pansysid = freadc(my_file, 9);
  a.anref = freadc(my_file, 3);
  a.prtcl = freadc(my_file, 25);
  a.start_ras = fread(my_file, 1, 'char');
  a.end_ras = fread(my_file, 1, 'char');
  for id = 1 : 32
    a.series_uid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.landmark_uid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.equipmnt_uid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.refsopcuids(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.refsopiuids(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.schacitval(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.schacitdesc(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 64
    a.schacitmea(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 65
    a.schprocstdesc(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.schprocstid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.reqprocstid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.perprocstid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 65
    a.perprocstdesc(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.reqprocstid2(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.reqprocstid3(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.schprocstid2(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.schprocstid3(id) = fread(my_file, 1, 'char');
  end
  for id1 = 1 : 4
    for id2 = 1 : 32
        a.refImgUID(id1,id2) = fread(my_file, 1, 'char');
    end
  end
  for id = 1 : 64
    a.PdgmStr(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 256
    a.PdgmDesc(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 64
    a.PdgmUID(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.ApplName(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.ApplVer(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 12
    a.asset_appl(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.scic_a(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.scic_s(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.scic_c(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 64
    a.pure_cfg_params(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 251
    a.se_padding(id) = fread(my_file, 1, 'char');
  end

end



function dat = mux_dftz(dat, dim, fov_shift, nz, enc_or_dec)
%
% function dat = mux_dftz(dat, [dim=ndims(dat)], [fov_shift=-size(dat,dim)], [nz=size(dat,dim)], [enc_or_dec='decode'])
%
% Encode or decode the DFTz encoding on simultaneous slices.
%
% Inputs
%   dat        - Original data. Slice indices are ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))).
%   dim        - The dimension for the simultaneous slices, along which the encoding or decoding will be conducted.
%   fov_shift  - CAIPI FOV shift. abs(fov_shift) is the length of the DFTz encoding. Negative for inverse DFTz encoding, positive for forward DFTz encoding.
%   nz         - Number of simultaneous slices.
%   enc_or_dec - 'encode' or 'e': Encode the input 'dat'.
%                'decode' or 'd': Decode the input 'dat'.
%
% Output
%   dat        - Input 'dat' with forward or inverse DFTz encoded or decoded along 'dim'. Slice indices are ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))).
%
% (c) Kangrong Zhu,     Stanford University     June 2013

if ~exist('dim', 'var') || isempty(dim)
    dim = ndims(dat);
end

if ~exist('fov_shift', 'var') || isempty(fov_shift)
    fov_shift = - size(dat, dim);
end

if ~exist('nz', 'var') || isempty(nz)
    nz = size(dat, dim);
end

if ~exist('enc_or_dec', 'var') || isempty(enc_or_dec)
    enc_or_dec = 'decode';
end

sz_orig = size(dat);
dims_all = 1 : length(sz_orig);
dims_non_slice = dims_all(dims_all ~= dim);                       % Dimensions except the slice dimension
dims_to_permute = [dim, dims_non_slice];                          % Matrix permutation order to make the slice dimension to be the 1st dimension
dims_to_permute_back = [2:dim, 1, dim+1:length(dims_to_permute)]; % Matrix permutation order to permute the data back to its orginal order

% Prepare data
dat = permute(dat, dims_to_permute);                              % Make the slice dimension to be the 1st dimension
sz_permuted = size(dat);
dat = reshape(dat, [sz_permuted(1), prod(sz_permuted(2:end))]);

% Calculate encoding matrix
encode_mtx = encode_dftz_mtx(fov_shift, nz);                         % Dim: [abs(fov_shift)(i.e. nkz), nz]

% Encode or decode
switch lower(enc_or_dec(1))
    case 'e'                                                      % Encode
        dat = encode_mtx * dat;
    case 'd'                                                      % Decode
        dat = encode_mtx \ dat;                                   % i.e. inv(encode_mtx)*dat
end
dat = reshape(dat, [size(dat,1), sz_permuted(2:end)]);
dat = permute(dat, dims_to_permute_back);



function a = read_image_header( my_file,rdbm_rev )
%read_image_header - Read GE image header
%
%  a = read_image_header( my_file, rdbm_rev );
%    my_file - string indicating file name to read
%    rdbm_rev - raw header (RDBM) revision number
%    a - structure with header values
%

% Copyright (c) 2012 by General Electric Company. All rights reserved.


% RDBM revision 24.000
if rdbm_rev == 24.000 
  a.autoSubParam.seriesUidToSubtract = fread(my_file, 32, 'char')';
  a.autoSubParam.imageNoToSubtract = fread(my_file, 1, 'int32');
  a.autoSubParam.destSeriesNo = fread(my_file, 1, 'int32');
  a.autoSubParam.destImageNo = fread(my_file, 1, 'int32');
  a.autoSubParam.dummy = fread(my_file, 1, 'int32');
  for id = 1 : 32
    a.double_padding(id) = fread(my_file, 1, 'float64');
  end
  a.dfov = fread(my_file, 1, 'float32');
  a.dfov_rect = fread(my_file, 1, 'float32');
  a.sctime = fread(my_file, 1, 'float32');
  a.slthick = fread(my_file, 1, 'float32');
  a.scanspacing = fread(my_file, 1, 'float32');
  a.loc = fread(my_file, 1, 'float32');
  a.tbldlta = fread(my_file, 1, 'float32');
  a.nex = fread(my_file, 1, 'float32');
  a.reptime = fread(my_file, 1, 'float32');
  a.saravg = fread(my_file, 1, 'float32');
  a.sarpeak = fread(my_file, 1, 'float32');
  a.pausetime = fread(my_file, 1, 'float32');
  a.vbw = fread(my_file, 1, 'float32');
  a.user0 = fread(my_file, 1, 'float32');
  a.user1 = fread(my_file, 1, 'float32');
  a.user2 = fread(my_file, 1, 'float32');
  a.user3 = fread(my_file, 1, 'float32');
  a.user4 = fread(my_file, 1, 'float32');
  a.user5 = fread(my_file, 1, 'float32');
  a.user6 = fread(my_file, 1, 'float32');
  a.user7 = fread(my_file, 1, 'float32');
  a.user8 = fread(my_file, 1, 'float32');
  a.user9 = fread(my_file, 1, 'float32');
  a.user10 = fread(my_file, 1, 'float32');
  a.user11 = fread(my_file, 1, 'float32');
  a.user12 = fread(my_file, 1, 'float32');
  a.user13 = fread(my_file, 1, 'float32');
  a.user14 = fread(my_file, 1, 'float32');
  a.user15 = fread(my_file, 1, 'float32');
  a.user16 = fread(my_file, 1, 'float32');
  a.user17 = fread(my_file, 1, 'float32');
  a.user18 = fread(my_file, 1, 'float32');
  a.user19 = fread(my_file, 1, 'float32');
  a.user20 = fread(my_file, 1, 'float32');
  a.user21 = fread(my_file, 1, 'float32');
  a.user22 = fread(my_file, 1, 'float32');
  a.proj_ang = fread(my_file, 1, 'float32');
  a.concat_sat = fread(my_file, 1, 'float32');
  a.user23 = fread(my_file, 1, 'float32');
  a.user24 = fread(my_file, 1, 'float32');
  a.x_axis_rot = fread(my_file, 1, 'float32');
  a.y_axis_rot = fread(my_file, 1, 'float32');
  a.z_axis_rot = fread(my_file, 1, 'float32');
  a.ihtagfa = fread(my_file, 1, 'float32');
  a.ihtagor = fread(my_file, 1, 'float32');
  a.ihbspti = fread(my_file, 1, 'float32');
  a.rtia_timer = fread(my_file, 1, 'float32');
  a.fps = fread(my_file, 1, 'float32');
  a.vencscale = fread(my_file, 1, 'float32');
  a.dbdt = fread(my_file, 1, 'float32');
  a.dbdtper = fread(my_file, 1, 'float32');
  a.estdbdtper = fread(my_file, 1, 'float32');
  a.estdbdtts = fread(my_file, 1, 'float32');
  a.saravghead = fread(my_file, 1, 'float32');
  a.neg_scanspacing = fread(my_file, 1, 'float32');
  a.user25 = fread(my_file, 1, 'float32');
  a.user26 = fread(my_file, 1, 'float32');
  a.user27 = fread(my_file, 1, 'float32');
  a.user28 = fread(my_file, 1, 'float32');
  a.user29 = fread(my_file, 1, 'float32');
  a.user30 = fread(my_file, 1, 'float32');
  a.user31 = fread(my_file, 1, 'float32');
  a.user32 = fread(my_file, 1, 'float32');
  a.user33 = fread(my_file, 1, 'float32');
  a.user34 = fread(my_file, 1, 'float32');
  a.user35 = fread(my_file, 1, 'float32');
  a.user36 = fread(my_file, 1, 'float32');
  a.user37 = fread(my_file, 1, 'float32');
  a.user38 = fread(my_file, 1, 'float32');
  a.user39 = fread(my_file, 1, 'float32');
  a.user40 = fread(my_file, 1, 'float32');
  a.user41 = fread(my_file, 1, 'float32');
  a.user42 = fread(my_file, 1, 'float32');
  a.user43 = fread(my_file, 1, 'float32');
  a.user44 = fread(my_file, 1, 'float32');
  a.user45 = fread(my_file, 1, 'float32');
  a.user46 = fread(my_file, 1, 'float32');
  a.user47 = fread(my_file, 1, 'float32');
  a.user48 = fread(my_file, 1, 'float32');
  a.RegressorVal = fread(my_file, 1, 'float32');
  a.SliceAsset = fread(my_file, 1, 'float32');
  a.PhaseAsset = fread(my_file, 1, 'float32');
  for id = 1 : 4
    a.sarValues(id) = fread(my_file, 1, 'float32');
  end
  for id = 1 : 2
    a.shim_fov(id) = fread(my_file, 1, 'float32');
  end
  for id = 1 : 2
    a.shim_ctr_R(id) = fread(my_file, 1, 'float32');
  end
  for id = 1 : 2
    a.shim_ctr_A(id) = fread(my_file, 1, 'float32');
  end
  for id = 1 : 2
    a.shim_ctr_S(id) = fread(my_file, 1, 'float32');
  end
  a.dim_X = fread(my_file, 1, 'float32');
  a.dim_Y = fread(my_file, 1, 'float32');
  a.pixsize_X = fread(my_file, 1, 'float32');
  a.pixsize_Y = fread(my_file, 1, 'float32');
  a.ctr_R = fread(my_file, 1, 'float32');
  a.ctr_A = fread(my_file, 1, 'float32');
  a.ctr_S = fread(my_file, 1, 'float32');
  a.norm_R = fread(my_file, 1, 'float32');
  a.norm_A = fread(my_file, 1, 'float32');
  a.norm_S = fread(my_file, 1, 'float32');
  a.tlhc_R = fread(my_file, 1, 'float32');
  a.tlhc_A = fread(my_file, 1, 'float32');
  a.tlhc_S = fread(my_file, 1, 'float32');
  a.trhc_R = fread(my_file, 1, 'float32');
  a.trhc_A = fread(my_file, 1, 'float32');
  a.trhc_S = fread(my_file, 1, 'float32');
  a.brhc_R = fread(my_file, 1, 'float32');
  a.brhc_A = fread(my_file, 1, 'float32');
  a.brhc_S = fread(my_file, 1, 'float32');
  a.menc = fread(my_file, 1, 'float32');
  a.normal_L = fread(my_file, 1, 'float32');
  a.normal_P = fread(my_file, 1, 'float32');
  a.normal_S = fread(my_file, 1, 'float32');
  a.osf = fread(my_file, 1, 'float32');
  a.fermi_radius = fread(my_file, 1, 'float32');
  a.fermi_width = fread(my_file, 1, 'float32');
  a.fermi_ecc = fread(my_file, 1, 'float32');
  for id = 1 : 25
    a.float_padding(id) = fread(my_file, 1, 'float32');
  end
  a.cal_fldstr = fread(my_file, 1, 'uint32');
  a.user_usage_tag = fread(my_file, 1, 'uint32');
  a.user_fill_mapMSW = fread(my_file, 1, 'uint32');
  a.user_fill_mapLSW = fread(my_file, 1, 'uint32');
  a.im_archived = fread(my_file, 1, 'int32');
  a.im_complete = fread(my_file, 1, 'int32');
  for id = 1 : 34
    a.int_padding1(id) = fread(my_file, 1, 'int32');
  end
  a.im_datetime = fread(my_file, 1, 'int32');
  a.im_actual_dt = fread(my_file, 1, 'int32');
  a.tr = fread(my_file, 1, 'int32');
  a.ti = fread(my_file, 1, 'int32');
  a.te = fread(my_file, 1, 'int32');
  a.te2 = fread(my_file, 1, 'int32');
  a.tdel = fread(my_file, 1, 'int32');
  a.mindat = fread(my_file, 1, 'int32');
  a.obplane = fread(my_file, 1, 'int32');
  a.slocfov = fread(my_file, 1, 'int32');
  a.obsolete1 = fread(my_file, 1, 'int32');
  a.obsolete2 = fread(my_file, 1, 'int32');
  a.user_bitmap = fread(my_file, 1, 'int32');
  a.iopt = fread(my_file, 1, 'int32');
  a.psd_datetime = fread(my_file, 1, 'int32');
  a.rawrunnum = fread(my_file, 1, 'int32');
  a.intr_del = fread(my_file, 1, 'int32');
  a.im_lastmod = fread(my_file, 1, 'int32');
  a.im_pds_a = fread(my_file, 1, 'int32');
  a.im_pds_c = fread(my_file, 1, 'int32');
  a.im_pds_u = fread(my_file, 1, 'int32');
  a.thresh_min1 = fread(my_file, 1, 'int32');
  a.thresh_max1 = fread(my_file, 1, 'int32');
  a.thresh_min2 = fread(my_file, 1, 'int32');
  a.thresh_max2 = fread(my_file, 1, 'int32');
  a.numslabs = fread(my_file, 1, 'int32');
  a.locsperslab = fread(my_file, 1, 'int32');
  a.overlaps = fread(my_file, 1, 'int32');
  a.slop_int_4 = fread(my_file, 1, 'int32');
  a.dfax = fread(my_file, 1, 'int32');
  a.fphase = fread(my_file, 1, 'int32');
  a.offsetfreq = fread(my_file, 1, 'int32');
  a.b_value = fread(my_file, 1, 'int32');
  a.iopt2 = fread(my_file, 1, 'int32');
  a.ihtagging = fread(my_file, 1, 'int32');
  a.ihtagspc = fread(my_file, 1, 'int32');
  a.ihfcineim = fread(my_file, 1, 'int32');
  a.ihfcinent = fread(my_file, 1, 'int32');
  a.num_seg = fread(my_file, 1, 'int32');
  a.oprtarr = fread(my_file, 1, 'int32');
  a.averages = fread(my_file, 1, 'int32');
  a.station_index = fread(my_file, 1, 'int32');
  a.station_total = fread(my_file, 1, 'int32');
  a.iopt3 = fread(my_file, 1, 'int32');
  a.delAcq = fread(my_file, 1, 'int32');
  a.rxmbloblen = fread(my_file, 1, 'int32');
  a.rxmblob_pad = fread(my_file, 1, 'int32');
  a.im_no = fread(my_file, 1, 'int32');
  a.imgrx = fread(my_file, 1, 'int32');
  a.temp_phases = fread(my_file, 1, 'int32');
  a.driver_freq = fread(my_file, 1, 'int32');
  a.driver_amp = fread(my_file, 1, 'int32');
  a.driverCyc_Trig = fread(my_file, 1, 'int32');
  a.MEG_dir = fread(my_file, 1, 'int32');
  a.rescan_time = fread(my_file, 1, 'int32');
  a.spokesPerSeg = fread(my_file, 1, 'int32');
  a.recoveryTime = fread(my_file, 1, 'int32');
  a.t2PrepTE = fread(my_file, 1, 'int32');
  a.hoecc = fread(my_file, 1, 'int32');
  a.user_bitmap2 = fread(my_file, 1, 'int32');
  for id = 1 : 20
    a.int_padding2(id) = fread(my_file, 1, 'int32');
  end
  a.imatrix_X = fread(my_file, 1, 'int16');
  a.imatrix_Y = fread(my_file, 1, 'int16');
  a.im_exno = fread(my_file, 1, 'uint16');
  a.img_window = fread(my_file, 1, 'uint16');
  a.img_level = fread(my_file, 1, 'int16');
  a.numecho = fread(my_file, 1, 'int16');
  a.echonum = fread(my_file, 1, 'int16');
  a.im_uniq = fread(my_file, 1, 'int16');
  a.im_seno = fread(my_file, 1, 'int16');
  a.contmode = fread(my_file, 1, 'int16');
  a.serrx = fread(my_file, 1, 'int16');
  a.screenformat = fread(my_file, 1, 'int16');
  a.plane = fread(my_file, 1, 'int16');
  a.im_compress = fread(my_file, 1, 'int16');
  a.im_scouttype = fread(my_file, 1, 'int16');
  a.contig = fread(my_file, 1, 'int16');
  a.hrtrate = fread(my_file, 1, 'int16');
  a.trgwindow = fread(my_file, 1, 'int16');
  a.imgpcyc = fread(my_file, 1, 'int16');
  a.obsolete3 = fread(my_file, 1, 'int16');
  a.obsolete4 = fread(my_file, 1, 'int16');
  a.obsolete5 = fread(my_file, 1, 'int16');
  a.mr_flip = fread(my_file, 1, 'int16');
  a.cphase = fread(my_file, 1, 'int16');
  a.swappf = fread(my_file, 1, 'int16');
  a.pauseint = fread(my_file, 1, 'int16');
  a.obsolete6 = fread(my_file, 1, 'int16');
  a.obsolete7 = fread(my_file, 1, 'int16');
  a.obsolete8 = fread(my_file, 1, 'int16');
  a.not_used_1 = fread(my_file, 1, 'int16');
  a.imode = fread(my_file, 1, 'int16');
  a.pseq = fread(my_file, 1, 'int16');
  a.pseqmode = fread(my_file, 1, 'int16');
  a.ctyp = fread(my_file, 1, 'int16');
  a.surfctyp = fread(my_file, 1, 'int16');
  a.surfcext = fread(my_file, 1, 'int16');
  a.supp_tech = fread(my_file, 1, 'int16');
  a.slquant = fread(my_file, 1, 'int16');
  a.gpre = fread(my_file, 1, 'int16');
  a.satbits = fread(my_file, 1, 'int16');
  a.scic = fread(my_file, 1, 'int16');
  a.satxloc1 = fread(my_file, 1, 'int16');
  a.satxloc2 = fread(my_file, 1, 'int16');
  a.satyloc1 = fread(my_file, 1, 'int16');
  a.satyloc2 = fread(my_file, 1, 'int16');
  a.satzloc1 = fread(my_file, 1, 'int16');
  a.satzloc2 = fread(my_file, 1, 'int16');
  a.satxthick = fread(my_file, 1, 'int16');
  a.satythick = fread(my_file, 1, 'int16');
  a.satzthick = fread(my_file, 1, 'int16');
  a.flax = fread(my_file, 1, 'int16');
  a.venc = fread(my_file, 1, 'int16');
  a.thk_disclmr = fread(my_file, 1, 'int16');
  a.obsolete9 = fread(my_file, 1, 'int16');
  a.obsolete10 = fread(my_file, 1, 'int16');
  a.image_type = fread(my_file, 1, 'int16');
  a.vas_collapse = fread(my_file, 1, 'int16');
  a.proj_alg = fread(my_file, 1, 'int16');
  a.echo_trn_len = fread(my_file, 1, 'int16');
  a.frac_echo = fread(my_file, 1, 'int16');
  a.prep_pulse = fread(my_file, 1, 'int16');
  a.cphasenum = fread(my_file, 1, 'int16');
  a.var_echo = fread(my_file, 1, 'int16');
  a.scanactno = fread(my_file, 1, 'int16');
  a.vasflags = fread(my_file, 1, 'int16');
  a.integrity = fread(my_file, 1, 'int16');
  a.freq_dir = fread(my_file, 1, 'int16');
  a.vas_mode = fread(my_file, 1, 'int16');
  a.pscopts = fread(my_file, 1, 'int16');
  a.obsolete11 = fread(my_file, 1, 'int16');
  a.obsolete12 = fread(my_file, 1, 'int16');
  a.obsolete13 = fread(my_file, 1, 'int16');
  a.unoriginal = fread(my_file, 1, 'int16');
  a.interleaves = fread(my_file, 1, 'int16');
  a.effechospace = fread(my_file, 1, 'int16');
  a.viewsperseg = fread(my_file, 1, 'int16');
  a.rbpm = fread(my_file, 1, 'int16');
  a.rtpoint = fread(my_file, 1, 'int16');
  a.rcvrtype = fread(my_file, 1, 'int16');
  a.sarMode = fread(my_file, 1, 'int16');
  a.dBdtMode = fread(my_file, 1, 'int16');
  a.govBody = fread(my_file, 1, 'int16');
  a.sarDefinition = fread(my_file, 1, 'int16');
  a.no_shimvol = fread(my_file, 1, 'int16');
  a.shim_vol_type = fread(my_file, 1, 'int16');
  a.current_phase = fread(my_file, 1, 'int16');
  a.art_level = fread(my_file, 1, 'int16');
  a.slice_group_number = fread(my_file, 1, 'int16');
  a.number_of_slice_groups = fread(my_file, 1, 'int16');
  a.show_in_autoview = fread(my_file, 1, 'int16');
  a.slice_number_inGroup = fread(my_file, 1, 'int16');
  a.specnuc = fread(my_file, 1, 'int16');
  a.label_duration = fread(my_file, 1, 'uint16');
  a.ihbsoffsetfreq = fread(my_file, 1, 'int16');
  a.scale_factor = fread(my_file, 1, 'int16');
  a.volume_prop = fread(my_file, 1, 'int16');
  a.excitation_mode = fread(my_file, 1, 'int16');
  for id = 1 : 35
    a.short_padding(id) = fread(my_file, 1, 'int16');
  end
  a.psdname = freadc(my_file, 33);
  a.proj_name = freadc(my_file, 13);
  a.psd_iname = freadc(my_file, 13);
  a.im_diskid = fread(my_file, 1, 'char');
  for id = 1 : 14
    a.pdid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 4
    a.im_suid(id) = fread(my_file, 1, 'char');
  end
  a.contrastIV = freadc(my_file, 17);
  a.contrastOral = freadc(my_file, 17);
  a.loc_ras = fread(my_file, 1, 'char');
  a.forimgrev = freadc(my_file, 4);
  a.cname = freadc(my_file, 17);
  for id = 1 : 2
    a.im_verscre(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 2
    a.im_verscur(id) = fread(my_file, 1, 'char');
  end
  a.im_alloc_key = freadc(my_file, 13);
  a.ref_img = fread(my_file, 1, 'char');
  a.sum_img = fread(my_file, 1, 'char');
  a.filter_mode = freadc(my_file, 16);
  a.slop_str_2 = freadc(my_file, 16);
  for id = 1 : 32
    a.image_uid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.sop_uid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 24
    a.GEcname(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 100
    a.usedCoilData(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.astcalseriesuid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.purecalseriesuid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.xml_psc_shm_vol(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 64
    a.rxmpath(id) = fread(my_file, 1, 'char');
  end
  a.psdnameannot = freadc(my_file, 33);
  for id = 1 : 250
    a.img_hdr_padding(id) = fread(my_file, 1, 'char');
  end

end



function a = read_exam_header( my_file,rdbm_rev )
%read_exam_header - Read GE exam header
%
%  a = read_exam_header( my_file, rdbm_rev );
%    my_file - string indicating file name to read
%    rdbm_rev - raw header (RDBM) revision number
%    a - structure with header values
%

% Copyright (c) 2012 by General Electric Company. All rights reserved.


% RDBM revision 24.000
if rdbm_rev == 24.000 
  a.firstaxtime = fread(my_file, 1, 'float64');
  for id = 1 : 31
    a.double_padding(id) = fread(my_file, 1, 'float64');
  end
  a.zerocell = fread(my_file, 1, 'float32');
  a.cellspace = fread(my_file, 1, 'float32');
  a.srctodet = fread(my_file, 1, 'float32');
  a.srctoiso = fread(my_file, 1, 'float32');
  for id = 1 : 32
    a.float_padding(id) = fread(my_file, 1, 'float32');
  end
  a.ex_delta_cnt = fread(my_file, 1, 'int32');
  a.ex_complete = fread(my_file, 1, 'int32');
  a.ex_seriesct = fread(my_file, 1, 'int32');
  a.ex_numarch = fread(my_file, 1, 'int32');
  a.ex_numseries = fread(my_file, 1, 'int32');
  a.ex_numunser = fread(my_file, 1, 'int32');
  a.ex_toarchcnt = fread(my_file, 1, 'int32');
  a.ex_prospcnt = fread(my_file, 1, 'int32');
  a.ex_modelnum = fread(my_file, 1, 'int32');
  a.ex_modelcnt = fread(my_file, 1, 'int32');
  a.patCheckSum = fread(my_file, 1, 'int32');
  for id = 1 : 31
    a.int_padding1(id) = fread(my_file, 1, 'int32');
  end
  a.numcells = fread(my_file, 1, 'int32');
  a.magstrength = fread(my_file, 1, 'int32');
  a.patweight = fread(my_file, 1, 'int32');
  a.ex_datetime = fread(my_file, 1, 'int32');
  a.ex_lastmod = fread(my_file, 1, 'int32');
  a.patChecksumType = fread(my_file, 1, 'int32');
  for id = 1 : 26
    a.int_padding2(id) = fread(my_file, 1, 'int32');
  end
  a.ex_no = fread(my_file, 1, 'uint16');
  a.ex_uniq = fread(my_file, 1, 'int16');
  a.detect = fread(my_file, 1, 'int16');
  a.tubetyp = fread(my_file, 1, 'int16');
  a.dastyp = fread(my_file, 1, 'int16');
  a.num_dcnk = fread(my_file, 1, 'int16');
  a.dcn_len = fread(my_file, 1, 'int16');
  a.dcn_density = fread(my_file, 1, 'int16');
  a.dcn_stepsize = fread(my_file, 1, 'int16');
  a.dcn_shiftcnt = fread(my_file, 1, 'int16');
  a.patage = fread(my_file, 1, 'int16');
  a.patian = fread(my_file, 1, 'int16');
  a.patsex = fread(my_file, 1, 'int16');
  a.ex_format = fread(my_file, 1, 'int16');
  a.trauma = fread(my_file, 1, 'int16');
  a.protocolflag = fread(my_file, 1, 'int16');
  a.study_status = fread(my_file, 1, 'int16');
  for id = 1 : 35
    a.short_padding(id) = fread(my_file, 1, 'int16');
  end
  a.hist = freadc(my_file, 257);
  a.refphy = freadc(my_file, 65);
  a.diagrad = freadc(my_file, 65);
  a.operator_new = freadc(my_file, 65);
  a.ex_desc = freadc(my_file, 65);
  a.ex_typ = freadc(my_file, 3);
  a.ex_sysid = freadc(my_file, 17);
  a.ex_alloc_key = freadc(my_file, 13);
  a.ex_diskid = fread(my_file, 1, 'char');
  a.hospname = freadc(my_file, 33);
  for id = 1 : 4
    a.ex_suid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 2
    a.ex_verscre(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 2
    a.ex_verscur(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.uniq_sys_id(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.service_id(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 4
    a.mobile_loc(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.study_uid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.refsopcuid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.refsopiuid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 65
    a.patnameff(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 65
    a.patidff(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 17
    a.reqnumff(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 9
    a.dateofbirth(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 32
    a.mwlstudyuid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 16
    a.mwlstudyid(id) = fread(my_file, 1, 'char');
  end
  for id = 1 : 232
    a.ex_padding(id) = fread(my_file, 1, 'char');
  end

end



function pha_flt = get_pha_flt(pha_coe, nx)
% function pha_flt = get_pha_flt(pha_coe, nx)
%
% Calculates the phase filter in the x-ky space for EPI ghosting correction.
%
% Inputs
%   pha_coe - Coefficients for x-ky phase correction, loaded from ref.dat 
%             file or calculated from reference scan p-file or reference
%             scan k-space data. Dim: [num_coe_per_ky_line(0th order, 1st order, 2nd order...),
%             other dimensions(e.g. ny, nsl, nc].
%   nx      - Matrix size in x.
%
% Output
%   pha_flt - X-ky phase filter. xKy_corrected = xKy_uncorrected .* pha_flt.
%             Dim: [nx, other dimensions(same as other dimensions in the input pha_coe)].
%
% (c) Kangrong Zhu,     Stanford University     Nov 2013

sz = size(pha_coe);
num_coe_per_ky_line = sz(1);                                        % Number of coefficients for each ky line.

xidx = (-nx/2 : 1 : (nx/2-1)).';
pix_idx = zeros(nx, num_coe_per_ky_line);                           % Dim: [X, PixelIndex(0thOrder-1stOrder-2ndOrder...)]
for order = 0 : 1 : num_coe_per_ky_line-1
    pix_idx(:, order+1) = xidx.^order;
end
pha_coe = reshape(pha_coe, [num_coe_per_ky_line, prod(sz(2:end))]); % Dim: [PhaseCoefficients(0thOrder-1stOrder-2ndOrder...), OtherDimensions(e.g. Ky-Slice-Coil)]
pha_flt = pix_idx * pha_coe;                                        % Dim: [X, OtherDimensions(e.g. Ky-Slice-Coil)]
pha_flt = exp(1i * pha_flt);
pha_flt = reshape(pha_flt, [nx, sz(2:end)]);                        % Dim: [X, OtherDimensions(e.g. Ky, Slice, Coil)]



function a = read_rdb_hdr( my_file,rdbm_rev )
%read_rdb_hdr - Read GE raw (RDB) header
%
%  a = read_rdb_hdr( my_file, rdbm_rev );
%    my_file - string indicating file name to read
%    rdbm_rev - raw header (RDBM) revision number
%    a - structure with header values
%

% Copyright (c) 2012 by General Electric Company. All rights reserved.

% RDBM revision 24.000
if rdbm_rev == 24.000 
  a.rdbm_rev = fread(my_file, 1, 'float32');
  a.run_int = fread(my_file, 1, 'int32');
  a.scan_seq = fread(my_file, 1, 'int16');
  a.run_char = freadc(my_file, 6);
  a.scan_date = freadc(my_file, 10);
  a.scan_time = freadc(my_file, 8);
  a.logo = freadc(my_file, 10);
  a.file_contents = fread(my_file, 1, 'int16');
  a.lock_mode = fread(my_file, 1, 'int16');
  a.dacq_ctrl = fread(my_file, 1, 'int16');
  a.recon_ctrl = fread(my_file, 1, 'int16');
  a.exec_ctrl = fread(my_file, 1, 'uint16');
  a.scan_type = fread(my_file, 1, 'int16');
  a.data_collect_type = fread(my_file, 1, 'int16');
  a.data_format = fread(my_file, 1, 'int16');
  a.recon = fread(my_file, 1, 'int16');
  a.datacq = fread(my_file, 1, 'int16');
  a.npasses = fread(my_file, 1, 'int16');
  a.npomp = fread(my_file, 1, 'int16');
  a.nslices = fread(my_file, 1, 'uint16');
  a.nechoes = fread(my_file, 1, 'int16');
  a.navs = fread(my_file, 1, 'int16');
  a.nframes = fread(my_file, 1, 'int16');
  a.baseline_views = fread(my_file, 1, 'int16');
  a.hnover = fread(my_file, 1, 'int16');
  a.frame_size = fread(my_file, 1, 'uint16');
  a.point_size = fread(my_file, 1, 'int16');
  a.vquant = fread(my_file, 1, 'int16');
  a.cheart = fread(my_file, 1, 'int16');
  a.ctr = fread(my_file, 1, 'float32');
  a.ctrr = fread(my_file, 1, 'float32');
  a.initpass = fread(my_file, 1, 'int16');
  a.incrpass = fread(my_file, 1, 'int16');
  a.method_ctrl = fread(my_file, 1, 'int16');
  a.da_xres = fread(my_file, 1, 'uint16');
  a.da_yres = fread(my_file, 1, 'int16');
  a.rc_xres = fread(my_file, 1, 'int16');
  a.rc_yres = fread(my_file, 1, 'int16');
  a.im_size = fread(my_file, 1, 'int16');
  a.rc_zres = fread(my_file, 1, 'int32');
  a.raw_pass_size_deprecated = fread(my_file, 1, 'int32');
  a.sspsave_deprecated = fread(my_file, 1, 'int32');
  a.udasave_deprecated = fread(my_file, 1, 'int32');
  a.fermi_radius = fread(my_file, 1, 'float32');
  a.fermi_width = fread(my_file, 1, 'float32');
  a.fermi_ecc = fread(my_file, 1, 'float32');
  a.clip_min = fread(my_file, 1, 'float32');
  a.clip_max = fread(my_file, 1, 'float32');
  a.default_offset = fread(my_file, 1, 'float32');
  a.xoff = fread(my_file, 1, 'float32');
  a.yoff = fread(my_file, 1, 'float32');
  a.nwin = fread(my_file, 1, 'float32');
  a.ntran = fread(my_file, 1, 'float32');
  a.scalei = fread(my_file, 1, 'float32');
  a.scaleq = fread(my_file, 1, 'float32');
  a.rotation = fread(my_file, 1, 'int16');
  a.transpose = fread(my_file, 1, 'int16');
  a.kissoff_views = fread(my_file, 1, 'int16');
  a.slblank = fread(my_file, 1, 'int16');
  a.gradcoil = fread(my_file, 1, 'int16');
  a.ddaover = fread(my_file, 1, 'int16');
  a.sarr = fread(my_file, 1, 'int16');
  a.fd_tr = fread(my_file, 1, 'int16');
  a.fd_te = fread(my_file, 1, 'int16');
  a.fd_ctrl = fread(my_file, 1, 'int16');
  a.algor_num = fread(my_file, 1, 'int16');
  a.fd_df_dec = fread(my_file, 1, 'int16');
  for id = 1 : 8
    a.dab(id) = fread(my_file, 1, 'int16');
  end
  a.user0 = fread(my_file, 1, 'float32');
  a.user1 = fread(my_file, 1, 'float32');
  a.user2 = fread(my_file, 1, 'float32');
  a.user3 = fread(my_file, 1, 'float32');
  a.user4 = fread(my_file, 1, 'float32');
  a.user5 = fread(my_file, 1, 'float32');
  a.user6 = fread(my_file, 1, 'float32');
  a.user7 = fread(my_file, 1, 'float32');
  a.user8 = fread(my_file, 1, 'float32');
  a.user9 = fread(my_file, 1, 'float32');
  a.user10 = fread(my_file, 1, 'float32');
  a.user11 = fread(my_file, 1, 'float32');
  a.user12 = fread(my_file, 1, 'float32');
  a.user13 = fread(my_file, 1, 'float32');
  a.user14 = fread(my_file, 1, 'float32');
  a.user15 = fread(my_file, 1, 'float32');
  a.user16 = fread(my_file, 1, 'float32');
  a.user17 = fread(my_file, 1, 'float32');
  a.user18 = fread(my_file, 1, 'float32');
  a.user19 = fread(my_file, 1, 'float32');
  a.v_type = fread(my_file, 1, 'int32');
  a.v_coefxa = fread(my_file, 1, 'float32');
  a.v_coefxb = fread(my_file, 1, 'float32');
  a.v_coefxc = fread(my_file, 1, 'float32');
  a.v_coefxd = fread(my_file, 1, 'float32');
  a.v_coefya = fread(my_file, 1, 'float32');
  a.v_coefyb = fread(my_file, 1, 'float32');
  a.v_coefyc = fread(my_file, 1, 'float32');
  a.v_coefyd = fread(my_file, 1, 'float32');
  a.v_coefza = fread(my_file, 1, 'float32');
  a.v_coefzb = fread(my_file, 1, 'float32');
  a.v_coefzc = fread(my_file, 1, 'float32');
  a.v_coefzd = fread(my_file, 1, 'float32');
  a.vm_coef1 = fread(my_file, 1, 'float32');
  a.vm_coef2 = fread(my_file, 1, 'float32');
  a.vm_coef3 = fread(my_file, 1, 'float32');
  a.vm_coef4 = fread(my_file, 1, 'float32');
  a.v_venc = fread(my_file, 1, 'float32');
  a.spectral_width = fread(my_file, 1, 'float32');
  a.csi_dims = fread(my_file, 1, 'int16');
  a.xcsi = fread(my_file, 1, 'int16');
  a.ycsi = fread(my_file, 1, 'int16');
  a.zcsi = fread(my_file, 1, 'int16');
  a.roilenx = fread(my_file, 1, 'float32');
  a.roileny = fread(my_file, 1, 'float32');
  a.roilenz = fread(my_file, 1, 'float32');
  a.roilocx = fread(my_file, 1, 'float32');
  a.roilocy = fread(my_file, 1, 'float32');
  a.roilocz = fread(my_file, 1, 'float32');
  a.numdwell = fread(my_file, 1, 'float32');
  a.ps_command = fread(my_file, 1, 'int32');
  a.ps_mps_r1 = fread(my_file, 1, 'int32');
  a.ps_mps_r2 = fread(my_file, 1, 'int32');
  a.ps_mps_tg = fread(my_file, 1, 'int32');
  a.ps_mps_freq = fread(my_file, 1, 'uint32');
  a.ps_aps_r1 = fread(my_file, 1, 'int32');
  a.ps_aps_r2 = fread(my_file, 1, 'int32');
  a.ps_aps_tg = fread(my_file, 1, 'int32');
  a.ps_aps_freq = fread(my_file, 1, 'uint32');
  a.ps_scalei = fread(my_file, 1, 'float32');
  a.ps_scaleq = fread(my_file, 1, 'float32');
  a.ps_snr_warning = fread(my_file, 1, 'int32');
  a.ps_aps_or_mps = fread(my_file, 1, 'int32');
  a.ps_mps_bitmap = fread(my_file, 1, 'int32');
  a.ps_powerspec = freadc(my_file, 256);
  a.ps_filler1 = fread(my_file, 1, 'int32');
  a.ps_filler2 = fread(my_file, 1, 'int32');
  for id = 1 : 16
    a.obsolete1(id) = fread(my_file, 1, 'float32');
  end
  for id = 1 : 16
    a.obsolete2(id) = fread(my_file, 1, 'float32');
  end
  a.halfecho = fread(my_file, 1, 'int16');
  a.im_size_y = fread(my_file, 1, 'int16');
  a.data_collect_type1 = fread(my_file, 1, 'int32');
  a.freq_scale = fread(my_file, 1, 'float32');
  a.phase_scale = fread(my_file, 1, 'float32');
  a.ovl = fread(my_file, 1, 'int16');
  a.pclin = fread(my_file, 1, 'int16');
  a.pclinnpts = fread(my_file, 1, 'int16');
  a.pclinorder = fread(my_file, 1, 'int16');
  a.pclinavg = fread(my_file, 1, 'int16');
  a.pccon = fread(my_file, 1, 'int16');
  a.pcconnpts = fread(my_file, 1, 'int16');
  a.pcconorder = fread(my_file, 1, 'int16');
  a.pcextcorr = fread(my_file, 1, 'int16');
  a.pcgraph = fread(my_file, 1, 'int16');
  a.pcileave = fread(my_file, 1, 'int16');
  a.hdbestky = fread(my_file, 1, 'int16');
  a.pcctrl = fread(my_file, 1, 'int16');
  a.pcthrespts = fread(my_file, 1, 'int16');
  a.pcdiscbeg = fread(my_file, 1, 'int16');
  a.pcdiscmid = fread(my_file, 1, 'int16');
  a.pcdiscend = fread(my_file, 1, 'int16');
  a.pcthrespct = fread(my_file, 1, 'int16');
  a.pcspacial = fread(my_file, 1, 'int16');
  a.pctemporal = fread(my_file, 1, 'int16');
  a.pcspare = fread(my_file, 1, 'int16');
  a.ileaves = fread(my_file, 1, 'int16');
  a.kydir = fread(my_file, 1, 'int16');
  a.alt = fread(my_file, 1, 'int16');
  a.reps = fread(my_file, 1, 'int16');
  a.ref = fread(my_file, 1, 'int16');
  a.pcconnorm = fread(my_file, 1, 'float32');
  a.pcconfitwt = fread(my_file, 1, 'float32');
  a.pclinnorm = fread(my_file, 1, 'float32');
  a.pclinfitwt = fread(my_file, 1, 'float32');
  a.pcbestky = fread(my_file, 1, 'float32');
  a.vrgf = fread(my_file, 1, 'int32');
  a.vrgfxres = fread(my_file, 1, 'int32');
  a.bp_corr = fread(my_file, 1, 'int32');
  a.recv_freq_s = fread(my_file, 1, 'float32');
  a.recv_freq_e = fread(my_file, 1, 'float32');
  a.hniter = fread(my_file, 1, 'int32');
  a.fast_rec = fread(my_file, 1, 'int32');
  a.refframes = fread(my_file, 1, 'int32');
  a.refframep = fread(my_file, 1, 'int32');
  a.scnframe = fread(my_file, 1, 'int32');
  a.pasframe = fread(my_file, 1, 'int32');
  a.user_usage_tag = fread(my_file, 1, 'uint32');
  a.user_fill_mapMSW = fread(my_file, 1, 'uint32');
  a.user_fill_mapLSW = fread(my_file, 1, 'uint32');
  a.user20 = fread(my_file, 1, 'float32');
  a.user21 = fread(my_file, 1, 'float32');
  a.user22 = fread(my_file, 1, 'float32');
  a.user23 = fread(my_file, 1, 'float32');
  a.user24 = fread(my_file, 1, 'float32');
  a.user25 = fread(my_file, 1, 'float32');
  a.user26 = fread(my_file, 1, 'float32');
  a.user27 = fread(my_file, 1, 'float32');
  a.user28 = fread(my_file, 1, 'float32');
  a.user29 = fread(my_file, 1, 'float32');
  a.user30 = fread(my_file, 1, 'float32');
  a.user31 = fread(my_file, 1, 'float32');
  a.user32 = fread(my_file, 1, 'float32');
  a.user33 = fread(my_file, 1, 'float32');
  a.user34 = fread(my_file, 1, 'float32');
  a.user35 = fread(my_file, 1, 'float32');
  a.user36 = fread(my_file, 1, 'float32');
  a.user37 = fread(my_file, 1, 'float32');
  a.user38 = fread(my_file, 1, 'float32');
  a.user39 = fread(my_file, 1, 'float32');
  a.user40 = fread(my_file, 1, 'float32');
  a.user41 = fread(my_file, 1, 'float32');
  a.user42 = fread(my_file, 1, 'float32');
  a.user43 = fread(my_file, 1, 'float32');
  a.user44 = fread(my_file, 1, 'float32');
  a.user45 = fread(my_file, 1, 'float32');
  a.user46 = fread(my_file, 1, 'float32');
  a.user47 = fread(my_file, 1, 'float32');
  a.user48 = fread(my_file, 1, 'float32');
  a.pcfitorig = fread(my_file, 1, 'int16');
  a.pcshotfirst = fread(my_file, 1, 'int16');
  a.pcshotlast = fread(my_file, 1, 'int16');
  a.pcmultegrp = fread(my_file, 1, 'int16');
  a.pclinfix = fread(my_file, 1, 'int16');
  a.pcconfix = fread(my_file, 1, 'int16');
  a.pclinslope = fread(my_file, 1, 'float32');
  a.pcconslope = fread(my_file, 1, 'float32');
  a.pccoil = fread(my_file, 1, 'int16');
  a.vvsmode = fread(my_file, 1, 'int16');
  a.vvsaimgs = fread(my_file, 1, 'int16');
  a.vvstr = fread(my_file, 1, 'int16');
  a.vvsgender = fread(my_file, 1, 'int16');
  a.zip_factor = fread(my_file, 1, 'int16');
  a.maxcoef1a = fread(my_file, 1, 'float32');
  a.maxcoef1b = fread(my_file, 1, 'float32');
  a.maxcoef1c = fread(my_file, 1, 'float32');
  a.maxcoef1d = fread(my_file, 1, 'float32');
  a.maxcoef2a = fread(my_file, 1, 'float32');
  a.maxcoef2b = fread(my_file, 1, 'float32');
  a.maxcoef2c = fread(my_file, 1, 'float32');
  a.maxcoef2d = fread(my_file, 1, 'float32');
  a.maxcoef3a = fread(my_file, 1, 'float32');
  a.maxcoef3b = fread(my_file, 1, 'float32');
  a.maxcoef3c = fread(my_file, 1, 'float32');
  a.maxcoef3d = fread(my_file, 1, 'float32');
  a.ut_ctrl = fread(my_file, 1, 'int32');
  a.dp_type = fread(my_file, 1, 'int16');
  a.arw = fread(my_file, 1, 'int16');
  a.vps = fread(my_file, 1, 'int16');
  a.mcReconEnable = fread(my_file, 1, 'int16');
  a.fov = fread(my_file, 1, 'float32');
  a.te = fread(my_file, 1, 'int32');
  a.te2 = fread(my_file, 1, 'int32');
  a.dfmrbw = fread(my_file, 1, 'float32');
  a.dfmctrl = fread(my_file, 1, 'int32');
  a.raw_nex = fread(my_file, 1, 'int32');
  a.navs_per_pass = fread(my_file, 1, 'int32');
  a.dfmxres = fread(my_file, 1, 'int32');
  a.dfmptsize = fread(my_file, 1, 'int32');
  a.navs_per_view = fread(my_file, 1, 'int32');
  a.dfmdebug = fread(my_file, 1, 'int32');
  a.dfmthreshold = fread(my_file, 1, 'float32');
  a.grid_control = fread(my_file, 1, 'int16');
  a.b0map = fread(my_file, 1, 'int16');
  a.grid_tediff = fread(my_file, 1, 'int16');
  a.grid_motion_comp = fread(my_file, 1, 'int16');
  a.grid_radius_a = fread(my_file, 1, 'float32');
  a.grid_radius_b = fread(my_file, 1, 'float32');
  a.grid_max_gradient = fread(my_file, 1, 'float32');
  a.grid_max_slew = fread(my_file, 1, 'float32');
  a.grid_scan_fov = fread(my_file, 1, 'float32');
  a.grid_a2d_time = fread(my_file, 1, 'float32');
  a.grid_density_factor = fread(my_file, 1, 'float32');
  a.grid_display_fov = fread(my_file, 1, 'float32');
  a.fatwater = fread(my_file, 1, 'int16');
  a.fiestamlf = fread(my_file, 1, 'int16');
  a.app = fread(my_file, 1, 'int16');
  a.rhncoilsel = fread(my_file, 1, 'int16');
  a.rhncoillimit = fread(my_file, 1, 'int16');
  a.app_option = fread(my_file, 1, 'int16');
  a.grad_mode = fread(my_file, 1, 'int16');
  a.pfile_passes = fread(my_file, 1, 'int16');
  a.asset = fread(my_file, 1, 'int32');
  a.asset_calthresh = fread(my_file, 1, 'int32');
  a.asset_R = fread(my_file, 1, 'float32');
  a.coilConfigUID = fread(my_file, 1, 'uint32');
  a.asset_phases = fread(my_file, 1, 'int32');
  a.scancent = fread(my_file, 1, 'float32');
  a.position = fread(my_file, 1, 'int32');
  a.entry = fread(my_file, 1, 'int32');
  a.lmhor = fread(my_file, 1, 'float32');
  a.last_slice_num = fread(my_file, 1, 'int32');
  a.asset_slice_R = fread(my_file, 1, 'float32');
  a.asset_slabwrap = fread(my_file, 1, 'float32');
  a.dwnav_coeff = fread(my_file, 1, 'float32');
  a.dwnav_cor = fread(my_file, 1, 'int16');
  a.dwnav_view = fread(my_file, 1, 'int16');
  a.dwnav_corecho = fread(my_file, 1, 'int16');
  a.dwnav_sview = fread(my_file, 1, 'int16');
  a.dwnav_eview = fread(my_file, 1, 'int16');
  a.dwnav_sshot = fread(my_file, 1, 'int16');
  a.dwnav_eshot = fread(my_file, 1, 'int16');
  a.a3dwin_type = fread(my_file, 1, 'int16');
  a.a3dwin_apod = fread(my_file, 1, 'float32');
  a.a3dwin_q = fread(my_file, 1, 'float32');
  a.ime_scic_enable = fread(my_file, 1, 'int16');
  a.clariview_type = fread(my_file, 1, 'int16');
  a.ime_scic_edge = fread(my_file, 1, 'float32');
  a.ime_scic_smooth = fread(my_file, 1, 'float32');
  a.ime_scic_focus = fread(my_file, 1, 'float32');
  a.clariview_edge = fread(my_file, 1, 'float32');
  a.clariview_smooth = fread(my_file, 1, 'float32');
  a.clariview_focus = fread(my_file, 1, 'float32');
  a.scic_reduction = fread(my_file, 1, 'float32');
  a.scic_gauss = fread(my_file, 1, 'float32');
  a.scic_threshold = fread(my_file, 1, 'float32');
  a.ectricks_no_regions = fread(my_file, 1, 'int32');
  a.ectricks_input_regions = fread(my_file, 1, 'int32');
  a.psc_reuse = fread(my_file, 1, 'int16');
  a.left_blank = fread(my_file, 1, 'int16');
  a.right_blank = fread(my_file, 1, 'int16');
  a.acquire_type = fread(my_file, 1, 'int16');
  a.retro_control = fread(my_file, 1, 'int16');
  a.etl = fread(my_file, 1, 'int16');
  a.pcref_start = fread(my_file, 1, 'int16');
  a.pcref_stop = fread(my_file, 1, 'int16');
  a.ref_skip = fread(my_file, 1, 'int16');
  a.extra_frames_top = fread(my_file, 1, 'int16');
  a.extra_frames_bot = fread(my_file, 1, 'int16');
  a.multiphase_type = fread(my_file, 1, 'int16');
  a.nphases = fread(my_file, 1, 'int16');
  a.pure = fread(my_file, 1, 'int16');
  a.pure_scale = fread(my_file, 1, 'float32');
  a.off_data = fread(my_file, 1, 'int32');
  a.off_per_pass = fread(my_file, 1, 'int32');
  a.off_unlock_raw = fread(my_file, 1, 'int32');
  a.off_data_acq_tab = fread(my_file, 1, 'int32');
  a.off_nex_tab = fread(my_file, 1, 'int32');
  a.off_nex_abort_tab = fread(my_file, 1, 'int32');
  a.off_tool = fread(my_file, 1, 'int32');
  a.off_exam = fread(my_file, 1, 'int32');
  a.off_series = fread(my_file, 1, 'int32');
  a.off_image = fread(my_file, 1, 'int32');
  a.off_ps = fread(my_file, 1, 'int32');
  a.off_spare_b = fread(my_file, 1, 'int32');
  a.new_wnd_level_flag = fread(my_file, 1, 'int32');
  a.wnd_image_hist_area = fread(my_file, 1, 'int32');
  a.wnd_high_hist = fread(my_file, 1, 'float32');
  a.wnd_lower_hist = fread(my_file, 1, 'float32');
  a.pure_filter = fread(my_file, 1, 'int16');
  a.cfg_pure_filter = fread(my_file, 1, 'int16');
  a.cfg_pure_fit_order = fread(my_file, 1, 'int16');
  a.cfg_pure_kernelsize_z = fread(my_file, 1, 'int16');
  a.cfg_pure_kernelsize_xy = fread(my_file, 1, 'int16');
  a.cfg_pure_weight_radius = fread(my_file, 1, 'int16');
  a.cfg_pure_intensity_scale = fread(my_file, 1, 'int16');
  a.cfg_pure_noise_threshold = fread(my_file, 1, 'int16');
  a.wienera = fread(my_file, 1, 'float32');
  a.wienerb = fread(my_file, 1, 'float32');
  a.wienert2 = fread(my_file, 1, 'float32');
  a.wieneresp = fread(my_file, 1, 'float32');
  a.wiener = fread(my_file, 1, 'int16');
  a.flipfilter = fread(my_file, 1, 'int16');
  a.dbgrecon = fread(my_file, 1, 'int16');
  a.ech2skip = fread(my_file, 1, 'int16');
  a.tricks_type = fread(my_file, 1, 'int32');
  a.lcfiesta_phase = fread(my_file, 1, 'float32');
  a.lcfiesta = fread(my_file, 1, 'int16');
  a.herawflt = fread(my_file, 1, 'int16');
  a.herawflt_befnwin = fread(my_file, 1, 'int16');
  a.herawflt_befntran = fread(my_file, 1, 'int16');
  a.herawflt_befamp = fread(my_file, 1, 'float32');
  a.herawflt_hpfamp = fread(my_file, 1, 'float32');
  a.heover = fread(my_file, 1, 'int16');
  a.pure_correction_threshold = fread(my_file, 1, 'int16');
  a.swiftenable = fread(my_file, 1, 'int32');
  a.numslabs = fread(my_file, 1, 'int16');
  a.numCoilConfigs = fread(my_file, 1, 'uint16');
  a.ps_autoshim_status = fread(my_file, 1, 'int32');
  a.dynaplan_numphases = fread(my_file, 1, 'int32');
  a.medal_cfg = fread(my_file, 1, 'int16');
  a.medal_nstack = fread(my_file, 1, 'int16');
  a.medal_echo_order = fread(my_file, 1, 'int16');
  a.medal_kernel_up = fread(my_file, 1, 'int16');
  a.medal_kernel_down = fread(my_file, 1, 'int16');
  a.medal_kernel_smooth = fread(my_file, 1, 'int16');
  a.medal_start = fread(my_file, 1, 'int16');
  a.medal_end = fread(my_file, 1, 'int16');
  a.rcideal = fread(my_file, 1, 'uint32');
  a.rcdixproc = fread(my_file, 1, 'uint32');
  a.df = fread(my_file, 1, 'float32');
  a.bw = fread(my_file, 1, 'float32');
  a.te1_deprecated = fread(my_file, 1, 'float32');
  a.esp_deprecated = fread(my_file, 1, 'float32');
  a.feextra = fread(my_file, 1, 'int32');
  a.raw_pass_size = fread(my_file, 1, 'uint64');
  a.sspsave = fread(my_file, 1, 'uint64');
  a.udasave = fread(my_file, 1, 'uint64');
  a.vibrant = fread(my_file, 1, 'int16');
  a.asset_torso = fread(my_file, 1, 'int16');
  a.asset_alt_cal = fread(my_file, 1, 'int32');
  a.kacq_uid = fread(my_file, 1, 'int32');
  for id = 1 : 4
    a.cttEntry(id).logicalCoilName = fread(my_file, 128, 'char')';
    a.cttEntry(id).clinicalCoilName = fread(my_file, 32, 'char')';
    a.cttEntry(id).configUID = fread(my_file, 1, 'uint32');
    a.cttEntry(id).coilConnector = fread(my_file, 1, 'int32');
    a.cttEntry(id).isActiveConfig = fread(my_file, 1, 'uint32');
    a.cttEntry(id).channelTranslationMap = fread(my_file, 32, 'int16')';
  end
  a.psc_ta = fread(my_file, 1, 'int32');
  a.disk_acq_ctrl = fread(my_file, 1, 'int32');
  a.asset_localTx = fread(my_file, 1, 'int32');
  a.a3dscale = fread(my_file, 1, 'float32');
  a.broad_band_select = fread(my_file, 1, 'int32');
  a.scanner_mode = fread(my_file, 1, 'int16');
  a.numbvals = fread(my_file, 1, 'int16');
  a.numdifdirs = fread(my_file, 1, 'int16');
  a.difnext2 = fread(my_file, 1, 'int16');
  for id = 1 : 100
    a.difnextab(id) = fread(my_file, 1, 'int16');
  end
  a.channel_combine_method = fread(my_file, 1, 'int16');
  a.nexForUnacquiredEncodes = fread(my_file, 1, 'int16');
  a.a2dscale = fread(my_file, 1, 'float32');
  a.dd_mode = fread(my_file, 1, 'int16');
  a.dd_q_ta_offset = fread(my_file, 1, 'int16');
  a.dd_q_phase_delay = fread(my_file, 1, 'float32');
  a.dacq_ctrl_chksum = fread(my_file, 1, 'uint32');
  a.patient_checksum = fread(my_file, 1, 'uint32');
  a.rcidealctrl = fread(my_file, 1, 'uint32');
  for id = 1 : 64
    a.echotimes(id) = fread(my_file, 1, 'float32');
  end
  a.asl_perf_weighted_scale = fread(my_file, 1, 'int16');
  a.echo_pc_extra_frames_bot = fread(my_file, 1, 'uint16');
  a.echo_pc_ctrl = fread(my_file, 1, 'uint32');
  a.echo_pc_ylines = fread(my_file, 1, 'uint16');
  a.echo_pc_primary_yfirst = fread(my_file, 1, 'uint16');
  a.echo_pc_reverse_yfirst = fread(my_file, 1, 'uint16');
  a.echo_pc_zlines = fread(my_file, 1, 'uint16');
  a.echo_pc_yxfitorder = fread(my_file, 1, 'uint16');
  a.channel_combine_filter_type = fread(my_file, 1, 'int16');
  a.mavric_control = fread(my_file, 1, 'uint32');
  a.mavric_ImageType = fread(my_file, 1, 'uint32');
  a.mavric_bin_separation = fread(my_file, 1, 'int32');
  for id = 1 : 40
    a.mavric_b0_offset(id) = fread(my_file, 1, 'float32');
  end
  a.channel_combine_filter_width = fread(my_file, 1, 'float32');
  a.channel_combine_filter_beta = fread(my_file, 1, 'float32');
  a.low_pass_nex_filter_width = fread(my_file, 1, 'float32');
  a.aps_tg_status = fread(my_file, 1, 'uint32');
  a.cal_pass_set_vector = fread(my_file, 1, 'int32');
  a.cal_nex_vector = fread(my_file, 1, 'int32');
  a.cal_weight_vector = fread(my_file, 1, 'int32');
  a.pure_filtering_mode = fread(my_file, 1, 'int32');
  a.pure_lambda = fread(my_file, 1, 'float32');
  a.pure_tuning_factor_surface = fread(my_file, 1, 'float32');
  a.pure_tuning_factor_body = fread(my_file, 1, 'float32');
  a.pure_derived_cal_fraction = fread(my_file, 1, 'float32');
  a.pure_derived_cal_reapodization = fread(my_file, 1, 'float32');
  a.noncart_grid_factor = fread(my_file, 1, 'float32');
  a.noncart_dual_traj = fread(my_file, 1, 'int32');
  a.noncart_traj_kmax_ratio = fread(my_file, 1, 'int32');
  a.noncart_traj_merge_start = fread(my_file, 1, 'int32');
  a.noncart_traj_merge_end = fread(my_file, 1, 'int32');
  a.oversamplingfactor = fread(my_file, 1, 'float32');
  a.nspokes_highres = fread(my_file, 1, 'int32');
  a.nspokes_lowres = fread(my_file, 1, 'int32');
  a.nrefslices = fread(my_file, 1, 'int32');
  a.psmde_pixel_offset = fread(my_file, 1, 'int32');
  a.off_grad_data = fread(my_file, 1, 'int32');
  a.hoecc = fread(my_file, 1, 'int32');
  a.hoec_fit_order = fread(my_file, 1, 'int32');
  a.esp = fread(my_file, 1, 'int32');
  a.excess = fread(my_file, 322, 'int16');

end



function res = zpad(x, varargin)
%
% res = zpad(x, s1, s2)
% Zero pads a 2D matrix around its center to a size of [s1, s2].
%
% res = zpad(x, [s1, s2])
% Same as the previous example.
%
% res = zpad(x, s1, s2, s3)
% Zero pads a 3D matrix around its center.
%
% ...
%
% res = zpad(x, s1, s2, s3, s4, s5, s6)
% Zero pads a 6D matrix around its center.
%
% (c) Michael Lustig,           Stanford University
% Modified by Kangrong Zhu,     Stanford University      2012

s = cat(2, varargin{1:end});

m = size(x);
if length(m) < length(s)
    m = [ m, ones(1, length(s)-length(m))];
end

if all(m == s)
    res = x;
    return;
end

res = zeros(s);%, 'single');
idx = cell(1, numel(s));
for n = 1:numel(s)
    idx{n} = floor(s(n)/2)+1+ceil(-m(n)/2) : floor(s(n)/2)+ceil(m(n)/2);
end

res(idx{:}) = x;



function C = freadc (fid, len)
C = deblank(fread(fid, [1, len], 'uchar=>char'));



function a = read_data_acq_tab( my_file, rdbm_rev)
%read_data_acq_tab - Read GE data acquisition table
%
%  a = read_data_acq_tab( my_file, rdbm_rev );
%    my_file - string indicating file name to read
%    rdbm_rev - raw header (RDBM) revision number
%    a - structure with header values
%

% Copyright (c) 2012 by General Electric Company. All rights reserved.

% RDBM revision 24.000
if rdbm_rev == 24.000
  a.pass_number = fread(my_file, 1, 'int16');
  a.slice_in_pass = fread(my_file, 1, 'int16');
  for id = 1 : 3
    a.gw_point1(id) = fread(my_file, 1, 'float32');
  end
  for id = 1 : 3
    a.gw_point2(id) = fread(my_file, 1, 'float32');
  end
  for id = 1 : 3
    a.gw_point3(id) = fread(my_file, 1, 'float32');
  end
  a.transpose = fread(my_file, 1, 'int16');
  a.rotate = fread(my_file, 1, 'int16');
  a.coilConfigUID = fread(my_file, 1, 'uint32');


end



% 
% read_psc_header.m
%
% Copyright (c) 2012 The General Electric Company
%
%
function  a = read_psc_header( my_file,rdbm_rev )
if rdbm_rev >= 20.006
 a.command = fread(my_file, 1, 'int32');
 a.mps_r1 = fread(my_file, 1, 'int32');
 a.mps_r2 = fread(my_file, 1, 'int32');
 a.mps_tg = fread(my_file, 1, 'int32');
 a.mps_freq = fread(my_file, 1, 'uint32');
 a.aps_r1 = fread(my_file, 1, 'int32');
 a.aps_r2 = fread(my_file, 1, 'int32');
 a.aps_tg = fread(my_file, 1, 'int32');
 a.aps_freq = fread(my_file, 1, 'uint32');
 a.scalei = fread(my_file, 1, 'float32');
 a.scaleq = fread(my_file, 1, 'float32');
 a.snr_warning = fread(my_file, 1, 'int32');
 a.aps_or_mps = fread(my_file, 1, 'int32');
 a.mps_bitmap = fread(my_file, 1, 'int32');
 for id = 1:256
  a.powerspec(id) = fread(my_file, 1, 'char');
 end
 a.filler1 = fread(my_file, 1, 'int32');
 a.filler2 = fread(my_file, 1, 'int32');
 a.xshim = fread(my_file, 1, 'int16');
 a.yshim = fread(my_file, 1, 'int16');
 a.zshim = fread(my_file, 1, 'int16');
 a.recon_enable = fread(my_file, 1, 'int16');
 a.autoshim_status = fread(my_file, 1, 'int32'); 
 for id = 1 : 128
  a.rec_std(id) = fread(my_file, 1, 'float32');
 end
 for id = 1 : 128
  a.rec_mean(id) = fread(my_file, 1, 'float32');
 end
 a.line_width = fread(my_file, 1, 'int32');
 a.ws_flip = fread(my_file, 1, 'int32');
 a.supp_lvl = fread(my_file, 1, 'int32');
 a.psc_reuse = fread(my_file, 1, 'int32');
 for id = 1:52
  a.psc_reuse_string(id) = fread(my_file, 1, 'char');
 end 
 a.psc_ta = fread(my_file, 1, 'int32');
 a.phase_correction_status = fread(my_file, 1, 'int32');
 a.broad_band_select = fread(my_file, 1, 'int32');
 a.dd_q_phase_delay_qd = fread(my_file, 1, 'float32');      
 a.dd_q_phase_delay_from_qd = fread(my_file, 1, 'float32');
 a.dd_q_ta_offset_qd = fread(my_file, 1, 'int16');      
 a.dd_q_ta_offset_from_qd = fread(my_file, 1, 'int16');
 a.dd_mode = fread(my_file, 1, 'int16');              
 a.dummy_for_32bit_align = fread(my_file, 1, 'int16');
 for id = 1:48
  a.buffer(id) = fread(my_file, 1, 'char');
 end
end



function dftz_mtx = encode_dftz_mtx(fov_shift, nz)
%
% function dftz_mtx = encode_dftz_mtx(fov_shift, nz)
%
% Calculate Discrete Fourier Transform (DFT) encoding matrix for the slice dimension.
%
% Inputs
%   fov_shift - CAIPI FOV shift. abs(fov_shift) = nkz is the length of the DFTz encoding. Negative for inverse DFTz encoding, positive for forward DFTz encoding.
%   nz        - Number of simultaneous slices.
%
% Output
%   dftz_mtx  - DFTz encoding matrix. Dim: [nkz, nz(z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))].
%
% (c) Kangrong Zhu      Stanford University     Mar 2014

nkz = abs(fov_shift);
msk.omegaz = encode_dftz_omegaz(fov_shift); % Dim: [1, nkz]
msk.kz = 1 : nkz;
dftz_mtx = encode_ftyz_mtx(msk, 1, nz, 0);  % Dim: [nkz, nz]



function [acs_pos, sz_acs] = grappa_acspos(sz_dat, sz_acs, sz_lim, debug)
%
% function [acs_pos, sz_acs] = grappa_acspos(sz_dat, sz_acs[, sz_lim, debug])
%
% Calculates the indices of the Autocalibration area for GRAPPA reconstruction. 
%
% Inputs:
%   sz_dat  - Size of the whole k-space matrix. A structure array with fields 'x' and 'y'.
%   sz_acs  - Size of the ACS area. A structure array with fields 'x' and 'y'.
%   sz_lim  - Limits on acs_pos. A structure array with fields 'x' and/or 'y'.
%             Will check the x/y indices if field x/y exists and is non-empty.
%   debug   - True: print debug messages. Only needed when sz_lim has a
%             non-empty 'x' and/or 'y' field.
%
% Output:
%   acs_pos - Indices for the ACS area. A structure array with fields 'x' and 'y'.
%   sz_acs  - Size of the ACS area, having the same fields as the input 'sz_acs'.
%             The values might have been changed in this function if partial k is used.
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

if ~exist('sz_lim', 'var')
    sz_lim = [];
end

if ~exist('debug', 'var') || isempty(debug)
    debug = false;
end

acs_pos = struct('x', {[]}, 'y', {[]});

acs_pos.x = get_acs_pos(sz_dat.x, sz_acs.x);
acs_pos.y = get_acs_pos(sz_dat.y, sz_acs.y);

% Remove the out-of-range k-space indices. Useful for partial k acquisitions.
if ~isempty(sz_lim)
    if isfield(sz_lim, 'x') && ~isempty(sz_lim.x)
        [acs_pos.x, sz_acs.x] = checkrange(acs_pos.x, sz_lim.x(1), sz_lim.x(2), 'ACS kx positions', debug); % acs_pos.x, sz_acs.x change if partial kx 
    end
    if isfield(sz_lim, 'y') && ~isempty(sz_lim.y)
        [acs_pos.y, sz_acs.y] = checkrange(acs_pos.y, sz_lim.y(1), sz_lim.y(2), 'ACS ky positions', debug); % acs_pos.y, sz_acs.y change if partial ky
    end
end

function pos = get_acs_pos(n_total, n_acs)

pos = floor(n_total/2) + ( -floor(n_acs/2)+1 : ceil(n_acs/2) );



function [arrout, numelout] = checkrange(arr, amin, amax, atype, debug)
%
% Check values in array are within range, and remove ones that are not.
% 
% Inputs:
%   arr      - The input array.
%   amin     - Minimum for the range(inclusive).
%   amax     - Maximum for the range(inclusive).
%   atype    - A string specifying the type of the array. Used for
%              displaying the out-of-range elements in the array.
%              Only needed when debug==true.
%   debug    - True: Print out some messages for debugging;
%              False(Default): No messages.
%
% Output:
%   arrout   - The output array with out-of-range elements removed.
%   numelout - The number of elements in the output array 'arrout';
%
% (c) Brian Hargreaves,         Stanford University     2008
% Modified by Kangrong Zhu,     Stanford Univeristy     2012

if ~exist('debug', 'var') || isempty('debug')
    debug = false;
end

arrout = arr(arr >= amin);
arrout = arrout(arrout <= amax);

numelout = numel(arrout);

if debug && (numelout < numel(arr))
    arr = arr(:);
    f1 = find(arr < amin);
    f2 = find(arr > amax);
    fprintf('   %d %s out of range: %s \n', ...
        length(f1)+length(f2), atype, num2str([arr(f1); arr(f2)].'));
    fprintf('   Removed those %s from the array.\n', atype);
end



function p = mux_epi_params_set_tpoints(p, nt_to_recon)
% function p = mux_epi_params_set_tpoints(p, nt_to_recon)
%
% Set up time-point-related fields in the parameter structure.
%
% Inputs
%   p           - Parameter structure. See mux_epi_params.m for details.
%                 Fields used in this function: num_passes, mux, num_mux_cycle.
%   nt_to_recon - Number of time points to reconstruct, excluding the first
%                 few mux phase cycling time points. Default: All time points.
%
% Output
%   p           - Parameter structure, with the following fields updated according to the input value of nt_to_recon:
%                 nt_to_recon - Number of accelerated time points(excluding the first few slice phase cycling time points) to reconstruct.
%                 tpoints_to_load - Time points to be loaded from the mux epi p-file.
%
% (c) Kangrong Zhu      Stanford University     Jan 2014

if ~exist('nt_to_recon', 'var') || isempty(nt_to_recon)
    nt_to_recon = p.num_passes - p.mux_excited*p.cap_get_ecc - p.mux_encoded*p.num_mux_cycle;
end
if numel(nt_to_recon) ~= 1
    error('''nt_to_recon'' must be a scalar.');
end

p.nt_to_recon = nt_to_recon;

end_time_indx = p.mux_excited*p.cap_get_ecc + p.mux_encoded*p.num_mux_cycle + p.nt_to_recon;
if end_time_indx > p.num_passes
    error('Number of time points to reconstruct exceeds total number of time points in data.');
end
if end_time_indx <= 0
    error('Number of time points to reconstruct requires data collectd at negative time.');
end
p.tpoints_to_load = 1 : end_time_indx;



function sz = get_dat_sz(dat, p)
%
% function sz = get_dat_sz(dat, p)
%
% Returns the size of the input matrix 'dat' in the structure array 'sz'.
%
% Inputs
%   dat - Input data. Dim: [FE, PE, Echo, Slice, Coil, Time, Kz].
%   p   - Parameter structure. See mux_epi_params.m for details. The following
%         fields are used in this function: FE_DIM, PE_DIM, EC_DIM, SL_DIM, C_DIM, T_DIM, KZ_DIM.
%
% Output
%   sz  - A structure array with the following fields: x(FE), y(PE),
%         ec(Echo), sl(Slice), c(Coil), t(Time), kz(Kz).
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

sz = struct('x', {size(dat{1}, p.FE_DIM)}, 'y', {size(dat{1},p.PE_DIM)}, ...
    'ec', {size(dat{1}, p.EC_DIM)}, 'sl', {size(dat{1}, p.SL_DIM)}, ...
    'c', {size(dat{1}, p.C_DIM)}, 't', {numel(dat)}, ...
    'kz', {size(dat{1}, p.KZ_DIM)});



function [pha_coe, uid] = rawload_ref(ny, nsl, nc, frames, slices, coils, fname)
%
% function [pha_coe, uid] = rawload_ref(ny, nsl, nc[, frames, slices, coils, fname])
%
% Load the GE ref.dat file for phase correction in EPI.
%
% Inputs
%   ny     - Total # of ky lines in the acquisition corresponding to the
%            ref.dat file to be read.
%   nsl    - Total # of slices in the acquisition.
%   nc     - Total # of coils in  the acquisition.
%   frames - Frames we want the phase coefficients of. Default: all frames.
%   slices - Slices we want the phase coefficients of. Default: all slices.
%   coils  - Coils we want the phase coefficients of. Default: all coils.
%   fname  - File name. Default: 'ref.dat'.
%
% Output
%   pha_coe - Coefficients for x-ky phase correction. Dim: [2(0th order, 1st order), ny, nsl, nc].
%   uid     - the uinique id from the corresponding p-file (if available)
%             (will be empty if the file does not have a uid header)
%
% (c) Kangrong Zhu,     Stanford University     July 2012

% -- Parse inputs.
if ~exist('frames', 'var');      frames = [];            end;
if ~exist('slices', 'var');      slices = [];            end;
if ~exist('coils', 'var');       coils  = [];            end;
if ~exist('fname', 'var');       fname = 'ref.dat';      end;

if (isempty(frames));            frames = 1:ny;          end;
if (isempty(slices));            slices = 1:nsl;         end;
if (isempty(coils));             coils  = 1:nc;          end;

% -- Specify some constants. (Guessed from the data loaded from the ref.dat file.)
MAX_NUM_FRAMES = 512; % Maximum # of frames which could be stored in the ref.dat file.
                      % Only the first ny points in a 512 point block are non-zero.
PHA_ORDER = 1;        % Order for the phase term.

% -- Open file.
fip = fopen(fname,'r','l');
if fip == -1
    error('Can not open file %s.\n',fname);
end

% Attempt to read the 32-byte UID header
uid = fread(fip, 32, 'unsigned char')';
% Check for the magic 7 bytes that the 32-byte UIDs always begin with
if ~all(uid(1:7) == [43 59 149 27 34 71 42])
    uid = [];
    frewind(fip);
end

% -- Read data
nexpected = MAX_NUM_FRAMES * (PHA_ORDER+1) * nc * nsl; % Expected # of points
[pha_coe, nread] = fread(fip, nexpected , 'float32');
fclose(fip);
if nread == nexpected
    pha_coe = reshape(pha_coe, MAX_NUM_FRAMES, PHA_ORDER+1, nc, nsl);
    pha_coe = permute(pha_coe, [2, 1, 4, 3]);
    pha_coe = pha_coe( :, frames, slices, coils );
else
    error('%d of %d expected data points are read.', nread, nexpected);
end



function [omegaz, cal_blips] = encode_dftz_omegaz(fov_shift)
%
% function [omegaz, cal_blips] = encode_dftz_omegaz(fov_shift)
%
% Calculate Discrete Fourier Transform (DFT) encoding frequencies for the slice dimension.
%
% Input
%   fov_shift - CAIPI FOV shift. abs(fov_shift) is the length of the DFTz encoding. Negative for inverse DFTz encoding, positive for forward DFTz encoding.
%
% Outputs
%   omegaz    - DFTz encoding frequencies, in range [-pi, pi]. The order of omegaz corresponds to the order used in the mux cycling calibration time points. Dim: [1, abs(fov_shift)].
%   cal_blips - Blip indices for the calibration time points if the calibration were acquired using a CAIPI FOV shift of 'fov_shift'. Blip indices 0~(abs(fov_shift)-1) correspond to -kzmax~kzmax. cal_blips is [1,2,0] for abs(fov_shift)=3, [2,3,0,1] for abs(fov_shift)=4.
%
% (c) Kangrong Zhu  Stanford University     Feb 2014

blip_polarity = sign(fov_shift);                        % -1 for inverse DFTz encoding, +1 for forward DFTz encoding.
fov_shift = abs(fov_shift);
cap_blip_start_cal = floor(fov_shift/2);                % Staring blip index of the calibration data when the calibration were acquired using a CAIPI FOV shift of 'fov_shift'. Indices 0~(abs(fov_shift)-1) correspond to -kzmax~kzmax.
cal_blips = mod( (cap_blip_start_cal : 1 : cap_blip_start_cal+fov_shift-1), fov_shift);
cal_kz = blip_polarity * (cal_blips - (fov_shift-1)/2); % kz values. When blip_polarity = -1, this is [0, -1, 1] for abs(fov_shift) = 3, and [-0.5, -1.5, 1.5, 0.5] for abs(fov_shift) = 4.
omegaz = - 2*pi * cal_kz / fov_shift;



function ft_mtx = encode_ftyz_mtx(msk, ny, nz, add_vcc)
%
% function ft_mtx = encode_ftyz_mtx(msk, ny, nz, add_vcc)
% 
% Calculate Fourier transform (FT) encoding matrix for the y-z plane.
% In y: DFT encoding.
% In z: DFT encoding(CAIPI) or DTFT encoding(MICA).
%
% Inputs
%   msk    - Structure for sampling mask on the ky-omegaz plane. Fields 'add_vcc', 'kz',
%            'omegaz', 'ftz_pha (if not exist, will be calculated using 'kz' and 'omegaz')'
%            and  'ky (used when ny > 1)' are used in this function. See get_ky_omegaz_us_msk.m for details.
%   ny     - Full matrix size in y. If ny <= 1, the DFTy part will not be included in the FT encoding matrix.
%   nz     - Number of simultaneous slices.
%   add_vcc- 0: No virtual coils, use only actual coils; 1: Use both actual and virtual coils; 2: Use only virtual coils.
%
% Output
%   ft_mtx - The FT encoding matix for the y-z plane.
%            Dim: [sample(nsamp*(1+extra_dat_vcc)), Y->Z(=ny*nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))].
%            if add_vcc == 0: extra_dat_vcc = 0, Order in 1st dimension is sample ordering in actual coils.
%            if add_vcc == 1: extra_dat_vcc = 1, Order in 1st dimension is sample ordering in actual coils -> sample ordering in virtual coils.
%            if add_vcc == 2: extra_dat_vcc = 0, Order in 1st dimension is sample ordering in virtual coils.
%
% (c) Kangrong Zhu,     Stanford University     Oct 2013

KEEP_ORIG_SZ             = 1;
NEG_OMEGAZ_ENC           = 1;
SPECIFIED_OMEGAZ_ENC     = 0;
ONLY_ACTUAL_COILS        = 0;
ACTUAL_AND_VIRTUAL_COILS = 1;
ONLY_VIRTUAL_COILS       = 2;

if ~exist('add_vcc', 'var') || isempty(add_vcc)
    add_vcc = ONLY_ACTUAL_COILS;
end

if (ny > 1) && (length(msk(1).ky) ~= length(msk(1).kz))
    error('Number of samples mismatch between msk.ky and msk.kz.');
end
nsamp = length(msk(1).kz);                           % Number of samples on the ky-omegaz plane (i.e. number of echos)

% Include the FTz(DFTz or randomly sampled DTFTz) part in the FT encoding matrix
if ~isfield(msk, 'ftz_pha')
    msk = encode_ftz_pha(msk, nz, SPECIFIED_OMEGAZ_ENC); % msk.ftz_pha is the encoding phase added to each individual slice for each sample. Dim: [nsamp, nz(z indices (-floor(nz/2):1:(ceil(nz/2)-1))])
end
switch add_vcc
    case ONLY_ACTUAL_COILS
        extra_dat_vcc = 0;
    case ACTUAL_AND_VIRTUAL_COILS
        extra_dat_vcc = 1;
        tmp_msk = encode_ftz_pha(msk, nz, NEG_OMEGAZ_ENC);
        msk.ftz_pha = cat(1, msk.ftz_pha, tmp_msk.ftz_pha);
    case ONLY_VIRTUAL_COILS
        extra_dat_vcc = 0;
        msk = encode_ftz_pha(msk, nz, NEG_OMEGAZ_ENC);
end
msk.ftz_pha = ifftshift(msk.ftz_pha, 2);             % z indices become ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))
ftz_mtx = zeros(nsamp*(1+extra_dat_vcc), ny*nz);     % The FTz part in the encoding matrix. Dim: [sample(=nsamp*(1+extra_dat_vcc)), Y->Z(=ny*nz)(z indices are ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))]
for slice = 1 : nz
    ftz_mtx(:, (slice-1)*ny+1 : slice*ny) = repmat( msk.ftz_pha(:, slice)./sqrt(nz), [KEEP_ORIG_SZ, ny]); % Divide by sqrt(nz) to scale for orthogonal DFTz, when DFTz is used.
end

ft_mtx = ftz_mtx;

% If ny > 1, include the DFTy part in the FT encoding matrix
if ny > 1
    ky_indices = -ny/2 : 1 : ny/2-1;                 % Raw data matrix always arranged top->bottom(-ky->+ky) in ky no matter whether it was acquired top->bottom or bottom->top.
    ky_indices = ky_indices(msk.ky(:)).';
    y_indices = fftshift(0 : ny-1);                  % fftshift is essential because the 0-th indexed PE pixel is in the center along y.
    dfty_mtx = exp(-1i*2*pi*ky_indices*y_indices/ny) ./ sqrt(ny); % Encoding using orthogonal DFT, corresponding to fftc.m
    switch add_vcc
        case ONLY_ACTUAL_COILS
            % DO NOTHING
        case ACTUAL_AND_VIRTUAL_COILS
            dfty_mtx = cat(1, dfty_mtx, exp(-1i*2*pi*(-ky_indices)*y_indices/ny) ./ sqrt(ny));
        case ONLY_VIRTUAL_COILS
            dfty_mtx = exp(-1i*2*pi*(-ky_indices)*y_indices/ny) ./ sqrt(ny);
    end
    dfty_mtx = repmat(dfty_mtx, [KEEP_ORIG_SZ, nz]); % Dim: [nsamp*(1+extra_dat_vcc), ny*nz]
    
    ft_mtx = ft_mtx .* dfty_mtx;
end



function msk = encode_ftz_pha(msk, nz, neg)
%
% function msk = encode_ftz_pha(msk, nz, neg)
%
% Calculate the encoding phase added to each individual slice for each sample on the ky-omegaz plane.
%
% Inputs
%   msk - Structure for sampling mask on the ky-omegaz plane. Fields 'omegaz' and 'kz' are used in this function. See get_ky_omegaz_us_msk.m for details. 
%   nz  - Number of simultaneous slices.
%   neg - True: use negative value of specified frequencies for slice encoding. False: use specified frequencies.
%
% Output
%   msk - The input 'msk' with the field 'ftz_pha' added.
%         ftz_pha: The FTz encoding phase added to each individual slice for each sample on the ky-omegaz plane. Dim: [Sample(=nsamp), SimultaneousSliceZ(=nz, z indices (-floor(nz/2):1:ceil(nz/2)-1)].
%
% (c) Kangrong Zhu      Stanford University     Feb 2014

if ~exist('neg', 'var') || isempty(neg)
    neg = 0;
end
if neg
    neg = 1;                                   % Ensure a known number is used for power of -1.
else
    neg = 0;
end

nmsk = length(msk);
nsamp = length(msk(1).kz);
z_indices = -floor(nz/2) : 1 : (ceil(nz/2)-1); % Slice indices must center around index 0. e.g. [-1,0,1] for nz=3, [-2,-1,0,1] for nz=4. ([-1,0,1] and [2,0,1] are the same for CAIPI(DFTz encoding) when abs(cap_fov_shift)==nz, but different for MICA(Randomly sampled DTFTz encoding) or for CAIPI when abs(cap_fov_shift)~=nz)
for msk_idx = 1 : nmsk
    msk(msk_idx).ftz_pha = zeros(nsamp, nz);
    for slice = 1 : nz
        z = z_indices(slice);                  % Slice index
        msk(msk_idx).ftz_pha(:, slice) = exp(1i * ((-1)^neg) * msk(msk_idx).omegaz(msk(msk_idx).kz(:)) * z);
    end
end



function [ramp_flt, uid] = rawload_vrgf(nx_ramp, nx_pres, fname)
%
% function [ramp_flt, uid] = rawload_vrgf(nx_ramp, [nx_pres, fname])
%
% Load the GE vrgf.dat file for the ramp sampling filter in EPI.
%
% Inputs:
%   nx_ramp  - Size in FE, with ramp sampling on.
%   nx_pres  - Prescribed size in FE.
%   fname    - File name. Default: 'vrgf.dat'
%
% Output:
%   ramp_flt - The filter for correcting the ramp sampling, size: [nx_pres, nx_ramp].
%              ksp_corrected(nx_pres-by-ny) = ramp_flt * ksp(nx_ramp-by-ny).
%   uid      - the uinique id from the corresponding p-file (if available)
%              (will be empty if the file does not have a uid header)
%
% (c) Kangrong Zhu,     Stanford University     June 2012

if ~exist('nx_pres', 'var')
    nx_pres = [];
end
if ~exist('fname', 'var')
    fname = 'vrgf.dat';
end

% Open file
fip = fopen(fname, 'r', 'l');
if fip == -1
  error('File %s not found\n', fname);
end

% Attempt to read the 32-byte UID header
uid = fread(fip, 32, 'unsigned char')';
% Check for the magic 7 bytes that the 32-byte UIDs always begin with
if ~all(uid(1:7) == [43 59 149 27 34 71 42])
    uid = [];
    frewind(fip);
end

% Read data
if isempty(nx_pres)
    [ramp_flt, nread] = fread(fip, Inf, 'float32');
    nx_pres = floor(nread/nx_ramp);
    nexpected = nx_ramp * nx_pres;
else
    nexpected = nx_ramp * nx_pres;
    [ramp_flt, nread] = fread(fip, nexpected, 'float32');
end
fclose(fip);
if nread == nexpected
    ramp_flt = permute( reshape(ramp_flt, nx_ramp, nx_pres), [2,1] );
else
    error('%d of %d expected data points are read.', nread, nexpected);
end



function ff = gen_fermi_filter(np, fr, fw, ty)
%
% function ff = gen_fermi_filter(np, fr, fw [, ty])
%
% Generates Fermi filter.
%
% Inputs
%   np - Number of pixels, [np(1), np(2)(optional)].
%   fr - Fermi radius, in pixels.
%   fw - Fermi width, in pixels.
%   ty - Filter type to generate when ( numel(np)==2 && np(2)>1 ).
%        'circ'(default): Circularly symmetric.
%        'sepa':          Separable.
%
% Output
%   ff - Fermi filter.
%        If (numel(np)==1 || (numel(np)==2 && np(2)==1)), ff is a 1D Fermi filter of size [np(1), 1].
%        If (numel(np)==2 && np(2)>1 && strcmp(ty, 'circ')), ff is a 2D circularly symmetric Fermi filter of size [np(1), np(2)].
%        If (numel(np)==2 && np(2)>1 && strcmp(ty, 'sepa')), ff is a 2D separable Fermi filter of size [np(1), np(2)].
%
% (c) Kangrong Zhu,     Stanford University     Sep 2013

if ~exist('ty', 'var') || isempty(ty)
    ty = 'circ';
end

if (numel(np)==1) || ( (numel(np)==2) && (np(2)==1) ) % 1D Fermi filter
    ff = fermf_1d(np(1), fr, fw);
end

if (numel(np) == 2) && (np(2) > 1)                    % 2D Fermi filter
    switch ty
        case 'circ'                                   % Circularly symmetric filter
            ff = zeros(np(1), np(2));
            or = np/2 + 1;                            % Origin
            for x = 1 : np(1)
                for y = 1 : np(2)
                    d = sqrt((x-or(1))^2 + (y-or(2))^2);
                    ff(x, y) = ( 1 + exp((d-fr)/fw) ) .^ (-1);
                end
            end
        case 'sepa'                                   % Separable filter
            ff = fermf_1d(np(1), fr, fw);
            ff2 = fermf_1d(np(2), fr, fw);
            ff = repmat(ff, [1, np(2)]) .* repmat(ff2.', [np(1), 1]);
    end
end


function f = fermf_1d(n, r, w)
%
%   n - Number of pixels.
%   r - Fermi radius, in pixels.
%   w - Fermi width, in pixels.
%
%   f - 1D Fermi filter.
%

d = abs( -floor(n/2) : 1 : ceil(n/2-1) );
d = d(:);
f = ( 1 + exp((d-r)/w) ) .^ (-1);
