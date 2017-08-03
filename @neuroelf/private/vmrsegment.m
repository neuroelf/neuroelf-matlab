function vmr = vmrsegment(vmr, opts)
% vmrsegment  - apply SPM's segmentation to a VMR
%
% FORMAT:       svmr = vmrsegment(vmr [, opts])
%
% Input fields:
%
%       vmr         VMR to segment
%       opts        optional settings
%        .c1minr    minimum c1 ./ (c1+c2) ratio to be considered c1 (2/3)
%        .c1minv    minimum c1 value to be considered c1 (0.5)
%        .c2minr    minimum c2 ./ (c1+c2) ratio to be considered c2 (0.6)
%        .c2minv    minimum c2 value to be considered c2 (0.4)
%        .cleanup   set output.cleanup flag in SPM segmentation job
%        .gm        gray matter value in target VMR (default: 236)
%        .hiinterp  hires interpolation method (default: 'cubic')
%        .hires     perform segmentation in hires (0.5mm) resolution
%        .ihcfirst  first apply inhomogeneity correction (false)
%        .wm        white matter value in target VMR (default: 240)
%
% Output fields:
%
%       svmr        segmentated VMR

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:27 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, 2016, Jochen Weber
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

% argument check
vmro = {[]};
if nargin > 0 && ...
    ischar(vmr) && ...
   ~isempty(vmr) && ...
    exist(vmr(:)', 'file') == 2
    try
        vmro{1} = xff(vmr(:)', 'vmr');
        if ~isxff(vmro{1}, 'vmr')
            error( ...
                'neuroelf:BadArgument', ...
                'Not a valid VMR filename.' ...
            );
        end
    catch ne_eo;
        clearxffobjects(vmro);
        rethrow(ne_eo);
    end
    vmr = vmro{1};
end
if nargin < 1 || ...
   ~isxff(vmr, 'vmr')
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing VMR object.' ...
    );
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'c1minr') || ...
   ~isa(opts.c1minr, 'double') || ...
    numel(opts.c1minr) ~= 1 || ...
    isinf(opts.c1minr) || ...
    isnan(opts.c1minr) || ...
    opts.c1minr < 0.5 || ...
    opts.c1minr > 1
    opts.c1minr = 2 / 3;
end
if ~isfield(opts, 'c1minv') || ...
   ~isa(opts.c1minv, 'double') || ...
    numel(opts.c1minv) ~= 1 || ...
    isinf(opts.c1minv) || ...
    isnan(opts.c1minv) || ...
    opts.c1minv < 0.25 || ...
    opts.c1minv > 1
    opts.c1minv = 0.5;
end
if ~isfield(opts, 'c2minr') || ...
   ~isa(opts.c2minr, 'double') || ...
    numel(opts.c2minr) ~= 1 || ...
    isinf(opts.c2minr) || ...
    isnan(opts.c2minr) || ...
    opts.c2minr < 0.5 || ...
    opts.c2minr > 1
    opts.c2minr = 0.6;
end
if ~isfield(opts, 'c2minv') || ...
   ~isa(opts.c2minv, 'double') || ...
    numel(opts.c2minv) ~= 1 || ...
    isinf(opts.c2minv) || ...
    isnan(opts.c2minv) || ...
    opts.c2minv < 0.25 || ...
    opts.c2minv > 1
    opts.c2minv = 0.4;
end
if ~isfield(opts, 'cleanup') || ...
   (~isa(opts.cleanup, 'double') && ...
    ~islogical(opts.cleanup)) || ...
    numel(opts.cleanup) ~= 1 || ...
    isinf(opts.cleanup) || ...
    isnan(opts.cleanup) || ...
    ~any(double(opts.cleanup) == [0, 1])
    opts.cleanup = 0;
else
    opts.cleanup = double(opts.cleanup);
end
if ~isfield(opts, 'gm') || ...
   ~isnumeric(opts.gm) || ...
    numel(opts.gm) ~= 1 || ...
    isinf(opts.gm) || ...
    isnan(opts.gm) || ...
    opts.gm < 0 || ...
    opts.gm > 255
    opts.gm = 236;
end
if ~isfield(opts, 'hiinterp') || ...
   ~ischar(opts.hiinterp) || ...
   ~any(strcmpi(opts.hiinterp(:)', ...
        {'cubic', 'lanczos3', 'lanczos8', 'linear', 'nearest'}))
    opts.hiinterp = 'linear';
else
    opts.hiinterp = lower(opts.hiinterp(:)');
end
if ~isfield(opts, 'hires') || ...
   ~islogical(opts.hires) || ...
    numel(opts.hires) ~= 1
    opts.hires = false;
end
if ~isfield(opts, 'ihcfirst') || ...
   ~islogical(opts.ihcfirst) || ...
    numel(opts.ihcfirst) ~= 1
    opts.ihcfirst = false;
end
if ~isfield(opts, 'wm') || ...
   ~isnumeric(opts.wm) || ...
    numel(opts.wm) ~= 1 || ...
    isinf(opts.wm) || ...
    isnan(opts.wm) || ...
    opts.wm < 0 || ...
    opts.wm > 255
    opts.wm = 240;
end

% keep track of clearing flag
clearvmr = false;

% is SPM available
try
    sver = lower(spm('ver'));
    if ~any(strcmp(sver, {'spm5', 'spm8', 'spm12b', 'spm12'}))
        error( ...
            'neuroelf:SPMError', ...
            'Required SPM version not on the path.' ...
        );
    end
    sver = str2double(regexprep(sver, '^spm(\d+).*$', '$1'));
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    clearxffobjects(vmro);
    error( ...
        'neuroelf:SPMError', ...
        'SPM not found on the path.' ...
    );
end

% load appropriate template
try
    if sver == 5
        jobs = neuroelf_file('p', 'spm5_segment');
    else
        jobs = neuroelf_file('p', 'spm8_segment');
        if sver > 8
            jobs.jobs{1}.spm.tools.oldseg = jobs.jobs{1}.spm.spatial.preproc;
            jobs.jobs{1}.spm = rmfield(jobs.jobs{1}.spm, 'spatial');
        end
    end
    jobs = jobs.jobs;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    clearxffobjects(vmro);
    error( ...
        'neuroelf:PackageError', ...
        'Required coregistration job not found.' ...
    );
end

% IHC first
if opts.ihcfirst
    disp('Performing inhomogeneity correction (auto-mode)...');
    vmr.InhomogeneityCorrect;
end

% hires
if opts.hires
    try
        disp('Resampling to 0.5mm resolution...');
        vmr = vmr.HiResRescale(0.5, opts.hiinterp);
        clearvmr = true;
    catch ne_eo;
        clearxffobjects(vmro);
        rethrow(ne_eo);
    end
end

% progress in console
disp('Writing VMRs -> Analyze...');

% write both as temporary files
vmrimg = [tempname(tempdir) '.nii'];
[vmrdir, vmrfile] = fileparts(vmrimg);
vmr.ExportNifti(vmrimg);

% clear temporary object
clearxffobjects(vmro);

% get TPM files
tpmdir = [spm('dir') '/tpm/'];
if sver <= 8
    tpm = {[tpmdir 'grey.nii']; [tpmdir 'white.nii']; [tpmdir 'csf.nii']};
else
    tpm = {[tpmdir 'TPM.nii,1']; [tpmdir 'TPM.nii,2']; [tpmdir 'TPM.nii,3']};
end

% put files into structure
if sver == 5
    jobs{1}.spatial{1}.preproc.data = {[vmrimg ',1']};
    jobs{1}.spatial{1}.preproc.output.cleanup = opts.cleanup;
    jobs{1}.spatial{1}.preproc.opts.tpm = tpm;
    spm_defaults;
elseif sver == 8
    jobs{1}.spm.spatial.preproc.data = {[vmrimg ',1']};
    jobs{1}.spm.spatial.preproc.output.cleanup = opts.cleanup;
    jobs{1}.spm.spatial.preproc.opts.tpm = tpm;
    spm('defaults', 'FMRI');
else
    jobs{1}.spm.tools.oldseg.data = {[vmrimg ',1']};
    jobs{1}.spm.tools.oldseg.output.cleanup = opts.cleanup;
    jobs{1}.spm.tools.oldseg.opts.tpm = tpm;
    spm('defaults', 'FMRI');
end

% try to run job
disp('Running SPM segmentation...');
spm_jobman('run', jobs);

% delete actual (temporary) file
try
    delete(vmrimg);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

% remove superfluous files (loading the content)
if exist([vmrdir '/' vmrfile '_seg_sn.mat'], 'file') > 0
    segsn = load([vmrdir '/' vmrfile '_seg_sn.mat']);
    try
        delete([vmrdir '/' vmrfile '_seg_sn.mat']);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
else
    segsn = [];
end
if exist([vmrdir '/' vmrfile '_seg_inv_sn.mat'], 'file') > 0
    segisn = load([vmrdir '/' vmrfile '_seg_inv_sn.mat']);
    try
        delete([vmrdir '/' vmrfile '_seg_inv_sn.mat']);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
else
    segisn = [];
end
if exist([vmrdir '/m' vmrfile '.nii'], 'file') > 0
    try
        delete([vmrdir '/m' vmrfile '.nii']);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

% read output images
clearo = cell(1, 2);
try
    c1hdr = xff([vmrdir '/c1' vmrfile '.nii']);
    clearo{1} = c1hdr;
    c2hdr = xff([vmrdir '/c2' vmrfile '.nii']);
    clearo{2} = c2hdr;
catch ne_eo;
    if clearvmr
        clearxffobjects({vmr});
    end
    neuroelf_lasterr(ne_eo);
    clearxffobjects(clearo);
    error( ...
        'neuroelf:FileError', ...
        'Error opening header files after segmentation.' ...
    );
end

% load data
c1hdr.VoxelData = c1hdr.VoxelData(:, :, :);
c2hdr.VoxelData = c2hdr.VoxelData(:, :, :);

% delete files
try
    delete(c1hdr.FilenameOnDisk);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end
try
    delete(c2hdr.FilenameOnDisk);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

% copy object now
if ~clearvmr
    vmr = vmr.CopyObject;
end

% assign data
vmr.RunTimeVars.SegSN = segsn;
vmr.RunTimeVars.SegInvSN = segisn;

% get object data
c1d = c1hdr.VoxelData;
c2d = c2hdr.VoxelData;

% clear objects
clearxffobjects(clearo);

% perform computation
vd = vmr.VMRData;
if opts.gm > 0
    vmr.VMRData((c1d > opts.c1minv) & ((c1d ./ (c1d + c2d)) > opts.c1minr)) = opts.gm;
end
if opts.wm > 0
    vmr.VMRData((c2d > opts.c2minv) & ((c2d ./ (c1d + c2d)) > opts.c2minr)) = opts.wm;
end

% fill holes (in white matter only)
if opts.wm > 0
    c2d = (vmr.VMRData == opts.wm);
    [cs, cv] = clustercoordsc(c2d, 1, 9);
    if isempty(cs)
        vmr.ClearObject;
        error( ...
            'neuroelf:SegmentationFailed', ...
            'Segmentation failed.' ...
        );
    end
    c2d = (cv == maxpos(cs));
    c1d = (vmr.VMRData == opts.wm) & ~c2d;
    vmr.VMRData(c1d) = vd(c1d);
    [cs, cv] = clustercoordsc(~c2d, 1, 9);
    c1d = (cv ~= maxpos(cs));
    vmr.VMRData(c1d) = opts.wm;
end
